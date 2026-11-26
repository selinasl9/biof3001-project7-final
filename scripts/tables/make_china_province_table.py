#!/usr/bin/env python3
"""
make_china_province_table.py

Find a province for every China run with a host-cleaned hit, trying three routes:

    metadata_table   a province column we already have
    geo_loc_name     the BioSample's own "China: Yunnan, Mengla"
    lat_lon          reverse geocoding the coordinates, one lookup per distinct pair

Blanks are written rather than dropped, so a gap stays visible in the table.

Inputs:
    --hits        host-cleaned hit tables, to pick the China runs out of (or --runs)
    --metadata    run and sample metadata tables, for the BioSample ids
    --cache       lat_lon -> province cache, read and rewritten

Outputs:
    china_srr_provinces.tsv   sample, province, source, what_biosample_says
    the coordinate cache, rewritten with any new lookups
"""

import argparse
import csv
import os
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter

sys.path.insert(0, __file__.rsplit("/", 2)[0])
from common.geography import (NOMINATIM_PAUSE_SECONDS, normalise_province,
                              parse_lat_lon, reverse_geocode, split_geo_loc_name)

# Set the NCBI request settings, the column names to try + the group-size rule
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
RETRIES = 4
TIMEOUT_SECONDS = 120
PAUSE_SECONDS = 0.4

RUN_COLUMNS = ["acc", "run", "Run", "sample"]
LAT_LON_COLUMNS = ["lat_lon", "lat lon"]
# Where a BioSample can put the country and province
REGION_KEYS = ["geo_loc_name", "geographic location", "country", "geo loc name"]

# 04_common.R drops a group with fewer than this many runs, so flag them here
MIN_GROUP_RUNS = 5

COLUMNS = ["sample", "province", "source", "what_biosample_says"]


def read_tsv(path):
    """Yield the rows of a tab-separated file as dicts, skipping a missing file."""
    try:
        handle = open(path)
    except FileNotFoundError:
        print("missing table:", path, file=sys.stderr)
        return
    with handle:
        yield from csv.DictReader(handle, delimiter="\t")


def china_runs(args):
    """Return the China runs to resolve, from --runs or the host-cleaned hit tables."""
    if args.runs:
        return {line.strip() for line in open(args.runs) if line.strip()} - {""}

    runs = set()
    for path in args.hits:
        for row in read_tsv(path):
            if (row.get("country") or "").strip().lower() == "china":
                runs.add((row.get("sample") or "").strip())
    return runs - {""}


def read_metadata(paths, runs):
    """Return {run: biosample} + any province or coordinates the tables already carry."""
    biosample, province, coordinates = {}, {}, {}
    for path in paths:
        for row in read_tsv(path):
            run = next((row[column].strip() for column in RUN_COLUMNS
                        if row.get(column)), "")
            if run not in runs:
                continue

            accession = (row.get("biosample") or row.get("BioSample") or "").strip()
            if accession.startswith("SAM"):
                biosample.setdefault(run, accession)
            if (row.get("province") or "").strip():
                province.setdefault(run, row["province"].strip())
            value = next((row[column].strip() for column in LAT_LON_COLUMNS
                          if row.get(column)), "")
            if value:
                coordinates.setdefault(run, value)
    return biosample, province, coordinates


def efetch_biosamples(accessions, email):
    """Fetch {accession: {attribute: value}} for one batch of BioSample records."""
    url = EUTILS + "?" + urllib.parse.urlencode(
        {"db": "biosample", "id": ",".join(accessions), "retmode": "xml",
         "tool": "biof3001", "email": email})

    raw = b""
    for attempt in range(RETRIES):
        try:
            with urllib.request.urlopen(url, timeout=TIMEOUT_SECONDS) as response:
                raw = response.read()
            break
        except Exception as error:
            print(f"  retry {attempt + 1}: {error}", file=sys.stderr)
            time.sleep(3 * (attempt + 1))
    if not raw:
        return {}

    try:
        root = ET.fromstring(raw)
    except ET.ParseError as error:
        print("  parse error:", error, file=sys.stderr)
        return {}

    records = {}
    for sample in root.iter("BioSample"):
        attributes = {}
        for attribute in sample.iter("Attribute"):
            # Key on the attribute name, not its position, as submitters reorder them
            key = (attribute.get("harmonized_name")
                   or attribute.get("attribute_name") or "").lower()
            if key and attribute.text:
                attributes.setdefault(key, attribute.text.strip())
        records[sample.get("accession")] = attributes
    return records


def fetch_attributes(accessions, args):
    """Fetch the BioSample attributes for every accession we still need, in batches."""
    needed = sorted(accessions)
    if not needed:
        return {}
    print(f"fetching {len(needed)} BioSample records "
          f"(about {len(needed) / args.batch * 2 / 60:.0f} min)")

    attributes = {}
    for start in range(0, len(needed), args.batch):
        attributes.update(efetch_biosamples(needed[start:start + args.batch], args.email))
        time.sleep(PAUSE_SECONDS)
        print(f"  {min(start + args.batch, len(needed))}/{len(needed)} biosamples")
    return attributes


def region_text(attributes):
    """Return whatever the BioSample says about where the sample came from."""
    return next((attributes[key] for key in REGION_KEYS if attributes.get(key)), "")


def read_cache(paths):
    """Return {lat_lon: province} from every cache file that exists."""
    cache = {}
    for path in paths:
        for row in read_tsv(path) if os.path.exists(path) else []:
            if row.get("lat_lon") and row.get("province"):
                cache.setdefault(row["lat_lon"], row["province"])
    return cache


def write_cache(path, cache):
    """Write the coordinate cache back out, sorted."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["lat_lon", "province"])
        for value, province in sorted(cache.items()):
            writer.writerow([value, province])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hits", nargs="*", default=[],
                        help="host-cleaned hit tables, to pick the China runs out of")
    parser.add_argument("--runs", help="one run accession per line, instead of --hits")
    parser.add_argument("--metadata", nargs="+", required=True)
    parser.add_argument("-o", "--out", required=True)
    parser.add_argument("--cache", required=True,
                        help="lat_lon -> province cache, read and rewritten")
    parser.add_argument("--extra-cache", nargs="*", default=[],
                        help="older cache files to read as well")
    parser.add_argument("--email", default="")
    parser.add_argument("--batch", type=int, default=100)
    args = parser.parse_args()

    # Step 1: list the China runs that need a province
    runs = china_runs(args)
    print(len(runs), "China runs with at least one host-cleaned hit")

    # Step 2: read what the metadata gives us, BioSample id, province + coordinates
    biosample, known, coordinates = read_metadata(args.metadata, runs)
    print(f"{len(biosample)} runs have a BioSample id, "
          f"{len(known)} already have a province")

    # Step 3: fetch the BioSample records for the runs with no province yet
    attributes = fetch_attributes(
        {biosample[run] for run in runs if run in biosample and run not in known}, args)

    # Step 4: take the first route, a province column we have or geo_loc_name
    province, source, evidence = {}, {}, {}
    for run in sorted(runs):
        sample = attributes.get(biosample.get(run, ""), {})
        region = region_text(sample)
        position = sample.get("lat_lon") or coordinates.get(run, "")
        evidence[run] = f"geo_loc_name={region or '-'} lat_lon={position or '-'}"

        if run in known:
            province[run], source[run] = normalise_province(known[run]), "metadata_table"
            continue
        name = normalise_province(split_geo_loc_name(region)[0])
        if name:
            province[run], source[run] = name, "geo_loc_name"

    # Step 5: take the second route, reverse geocoding one lookup per coordinate pair
    cache = read_cache([args.cache] + args.extra_cache)
    print(len(cache), "coordinates already in the cache")

    positions = {}
    for run in sorted(runs):
        if run in province:
            continue
        value = attributes.get(biosample.get(run, ""), {}).get("lat_lon") or coordinates.get(run, "")
        if value and parse_lat_lon(value):
            positions[run] = value

    wanted = sorted(set(positions.values()) - set(cache))
    if wanted:
        print(f"looking up {len(wanted)} coordinates "
              f"(about {len(wanted) * NOMINATIM_PAUSE_SECONDS / 60:.0f} min)")
        for index, value in enumerate(wanted, 1):
            latitude, longitude = parse_lat_lon(value)
            name = reverse_geocode(latitude, longitude).get("province", "")
            if name:
                cache[value] = name
            time.sleep(NOMINATIM_PAUSE_SECONDS)
            if index % 10 == 0:
                print(f"  {index}/{len(wanted)}")
        write_cache(args.cache, cache)

    for run, value in positions.items():
        name = normalise_province(cache.get(value, ""))
        if name:
            province[run], source[run] = name, "lat_lon"

    # Step 6: write every run, blanks included, so the gaps stay visible
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(COLUMNS)
        for run in sorted(runs):
            writer.writerow([run, province.get(run, ""), source.get(run, ""), evidence[run]])

    # Step 7: return summary in the terminal, flagging the provinces 04_common.R drops
    counts = Counter(province.values())
    print(f"\n{len(province)}/{len(runs)} China runs have a province")
    print("by source:", Counter(source.values()).most_common())
    print("runs per province:")
    for name, count in counts.most_common():
        print(f"  {name:<20} {count}")

    small = sorted(name for name, count in counts.items() if count < MIN_GROUP_RUNS)
    if small:
        print(f"\n{len(small)} provinces have fewer than {MIN_GROUP_RUNS} runs and will be "
              f"dropped by 04_common.R (MIN_GROUP_N): " + ", ".join(small))
    print("\nwritten to", args.out)


if __name__ == "__main__":
    main()

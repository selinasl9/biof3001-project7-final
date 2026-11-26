#!/usr/bin/env python3
"""
summarise_bioproject.py

Build the summary tables for one BioProject run on its own, from its DIAMOND output.

Inputs:
    data/runs/logan_matches_rna_stats.tsv       run metadata
    results/<project>/runs.txt                  the runs that were screened
    results/<project>/diamond/<run>.blastx.tsv  DIAMOND output, one file per run
    results/<project>/biosample_attributes.tsv  cached, fetched from NCBI if absent
    results/<project>/latlon_locations.tsv      cached, geocoded if absent
    data/reference/taxonomy/{nodes,names}.dmp   NCBI taxonomy

Outputs, all under results/<project>/summary/ with a <CC>_<project>_ prefix:
    all_best_hits_with_metadata.tsv   sample_metadata.tsv
    sample_virus_contig_counts.tsv    sample_virus_presence_absence.tsv
    sample_family_presence_absence.tsv  sample_genus_presence_absence.tsv
    zero_hit_samples.tsv              no_logan_contigs_samples.tsv
"""

import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd

sys.path.append(str(Path(__file__).resolve().parents[1]))

from common.geography import parse_lat_lon, reverse_geocode, split_geo_loc_name

# Set the input paths, the request settings + the output column order
RUN_METADATA = "data/runs/logan_matches_rna_stats.tsv"
TAXONOMY_DIR = Path("data/reference/taxonomy")

BIOSAMPLE_BATCH = 200
BIOSAMPLE_RETRIES = 3
NCBI_PAUSE_SECONDS = 0.4
NOMINATIM_PAUSE_SECONDS = 1.1

# DIAMOND was run with these columns, in this order
DIAMOND_COLUMNS = ["qseqid", "sseqid", "pident", "length", "qlen", "slen",
                   "qstart", "qend", "sstart", "send", "evalue", "bitscore",
                   "qcovhsp", "scovhsp", "staxids", "sscinames", "stitle"]

NUMERIC_COLUMNS = ["pident", "length", "qlen", "slen", "evalue",
                   "bitscore", "qcovhsp", "scovhsp"]

METADATA_COLUMNS = [
    "sample", "bioproject", "country", "province", "city_latlon", "district_latlon",
    "site", "province_source", "geo_loc_name", "lat_lon", "host_organism",
    "host_biosample", "biosample", "collection_date_sam", "collection_date_biosample",
    "seqstats_contigs_n50", "seqstats_contigs_sumlen",
    "seqstats_contigs_nbseq", "seqstats_contigs_maxlen",
]

# Location first, then the hit, then the numbers, then the raw fields
COLUMN_ORDER = METADATA_COLUMNS[:15] + [
    "qseqid", "sseqid",
    "sscinames", "family", "genus", "species", "full_lineage",
    "pident", "length", "qlen", "slen", "qcovhsp", "scovhsp", "evalue", "bitscore",
    "staxids", "stitle",
    "seqstats_contigs_n50", "seqstats_contigs_sumlen",
    "seqstats_contigs_nbseq", "seqstats_contigs_maxlen",
    "qstart", "qend", "sstart", "send",
]

# These three describe where the province came from, so they stay in sample_metadata
METADATA_ONLY_COLUMNS = ["geo_loc_name", "site", "province_source"]


def fetch_biosamples(accessions):
    """Fetch every attribute of each BioSample, straight from NCBI."""
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    rows = []

    for start in range(0, len(accessions), BIOSAMPLE_BATCH):
        batch = accessions[start:start + BIOSAMPLE_BATCH]
        data = urllib.parse.urlencode(
            {"db": "biosample", "id": ",".join(batch), "retmode": "xml"}).encode()

        for attempt in range(BIOSAMPLE_RETRIES):
            try:
                with urllib.request.urlopen(url, data=data, timeout=120) as response:
                    xml = response.read()
                break
            except Exception as error:
                print(f"  NCBI request failed ({error}), retrying...")
                time.sleep(5)
        else:
            sys.exit("Could not reach NCBI; check the server's internet access")

        for biosample in ET.fromstring(xml).iter("BioSample"):
            row = {"biosample": biosample.get("accession")}
            for attribute in biosample.iter("Attribute"):
                name = attribute.get("harmonized_name") or attribute.get("attribute_name")
                row[name] = (attribute.text or "").strip()
            rows.append(row)

        print(f"  BioSamples fetched: {min(start + BIOSAMPLE_BATCH, len(accessions))} "
              f"/ {len(accessions)}")
        time.sleep(NCBI_PAUSE_SECONDS)

    return pd.DataFrame(rows).fillna("")


def load_biosample_attributes(path, accessions):
    """Return the cached BioSample table, fetching and saving it the first time."""
    if path.exists():
        print("BioSample attributes loaded from", path)
        return pd.read_csv(path, sep="\t", dtype=str).fillna("")

    print("Fetching BioSample attributes from NCBI...")
    biosample = fetch_biosamples(accessions)
    biosample.to_csv(path, sep="\t", index=False)
    return biosample


def geocode_missing(latlon_path, wanted):
    """Return the cached coordinate table, looking up anything new one at a time."""
    if latlon_path.exists():
        places = pd.read_csv(latlon_path, sep="\t", dtype=str).fillna("")
    else:
        places = pd.DataFrame(columns=["lat_lon", "province", "city", "district"])

    todo = sorted(set(wanted) - set(places["lat_lon"]))
    if not todo:
        return places

    print(f"Looking up {len(todo)} unique coordinates "
          f"(about {len(todo) * NOMINATIM_PAUSE_SECONDS / 60:.0f} min)...")

    new_rows = []
    for i, value in enumerate(todo, start=1):
        latitude, longitude = parse_lat_lon(value)
        place = reverse_geocode(latitude, longitude)
        # Save a successful lookup only, so a failure is retried on the next run
        if place:
            new_rows.append({"lat_lon": value, **place})
        time.sleep(NOMINATIM_PAUSE_SECONDS)
        if i % 20 == 0:
            print(f"  {i} / {len(todo)}")

    places = pd.concat([places, pd.DataFrame(new_rows)], ignore_index=True).fillna("")
    places.to_csv(latlon_path, sep="\t", index=False)

    if len(new_rows) < len(todo):
        print(f"WARNING: {len(todo) - len(new_rows)} coordinates could not be "
              f"looked up; re-run to retry")

    return places


def add_location(biosample, latlon_path):
    """Add the province, site + finer location, from geo_loc_name then the coordinates."""
    for column in ["geo_loc_name", "lat_lon", "host", "collection_date"]:
        if column not in biosample.columns:
            biosample[column] = ""

    # Take the first choice, where the submitter wrote the province into geo_loc_name
    split = biosample["geo_loc_name"].apply(split_geo_loc_name)
    biosample["province"] = [province for province, _ in split]
    biosample["site"] = [site for _, site in split]
    biosample["province_source"] = biosample["province"].apply(
        lambda province: "geo_loc_name" if province else "")

    # Take the second choice, where only "China" was given but there are coordinates
    coordinates = [value for value in biosample["lat_lon"] if parse_lat_lon(value)]
    places = geocode_missing(latlon_path, coordinates).rename(
        columns={"province": "province_latlon", "city": "city_latlon",
                 "district": "district_latlon"})

    biosample = biosample.merge(places, on="lat_lon", how="left").fillna("")

    from_latlon = (biosample["province"] == "") & (biosample["province_latlon"] != "")
    biosample.loc[from_latlon, "province"] = biosample.loc[from_latlon, "province_latlon"]
    biosample.loc[from_latlon, "province_source"] = "lat_lon"

    return biosample.rename(columns={"host": "host_biosample",
                                     "collection_date": "collection_date_biosample"})


def load_taxonomy():
    """Read {taxid: parent, rank, scientific name} from the NCBI taxonomy dump."""
    print("Loading NCBI taxonomy...")

    names = {}
    with open(TAXONOMY_DIR / "names.dmp") as handle:
        for line in handle:
            parts = [field.strip() for field in line.split("|")]
            if len(parts) >= 4 and parts[3] == "scientific name":
                names[parts[0]] = parts[1]

    parent, rank = {}, {}
    with open(TAXONOMY_DIR / "nodes.dmp") as handle:
        for line in handle:
            parts = [field.strip() for field in line.split("|")]
            if len(parts) >= 3:
                parent[parts[0]], rank[parts[0]] = parts[1], parts[2]

    return names, parent, rank


def lineage_of(taxid_value, names, parent, rank):
    """Return the family, genus, species + full lineage for one DIAMOND staxids value."""
    if pd.isna(taxid_value):
        return "", "", "", ""

    # Take the first taxid we know, as DIAMOND sometimes gives several
    taxid = next((candidate.strip() for candidate in str(taxid_value).split(";")
                  if candidate.strip() in parent), None)
    if taxid is None:
        return "", "", "", ""

    lineage, ranked, seen = [], {}, set()
    while taxid in parent and taxid not in seen:
        seen.add(taxid)
        lineage.append(names.get(taxid, taxid))
        if rank.get(taxid, "no rank") != "no rank":
            ranked[rank[taxid]] = names.get(taxid, taxid)
        if taxid == "1" or parent[taxid] == taxid:
            break
        taxid = parent[taxid]

    lineage.reverse()
    return (ranked.get("family", ""), ranked.get("genus", ""),
            ranked.get("species", ""), "; ".join(lineage))


def read_best_hits(diamond_dir, samples):
    """Read one best hit per contig for every run that has DIAMOND output."""
    best_hits, zero_hit, missing = [], [], []

    print("Reading DIAMOND outputs...")
    for i, sample in enumerate(samples, start=1):
        path = diamond_dir / f"{sample}.blastx.tsv"

        # Treat a missing file as Logan having no contigs for this run
        if not path.exists():
            missing.append(sample)
            continue

        # Treat an empty file as DIAMOND having finished and found nothing
        if path.stat().st_size == 0:
            zero_hit.append(sample)
            continue

        try:
            hits = pd.read_csv(path, sep="\t", dtype={"staxids": str})
            hits.columns = [str(column).lstrip("#") for column in hits.columns]
            if "qseqid" not in hits.columns:
                hits = pd.read_csv(path, sep="\t", comment="#", header=None,
                                   names=DIAMOND_COLUMNS, dtype={"staxids": str})
        except pd.errors.EmptyDataError:
            zero_hit.append(sample)
            continue

        if len(hits) == 0:
            zero_hit.append(sample)
            continue

        for column in NUMERIC_COLUMNS:
            hits[column] = pd.to_numeric(hits[column], errors="coerce")

        # Take the best hit, highest bitscore then lowest E-value then highest coverage
        hits = hits.sort_values(["qseqid", "bitscore", "evalue", "qcovhsp"],
                                ascending=[True, False, True, False])
        best = hits.drop_duplicates("qseqid").copy()
        best["sample"] = sample
        best_hits.append(best)

        if i % 200 == 0:
            print(f"Processed {i} / {len(samples)} samples")

    return pd.concat(best_hits, ignore_index=True), zero_hit, missing


def presence_matrix(best_all, samples, column):
    """Build a runs x values count matrix, zero-hit runs added back as rows of zeros."""
    hits = best_all[best_all[column].notna() & (best_all[column] != "")]
    counts = pd.crosstab(hits["sample"], hits[column]).reindex(samples, fill_value=0)
    counts.index.name = "sample"
    return counts


def main():
    project = sys.argv[1] if len(sys.argv) > 1 else "PRJNA1003643"
    country_code = sys.argv[2] if len(sys.argv) > 2 else "CN"
    prefix = f"{country_code}_{project}_"

    project_dir = Path("results") / project
    outdir = project_dir / "summary"
    outdir.mkdir(parents=True, exist_ok=True)

    # Step 1: read the runs that were screened + their Logan metadata
    metadata = pd.read_csv(RUN_METADATA, sep="\t", dtype=str).rename(
        columns={"acc": "sample", "geo_loc_name_country_calc": "country",
                 "organism": "host_organism"})

    samples = [line.strip() for line in open(project_dir / "runs.txt") if line.strip()]
    metadata = metadata[metadata["sample"].isin(samples)].drop_duplicates("sample").copy()
    print("Analysed samples in metadata:", len(metadata))

    # Step 2: work out where each run came from, finer than Logan's country
    biosample = load_biosample_attributes(
        project_dir / "biosample_attributes.tsv",
        sorted(metadata["biosample"].dropna().unique()))
    biosample = add_location(biosample, project_dir / "latlon_locations.tsv")

    metadata = metadata.merge(
        biosample[["biosample", "province", "city_latlon", "district_latlon", "site",
                   "province_source", "geo_loc_name", "lat_lon",
                   "host_biosample", "collection_date_biosample"]],
        on="biosample", how="left")
    print("Samples with a province:", (metadata["province"].fillna("") != "").sum())

    # Step 3: read the best hit per contig
    best_all, zero_hit, missing = read_best_hits(project_dir / "diamond", samples)
    samples = [sample for sample in samples if sample not in set(missing)]

    print()
    print("Runs without Logan contigs:", len(missing))
    print("Samples with zero hits:", len(zero_hit))
    print("Best-hit rows:", len(best_all))

    # Step 4: add the taxonomy for each hit
    print("Adding taxonomy...")
    names, parent, rank = load_taxonomy()
    taxonomy = best_all["staxids"].apply(lambda value: lineage_of(value, names, parent, rank))
    best_all[["family", "genus", "species", "full_lineage"]] = pd.DataFrame(
        taxonomy.tolist(), index=best_all.index)

    # Step 5: join the metadata on + put the columns in reading order
    metadata_small = metadata[[c for c in METADATA_COLUMNS if c in metadata.columns]].copy()
    best_all = best_all.merge(metadata_small, on="sample", how="left")

    ordered = [c for c in COLUMN_ORDER if c in best_all.columns]
    best_all = best_all[ordered + [c for c in best_all.columns if c not in ordered]]

    # Step 6: write the tables
    def save(frame, name, index=False):
        frame.to_csv(outdir / f"{prefix}{name}", sep="\t", index=index)

    save(best_all.drop(columns=METADATA_ONLY_COLUMNS, errors="ignore"),
         "all_best_hits_with_metadata.tsv")

    sample_metadata = metadata_small[metadata_small["sample"].isin(samples)]
    save(sample_metadata.sort_values("sample"), "sample_metadata.tsv")

    counts = presence_matrix(best_all, samples, "sscinames")
    save(counts, "sample_virus_contig_counts.tsv", index=True)
    save((counts > 0).astype(int), "sample_virus_presence_absence.tsv", index=True)
    save((presence_matrix(best_all, samples, "family") > 0).astype(int),
         "sample_family_presence_absence.tsv", index=True)
    save((presence_matrix(best_all, samples, "genus") > 0).astype(int),
         "sample_genus_presence_absence.tsv", index=True)

    save(pd.DataFrame({"sample": zero_hit}), "zero_hit_samples.tsv")
    save(pd.DataFrame({"sample": missing}), "no_logan_contigs_samples.tsv")

    print()
    print("DIAMOND samples analysed:", len(samples))
    print("Samples with zero hits:", len(zero_hit))
    print("Best-hit contigs:", len(best_all))
    print("Unique scientific names:", best_all["sscinames"].nunique())
    print("Unique families:", best_all["family"].replace("", pd.NA).nunique())
    print("Unique genera:", best_all["genus"].replace("", pd.NA).nunique())

    print()
    print("Samples per province:")
    print(sample_metadata["province"].replace("", "(not given)").value_counts().to_string())

    print()
    print("Files created in", outdir)
    for path in sorted(outdir.glob(f"{prefix}*")):
        print(" ", path.name)


if __name__ == "__main__":
    main()

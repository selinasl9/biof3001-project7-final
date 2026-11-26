#!/usr/bin/env python3
"""
lookup_host_species.py

Find which mosquito each run came from, trying four sources in order of trust:

    lib-meta    the source study's own COI-confirmed species table (PRJNA1003643)
    metadata    any column of the metadata tables we already downloaded
    sra         the SRA record: scientific name, title, design, sample attributes
    bioproject  the BioProject blurb, or a two-letter code in a sample name

The last is a guess about a whole BioProject, so those rows are marked "check by hand".

Inputs:
    --annotation  contig annotation table, read for its run list (or --runs)
    --metadata    run and sample metadata tables we already downloaded
    --lib-meta    the source study's library table (lib_id, species, genus, province)

Outputs:
    host_lookup.tsv   run, bioproject, host_species, host_genus, host_group,
                      province, confidence, source_field, evidence
"""

import argparse
import csv
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter

sys.path.insert(0, __file__.rsplit("/", 2)[0])
from common.geography import normalise_province
from common.mosquito import find_species, genus_from_code, genus_and_group

# Set the source order, the NCBI request settings + the output columns
# NCBI allows three requests a second without a key, so stay well under it
SOURCES = ["lib-meta", "metadata", "sra", "bioproject"]
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
RETRIES = 4
TIMEOUT_SECONDS = 120
PAUSE_SECONDS = 0.4

# The SRA fields worth reading, best first
SRA_FIELDS = ["SCIENTIFIC_NAME", "TITLE", "DESIGN_DESCRIPTION", "LIBRARY_NAME", "ALIAS"]

# A run column is called something different in each metadata table
RUN_COLUMNS = ["acc", "sample", "run", "Run"]

# The source study's library ids look like 19BJKW1832
LIBRARY_ID = re.compile(r"\b(\d{2}[A-Z]{2,6}[A-Z0-9]*\d{2,6})\b")

COLUMNS = ["run", "bioproject", "host_species", "host_genus", "host_group",
           "province", "confidence", "source_field", "evidence"]
EVIDENCE_CHARACTERS = 120


def read_tsv(path):
    """Yield the rows of a tab-separated file as dicts."""
    with open(path) as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def runs_to_look_up(args):
    """Return the runs we need a host for + the BioProject each belongs to."""
    if args.runs:
        return [line.strip() for line in open(args.runs) if line.strip()], {}

    runs, projects = [], {}
    for row in read_tsv(args.annotation):
        run = row.get("run") or row["contig"].split("_")[0]
        if run not in projects:
            runs.append(run)
        projects.setdefault(run, row.get("bioproject", ""))
    return runs, projects


def eutil(cgi, params, email):
    """Make one Entrez request, retried with a growing pause, empty bytes if it fails."""
    url = f"{EUTILS}/{cgi}?" + urllib.parse.urlencode(
        dict(params, tool="biof3001", email=email))

    for attempt in range(RETRIES):
        try:
            with urllib.request.urlopen(url, timeout=TIMEOUT_SECONDS) as response:
                return response.read()
        except Exception as error:
            print(f"  retry {attempt + 1}: {error}", file=sys.stderr)
            time.sleep(3 * (attempt + 1))
    return b""


def parse_xml(raw):
    """Return the root element, or None if NCBI sent something that will not parse."""
    if not raw:
        return None
    try:
        return ET.fromstring(raw)
    except ET.ParseError as error:
        print(f"  could not parse the reply: {error}", file=sys.stderr)
        return None


def from_lib_meta(runs, found, args):
    """Take the source study's COI-confirmed species, joined on its library id."""
    if not args.lib_meta:
        return

    try:
        libraries = {row["lib_id"]: row for row in csv.DictReader(open(args.lib_meta))}
    except FileNotFoundError:
        print(f"no lib_meta table at {args.lib_meta}, skipping that source")
        return
    print(f"{len(libraries)} libraries in {args.lib_meta}")

    for row in metadata_rows(runs, found, args):
        run, column, value = row
        library = value if value in libraries else next(
            (m.group(1) for m in LIBRARY_ID.finditer(value) if m.group(1) in libraries), None)
        if not library:
            continue

        record = libraries[library]
        # Take the province too, as this table is the only place some runs have one
        found[run] = {"species": record["species"], "confidence": "good",
                      "province": normalise_province(record.get("province", "")),
                      "field": f"lib_meta:{column}",
                      "evidence": f"{library}, COI-confirmed in the source study"}


def metadata_rows(runs, found, args):
    """Yield (run, column, value) for every non-empty cell of a run we still need."""
    wanted = set(runs)
    for path in args.metadata:
        try:
            handle = open(path)
        except FileNotFoundError:
            continue
        with handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                run = next((row[column].strip() for column in RUN_COLUMNS
                            if row.get(column)), "")
                if run not in wanted or run in found:
                    continue
                for column, value in row.items():
                    if value and value.strip():
                        yield run, column, value.strip()


def from_metadata(runs, found, args):
    """Look for a mosquito name in any column of the metadata tables we already have."""
    before = len(found)
    for run, column, value in metadata_rows(runs, found, args):
        species = find_species(value)
        if species:
            found[run] = {"species": species, "confidence": "good",
                          "field": f"metadata:{column}",
                          "evidence": value[:EVIDENCE_CHARACTERS]}
    print(f"{len(found) - before} runs resolved from the local metadata tables")


def from_sra(runs, found, args):
    """Read the SRA record for each remaining run, in batches."""
    todo = [run for run in runs if run not in found]
    if not todo:
        return
    print(f"asking SRA about {len(todo)} runs, {args.batch} at a time")

    for start in range(0, len(todo), args.batch):
        chunk = todo[start:start + args.batch]
        root = parse_xml(eutil("efetch.fcgi",
                               {"db": "sra", "id": ",".join(chunk),
                                "rettype": "xml", "retmode": "xml"}, args.email))
        if root is not None:
            for package in root.iter("EXPERIMENT_PACKAGE"):
                read_package(package, chunk, found)
        time.sleep(PAUSE_SECONDS)
        print(f"  {min(start + args.batch, len(todo))}/{len(todo)} done, "
              f"{len(found)} runs with a host")


def read_package(package, chunk, found):
    """Take the first field of one SRA experiment package that names a mosquito."""
    mine = [element.get("accession") for element in package.iter("RUN")
            if element.get("accession") in chunk]
    if not mine:
        return

    fields = [(tag, element.text) for tag in SRA_FIELDS
              for element in package.iter(tag) if element.text]
    for attribute in package.iter("SAMPLE_ATTRIBUTE"):
        fields.append((attribute.findtext("TAG") or "attribute",
                       attribute.findtext("VALUE")))
    for sample in package.iter("SAMPLE"):
        if sample.get("alias"):
            fields.append(("alias", sample.get("alias")))

    for tag, value in fields:
        species = find_species(value)
        if species:
            for run in mine:
                found.setdefault(run, {"species": species, "confidence": "good",
                                       "field": f"sra:{tag}",
                                       "evidence": (value or "")[:EVIDENCE_CHARACTERS]})
            return

    # Fall back to a code at the front of a name, as no genus is spelled out anywhere
    for tag, value in fields:
        genus = genus_from_code(value)
        if genus:
            for run in mine:
                found.setdefault(run, {"species": genus, "confidence": "check by hand",
                                       "field": f"sra:{tag} (name code)",
                                       "evidence": (value or "")[:EVIDENCE_CHARACTERS]})
            return


def from_bioproject(runs, found, projects, args):
    """Take a mosquito named in the BioProject blurb, for every run of that project."""
    todo = [run for run in runs if run not in found]
    accessions = sorted({projects.get(run, "") for run in todo} - {""})
    if not accessions:
        return
    print(f"checking {len(accessions)} BioProject records")

    for accession in accessions:
        root = parse_xml(eutil("esearch.fcgi",
                               {"db": "bioproject", "term": accession}, args.email))
        uid = root.findtext(".//Id") if root is not None else None
        if not uid:
            continue

        raw = eutil("efetch.fcgi",
                    {"db": "bioproject", "id": uid, "retmode": "xml"}, args.email)
        # Join every text field, as the useful text is spread over title, description + summary
        text = " ".join(match.strip() for match in
                        re.findall(r">([^<>]{3,})<", raw.decode("utf8", "replace")))
        species = find_species(text)
        print(f"  {accession}: {species or 'no mosquito named'}")

        if species:
            for run in todo:
                if projects.get(run) == accession:
                    found.setdefault(run, {
                        "species": species, "confidence": "check by hand",
                        "field": f"bioproject:{accession}",
                        "evidence": text[:EVIDENCE_CHARACTERS]})
        time.sleep(PAUSE_SECONDS)


def write_table(path, runs, projects, found):
    """Write one row per run, including the runs nothing could name."""
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(COLUMNS)
        for run in runs:
            hit = found.get(run, {})
            species = hit.get("species", "")
            genus, group = genus_and_group(species)
            writer.writerow([run, projects.get(run, ""), species, genus, group,
                             hit.get("province", ""), hit.get("confidence", ""),
                             hit.get("field", ""), hit.get("evidence", "")])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotation", help="contig annotation table, for the run list")
    parser.add_argument("--runs", help="one run accession per line, instead of --annotation")
    parser.add_argument("--metadata", nargs="*", default=[],
                        help="run and sample metadata tables")
    parser.add_argument("--lib-meta", default="",
                        help="the source study's library table (lib_id, species, genus)")
    parser.add_argument("-o", "--out", required=True)
    parser.add_argument("--source", nargs="+", choices=SOURCES + ["all"], default=["all"],
                        help="which sources to try, in the order given")
    parser.add_argument("--email", default="", help="NCBI asks for one, it raises the rate limit")
    parser.add_argument("--batch", type=int, default=50)
    args = parser.parse_args()

    if not (args.annotation or args.runs):
        parser.error("give either --annotation or --runs")
    sources = SOURCES if "all" in args.source else args.source

    # Step 1: list the runs we need a host for + the BioProject each belongs to
    runs, projects = runs_to_look_up(args)
    print(f"{len(runs)} runs to look up, trying {', '.join(sources)}")

    # Step 2: try each source in turn, each one filling only the runs still missing
    found = {}
    for source in sources:
        if source == "lib-meta":
            from_lib_meta(runs, found, args)
        elif source == "metadata":
            from_metadata(runs, found, args)
        elif source == "sra":
            from_sra(runs, found, args)
        elif source == "bioproject":
            from_bioproject(runs, found, projects, args)

    # Step 3: write every run, including the ones nothing could name
    write_table(args.out, runs, projects, found)

    # Step 4: return summary in the terminal, records and guesses counted apart
    good = sum(1 for hit in found.values() if hit["confidence"] == "good")
    print(f"\n{good} runs have a host from a record for that run, "
          f"{len(found) - good} are guesses to check by hand, "
          f"{len(runs) - len(found)} still unknown")
    print("genera:", ", ".join(f"{genus} {count}" for genus, count in
                               Counter(genus_and_group(hit["species"])[0]
                                       for hit in found.values()).most_common()))
    print("wrote", args.out)


if __name__ == "__main__":
    main()

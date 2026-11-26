#!/usr/bin/env python3
"""
lookup_reference_country.py

Read the country of origin for every GenBank reference in the tree, off the records.

Inputs:
    scratch_gt08.treefile        the tree, for the tip names
    contig_annotation.tsv        our own contigs, so they can be told from references

Outputs:
    reference_country.tsv        accession, country, host, organism
"""

import argparse
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

# Set the NCBI request settings, the tip pattern + the GenBank fields to read
# NCBI asks for no more than 3 requests a second without an API key
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
TOOL_NAME = "biof3001"

RETRIES = 4
TIMEOUT_SECONDS = 180
PAUSE_SECONDS = 0.4

# Tip names look like "(MH620810.1:0.012," so take the name before the colon
TIP_NAME = re.compile(r"[(,]([A-Za-z0-9_.|]+):")

COUNTRY_QUALIFIERS = ("country", "geo_loc_name")


def parse_arguments():
    """Read the command-line options."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", default="results/msa/tree_v2/scratch_gt08.treefile")
    parser.add_argument("--annotation",
                        default="results/msa/complete/tree_gt08/contig_annotation.tsv")
    parser.add_argument("-o", "--out", default="results/msa/tree_v2/reference_country.tsv")
    parser.add_argument("--email", default="")
    parser.add_argument("--batch", type=int, default=100)
    return parser.parse_args()


def reference_accessions(tree_path, annotation_path):
    """Return the tip names that are GenBank references, i.e. not one of ours."""
    ours = {line.split("\t")[0].replace("_R_", "") for line in open(annotation_path)}

    names = set(TIP_NAME.findall(open(tree_path).read()))
    return sorted({name.replace("_R_", "") for name in names
                   if name.replace("_R_", "") not in ours and re.match(r"[A-Z]", name)})


def efetch(accessions, email):
    """Fetch GenBank XML for one batch of accessions, retrying on a network error."""
    query = urllib.parse.urlencode({"db": "nuccore", "id": ",".join(accessions),
                                    "rettype": "gb", "retmode": "xml",
                                    "tool": TOOL_NAME, "email": email})

    for attempt in range(RETRIES):
        try:
            with urllib.request.urlopen(f"{EUTILS}?{query}", timeout=TIMEOUT_SECONDS) as response:
                return response.read()
        except Exception as error:
            print("  retry", attempt + 1, error, file=sys.stderr)
            time.sleep(3 * (attempt + 1))

    return b""


def read_records(xml):
    """Pull accession -> (country, host, organism) out of one batch of GenBank XML."""
    try:
        root = ET.fromstring(xml)
    except ET.ParseError as error:
        print("  parse error:", error, file=sys.stderr)
        return {}

    records = {}
    for sequence in root.iter("GBSeq"):
        accession = (sequence.findtext("GBSeq_accession-version")
                     or sequence.findtext("GBSeq_primary-accession"))
        organism = sequence.findtext("GBSeq_organism") or ""

        # Read country and host from the source qualifiers, not the top-level fields
        country = host = ""
        for qualifier in sequence.iter("GBQualifier"):
            name = qualifier.findtext("GBQualifier_name")
            value = qualifier.findtext("GBQualifier_value")
            if name in COUNTRY_QUALIFIERS and value and not country:
                country = value.split(":")[0].strip()
            if name == "host" and value and not host:
                host = value.strip()

        if accession:
            records[accession] = (country, host, organism)

    return records


def main():
    args = parse_arguments()

    # Step 1: list the reference tips in the tree, which is everything not our contig
    accessions = reference_accessions(args.tree, args.annotation)
    print(f"{len(accessions)} reference accessions in the tree")

    # Step 2: fetch the GenBank records in batches, pausing between requests
    found = {}
    for start in range(0, len(accessions), args.batch):
        batch = accessions[start:start + args.batch]
        found.update(read_records(efetch(batch, args.email)))
        time.sleep(PAUSE_SECONDS)
        print(f"  {min(start + args.batch, len(accessions))}/{len(accessions)}")

    # Step 3: write a row per accession, blank where NCBI had nothing, to match the tree
    with open(args.out, "w") as out:
        out.write("accession\tcountry\thost\torganism\n")
        for accession in accessions:
            country, host, organism = found.get(accession, ("", "", ""))
            out.write(f"{accession}\t{country}\t{host}\t{organism}\n")

    # Step 4: return summary in the terminal, how many had a country recorded
    named = sum(1 for a in accessions if found.get(a, ("",))[0])
    print(f"{named} of {len(accessions)} references have a country")
    print("wrote", args.out)


if __name__ == "__main__":
    main()

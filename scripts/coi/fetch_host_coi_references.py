#!/usr/bin/env python3
"""
fetch_host_coi_references.py

Fetch one COI sequence per mosquito species in the tree + a fly outgroup, for the host tree.

Inputs:
    results/msa/complete/tree_gt08/contig_annotation.tsv   host_species column

Outputs:
    results/host_tree/host_coi.fasta        one sequence per species
    results/host_tree/host_coi_sources.tsv  species, accession, title
"""

import csv
import os
import re
import sys
import time
import urllib.parse
import urllib.request

# Set file paths, NCBI eutils settings, the outgroup + the COI length window
# 500-2000 bp keeps full barcodes, leaves out whole mitochondrial genomes
ANNOTATION = "results/msa/complete/tree_gt08/contig_annotation.tsv"
OUTDIR = "results/host_tree"
OUT_FASTA = f"{OUTDIR}/host_coi.fasta"
OUT_LOG = f"{OUTDIR}/host_coi_sources.tsv"

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOL_NAME = "biof3001"

RETRIES = 4
TIMEOUT_SECONDS = 120
PAUSE_SECONDS = 0.4

OUTGROUP = "Drosophila melanogaster"

MIN_LENGTH = 500
MAX_LENGTH = 2000

FASTA_WIDTH = 60


def eutil(tool, email, **params):
    """Make one eutils call, retried a few times, and return the raw bytes."""
    params.update(tool=TOOL_NAME, email=email)
    url = f"{EUTILS}/{tool}.fcgi?" + urllib.parse.urlencode(params)

    for attempt in range(RETRIES):
        try:
            with urllib.request.urlopen(url, timeout=TIMEOUT_SECONDS) as response:
                return response.read()
        except Exception as error:
            print("  retry", attempt + 1, error, file=sys.stderr)
            time.sleep(3 * (attempt + 1))

    return b""


def host_species(path):
    """Return the binomial host species in first-seen order, + the outgroup."""
    species = []
    for row in csv.DictReader(open(path), delimiter="\t"):
        name = (row.get("host_species") or "").strip()
        # Keep binomials only, as a bare genus has no COI of its own
        if name and len(name.split()) >= 2 and name not in species:
            species.append(name)

    species.append(OUTGROUP)
    return species


def fetch_one_coi(species, email):
    """Return the longest COI record in the length window for one species, or None."""
    query = (f'"{species}"[Organism] AND (COI[Gene] OR COX1[Gene] OR '
             f'"cytochrome c oxidase subunit 1"[All Fields]) AND '
             f'{MIN_LENGTH}:{MAX_LENGTH}[SLEN]')

    xml = eutil("esearch", email, db="nuccore", term=query, retmax="20", sort="slen")
    ids = re.findall(r"<Id>(\d+)</Id>", xml.decode())
    if not ids:
        return None

    fasta = eutil("efetch", email, db="nuccore", id=ids[0],
                  rettype="fasta", retmode="text").decode()
    header, *body = fasta.strip().split("\n")
    return header.split()[0].lstrip(">"), "".join(body), header[1:120]


def main():
    email = sys.argv[1] if len(sys.argv) > 1 else ""
    os.makedirs(OUTDIR, exist_ok=True)

    # Step 1: list the species that need a reference
    species = host_species(ANNOTATION)
    print(len(species), "species:", ", ".join(species))

    # Step 2: fetch one record per species, in the same order
    sequences, log = [], []
    for name in species:
        record = fetch_one_coi(name, email)
        if record is None:
            print("  no COI found for", name)
            log.append((name, "", "not found"))
            continue

        accession, sequence, title = record
        sequences.append((name, accession, sequence))
        log.append((name, accession, title))
        print(f"  {name}: {accession} ({len(sequence)} bp)")
        time.sleep(PAUSE_SECONDS)

    # Step 3: write the FASTA named species|accession, so each tip carries both
    with open(OUT_FASTA, "w") as out:
        for name, accession, sequence in sequences:
            out.write(f">{name.replace(' ', '_')}|{accession}\n")
            for start in range(0, len(sequence), FASTA_WIDTH):
                out.write(sequence[start:start + FASTA_WIDTH] + "\n")

    # Step 4: write the sources table
    with open(OUT_LOG, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["species", "accession", "title"])
        writer.writerows(log)

    # Step 5: return summary in the terminal + the next commands
    print(len(sequences), "sequences written to", OUT_FASTA)
    print("next:")
    print(f"  mafft --auto {OUT_FASTA} > {OUTDIR}/host_coi_msa.fasta")
    print(f"  iqtree2 -s {OUTDIR}/host_coi_msa.fasta -m MFP -B 1000 -T 4 \\")
    print(f"          --prefix {OUTDIR}/host_coi")


if __name__ == "__main__":
    main()

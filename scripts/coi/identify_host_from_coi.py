#!/usr/bin/env python3
"""
identify_host_from_coi.py

Identify each run's mosquito species by BLASTing its assembly against the Culicidae COI database.

Inputs:
    --assemblies                   directory of per-run assemblies
    reference/coi/culicidae_coi    the BLAST database

Outputs:
    results/msa/tree_v2/host_from_coi.tsv   run, species, pident, bitscore, evidence
"""

import argparse
import csv
import glob
import os
import re
import subprocess
import sys

# Set the blastn call, its cut-offs + the patterns for reading names out of text
OUTFMT = "6 qseqid sseqid stitle pident length bitscore"
PIDENT_FIELD, BITSCORE_FIELD, TITLE_FIELD = 3, 5, 2

MAX_TARGET_SEQS = 5
EVALUE = "1e-20"
BLAST_TIMEOUT_SECONDS = 1800

BINOMIAL = re.compile(r"\b([A-Z][a-z]+)\s+([a-z]{3,})")
RUN_ACCESSION = re.compile(r"([SED]RR\d+)")


def parse_arguments():
    """Read the command-line options."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--assemblies", required=True,
                        help="directory of per-run assemblies")
    parser.add_argument("--pattern", default="*/contigs.fasta",
                        help="glob inside that directory, the run accession comes from the path")
    parser.add_argument("--db", default="reference/coi/culicidae_coi",
                        help="makeblastdb prefix built from the Culicidae COI sequences")
    parser.add_argument("--runs", help="file with the run accessions to do, one per line")
    parser.add_argument("-o", "--out", default="results/msa/tree_v2/host_from_coi.tsv")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--min-bits", type=float, default=200.0)
    parser.add_argument("--min-pident", type=float, default=95.0)
    return parser.parse_args()


def run_accession(path):
    """Return the run this assembly belongs to, from the path."""
    match = RUN_ACCESSION.search(path)
    return match.group(1) if match else os.path.basename(os.path.dirname(path))


def find_assemblies(args):
    """Return the (run, path) pairs to search, narrowed to --runs if given."""
    paths = sorted(glob.glob(os.path.join(args.assemblies, args.pattern)))
    if not paths:
        sys.exit(f"nothing matched {os.path.join(args.assemblies, args.pattern)}")

    jobs = [(run_accession(path), path) for path in paths]

    if args.runs:
        wanted = {line.strip() for line in open(args.runs) if line.strip()}
        jobs = [job for job in jobs if job[0] in wanted]

    return jobs


def best_coi_hit(path, args):
    """Return the highest-scoring COI hit for one assembly, or None if nothing passes."""
    command = ["blastn", "-query", path, "-db", args.db, "-outfmt", OUTFMT,
               "-max_target_seqs", str(MAX_TARGET_SEQS), "-evalue", EVALUE,
               "-num_threads", str(args.threads)]
    try:
        output = subprocess.run(command, capture_output=True, text=True,
                                timeout=BLAST_TIMEOUT_SECONDS).stdout
    except Exception as error:
        print(f"  blastn failed on {path}, {error}", file=sys.stderr)
        return None

    best = None
    for line in output.splitlines():
        fields = line.split("\t")
        if len(fields) < 6:
            continue

        pident, bits = float(fields[PIDENT_FIELD]), float(fields[BITSCORE_FIELD])
        if bits < args.min_bits or pident < args.min_pident:
            continue
        if best is None or bits > best[0]:
            best = (bits, pident, fields[TITLE_FIELD])

    return best


def species_from_title(title):
    """Return the binomial in a subject title, or the title itself if there isn't one."""
    match = BINOMIAL.search(title)
    return f"{match.group(1)} {match.group(2)}" if match else title[:60]


def main():
    args = parse_arguments()

    # Step 1: find the assemblies to search, one per run
    jobs = find_assemblies(args)
    print(f"{len(jobs)} assemblies to search")

    # Step 2: blastn each one against the COI database + keep its best hit
    rows = []
    for i, (run, path) in enumerate(jobs, 1):
        best = best_coi_hit(path, args)

        # Keep a row for runs with no hit, so the table matches the run list
        if best is None:
            rows.append((run, "", "", "", "no COI hit above the thresholds"))
            print(f"  {i}/{len(jobs)} {run}: nothing")
            continue

        bits, pident, title = best
        species = species_from_title(title)
        rows.append((run, species, f"{pident:.1f}", f"{bits:.0f}", title[:120]))
        print(f"  {i}/{len(jobs)} {run}: {species} ({pident:.1f}% id)")

    # Step 3: write one row per run, including runs with no hit
    with open(args.out, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["run", "species", "pident", "bitscore", "evidence"])
        writer.writerows(rows)

    # Step 4: return summary in the terminal
    named = sum(1 for row in rows if row[1])
    print(f"{named} of {len(rows)} runs got a species from their own reads")
    print("wrote", args.out)


if __name__ == "__main__":
    main()

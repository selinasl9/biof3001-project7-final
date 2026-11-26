#!/usr/bin/env python3
"""
merge_host_lookup.py

Put the hosts found by lookup_host_species.py back into the contig annotation table.

Only the "good" rows are merged by default, as the rest are guesses to check by hand.

Inputs:
    --annotation  contig_annotation.tsv from build_contig_annotation.py
    --lookup      host_lookup.tsv from lookup_host_species.py

Outputs:
    the annotation table with host_species, host_genus, host_group + province filled in
    contig_annotation.tsv.bak   the original, when rewriting in place
"""

import argparse
import csv
import shutil
from collections import Counter

# Set the columns the lookup fills in
FILLED = ["host_species", "host_genus", "host_group", "province"]


def read_lookup(path, confidences):
    """Return {run: lookup row} for the rows that name a host at an accepted confidence."""
    found = {}
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("host_species") and row.get("confidence", "") in confidences:
                found[row["run"]] = row
    return found


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotation", required=True)
    parser.add_argument("--lookup", required=True)
    parser.add_argument("-o", "--out", default="",
                        help="default is to rewrite the annotation table, keeping a .bak")
    parser.add_argument("--overwrite", action="store_true",
                        help="replace hosts already in the table, not just the unknown ones")
    parser.add_argument("--include-unchecked", action="store_true",
                        help="merge the \"check by hand\" rows as well")
    args = parser.parse_args()
    out = args.out or args.annotation

    # Step 1: read the lookup rows we are willing to merge
    confidences = {"good", "check by hand"} if args.include_unchecked else {"good"}
    found = read_lookup(args.lookup, confidences)

    # Step 2: read the annotation table, keeping its column order + adding missing columns
    with open(args.annotation) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    columns = list(rows[0].keys()) + [c for c in FILLED if c not in rows[0]]

    # Step 3: fill the unknown hosts, leaving the ones already named alone
    changed = 0
    for row in rows:
        hit = found.get(row.get("run") or row["contig"].split("_")[0])
        if not hit:
            continue
        if (row.get("host_group") or "unknown") != "unknown" and not args.overwrite:
            continue

        row["host_species"] = hit["host_species"]
        row["host_genus"] = hit["host_genus"]
        row["host_group"] = hit["host_group"]
        # Fill a blank province only, as the lookup never argues with the metadata
        if hit.get("province") and not (row.get("province") or "").strip():
            row["province"] = hit["province"]
        changed += 1

    # Step 4: back the old table up before overwriting it in place
    if out == args.annotation:
        shutil.copy(args.annotation, args.annotation + ".bak")
        print("kept the old table as", args.annotation + ".bak")

    # Step 5: write the table back out
    with open(out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t",
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})

    # Step 6: return summary in the terminal, with the host_group counts now
    groups = Counter(row.get("host_group") or "unknown" for row in rows)
    print(f"{changed} contigs given a host from the lookup")
    print("host_group now:", ", ".join(f"{k} {v}" for k, v in sorted(groups.items())))
    print("wrote", out)


if __name__ == "__main__":
    main()

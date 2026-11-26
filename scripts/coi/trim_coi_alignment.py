#!/usr/bin/env python3
"""
trim_coi_alignment.py

Trim the host COI alignment to the region every sequence covers.

Inputs:
    results/host_tree/host_coi_msa.fasta      the aligned COI sequences

Outputs:
    results/host_tree/host_coi_barcode.fasta  the alignment, trimmed
"""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))

from common.fasta_io import read_fasta_dict, write_fasta

# Set file paths, the column occupancy cutoff, the gap characters + the usable width
IN_MSA = "results/host_tree/host_coi_msa.fasta"
OUT_MSA = "results/host_tree/host_coi_barcode.fasta"

MIN_OCCUPANCY = 1.0

GAP_CHARACTERS = "-nN?"

MIN_USEFUL_COLUMNS = 300


def main():
    in_msa = sys.argv[1] if len(sys.argv) > 1 else IN_MSA
    out_msa = sys.argv[2] if len(sys.argv) > 2 else OUT_MSA

    # Step 1: read the alignment
    sequences = read_fasta_dict(in_msa)
    names = list(sequences)
    n_columns = len(sequences[names[0]])
    print(len(names), "sequences,", n_columns, "columns")

    # Step 2: keep a column only if enough sequences have a real base in it
    needed = MIN_OCCUPANCY * len(names)
    keep = [i for i in range(n_columns)
            if sum(sequences[name][i] not in GAP_CHARACTERS for name in names) >= needed]
    print(len(keep), "columns kept (every sequence has a base)")

    # Step 3: warn if too few columns are left for a trustworthy tree
    if len(keep) < MIN_USEFUL_COLUMNS:
        print("WARNING: that is a very short region; the tree will be weakly supported")

    # Step 4: write the same columns out for every sequence
    write_fasta(out_msa, ((name, "".join(sequences[name][i] for i in keep))
                          for name in names))
    print("wrote", out_msa)

    # Step 5: return the next command in the terminal
    print("next:")
    print(f"  iqtree2 -s {out_msa} -m MFP -B 1000 -T 4 "
          f"--prefix results/host_tree/host_coi_barcode")


if __name__ == "__main__":
    main()

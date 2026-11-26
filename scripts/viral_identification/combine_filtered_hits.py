#!/usr/bin/env python3
"""
combine_filtered_hits.py

Apply the three hit thresholds to the PRJNA1003643 best hits, then combine each
filtered table with Cammy's table for the other BioProjects.

Inputs:
    CN_PRJNA1003643_all_best_hits_with_metadata.tsv  China best hits
    {mode}_filtered_hits.tsv                         Cammy's half, one per mode

Outputs:
    CN_PRJNA1003643_{mode}_filtered_hits.tsv         China half, one per mode
    combined_{mode}_filtered_hits.tsv                both halves, one per mode
"""

from pathlib import Path

import pandas as pd

CN_INPUT = Path("results/viral_identification/PRJNA1003643/summary/CN_PRJNA1003643_all_best_hits_with_metadata.tsv")
CN_OUTDIR = Path("results/viral_identification/PRJNA1003643/filtered")
CN_PREFIX = "CN_PRJNA1003643_"

CAMMY_DIR = Path("results/viral_identification/filtered")
COMBINED_PREFIX = "combined_"

# Set the per-mode thresholds, the same values as Cammy's filter_hits_by_threshold.py
# evalue has no identity or coverage cutoff, it is E-value < 1e-5 on its own
THRESHOLDS = {"strict": (70, 80), "loose": (30, 50), "evalue": None}
MAX_EVALUE = 1e-5

NUMERIC_COLUMNS = ["pident", "qstart", "qend", "sstart", "send", "qlen", "slen", "evalue"]


def load_china_best_hits():
    """Read the China best hits + add the recalculated coverage columns."""
    df = pd.read_csv(CN_INPUT, sep="\t", dtype={"staxids": str})
    print("China best-hit rows:", len(df))

    for column in NUMERIC_COLUMNS:
        df[column] = pd.to_numeric(df[column], errors="coerce")

    # Take coverage from the alignment spans, so both halves are filtered the same way
    df["query_aligned_span"] = (df["qend"] - df["qstart"]).abs() + 1
    df["subject_aligned_span"] = (df["send"] - df["sstart"]).abs() + 1
    df["query_coverage_calc"] = df["query_aligned_span"] / df["qlen"] * 100
    df["subject_coverage_calc"] = df["subject_aligned_span"] / df["slen"] * 100

    return df


def apply_threshold(df, mode):
    """Return the rows of df that pass one of the three thresholds."""
    cutoffs = THRESHOLDS[mode]
    if cutoffs is None:
        return df[df["evalue"] < MAX_EVALUE]

    min_pident, min_coverage = cutoffs
    return df[(df["pident"] >= min_pident)
              & ((df["query_coverage_calc"] >= min_coverage)
                 | (df["subject_coverage_calc"] >= min_coverage))]


def combine_with_cammy(cn, mode):
    """Concatenate the China table with Cammy's table for the same mode."""
    cammy = pd.read_csv(CAMMY_DIR / f"{mode}_filtered_hits.tsv", sep="\t", dtype={"staxids": str})

    # Drop China rows from her table, in case the China runs were ever added there
    cammy = cammy[cammy["bioproject"] != "PRJNA1003643"]

    combined = pd.concat([cammy, cn], ignore_index=True)

    # Keep her column order, with the China-only location columns after country
    extra = [c for c in cn.columns if c not in cammy.columns]
    cols = list(cammy.columns)
    cut = cols.index("country") + 1

    return cammy, combined[cols[:cut] + extra + cols[cut:]]


def main():
    CN_OUTDIR.mkdir(parents=True, exist_ok=True)

    # Step 1: read the China best hits + recalculate coverage
    df = load_china_best_hits()

    written = []

    for mode in THRESHOLDS:
        # Step 2: filter the China half on this mode's thresholds
        cn = apply_threshold(df, mode).copy()
        cn_path = CN_OUTDIR / f"{CN_PREFIX}{mode}_filtered_hits.tsv"
        cn.to_csv(cn_path, sep="\t", index=False)

        # Step 3: add Cammy's half + keep her column order
        cammy, combined = combine_with_cammy(cn, mode)
        combined_path = CAMMY_DIR / f"{COMBINED_PREFIX}{mode}_filtered_hits.tsv"
        combined.to_csv(combined_path, sep="\t", index=False)

        written += [cn_path, combined_path]

        print(f"\n{mode.upper()}")
        print(f"  China:    {len(cn)} hits, {cn['sample'].nunique()} samples, "
              f"{cn['sscinames'].nunique()} viral names")
        print(f"  Cammy:    {len(cammy)} hits, {cammy['sample'].nunique()} samples")
        print(f"  Combined: {len(combined)} hits, {combined['sample'].nunique()} samples, "
              f"{combined['sscinames'].nunique()} viral names")

    print("\nCreated:")
    for path in written:
        print(f"  {path}")


if __name__ == "__main__":
    main()

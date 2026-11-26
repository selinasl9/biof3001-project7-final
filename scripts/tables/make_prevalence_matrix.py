#!/usr/bin/env python3
"""
make_prevalence_matrix.py

Build a species-level prevalence matrix from the host-cleaned hits, one per threshold.

    rows    = datasets (BioProjects)
    columns = virus species
    cells   = percentage of that dataset's runs in which the species was detected

Inputs:
    results/viral_identification/host_filter/{mode}_host_cleaned_hits.tsv
    results/viral_identification/PRJNA1003643/host_filter/CN_PRJNA1003643_{mode}_host_cleaned_hits.tsv
    results/summary/sample_metadata.tsv
    results/viral_identification/PRJNA1003643/summary/CN_PRJNA1003643_sample_metadata.tsv

Outputs:
    results/prevalence/{mode}_prevalence_matrix_species.tsv
    results/prevalence/prevalence_matrix_dimensions.tsv
"""

import argparse
from pathlib import Path

import pandas as pd

# Set the thresholds, the input templates + the output folder
MODES = ["strict", "loose", "evalue"]

CAMMY_HITS = "results/viral_identification/host_filter/{mode}_host_cleaned_hits.tsv"
CHINA_HITS = ("results/viral_identification/PRJNA1003643/host_filter/"
              "CN_PRJNA1003643_{mode}_host_cleaned_hits.tsv")
META = ["results/summary/sample_metadata.tsv",
        "results/viral_identification/PRJNA1003643/summary/CN_PRJNA1003643_sample_metadata.tsv"]

OUTDIR = Path("results/prevalence")

UNCLASSIFIED = "unclassified virus"


def load_runs():
    """Return every run screened, with the BioProject it belongs to."""
    frames = [pd.read_csv(path, sep="\t", low_memory=False)[["sample", "bioproject"]]
              for path in META]
    return pd.concat(frames, ignore_index=True).drop_duplicates("sample")


def load_hits(mode):
    """Read the host-cleaned hits for one threshold, both halves of the dataset."""
    keep = ["sample", "species", "sscinames"]

    frames = []
    for template in (CAMMY_HITS, CHINA_HITS):
        path = Path(template.format(mode=mode))
        if not path.exists():
            raise SystemExit(f"Missing input file: {path}")
        frames.append(pd.read_csv(path, sep="\t", low_memory=False)[keep])

    hits = pd.concat(frames, ignore_index=True)

    # Take the species name, then the subject name, then a placeholder
    # Blanks have to become NA first or fillna skips them
    hits["species"] = (hits["species"].replace("", pd.NA)
                       .fillna(hits["sscinames"].replace("", pd.NA))
                       .fillna(UNCLASSIFIED))

    return hits[["sample", "species"]]


def build_matrix(hits, runs, as_counts=False):
    """Build a dataset x virus species matrix of prevalence percentages."""
    # Step 1: count the runs screened per dataset, the denominator for each row
    screened = runs.groupby("bioproject")["sample"].nunique()

    # Step 2: count the runs positive per dataset and species, the numerators
    detections = (hits.merge(runs, on="sample", how="inner")
                  .drop_duplicates(["bioproject", "species", "sample"]))
    positives = (detections.groupby(["bioproject", "species"])["sample"]
                 .nunique().unstack(fill_value=0))

    # Step 3: give every screened dataset a row, even with no detections
    positives = positives.reindex(index=screened.index, fill_value=0)

    matrix = positives if as_counts else (positives.div(screened, axis=0) * 100).round(2)

    # Step 4: order the species, most prevalent first over all datasets
    order = positives.sum().sort_values(ascending=False).index
    matrix = matrix[order]

    matrix.insert(0, "runs_screened", screened)
    matrix.index.name = "dataset"
    return matrix


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--counts", action="store_true",
                        help="write positive run counts instead of percentages")
    args = parser.parse_args()

    OUTDIR.mkdir(parents=True, exist_ok=True)

    runs = load_runs()
    print(f"{len(runs)} runs screened across {runs['bioproject'].nunique()} datasets\n")

    dimensions = []

    for mode in MODES:
        hits = load_hits(mode)

        # Drop a hit whose run has no metadata row, as it has no denominator
        unknown = set(hits["sample"]) - set(runs["sample"])
        if unknown:
            print(f"{mode}: {len(unknown)} runs have hits but no metadata row, dropped")

        matrix = build_matrix(hits, runs, as_counts=args.counts)

        path = OUTDIR / f"{mode}_prevalence_matrix_species.tsv"
        matrix.to_csv(path, sep="\t")

        n_datasets, n_species = matrix.shape[0], matrix.shape[1] - 1
        print(f"{mode:7s} {n_datasets} datasets x {n_species} virus species  ->  {path}")

        dimensions.append({
            "mode": mode,
            "datasets": n_datasets,
            "virus_species": n_species,
            "runs_screened": int(matrix["runs_screened"].sum()),
        })

    pd.DataFrame(dimensions).to_csv(
        OUTDIR / "prevalence_matrix_dimensions.tsv", sep="\t", index=False)
    print(f"\nWrote {OUTDIR}/prevalence_matrix_dimensions.tsv")


if __name__ == "__main__":
    main()

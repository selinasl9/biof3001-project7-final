#!/usr/bin/env python3
"""BIOF3001 Project 7 — China PRJNA1003643 host-match filtering.

Objective
---------
Remove candidate viral contigs with strong Culicidae nucleotide matches
(>=80% identity AND >=80% query coverage) from each of three China-specific
DIAMOND candidate-hit tables. Preserve the remaining hit annotations.

Inputs (relative to repository root)
------------------------------------
results/viral_identification/PRJNA1003643/filtered/
    CN_PRJNA1003643_{strict,loose,evalue}_filtered_hits.tsv
results/viral_identification/PRJNA1003643/host_filter/blastn/
    CN_PRJNA1003643_union_vs_culicidae.tsv

Outputs in results/viral_identification/PRJNA1003643/host_filter/
----------------------------------------------------------------
    CN_PRJNA1003643_strong_host_alignments.tsv
    CN_PRJNA1003643_{strict,loose,evalue}_removed_mosquito_contigs.tsv
    CN_PRJNA1003643_{strict,loose,evalue}_host_cleaned_hits.tsv

Usage
-----
python scripts/viral_identification/china_supplement/clean_CN_PRJNA1003643_host_hits.py

Dependencies: Python 3, pandas.


Identification convention:
The extraction step uses FASTA identifiers sample|qseqid. Using this
composite key avoids incorrectly removing same-named contigs from different
sequencing runs. Similarity to Culicidae is not definitive proof of origin.
"""

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "results/viral_identification/PRJNA1003643"
FILTER = BASE / "filtered"
HOST = BASE / "host_filter"
BLAST = HOST / "blastn/CN_PRJNA1003643_union_vs_culicidae.tsv"
PREFIX = "CN_PRJNA1003643"
MODES = ("strict", "loose", "evalue")
BLAST_COLS = (
    "qseqid", "sseqid", "pident", "length", "qlen", "slen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
)


def main() -> int:
    # Validate all inputs before writing any outputs.
    if not BLAST.is_file():
        raise FileNotFoundError(f"Host BLASTN table not found: {BLAST}")
    tables = {mode: FILTER / f"{PREFIX}_{mode}_filtered_hits.tsv" for mode in MODES}
    for mode, file in tables.items():
        if not file.is_file():
            raise FileNotFoundError(f"{mode} candidate table missing: {file}")

    # An empty file is a valid BLASTN result: no host matches found.
    if BLAST.stat().st_size == 0:
        blast = pd.DataFrame(columns=BLAST_COLS)
    else:
        blast = pd.read_csv(BLAST, sep="\t", header=None, names=list(BLAST_COLS))

    for col in ("pident", "qstart", "qend", "qlen"):
        blast[col] = pd.to_numeric(blast[col], errors="coerce")
    if not blast.empty and (blast["qlen"] <= 0).any():
        raise ValueError("Host BLASTN data contain nonpositive query lengths")

    blast["query_coverage"] = (
        (blast["qend"] - blast["qstart"]).abs() + 1
    ) / blast["qlen"] * 100

    strong = blast.loc[
        (blast["pident"] >= 80) & (blast["query_coverage"] >= 80)
    ].copy()
    host_ids = set(strong["qseqid"].astype(str))

    HOST.mkdir(parents=True, exist_ok=True)
    strong.to_csv(HOST / f"{PREFIX}_strong_host_alignments.tsv", sep="\t", index=False)
    print(f"BLASTN alignments: {len(blast):,}; strong alignments: {len(strong):,}")
    print(f"Distinct sample-specific host-like contigs: {len(host_ids):,}")

    for mode, infile in tables.items():
        df = pd.read_csv(infile, sep="\t")
        if not {"sample", "qseqid"}.issubset(df.columns):
            raise ValueError(f"Missing sample/qseqid column in {infile}")
        # Must match the exact FASTA header created by extraction.
        ids = df["sample"].astype(str) + "|" + df["qseqid"].astype(str)
        flagged = ids.isin(host_ids)
        removed = df.loc[flagged].copy()
        cleaned = df.loc[~flagged].copy()

        removed.to_csv(HOST / f"{PREFIX}_{mode}_removed_mosquito_contigs.tsv", sep="\t", index=False)
        cleaned.to_csv(HOST / f"{PREFIX}_{mode}_host_cleaned_hits.tsv", sep="\t", index=False)
        print(
            f"{mode}: before={len(df):,}, removed={len(removed):,}, "
            f"after={len(cleaned):,}, runs remaining={cleaned['sample'].nunique():,}"
        )
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

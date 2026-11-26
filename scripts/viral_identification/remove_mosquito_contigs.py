#!/usr/bin/env python3
"""BIOF3001 Project 7 — Mosquito RNA Virome.

Script: remove_mosquito_contigs.py

Objective
---------
Remove candidate viral-hit rows for contigs with strong alignments to a
Culicidae nucleotide reference database, based on host BLASTN output.

Strong host-match rule (unchanged from original):
    nucleotide identity >= 80% AND query alignment coverage >= 80%.

Inputs
------
results/viral_identification/filtered/{mode}_filtered_hits.tsv
results/viral_identification/host_filter/blastn/{mode}_vs_culicidae.tsv

Outputs (under results/viral_identification/host_filter/)
--------------------------------------------------------
{mode}_mosquito_blast_hits.tsv
{mode}_removed_mosquito_contigs.tsv
{mode}_host_cleaned_hits.tsv

Modes: strict, loose, evalue.
Dependencies: Python 3 and pandas.
Run: python scripts/viral_identification/remove_mosquito_contigs.py

Interpretation: A strong Culicidae reference match is a screening criterion,
not definitive biological proof that a contig is host-derived or an EVE.
"""

from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
FILTERED_DIR = ROOT / "results/viral_identification/filtered"
OUTDIR = ROOT / "results/viral_identification/host_filter"
BLAST_DIR = OUTDIR / "blastn"
MODES = ("strict", "loose", "evalue")
BLAST_COLS = [
    "qseqid", "sseqid", "pident", "length", "qlen", "slen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
]


def main():
    # Require all BLAST outputs. An empty file is acceptable (successful zero hits),
    # but a missing file must not be interpreted as no host matches.
    for mode in MODES:
        for file in (
            FILTERED_DIR / f"{mode}_filtered_hits.tsv",
            BLAST_DIR / f"{mode}_vs_culicidae.tsv",
        ):
            if not file.is_file():
                raise FileNotFoundError(f"Required input missing: {file}")

    OUTDIR.mkdir(parents=True, exist_ok=True)
    for mode in MODES:
        blast_file = BLAST_DIR / f"{mode}_vs_culicidae.tsv"
        filtered_file = FILTERED_DIR / f"{mode}_filtered_hits.tsv"
        if blast_file.stat().st_size == 0:
            blast = pd.DataFrame(columns=BLAST_COLS)
        else:
            blast = pd.read_csv(blast_file, sep="\t", names=BLAST_COLS)

        for col in ("pident", "qstart", "qend", "qlen"):
            blast[col] = pd.to_numeric(blast[col], errors="coerce")
        blast["query_aligned_span"] = (blast["qend"] - blast["qstart"]).abs() + 1
        blast["query_coverage_calc"] = (
            blast["query_aligned_span"] / blast["qlen"].where(blast["qlen"] > 0) * 100
        )
        strong = blast[
            (blast["pident"] >= 80) & (blast["query_coverage_calc"] >= 80)
        ].copy()
        host_ids = set(strong["qseqid"].dropna().astype(str))

        original = pd.read_csv(filtered_file, sep="\t", dtype={"sample": str, "qseqid": str})
        if not {"sample", "qseqid"}.issubset(original.columns):
            raise ValueError(f"Missing sample/qseqid columns in {filtered_file}")

        # Must match the header format created by extract_filtered_contigs.py.
        identifiers = original["sample"].astype(str) + "|" + original["qseqid"].astype(str)
        is_host = identifiers.isin(host_ids)
        cleaned = original.loc[~is_host].copy()
        removed = original.loc[is_host].copy()

        strong.to_csv(OUTDIR / f"{mode}_mosquito_blast_hits.tsv", sep="\t", index=False)
        removed.to_csv(OUTDIR / f"{mode}_removed_mosquito_contigs.tsv", sep="\t", index=False)
        cleaned.to_csv(OUTDIR / f"{mode}_host_cleaned_hits.tsv", sep="\t", index=False)

        print(
            f"{mode.upper()}: original hits={len(original)}, "
            f"removed contigs={identifiers[is_host].nunique()}, "
            f"remaining hits={len(cleaned)}, remaining runs={cleaned['sample'].nunique()}"
        )
    print("Host reference-based filtering complete.")


if __name__ == "__main__":
    main()

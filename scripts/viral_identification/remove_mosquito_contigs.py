
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

# STEP 1: Import required libraries
# pathlib handles file paths; pandas processes tables.
from pathlib import Path
import pandas as pd


# STEP 2: Define input/output directories and filtering modes

# Find project root (two directories above this script).
ROOT = Path(__file__).resolve().parents[2]

# Directory containing DIAMOND-filtered viral hits.
FILTERED_DIR = ROOT / "results/viral_identification/filtered"

# Directory for host-filtering results.
OUTDIR = ROOT / "results/viral_identification/host_filter"

# Directory containing Culicidae BLASTN results.
BLAST_DIR = OUTDIR / "blastn"

# Process three filtering branches independently.
MODES = ("strict", "loose", "evalue")

# Define the 12 columns in the BLASTN output.
BLAST_COLS = [
    "qseqid", "sseqid", "pident", "length", "qlen", "slen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
]


def main():

    # STEP 3: Check all required input files
    # Missing files cause an error; empty BLASTN files are allowed.
    for mode in MODES:
        for file in (
            FILTERED_DIR / f"{mode}_filtered_hits.tsv",
            BLAST_DIR / f"{mode}_vs_culicidae.tsv",
        ):
            if not file.is_file():
                raise FileNotFoundError(f"Required input missing: {file}")


    # STEP 4: Create output directory and process each mode
    OUTDIR.mkdir(parents=True, exist_ok=True)

    for mode in MODES:

        # Define BLASTN and filtered-hit input files.
        blast_file = BLAST_DIR / f"{mode}_vs_culicidae.tsv"
        filtered_file = FILTERED_DIR / f"{mode}_filtered_hits.tsv"


        # STEP 5: Read BLASTN alignment results
        # Create an empty table if BLASTN reported zero hits.
        if blast_file.stat().st_size == 0:
            blast = pd.DataFrame(columns=BLAST_COLS)
        else:
            blast = pd.read_csv(blast_file, sep="\t", names=BLAST_COLS)


        # STEP 6: Calculate query alignment coverage

        # Convert identity and alignment positions to numbers.
        for col in ("pident", "qstart", "qend", "qlen"):
            blast[col] = pd.to_numeric(blast[col], errors="coerce")

        # Calculate the nucleotide span covered by each alignment.
        blast["query_aligned_span"] = (blast["qend"] - blast["qstart"]).abs() + 1

        # Calculate percentage of the full query covered by alignment.
        # Query coverage = aligned span / query length × 100.
        blast["query_coverage_calc"] = (
            blast["query_aligned_span"] / blast["qlen"].where(blast["qlen"] > 0) * 100
        )


        # STEP 7: Identify strong mosquito-reference matches

        # Keep alignments with >=80% identity AND >=80% query coverage.
        strong = blast[
            (blast["pident"] >= 80) & (blast["query_coverage_calc"] >= 80)
        ].copy()

        # Collect unique query contig IDs with strong host matches.
        host_ids = set(strong["qseqid"].dropna().astype(str))


        # STEP 8: Read original candidate viral hits

        # Load filtered hits, preserving sample and contig IDs as strings.
        original = pd.read_csv(filtered_file, sep="\t", dtype={"sample": str, "qseqid": str})

        # Check that the required identifier columns exist.
        if not {"sample", "qseqid"}.issubset(original.columns):
            raise ValueError(f"Missing sample/qseqid columns in {filtered_file}")


        # STEP 9: Remove contigs with strong host matches

        # Combine sample ID and contig ID to match BLASTN FASTA headers.
        # Example: SRR123456|contig_1
        identifiers = original["sample"].astype(str) + "|" + original["qseqid"].astype(str)

        # Identify candidate-hit rows belonging to host-matched contigs.
        is_host = identifiers.isin(host_ids)

        # Keep rows without strong mosquito-reference matches.
        cleaned = original.loc[~is_host].copy()

        # Collect rows removed by host filtering.
        removed = original.loc[is_host].copy()


        # STEP 10: Save filtering results

        # Save strong mosquito-reference BLASTN alignments.
        strong.to_csv(OUTDIR / f"{mode}_mosquito_blast_hits.tsv", sep="\t", index=False)

        # Save candidate-hit rows removed by host filtering.
        removed.to_csv(OUTDIR / f"{mode}_removed_mosquito_contigs.tsv", sep="\t", index=False)

        # Save remaining candidate viral hits for downstream analysis.
        cleaned.to_csv(OUTDIR / f"{mode}_host_cleaned_hits.tsv", sep="\t", index=False)


        # STEP 11: Print summary for each filtering mode

        # Report original hits, removed contigs, remaining hits and runs.
        print(
            f"{mode.upper()}: original hits={len(original)}, "
            f"removed contigs={identifiers[is_host].nunique()}, "
            f"remaining hits={len(cleaned)}, remaining runs={cleaned['sample'].nunique()}"
        )

    # Confirm completion of all three filtering branches.
    print("Host reference-based filtering complete.")


# STEP 12: Run the main function
# Execute main() only when this script is run directly.
if __name__ == "__main__":
    main()

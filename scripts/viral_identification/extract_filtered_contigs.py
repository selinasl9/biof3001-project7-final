#!/usr/bin/env python3
"""BIOF3001 Project 7 — Mosquito RNA Virome.

Script: extract_filtered_contigs.py
Responsible contributor: [Confirm name before submission]

Objective
---------
Extract the original nucleotide contigs for candidate DIAMOND hits under each
filtering condition, preparing FASTA inputs for Culicidae host BLASTN screening.
Each FASTA header is rewritten as >sample|qseqid so identical contig IDs from
different sequencing runs cannot collide.

Inputs
------
results/viral_identification/filtered/{strict,loose,evalue}_filtered_hits.tsv
data/raw/logan_contigs/<sample>.contigs.fa

Outputs
-------
results/viral_identification/host_filter/{mode}_candidate_contigs.fa
results/viral_identification/host_filter/{mode}_missing_contigs.tsv
    (written only when extraction is incomplete; script then exits nonzero)

Dependencies: Python 3 and pandas.
Run: python scripts/viral_identification/extract_filtered_contigs.py

Note: This script does not identify viruses or run sequence similarity searches;
it extracts previously annotated contigs for a separate host-screening step.
"""

from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
FILTERED_DIR = ROOT / "results/viral_identification/filtered"
FASTA_DIR = ROOT / "data/raw/logan_contigs"
OUTDIR = ROOT / "results/viral_identification/host_filter"
MODES = ("strict", "loose", "evalue")


def write_selected_fasta(source: Path, sample: str, wanted: set[str], output) -> set[str]:
    """Stream a FASTA file once; write selected records and return found IDs."""
    found = set()
    write_record = False
    current_id = None
    with source.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                current_id = line[1:].strip().split()[0]
                # Prevent duplicate sequence records from being written twice.
                write_record = current_id in wanted and current_id not in found
                if write_record:
                    output.write(f">{sample}|{current_id}\n")
                    found.add(current_id)
            elif write_record:
                output.write(line)
    return found


def main():
    # Check all required tables before generating any outputs.
    tables = {mode: FILTERED_DIR / f"{mode}_filtered_hits.tsv" for mode in MODES}
    for mode, table in tables.items():
        if not table.is_file():
            raise FileNotFoundError(f"Missing {mode} table: {table}")

    OUTDIR.mkdir(parents=True, exist_ok=True)
    errors = []
    for mode, table in tables.items():
        df = pd.read_csv(table, sep="\t", dtype={"sample": str, "qseqid": str})
        if not {"sample", "qseqid"}.issubset(df.columns):
            raise ValueError(f"Missing 'sample' or 'qseqid' column in {table}")
        candidates = df[["sample", "qseqid"]].dropna().drop_duplicates()
        final = OUTDIR / f"{mode}_candidate_contigs.fa"
        temp = OUTDIR / f"{mode}_candidate_contigs.fa.tmp"
        missing_file = OUTDIR / f"{mode}_missing_contigs.tsv"
        missing = []
        found_count = 0
        try:
            with temp.open("w", encoding="utf-8") as out:
                for sample, group in candidates.groupby("sample", sort=True):
                    wanted = set(group["qseqid"])
                    fasta = FASTA_DIR / f"{sample}.contigs.fa"
                    if not fasta.is_file():
                        print(f"WARNING: Missing contig FASTA: {fasta}")
                        missing.extend((sample, q) for q in sorted(wanted))
                        continue
                    found = write_selected_fasta(fasta, sample, wanted, out)
                    found_count += len(found)
                    missing.extend((sample, q) for q in sorted(wanted - found))
            if missing:
                # Never publish an incomplete candidate FASTA as a finished output.
                pd.DataFrame(missing, columns=["sample", "qseqid"]).to_csv(
                    missing_file, sep="\t", index=False
                )
                temp.unlink(missing_ok=True)
                final.unlink(missing_ok=True)
                errors.append(f"{mode}: {len(missing)} missing contigs; see {missing_file}")
            else:
                temp.replace(final)
                missing_file.unlink(missing_ok=True)  # remove stale reports
            print(f"{mode.upper()}: requested {len(candidates)}, extracted {found_count}, missing {len(missing)}")
        finally:
            temp.unlink(missing_ok=True)

    if errors:
        raise RuntimeError("Incomplete extraction: " + "; ".join(errors))
    print("All candidate FASTA files completed successfully.")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""BIOF3001 Project 7 — Mosquito RNA Virome.

Script: extract_filtered_contigs.py

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

# STEP 1 — Resolve paths relative to the repository, not the current working
# directory. This script lives two directories beneath the project root:
# scripts/viral_identification/extract_filtered_contigs.py.
ROOT = Path(__file__).resolve().parents[2]
FILTERED_DIR = ROOT / "results/viral_identification/filtered"
FASTA_DIR = ROOT / "data/raw/logan_contigs"
OUTDIR = ROOT / "results/viral_identification/host_filter"

# Process each DIAMOND candidate set separately so downstream BLASTN screening
# can compare results from the strict, loose and E-value-only filters.
MODES = ("strict", "loose", "evalue")


def write_selected_fasta(source: Path, sample: str, wanted: set[str], output) -> set[str]:
    """Stream a FASTA file once; write selected records and return found IDs."""
    # STEP 2 — Read each source FASTA sequentially, avoiding loading potentially
    # large Logan assemblies into memory. Keep track of IDs actually extracted.
    found = set()
    write_record = False
    current_id = None
    with source.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                # The first whitespace-delimited field of the FASTA header is
                # the contig ID used by the corresponding DIAMOND hit table.
                current_id = line[1:].strip().split()[0]

                # STEP 3 — Select a contig only when its ID is requested, and
                # prevent duplicate records with the same ID from being output.
                write_record = current_id in wanted and current_id not in found
                if write_record:
                    # Prefix the run accession because different SRRs may reuse
                    # contig IDs; BLASTN output must identify the source sample.
                    output.write(f">{sample}|{current_id}\n")
                    found.add(current_id)
            elif write_record:
                # Write the complete sequence, including multi-line sequences,
                # until the next FASTA header changes the selection state.
                output.write(line)
    return found


def main():
    # STEP 4 — Confirm that all three filtered-hit tables exist before starting.
    # A missing condition should not silently produce an incomplete analysis.
    tables = {mode: FILTERED_DIR / f"{mode}_filtered_hits.tsv" for mode in MODES}
    for mode, table in tables.items():
        if not table.is_file():
            raise FileNotFoundError(f"Missing {mode} table: {table}")

    OUTDIR.mkdir(parents=True, exist_ok=True)
    errors = []

    # STEP 5 — Extract candidates independently for each filtering condition.
    for mode, table in tables.items():
        # Keep IDs as strings and retain only the columns needed for extraction.
        # Multiple DIAMOND rows for the same sample/contig must count once.
        df = pd.read_csv(table, sep="\t", dtype={"sample": str, "qseqid": str})
        if not {"sample", "qseqid"}.issubset(df.columns):
            raise ValueError(f"Missing 'sample' or 'qseqid' column in {table}")
        candidates = df[["sample", "qseqid"]].dropna().drop_duplicates()

        # STEP 6 — Write to a temporary FASTA first. Only promote it to the
        # final filename after verifying that every requested contig was found.
        final = OUTDIR / f"{mode}_candidate_contigs.fa"
        temp = OUTDIR / f"{mode}_candidate_contigs.fa.tmp"
        missing_file = OUTDIR / f"{mode}_missing_contigs.tsv"
        missing = []
        found_count = 0
        try:
            with temp.open("w", encoding="utf-8") as out:
                # Group by SRR so its original assembly FASTA is scanned once,
                # even if several candidate contigs come from the same sample.
                for sample, group in candidates.groupby("sample", sort=True):
                    wanted = set(group["qseqid"])
                    fasta = FASTA_DIR / f"{sample}.contigs.fa"
                    if not fasta.is_file():
                        # STEP 7 — Audit missing source FASTAs: all contigs from
                        # that SRR are considered missing, not silently omitted.
                        print(f"WARNING: Missing contig FASTA: {fasta}")
                        missing.extend((sample, q) for q in sorted(wanted))
                        continue

                    found = write_selected_fasta(fasta, sample, wanted, out)
                    found_count += len(found)
                    # Report requested contig IDs absent from an existing FASTA.
                    missing.extend((sample, q) for q in sorted(wanted - found))

            if missing:
                # STEP 8 — Reject incomplete extraction and save an audit table.
                # Never expose a partial candidate FASTA as a completed input
                # to the subsequent host BLASTN workflow.
                pd.DataFrame(missing, columns=["sample", "qseqid"]).to_csv(
                    missing_file, sep="\t", index=False
                )
                temp.unlink(missing_ok=True)
                final.unlink(missing_ok=True)
                errors.append(f"{mode}: {len(missing)} missing contigs; see {missing_file}")
            else:
                # STEP 9 — Atomic rename publishes a complete FASTA; remove an
                # outdated missing-contigs report from an earlier failed run.
                temp.replace(final)
                missing_file.unlink(missing_ok=True)

            print(f"{mode.upper()}: requested {len(candidates)}, extracted {found_count}, missing {len(missing)}")
        finally:
            # Also remove temporary data if an unexpected exception occurs.
            temp.unlink(missing_ok=True)

    # STEP 10 — Report all incomplete modes after processing all three. A
    # nonzero exit prevents downstream automation treating them as successful.
    if errors:
        raise RuntimeError("Incomplete extraction: " + "; ".join(errors))
    print("All candidate FASTA files completed successfully.")


if __name__ == "__main__":
    main()

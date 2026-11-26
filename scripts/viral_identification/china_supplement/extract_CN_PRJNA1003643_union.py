#!/usr/bin/env python3
"""BIOF3001 Project 7 — China PRJNA1003643 candidate extraction.

Objective
---------
Extract the union of candidate viral contigs selected by the strict, loose,
and E-value-only DIAMOND filters. Each requested sample–contig pair is
written once, with a collision-safe FASTA identifier: ``sample|qseqid``.

Inputs (relative to repository root)
------------------------------------
results/viral_identification/PRJNA1003643/filtered/
    CN_PRJNA1003643_{strict,loose,evalue}_filtered_hits.tsv
    (required columns: sample, qseqid)
data/raw/logan_contigs/<sample>.contigs.fa

Outputs
-------
results/viral_identification/PRJNA1003643/host_filter/
    CN_PRJNA1003643_union_candidate_contigs.fa
    CN_PRJNA1003643_missing_candidate_contigs.tsv

Usage
-----
python scripts/viral_identification/china_supplement/extract_CN_PRJNA1003643_union.py

Dependencies: Python 3 standard library.


The extraction fails if any requested contigs cannot be found, to avoid
silently passing incomplete data to host screening. Edit LOGAN_DIR if this
China subset is stored in a different directory.
"""

import csv
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "results/viral_identification/PRJNA1003643"
FILTER_DIR = BASE / "filtered"
LOGAN_DIR = ROOT / "data/raw/logan_contigs"
OUT_DIR = BASE / "host_filter"
PREFIX = "CN_PRJNA1003643"
MODES = ("strict", "loose", "evalue")


def main() -> int:
    needed: dict[str, set[str]] = defaultdict(set)
    per_mode: dict[str, int] = {}

    for mode in MODES:
        table = FILTER_DIR / f"{PREFIX}_{mode}_filtered_hits.tsv"
        if not table.is_file():
            raise FileNotFoundError(f"Required filtered-hit table missing: {table}")
        pairs = set()
        with table.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not {"sample", "qseqid"}.issubset(reader.fieldnames or []):
                raise ValueError(f"Missing sample/qseqid column in {table}")
            for row in reader:
                sample, qid = row["sample"].strip(), row["qseqid"].strip()
                if not sample or not qid:
                    raise ValueError(f"Blank sample or qseqid in {table}")
                if "|" in sample:
                    raise ValueError(f"'|' in sample ID prevents unambiguous FASTA header: {sample}")
                pairs.add((sample, qid))
                needed[sample].add(qid)
        per_mode[mode] = len(pairs)
        print(f"{mode}: {len(pairs):,} unique sample-contig pairs")

    if not needed:
        raise ValueError("No candidate contigs found in the three filtered tables")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_fasta = OUT_DIR / f"{PREFIX}_union_candidate_contigs.fa"
    missing_file = OUT_DIR / f"{PREFIX}_missing_candidate_contigs.tsv"
    tmp_fasta = out_fasta.with_suffix(out_fasta.suffix + ".tmp")

    requested = sum(map(len, needed.values()))
    found: set[tuple[str, str]] = set()
    missing: list[tuple[str, str]] = []
    try:
        with tmp_fasta.open("w", encoding="utf-8") as output:
            for i, (sample, wanted) in enumerate(sorted(needed.items()), 1):
                fasta = LOGAN_DIR / f"{sample}.contigs.fa"
                if not fasta.is_file():
                    print(f"MISSING FASTA: {fasta}", file=sys.stderr)
                    missing.extend((sample, qid) for qid in sorted(wanted))
                    continue

                current_id = None
                capture = False
                with fasta.open(encoding="utf-8") as handle:
                    for line in handle:
                        if line.startswith(">"):
                            current_id = line[1:].split()[0]
                            capture = current_id in wanted and (sample, current_id) not in found
                            if capture:
                                output.write(f">{sample}|{current_id}\n")
                                found.add((sample, current_id))
                        elif capture:
                            output.write(line)
                missing.extend((sample, qid) for qid in sorted(wanted) if (sample, qid) not in found)
                if i % 100 == 0:
                    print(f"Processed {i:,}/{len(needed):,} samples")

        with missing_file.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(("sample", "qseqid"))
            writer.writerows(sorted(missing))

        print(f"Requested: {requested:,} | Extracted: {len(found):,} | Missing: {len(missing):,}")
        if missing:
            # Never keep a partial FASTA under the final output name.
            raise RuntimeError(f"Extraction incomplete; see {missing_file}")
        tmp_fasta.replace(out_fasta)
        print(f"Wrote: {out_fasta}")
        return 0
    finally:
        tmp_fasta.unlink(missing_ok=True)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

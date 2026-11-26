#!/usr/bin/env bash
# BIOF3001 Project 7 — Mosquito RNA Virome
# Script: download_all_logan_contigs.sh
# Responsible contributor: [confirm name]
# Objective: Download and decompress preassembled Logan contigs for all
#            accessions listed in results/run_selection/logan_available_runs.txt.
# Input: results/run_selection/logan_available_runs.txt (one SRA accession per line)
# Output: data/raw/logan_contigs/<accession>.contigs.fa
# Dependencies: bash, wget, zstd, grep
# Usage: bash scripts/viral_identification/download_all_logan_contigs.sh
# Notes: Completed nonempty FASTA files are skipped. Network failures are
#        reported and do not stop processing the remaining accessions.

set -uo pipefail
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
RUN_LIST="$ROOT_DIR/results/run_selection/logan_available_runs.txt"
OUTDIR="$ROOT_DIR/data/raw/logan_contigs"

for cmd in wget zstd grep; do
    command -v "$cmd" >/dev/null 2>&1 || { echo "ERROR: missing command $cmd" >&2; exit 1; }
done
[[ -f "$RUN_LIST" ]] || { echo "ERROR: no accession list: $RUN_LIST" >&2; exit 1; }
mkdir -p "$OUTDIR"

success=0; skipped=0; failed=0
while IFS= read -r acc || [[ -n "$acc" ]]; do
    acc="${acc%$'\r'}"
    [[ -z "$acc" || "$acc" == \#* ]] && continue
    if [[ ! "$acc" =~ ^(SRR|ERR|DRR)[0-9]+$ ]]; then
        echo "WARNING: invalid accession: $acc" >&2
        ((failed+=1)); continue
    fi

    archive="$OUTDIR/${acc}.contigs.fa.zst"
    fasta="$OUTDIR/${acc}.contigs.fa"
    temporary="$OUTDIR/${acc}.contigs.fa.tmp"
    url="https://s3.amazonaws.com/logan-pub/c/${acc}/${acc}.contigs.fa.zst"
    echo "Processing $acc"

    # Skip previous successful downloads, avoiding repeated large transfers.
    if [[ -s "$fasta" ]]; then
        echo "  Existing nonempty FASTA; skipping"
        ((skipped+=1)); continue
    fi
    rm -f "$temporary"

    # Continue partial downloads when possible; decompress only after a valid archive.
    if ! wget -c "$url" -O "$archive"; then
        echo "WARNING: download failed: $acc" >&2
        ((failed+=1)); continue
    fi
    if ! zstd -t "$archive" >/dev/null 2>&1; then
        echo "WARNING: invalid compressed archive: $acc" >&2
        rm -f "$archive"
        ((failed+=1)); continue
    fi
    if ! zstd -dc "$archive" > "$temporary"; then
        echo "WARNING: decompression failed: $acc" >&2
        rm -f "$temporary"
        ((failed+=1)); continue
    fi

    # A FASTA header begins with '>'; check presence, not full FASTA validity.
    contigs=$(grep -c '^>' "$temporary" || true)
    if [[ "$contigs" -eq 0 ]]; then
        echo "WARNING: no FASTA headers: $acc" >&2
        rm -f "$temporary"
        ((failed+=1)); continue
    fi
    mv -f "$temporary" "$fasta"
    echo "  Downloaded $acc: $contigs contig records"
    ((success+=1))
done < "$RUN_LIST"

echo "Summary: downloaded=$success, skipped=$skipped, failed=$failed"
if (( failed > 0 )); then
    echo "WARNING: some runs failed; review messages above" >&2
    exit 1
fi

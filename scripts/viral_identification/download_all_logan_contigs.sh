#!/usr/bin/env bash
# BIOF3001 Project 7 — Mosquito RNA Virome
# Script: download_all_logan_contigs.sh
# Objective: Download and decompress Logan-assembled nucleotide contigs for
#            SRA accessions selected by the run-selection stage.
# Input:  results/run_selection/logan_available_runs.txt (one accession per line)
# Output: data/raw/logan_contigs/<accession>.contigs.fa
# Requirements: bash, wget, zstd, grep; network access to the Logan public S3 bucket.
# Usage (from any directory): bash scripts/viral_identification/download_all_logan_contigs.sh
# Notes: Existing nonempty FASTA files are skipped. Download errors are logged,
#        remaining runs are attempted, and a nonzero status signals any failure.

set -uo pipefail

# STEP 1 — Resolve the repository root from this script's own location.
# This is independent of the directory where the user runs the script.
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
RUN_LIST="$ROOT_DIR/results/run_selection/logan_available_runs.txt"
OUTDIR="$ROOT_DIR/data/raw/logan_contigs"

# STEP 2 — Verify software and the run list before creating/downloading data.
for cmd in wget zstd grep; do
    command -v "$cmd" >/dev/null 2>&1 || {
        echo "ERROR: missing command $cmd" >&2
        exit 1
    }
done
[[ -f "$RUN_LIST" ]] || { echo "ERROR: missing accession list: $RUN_LIST" >&2; exit 1; }
mkdir -p "$OUTDIR"

success=0
skipped=0
failed=0

# STEP 3 — Read the run-selection accession list one line at a time.
# Support a final line with no newline, Windows CRLF endings, and comments.
while IFS= read -r acc || [[ -n "$acc" ]]; do
    acc="${acc%$'\r'}"
    [[ -z "$acc" || "$acc" == \#* ]] && continue

    # Prevent malformed accessions from producing unintended paths/URLs.
    if [[ ! "$acc" =~ ^(SRR|ERR|DRR)[0-9]+$ ]]; then
        echo "WARNING: invalid accession: $acc" >&2
        ((failed+=1))
        continue
    fi

    archive="$OUTDIR/${acc}.contigs.fa.zst"
    fasta="$OUTDIR/${acc}.contigs.fa"
    temporary="$OUTDIR/${acc}.contigs.fa.tmp"
    url="https://s3.amazonaws.com/logan-pub/c/${acc}/${acc}.contigs.fa.zst"
    echo "Processing $acc"

    # STEP 4 — Avoid repeated transfers for already-decompressed runs.
    # A nonempty file is a practical completion check, not exhaustive validation.
    if [[ -s "$fasta" ]]; then
        echo "  Existing nonempty FASTA; skipping"
        ((skipped+=1))
        continue
    fi
    rm -f "$temporary"

    # STEP 5 — Download the compressed Logan FASTA.
    # -c asks wget to resume interrupted downloads when supported by server.
    # Reuse a complete local archive rather than downloading it again.
    if [[ -s "$archive" ]] && zstd -t "$archive" >/dev/null 2>&1; then
        echo "  Reusing verified compressed archive"
    else
        if ! wget -c "$url" -O "$archive"; then
            echo "WARNING: download failed: $acc" >&2
            ((failed+=1))
            continue
        fi
    fi

    # STEP 6 — Check archive integrity before attempting decompression.
    if ! zstd -t "$archive" >/dev/null 2>&1; then
        echo "WARNING: invalid compressed archive: $acc" >&2
        rm -f "$archive"
        ((failed+=1))
        continue
    fi

    # STEP 7 — Decompress to a temporary file to avoid exposing partial FASTA.
    if ! zstd -dc "$archive" > "$temporary"; then
        echo "WARNING: decompression failed: $acc" >&2
        rm -f "$temporary"
        ((failed+=1))
        continue
    fi

    # STEP 8 — Perform a minimal FASTA sanity check and publish atomically.
    # Checking '>' headers is not equivalent to validating every sequence line.
    contigs=$(grep -c '^>' "$temporary" || true)
    if [[ "$contigs" -eq 0 ]]; then
        echo "WARNING: no FASTA headers: $acc" >&2
        rm -f "$temporary"
        ((failed+=1))
        continue
    fi
    mv -f "$temporary" "$fasta"

    # STEP 9 — Remove the redundant compressed copy after successful extraction.
    # Delete this line if your project policy requires keeping source archives.
    rm -f "$archive"
    echo "  Downloaded $acc: $contigs contig records"
    ((success+=1))
done < "$RUN_LIST"

# STEP 10 — Report outcomes; nonzero exit prevents treating partial downloads
# as a fully successful reproducible pipeline stage.
echo "Summary: downloaded=$success, skipped=$skipped, failed=$failed"
if (( failed > 0 )); then
    echo "WARNING: some runs failed; review messages above" >&2
    exit 1
fi

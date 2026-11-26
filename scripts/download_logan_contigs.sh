
#!/usr/bin/env bash

# script: download_logan_contigs.sh
# contributor: Cammy
#
# Objective:
# Download and decompress preassembled Logan contigs for
# selected mosquito-associated SRA sequencing runs.
#
# Input:
#   runs/pilot_runs.txt
#   A text file containing one SRA run accession per line.
#
# Output:
#   data/logan_contigs/<accession>.contigs.fa
#   Decompressed FASTA files containing assembled contigs.
#
# Requirements:
#   Bash, wget, zstd, grep
#
# Notes:
#   - Downloads contigs from the public Logan S3 repository.
#   - Skips existing nonempty FASTA files.
#   - Reports unsuccessful downloads and decompression errors.
#   - Uses temporary files to avoid retaining partial FASTA files.
#   - This version processes the pilot accession list.
# ============================================================

set -uo pipefail

# Determine the repository root, regardless of the current
# working directory.
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

RUN_LIST="$ROOT_DIR/runs/pilot_runs.txt"
OUTDIR="$ROOT_DIR/data/logan_contigs"

# Check that the required command-line tools are installed.
for cmd in wget zstd grep; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: Required program '$cmd' not found." >&2
        exit 1
    fi
done

# Verify that the accession list exists.
if [ ! -f "$RUN_LIST" ]; then
    echo "ERROR: Input file not found: $RUN_LIST" >&2
    exit 1
fi

# Create the output directory if necessary.
mkdir -p "$OUTDIR"

# Read and process one SRA run accession at a time.
while IFS= read -r acc || [ -n "$acc" ]; do

    # Remove carriage returns from Windows-formatted files.
    acc="${acc%$'\r'}"

    # Ignore blank lines and comments.
    [[ -z "$acc" || "$acc" == \#* ]] && continue

    # Validate the SRA run accession format.
    if [[ ! "$acc" =~ ^(SRR|ERR|DRR)[0-9]+$ ]]; then
        echo "WARNING: Invalid accession '$acc'. Skipping." >&2
        continue
    fi

    echo "======================================"
    echo "Processing $acc"
    echo "======================================"

    # Define filenames for the compressed and decompressed data.
    ZST="$OUTDIR/${acc}.contigs.fa.zst"
    FA="$OUTDIR/${acc}.contigs.fa"
    TMP_FA="$OUTDIR/${acc}.contigs.fa.tmp"

    # Avoid downloading data already available locally.
    if [ -s "$FA" ]; then
        echo "Existing FASTA file found. Skipping $acc."
        continue
    fi

    # Construct the Logan S3 download URL.
    URL="https://s3.amazonaws.com/logan-pub/c/${acc}/${acc}.contigs.fa.zst"

    echo "Downloading Logan contigs..."

    # Download with resume support for interrupted transfers.
    if ! wget -c "$URL" -O "$ZST"; then
        echo "WARNING: Download failed for $acc." >&2
        continue
    fi

    # Verify the integrity of the compressed Zstandard file.
    if ! zstd -t "$ZST" >/dev/null 2>&1; then
        echo "WARNING: Invalid or incomplete archive for $acc." >&2
        rm -f "$ZST"
        continue
    fi

    echo "Decompressing contigs..."

    # Decompress into a temporary file to prevent incomplete
    # FASTA files from being mistaken for completed outputs.
    if ! zstd -dc "$ZST" > "$TMP_FA"; then
        echo "WARNING: Decompression failed for $acc." >&2
        rm -f "$TMP_FA"
        continue
    fi

    # Count FASTA records by identifying header lines.
    CONTIG_COUNT=$(grep -c '^>' "$TMP_FA" || true)

    # Reject empty files or files without FASTA headers.
    if [ "$CONTIG_COUNT" -eq 0 ]; then
        echo "WARNING: No contigs found for $acc." >&2
        rm -f "$TMP_FA"
        continue
    fi

    # Move the validated output to its final filename.
    mv "$TMP_FA" "$FA"

    echo "Number of contigs: $CONTIG_COUNT"
    echo "Successfully processed $acc"

done < "$RUN_LIST"

echo "======================================"
echo "Logan contig download process completed."
echo "Output directory: $OUTDIR"

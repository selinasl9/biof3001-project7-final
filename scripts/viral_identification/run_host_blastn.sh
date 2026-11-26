
#!/usr/bin/env bash
# BIOF3001 Project 7 — Mosquito RNA Virome
# Script: run_host_blastn.sh
#
# Objective:
#   Align candidate nucleotide contigs from each DIAMOND filtering branch
#   against a Culicidae nucleotide BLAST database to identify strong
#   mosquito-reference matches for subsequent filtering.
#
# Inputs:
#   results/viral_identification/host_filter/{mode}_candidate_contigs.fa
#   data/reference/culicidae_blast/culicidae.* (BLAST database files)
#
# Outputs:
#   results/viral_identification/host_filter/blastn/{mode}_vs_culicidae.tsv
#   results/viral_identification/host_filter/logs/{mode}_blastn.log
#
# Requirements: Bash, NCBI blastn, blastdbcmd
# Run: bash scripts/viral_identification/run_host_blastn.sh
#
# A zero-byte BLAST result is valid when the search succeeds but has no hits.
# The downstream removal script must still check that the result exists.


# STEP 1: Enable strict error checking
# -e: Stop if a command fails.
# -u: Stop for undefined variables.
# pipefail: Detect errors in pipelines.
set -euo pipefail


# STEP 2: Set file paths and CPU threads

# Find project root (two directories above this script).
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

# Culicidae nucleotide reference database.
DB="$ROOT_DIR/data/reference/culicidae_blast/culicidae"

# Directory containing candidate contig FASTA files.
HOST_DIR="$ROOT_DIR/results/viral_identification/host_filter"

# Directories for BLASTN results and logs.
OUTDIR="$HOST_DIR/blastn"
LOGDIR="$HOST_DIR/logs"

# Use 8 CPU threads.
THREADS=8


# STEP 3: Check required software

# Check whether blastn and blastdbcmd are installed.
# Stop if either command is unavailable.
for cmd in blastn blastdbcmd; do
    command -v "$cmd" >/dev/null 2>&1 || { echo "ERROR: $cmd not installed" >&2; exit 1; }
done


# STEP 4: Check the Culicidae reference database

# Verify that the nucleotide BLAST database is readable.
# Stop if the database is missing or cannot be opened.
blastdbcmd -db "$DB" -dbtype nucl -info >/dev/null 2>&1 || {
    echo "ERROR: Culicidae nucleotide BLAST database missing or unreadable: $DB" >&2
    exit 1
}

# Create output directories if needed.
mkdir -p "$OUTDIR" "$LOGDIR"


# STEP 5: Process three filtering branches

# Run host BLASTN separately for strict, loose and evalue.
for mode in strict loose evalue; do

    # Define input FASTA for the current branch.
    QUERY="$HOST_DIR/${mode}_candidate_contigs.fa"

    # Define final BLASTN result file.
    OUT="$OUTDIR/${mode}_vs_culicidae.tsv"

    # Temporary output prevents saving incomplete results.
    TMP="$OUT.tmp"

    # Log file records BLASTN messages and errors.
    LOG="$LOGDIR/${mode}_blastn.log"

    echo "Host BLASTN: $mode"


    # STEP 6: Check input FASTA

    # Stop if the candidate FASTA file is missing.
    if [[ ! -f "$QUERY" ]]; then
        echo "ERROR: Missing extracted candidate FASTA: $QUERY" >&2
        exit 1
    fi

    # Remove any old temporary output.
    rm -f "$TMP"


    # STEP 7: Handle empty candidate files

    # If no candidate sequences exist, create an empty result.
    # Skip BLASTN and move to the next branch.
    if [[ ! -s "$QUERY" ]]; then
        : > "$OUT"
        printf 'No candidate contigs for %s; created empty result.\n' "$mode" > "$LOG"
        continue
    fi


    # STEP 8: Run BLASTN against the Culicidae database

    # Compare candidate nucleotide contigs with mosquito references.
    # Require at least 80% nucleotide identity for reported matches.
    # Save 12 alignment fields in a tab-separated file.
    if blastn \
        -task blastn \
        -query "$QUERY" \
        -db "$DB" \
        -num_threads "$THREADS" \
        -perc_identity 80 \
        -outfmt '6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore' \
        -out "$TMP" > "$LOG" 2>&1; then


        # STEP 9: Save successful BLASTN results

        # Rename temporary file to final output.
        # An empty file is valid if BLASTN finds no matches.
        mv "$TMP" "$OUT"

        # Print success message.
        echo "SUCCESS: $mode -> $OUT"

    else

        # STEP 10: Handle failed BLASTN searches

        # Remove incomplete results and stop execution.
        rm -f "$TMP" "$OUT"
        echo "ERROR: Host BLASTN failed for $mode; see $LOG" >&2
        exit 1
    fi

# Continue until all three branches have been processed.
done


# STEP 11: Confirm completion

# The resulting TSV files are used by remove_mosquito_contigs.py
# to remove contigs with strong Culicidae reference matches
# (identity >= 80% AND query alignment coverage >= 80%).
echo 'All three host BLASTN searches completed.'

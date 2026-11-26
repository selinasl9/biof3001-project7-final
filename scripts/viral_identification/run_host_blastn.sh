#!/usr/bin/env bash
# BIOF3001 Project 7 — Mosquito RNA Virome
# Script: run_host_blastn.sh
# Responsible contributor: [Confirm name before submission]
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

set -euo pipefail
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
DB="$ROOT_DIR/data/reference/culicidae_blast/culicidae"
HOST_DIR="$ROOT_DIR/results/viral_identification/host_filter"
OUTDIR="$HOST_DIR/blastn"
LOGDIR="$HOST_DIR/logs"
THREADS=8

for cmd in blastn blastdbcmd; do
    command -v "$cmd" >/dev/null 2>&1 || { echo "ERROR: $cmd not installed" >&2; exit 1; }
done
blastdbcmd -db "$DB" -dbtype nucl -info >/dev/null 2>&1 || {
    echo "ERROR: Culicidae nucleotide BLAST database missing or unreadable: $DB" >&2
    exit 1
}
mkdir -p "$OUTDIR" "$LOGDIR"

for mode in strict loose evalue; do
    QUERY="$HOST_DIR/${mode}_candidate_contigs.fa"
    OUT="$OUTDIR/${mode}_vs_culicidae.tsv"
    TMP="$OUT.tmp"
    LOG="$LOGDIR/${mode}_blastn.log"
    echo "Host BLASTN: $mode"

    if [[ ! -f "$QUERY" ]]; then
        echo "ERROR: Missing extracted candidate FASTA: $QUERY" >&2
        exit 1
    fi
    rm -f "$TMP"

    # If there are no candidates, there are no host alignments to report.
    if [[ ! -s "$QUERY" ]]; then
        : > "$OUT"
        printf 'No candidate contigs for %s; created empty result.\n' "$mode" > "$LOG"
        continue
    fi

    if blastn \
        -task blastn \
        -query "$QUERY" \
        -db "$DB" \
        -num_threads "$THREADS" \
        -perc_identity 80 \
        -outfmt '6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore' \
        -out "$TMP" > "$LOG" 2>&1; then
        mv "$TMP" "$OUT"
        echo "SUCCESS: $mode -> $OUT"
    else
        rm -f "$TMP" "$OUT"
        echo "ERROR: Host BLASTN failed for $mode; see $LOG" >&2
        exit 1
    fi
done

echo 'All three host BLASTN searches completed.'

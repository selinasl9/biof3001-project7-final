#!/usr/bin/env bash
# BIOF3001 Project 7 — Mosquito RNA Virome
# Script: run_diamond_all.sh

# Objective: Search assembled Logan contigs against an arbovirus protein
#            reference database using DIAMOND BLASTx, one run at a time.
# Input: results/viral_identification/diamond_available_runs.txt; data/raw/logan_contigs/<run>.contigs.fa
# Database: data/reference/arbovirus_refseq (DIAMOND database prefix)
# Output: results/viral_identification/diamond/<run>.blastx.tsv (17 columns, simple header)
#         results/viral_identification/timing/<run>.time.txt; results/viral_identification/logs/<run>.diamond.log
#         results/viral_identification/diamond/<run>.blastx.done (success marker, including no-hit)
# Dependencies: bash, diamond, GNU /usr/bin/time
# Usage: bash scripts/viral_identification/run_diamond_all.sh
# Note: Legacy nonempty output files are skipped as in the original script.
#       Their completion cannot be independently verified without logs.

set -uo pipefail
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
RUNLIST="$ROOT_DIR/results/viral_identification/diamond_available_runs.txt"
DB="$ROOT_DIR/data/reference/arbovirus_refseq"
THREADS="${THREADS:-8}"
OUTDIR="$ROOT_DIR/results/viral_identification/diamond"
TIMEDIR="$ROOT_DIR/results/viral_identification/timing"
LOGDIR="$ROOT_DIR/results/viral_identification/logs"

command -v diamond >/dev/null 2>&1 || { echo "ERROR: diamond not found" >&2; exit 1; }
[[ -x /usr/bin/time ]] || { echo "ERROR: GNU /usr/bin/time not found" >&2; exit 1; }
[[ -f "$RUNLIST" ]] || { echo "ERROR: missing run list: $RUNLIST" >&2; exit 1; }
[[ -f "$DB.dmnd" ]] || { echo "ERROR: DIAMOND database not found: $DB.dmnd" >&2; exit 1; }
mkdir -p "$OUTDIR" "$TIMEDIR" "$LOGDIR"

success=0; nohits=0; skipped=0; missing=0; failed=0
while IFS= read -r acc || [[ -n "$acc" ]]; do
    acc="${acc%$'\r'}"
    [[ -z "$acc" || "$acc" == \#* ]] && continue
    if [[ ! "$acc" =~ ^(SRR|ERR|DRR)[0-9]+$ ]]; then
        echo "WARNING: invalid accession: $acc" >&2
        ((failed+=1)); continue
    fi

    query="$ROOT_DIR/data/raw/logan_contigs/${acc}.contigs.fa"
    out="$OUTDIR/${acc}.blastx.tsv"
    temporary="$OUTDIR/${acc}.blastx.tsv.tmp"
    done_marker="$OUTDIR/${acc}.blastx.done"
    timefile="$TIMEDIR/${acc}.time.txt"
    logfile="$LOGDIR/${acc}.diamond.log"
    echo "Processing $acc"

    if [[ ! -s "$query" ]]; then
        echo "WARNING: missing/empty input FASTA: $acc" >&2
        ((missing+=1)); continue
    fi
    if [[ -f "$done_marker" && -f "$out" ]]; then
        echo "  Already completed (status marker); skipping"
        ((skipped+=1)); continue
    fi
    if [[ -s "$out" && ! -e "$done_marker" ]]; then
        echo "  Legacy nonempty output present; skipping (check old log for completion)"
        ((skipped+=1)); continue
    fi

    # Write to a temporary file so interrupted runs cannot be mistaken for completed hits.
    rm -f "$temporary" "$timefile"
    if /usr/bin/time -v -o "$timefile" \
        diamond blastx \
        --query "$query" \
        --db "$DB" \
        --out "$temporary" \
        --threads "$THREADS" \
        --sensitive \
        --evalue 1e-5 \
        --max-target-seqs 5 \
        --header simple \
        --outfmt 6 \
        qseqid sseqid pident length qlen slen \
        qstart qend sstart send evalue bitscore \
        qcovhsp scovhsp staxids sscinames stitle \
        > "$logfile" 2>&1; then
        # An empty table is a valid successful zero-hit search.
        [[ -f "$temporary" ]] || : > "$temporary"
        mv -f "$temporary" "$out"
        printf 'success\n' > "$done_marker"
        if [[ -s "$out" ]]; then
            echo "  SUCCESS: matches found"
            ((success+=1))
        else
            echo "  SUCCESS: zero hits"
            ((nohits+=1))
        fi
    else
        echo "WARNING: DIAMOND failed for $acc; see $logfile" >&2
        rm -f "$temporary" "$done_marker"
        ((failed+=1))
    fi
done < "$RUNLIST"

echo "Summary: hits=$success, zero_hits=$nohits, skipped=$skipped, missing_queries=$missing, failed=$failed"
if (( failed > 0 || missing > 0 )); then
    exit 1
fi

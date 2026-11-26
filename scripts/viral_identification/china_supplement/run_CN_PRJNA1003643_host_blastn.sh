#!/usr/bin/env bash
# BIOF3001 Project 7 — China PRJNA1003643 host BLASTN screening

#
# Objective:
# Compare the union of China candidate viral contigs against a Culicidae
# nucleotide BLAST database to identify potential host-associated matches.
#
# Input:
# results/viral_identification/PRJNA1003643/host_filter/
#   CN_PRJNA1003643_union_candidate_contigs.fa
# data/reference/culicidae_blast/culicidae.* (BLAST nucleotide DB)
#
# Output:
# results/viral_identification/PRJNA1003643/host_filter/blastn/
#   CN_PRJNA1003643_union_vs_culicidae.tsv
#   CN_PRJNA1003643_union.blastn.log
#   CN_PRJNA1003643_union.time.txt
#
# Usage from any working directory:
# bash scripts/viral_identification/china_supplement/run_CN_PRJNA1003643_host_blastn.sh
#
# Dependencies: bash, BLAST+ (blastn and blastdbcmd), /usr/bin/time optional.
# Note: 80% nucleotide identity is a search constraint; final host exclusion
# also requires >=80% query coverage in the downstream cleaning script.

set -euo pipefail
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
BASE="$ROOT_DIR/results/viral_identification/PRJNA1003643/host_filter"
QUERY="$BASE/CN_PRJNA1003643_union_candidate_contigs.fa"
DB="$ROOT_DIR/data/reference/culicidae_blast/culicidae"
OUTDIR="$BASE/blastn"
OUT="$OUTDIR/CN_PRJNA1003643_union_vs_culicidae.tsv"
LOG="$OUTDIR/CN_PRJNA1003643_union.blastn.log"
TIMEFILE="$OUTDIR/CN_PRJNA1003643_union.time.txt"
THREADS="${THREADS:-8}"

for cmd in blastn blastdbcmd; do
    command -v "$cmd" >/dev/null 2>&1 || { echo "ERROR: $cmd not found" >&2; exit 1; }
done
[[ -s "$QUERY" ]] || { echo "ERROR: candidate FASTA missing/empty: $QUERY" >&2; exit 1; }
blastdbcmd -db "$DB" -dbtype nucl -info >/dev/null || {
    echo "ERROR: Culicidae nucleotide BLAST database unavailable: $DB" >&2
    exit 1
}
mkdir -p "$OUTDIR"
TMP_OUT="${OUT}.tmp"
trap 'rm -f "$TMP_OUT"' EXIT

echo "Starting China Culicidae BLASTN: $(date)"
BLAST_ARGS=(
    -task blastn -query "$QUERY" -db "$DB"
    -num_threads "$THREADS" -perc_identity 80
    -outfmt '6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore'
    -out "$TMP_OUT"
)
if [[ -x /usr/bin/time ]]; then
    /usr/bin/time -v -o "$TIMEFILE" blastn "${BLAST_ARGS[@]}" >"$LOG" 2>&1
else
    blastn "${BLAST_ARGS[@]}" >"$LOG" 2>&1
fi
# Empty BLASTN output is valid: no reference matches at this threshold.
mv "$TMP_OUT" "$OUT"
echo "Host BLASTN complete: $(date)"
echo "Output: $OUT"

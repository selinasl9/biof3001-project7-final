#!/usr/bin/env bash
# BIOF3001 Project 7 — 70% VSEARCH clustering
# Objective: group strict host-cleaned viral contigs at 70% nucleotide identity
# across four query/target coverage settings; compare with 90% CD-HIT analysis.
# Input: results/viral_clustering/strict_viral_contigs.fa (generated upstream)
# Outputs: results/viral_clustering/vsearch/strict_vsearch70*.fa and *.uc
# Requirements: vsearch, Bash. Usage: bash scripts/viral_clustering/run_vsearch_70.sh
# Environment: THREADS (default 8), FORCE=1 to rerun.
# VSEARCH query/target coverage differs from CD-HIT shorter/longer definitions.
# Responsible contributor: [confirm before submission]
# BIOF3001 Project 7: 70% identity clustering using VSEARCH, not CD-HIT-EST.
# Different clustering algorithms can introduce sensitivity unrelated to identity.
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT_DIR"
INPUT="${1:-results/viral_clustering/strict_viral_contigs.fa}"
OUT="${2:-results/viral_clustering/vsearch}"
THREADS="${THREADS:-8}"
[[ -s "$INPUT" ]] || { echo "Missing FASTA: $INPUT" >&2; exit 1; }
command -v vsearch >/dev/null || { echo 'vsearch is not installed' >&2; exit 1; }
mkdir -p "$OUT"
for cov in identity aS60 aL60 aS60_aL60; do
  suffix='';opts=()
  case "$cov" in
    identity) ;;
    aS60) suffix='_aS60';opts=(--query_cov 0.60);;
    aL60) suffix='_aL60';opts=(--target_cov 0.60);;
    aS60_aL60) suffix='_aS60_aL60';opts=(--query_cov 0.60 --target_cov 0.60);;
  esac
  stem="$OUT/strict_vsearch70${suffix}"
  if [[ -s "$stem.uc" && "${FORCE:-0}" != "1" ]]; then
    echo "Skipping existing $stem.uc (FORCE=1 to rerun)";continue
  fi
  echo "Clustering 70% VSEARCH; $cov"
  # NOTE: query/target coverage in VSEARCH is not guaranteed to map 1-to-1
  # onto CD-HIT-EST's shorter/longer sequence definitions.
  vsearch --cluster_fast "$INPUT" --id 0.70 --strand both \
    --threads "$THREADS" --centroids "$stem.fa" --uc "$stem.uc" "${opts[@]}"
done

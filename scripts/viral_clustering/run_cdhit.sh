#!/usr/bin/env bash
# BIOF3001 Project 7 — 90% CD-HIT-EST clustering
# Objective: cluster strict host-cleaned viral contigs using 90% nucleotide identity
# under four coverage settings (identity, aS60, aL60, aS60_aL60).
# Input: results/viral_clustering/strict_viral_contigs.fa (must be generated upstream)
# Outputs: results/viral_clustering/cdhit/strict_cdhit90*.fa and *.fa.clstr
# Requirements: cd-hit-est, Bash. Usage: bash scripts/viral_clustering/run_cdhit.sh
# Environment: THREADS (default 8), MEM_MB (default 64000), FORCE=1 to rerun.
# Run location: any directory; paths are relative to the repository root.
# Responsible contributor: [confirm before submission]
# BIOF3001 Project 7: 90% nucleotide sequence clustering (CD-HIT-EST)
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT_DIR"

INPUT="${1:-results/viral_clustering/strict_viral_contigs.fa}"
OUT="${2:-results/viral_clustering/cdhit}"
THREADS="${THREADS:-8}"
MEM_MB="${MEM_MB:-64000}"
[[ -s "$INPUT" ]] || { echo "Missing FASTA: $INPUT" >&2; exit 1; }
command -v cd-hit-est >/dev/null || { echo 'cd-hit-est is not installed' >&2; exit 1; }
mkdir -p "$OUT"
# CD-HIT-EST nucleotide clustering does NOT reliably support 70% identity.
# See separate run_vsearch_70.sh for 70%.
# Coverage conditions: baseline; shorter sequence >=60%; longer sequence >=60%; both.
for id in 90; do
  cutoff=0.90; word=8
  for cov in identity aS60 aL60 aS60_aL60; do
    name="strict_cdhit${id}"
    args=()
    case "$cov" in
      identity) ;;
      aS60) name+="_aS60"; args=(-aS 0.60);;
      aL60) name+="_aL60"; args=(-aL 0.60);;
      aS60_aL60) name+="_aS60_aL60"; args=(-aS 0.60 -aL 0.60);;
    esac
    dest="$OUT/${name}.fa"
    if [[ -s "$dest" && -s "$dest.clstr" && "${FORCE:-0}" != "1" ]]; then
      echo "Skipping existing $dest (FORCE=1 to rerun)"; continue
    fi
    printf '\nClustering %s: identity=%s, coverage=%s\n' "$name" "$cutoff" "$cov"
    cd-hit-est -i "$INPUT" -o "$dest" -c "$cutoff" -n "$word" \
      -g 1 -r 1 -p 1 -d 0 -T "$THREADS" -M "$MEM_MB" "${args[@]}"
    test -s "$dest.clstr" || { echo "Missing .clstr after run" >&2; exit 1; }
  done
done

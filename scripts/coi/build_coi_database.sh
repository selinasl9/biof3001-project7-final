#!/usr/bin/env bash
#
# build_coi_database.sh
#
# Build the Culicidae COI BLAST database that identify_host_from_coi.py matches runs against.
#
# Inputs:
#   NCBI nuccore                     Culicidae COI sequences, via EDirect
#
# Outputs:
#   reference/coi/culicidae_coi.fa   the COI sequences
#   reference/coi/culicidae_coi.n*   the BLAST database
#
set -euo pipefail

# Set output paths, the COI length window + the NCBI query
# 400-2000 bp keeps the barcode + full-length COI, leaves out whole mitochondrial genomes
OUTDIR=reference/coi
FASTA="$OUTDIR/culicidae_coi.fa"
DB="$OUTDIR/culicidae_coi"

MIN_LENGTH=400
MAX_LENGTH=2000

QUERY="Culicidae[Organism] AND (COI[Gene] OR COX1[Gene] OR \"cytochrome c oxidase subunit 1\"[All Fields]) AND ${MIN_LENGTH}:${MAX_LENGTH}[SLEN]"

# Step 1: fetch the COI sequences from NCBI
mkdir -p "$OUTDIR"

esearch -db nuccore -query "$QUERY" | efetch -format fasta > "$FASTA"

# Step 2: count the sequences, stop if NCBI returned none
count=$(grep -c '^>' "$FASTA" || true)
if [ "$count" -eq 0 ]; then
    echo "no sequences came back from NCBI, not building a database" >&2
    exit 1
fi
echo "$count COI sequences"

# Step 3: build the nucleotide BLAST database for identify_host_from_coi.py
makeblastdb -in "$FASTA" -dbtype nucl -parse_seqids -out "$DB"
echo "database written to $DB"

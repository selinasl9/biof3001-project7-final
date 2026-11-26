#!/usr/bin/env bash
#
# run_alignment.sh
#
# Run every MAFFT stage of the project, as one script with a stage argument.
#
# Inputs:
#   STAGE      which stage to run, see below
#   THREADS    threads for mafft and cd-hit-est, default 8
#
# Outputs:
#   results/msa/SET/   the alignment and logs for that stage, printed as it is written
#
# Stages, in the order they were run:
#   refs       align the reference genomes on their own
#   contigs    cluster the contigs at 95% then add them as fragments  (strict set)
#   complete   the 5-11.5 kb window: references, then cluster members added
#   scratch    everything aligned from scratch, refs + members + BLAST downloads
#              -> this is the alignment the v2/v3 trees were built from
#   add        add the BLAST downloads onto the existing alignment, --keeplength
#   compare    try every trimming cutoff and report the columns each one leaves
#   trim       trimAl -gt 0.8, and print the IQ-TREE command to run next
#   china      keep the references and the China contigs only
#
set -euo pipefail

STAGE="${1:?usage: run_alignment.sh refs|contigs|complete|scratch|add|compare|trim|china [threads]}"
THREADS="${2:-8}"

# Set the paths + the fixed parameters
MODE="${MODE:-strict}"
REFS_FASTA="${REFS_FASTA:-reference/flaviviridae_refseq_nt.fasta}"
CONTIG_DIR="results/viral_identification/viral_contigs"
MSA_DIR="results/msa/${MODE}"
COMPLETE_DIR="results/msa/complete"
TREE_DIR="results/msa/tree_v2"
ANNOTATION="${ANNOTATION:-$COMPLETE_DIR/tree_gt08/contig_annotation.tsv}"

CDHIT_IDENTITY=0.95
TRIM_OCCUPANCY=0.8
# Cutoffs the compare stage tries, -gt N keeping a column N of the sequences cover,
# -gappyout and -automated1 picking the cutoff from the gap distribution themselves
TRIM_METHODS=("-gappyout" "-automated1" "-gt 0.3" "-gt 0.5" "-gt 0.8")
# The three references on very long branches, which drag the midpoint root
DROP="${DROP-NC_040589.1 NC_040555.1 NC_003996.1}"

# Stop with a clear message if a stage's tools are not installed
need() {
    for tool in "$@"; do
        command -v "$tool" >/dev/null || { echo "$tool not found on PATH" >&2; exit 1; }
    done
}

# Report the sequences and columns of a finished alignment, so a stage says what it made
report() {
    python3 - "$1" <<'PY'
import sys
path = sys.argv[1]
lengths, width = 0, 0
for line in open(path):
    if line.startswith(">"):
        lengths += 1
        if lengths == 2:
            break
    else:
        width += len(line.strip())
print(f"  {path}: {sum(1 for l in open(path) if l.startswith('>'))} sequences "
      f"x {width} columns")
PY
}

case "$STAGE" in

# Step 1: align the references on their own, which every --add stage builds on
refs)
    need mafft
    mkdir -p "$MSA_DIR"
    mafft --auto --adjustdirectionaccurately --reorder --thread "$THREADS" \
        "$REFS_FASTA" > "$MSA_DIR/refs_aligned.fasta" 2> "$MSA_DIR/mafft_refs.log"
    tail -3 "$MSA_DIR/mafft_refs.log" || true
    report "$MSA_DIR/refs_aligned.fasta"
    ;;

# Step 2: cluster the contig set, as 10,388 contigs is too many to align, then add
# them as fragments, as most are far shorter than a genome
contigs)
    need mafft cd-hit-est
    mkdir -p "$MSA_DIR"
    combined="$CONTIG_DIR/${MODE}_viral_contigs.fa"
    clustered="$CONTIG_DIR/${MODE}_cdhit95.fa"

    cd-hit-est -i "$combined" -o "$clustered" -c "$CDHIT_IDENTITY" -n 10 \
        -d 0 -M 0 -T "$THREADS" > "$MSA_DIR/cdhit.log" 2>&1
    echo "  $(grep -c '^>' "$clustered") representatives at ${CDHIT_IDENTITY}"

    [ -s "$MSA_DIR/refs_aligned.fasta" ] || \
        mafft --auto --thread "$THREADS" "$REFS_FASTA" \
            > "$MSA_DIR/refs_aligned.fasta" 2> "$MSA_DIR/mafft_refs.log"

    mafft --addfragments "$clustered" --auto --thread "$THREADS" \
        "$MSA_DIR/refs_aligned.fasta" > "$MSA_DIR/refs_plus_contigs.fasta" \
        2> "$MSA_DIR/mafft_add.log"
    tail -3 "$MSA_DIR/mafft_add.log" || true
    report "$MSA_DIR/refs_plus_contigs.fasta"
    ;;

# Step 3: align the complete-genome references, then add the cluster members onto them
complete)
    need mafft
    mkdir -p "$COMPLETE_DIR"
    mafft --auto --thread "$THREADS" "$COMPLETE_DIR/length_split/refs_moderate.fa" \
        > "$COMPLETE_DIR/refs_moderate_only_msa.fasta" 2> "$COMPLETE_DIR/mafft_refs.log"
    report "$COMPLETE_DIR/refs_moderate_only_msa.fasta"

    mafft --add "$COMPLETE_DIR/allmembers_c0-27.fa" --auto --adjustdirectionaccurately \
        --thread "$THREADS" "$COMPLETE_DIR/refs_moderate_only_msa.fasta" \
        > "$COMPLETE_DIR/allmembers_with_refs_msa.fasta" \
        2> "$COMPLETE_DIR/mafft_allmembers.log"
    tail -2 "$COMPLETE_DIR/mafft_allmembers.log" || true
    report "$COMPLETE_DIR/allmembers_with_refs_msa.fasta"
    ;;

# Step 4: align everything from scratch, references + members + BLAST downloads,
# which is the alignment behind the v2 and v3 trees
scratch)
    need mafft
    mkdir -p "$TREE_DIR"
    cat "$COMPLETE_DIR/length_split/refs_moderate.fa" \
        "$COMPLETE_DIR/allmembers_c0-27.fa" \
        "$TREE_DIR/blast_dedup.fa" > "$TREE_DIR/all_unaligned.fa"
    echo "  $(grep -c '^>' "$TREE_DIR/all_unaligned.fa") sequences to align"

    mafft --auto --adjustdirectionaccurately --thread "$THREADS" \
        "$TREE_DIR/all_unaligned.fa" > "$TREE_DIR/scratch_msa.fasta" \
        2> "$TREE_DIR/mafft_scratch.log"
    tail -3 "$TREE_DIR/mafft_scratch.log" || true
    report "$TREE_DIR/scratch_msa.fasta"
    ;;

# Step 5: add the downloads onto the existing alignment, --keeplength holding the
# reference columns fixed so the two are comparable
add)
    need mafft
    mkdir -p "$TREE_DIR"
    mafft --add "$TREE_DIR/blast_dedup.fa" --keeplength --adjustdirectionaccurately \
        --thread "$THREADS" "$COMPLETE_DIR/allmembers_with_refs_msa.fasta" \
        > "$TREE_DIR/add_msa.fasta" 2> "$TREE_DIR/mafft_add.log"
    tail -3 "$TREE_DIR/mafft_add.log" || true
    report "$TREE_DIR/add_msa.fasta"
    ;;

# Step 6: try every cutoff + report what each one leaves, which is how -gt 0.8 was chosen
compare)
    need trimal
    src="${SRC:-$TREE_DIR/scratch_msa.fasta}"
    outdir="${OUTDIR:-$(dirname "$src")/trim}"
    mkdir -p "$outdir"
    report "$src"

    for method in "${TRIM_METHODS[@]}"; do
        name=$(echo "$method" | tr -d ' .-')
        out="$outdir/trim_${name}.fasta"
        if ! trimal -in "$src" -out "$out" $method 2> "$outdir/trim_${name}.log"; then
            echo "  $method failed, see $outdir/trim_${name}.log"
            continue
        fi
        echo -n "  $method ->"
        report "$out"
    done

    echo "pick one, then either:"
    echo "  FastTree -nt -gtr <chosen fasta> > $outdir/tree.nwk 2> $outdir/fasttree.log"
    echo "  iqtree2 -s <chosen fasta> -m MFP -T AUTO --prefix $outdir/tree"
    ;;

# Step 7: drop the long-branch references + trim at the chosen cutoff, then IQ-TREE
trim)
    need trimal
    src="${SRC:-$TREE_DIR/scratch_msa.fasta}"
    out="${src%.fasta}_gt08.fasta"

    if [ -n "$DROP" ]; then
        python3 scripts/alignment/fasta_tools.py drop --alignment "$src" \
            --out "${src%.fasta}_dropped.fasta" --drop $DROP
        src="${src%.fasta}_dropped.fasta"
    fi

    trimal -in "$src" -out "$out" -gt "$TRIM_OCCUPANCY" -fasta
    report "$out"
    echo "then: iqtree2 -s $out -m MFP -T AUTO --prefix ${out%.fasta}"
    echo "      add -B 1000 -alrt 1000 for support values"
    ;;

# Step 8: build the China-only alignment, references kept + other countries dropped
china)
    need trimal
    src="${SRC:-$TREE_DIR/scratch_msa.fasta}"
    away=$(awk -F'\t' 'NR>1 && $6 != "China" {print $1}' "$ANNOTATION" | tr '\n' ' ')
    echo "  dropping $(echo "$away" | wc -w) contigs from outside China"

    python3 scripts/alignment/fasta_tools.py drop --alignment "$src" \
        --out "$TREE_DIR/china_only.fasta" --drop $away $DROP
    report "$TREE_DIR/china_only.fasta"
    echo "then: iqtree2 -s $TREE_DIR/china_only.fasta -T AUTO --prefix $TREE_DIR/china_only"
    ;;

*)
    echo "unknown stage: $STAGE" >&2
    exit 1
    ;;
esac

echo "stage $STAGE done"

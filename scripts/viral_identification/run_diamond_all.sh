
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


# STEP 1: Enable error checking
# -u: Stop for undefined variables.
# pipefail: Detect errors in pipelines.
set -uo pipefail


# STEP 2: Set file paths and CPU threads
# Find project root (two directories above this script).
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

# List of sequencing runs to process.
RUNLIST="$ROOT_DIR/results/viral_identification/diamond_available_runs.txt"

# Arbovirus protein reference database.
DB="$ROOT_DIR/data/reference/arbovirus_refseq"

# Use 8 threads unless THREADS is already defined.
THREADS="${THREADS:-8}"

# Directories for results, timing and logs.
OUTDIR="$ROOT_DIR/results/viral_identification/diamond"
TIMEDIR="$ROOT_DIR/results/viral_identification/timing"
LOGDIR="$ROOT_DIR/results/viral_identification/logs"


# STEP 3: Check required tools and input files
# Stop if DIAMOND, GNU time, run list or database is missing.
command -v diamond >/dev/null 2>&1 || { echo "ERROR: diamond not found" >&2; exit 1; }
[[ -x /usr/bin/time ]] || { echo "ERROR: GNU /usr/bin/time not found" >&2; exit 1; }
[[ -f "$RUNLIST" ]] || { echo "ERROR: missing run list: $RUNLIST" >&2; exit 1; }
[[ -f "$DB.dmnd" ]] || { echo "ERROR: DIAMOND database not found: $DB.dmnd" >&2; exit 1; }

# Create output directories if needed.
mkdir -p "$OUTDIR" "$TIMEDIR" "$LOGDIR"


# STEP 4: Initialise counters
# Track successful searches, zero hits, skipped runs and errors.
success=0; nohits=0; skipped=0; missing=0; failed=0


# STEP 5: Read and process each sequencing run
# Read accession IDs one by one, including a final line without a newline.
while IFS= read -r acc || [[ -n "$acc" ]]; do

    # Remove Windows carriage return characters.
    acc="${acc%$'\r'}"

    # Skip empty lines and comments.
    [[ -z "$acc" || "$acc" == \#* ]] && continue

    # Accept only valid SRR, ERR or DRR accession formats.
    if [[ ! "$acc" =~ ^(SRR|ERR|DRR)[0-9]+$ ]]; then
        echo "WARNING: invalid accession: $acc" >&2
        ((failed+=1)); continue
    fi


    # STEP 6: Define files for the current run
    # Query: assembled nucleotide contigs.
    query="$ROOT_DIR/data/raw/logan_contigs/${acc}.contigs.fa"

    # Output: DIAMOND BLASTx alignment results.
    out="$OUTDIR/${acc}.blastx.tsv"

    # Temporary output prevents incomplete results being saved.
    temporary="$OUTDIR/${acc}.blastx.tsv.tmp"

    # Marker indicating successful completion.
    done_marker="$OUTDIR/${acc}.blastx.done"

    # Files for recording runtime and error messages.
    timefile="$TIMEDIR/${acc}.time.txt"
    logfile="$LOGDIR/${acc}.diamond.log"

    echo "Processing $acc"


    # STEP 7: Check input and avoid repeated searches

    # Skip missing or empty contig FASTA files.
    if [[ ! -s "$query" ]]; then
        echo "WARNING: missing/empty input FASTA: $acc" >&2
        ((missing+=1)); continue
    fi

    # Skip runs already completed successfully.
    if [[ -f "$done_marker" && -f "$out" ]]; then
        echo "  Already completed (status marker); skipping"
        ((skipped+=1)); continue
    fi

    # Skip existing nonempty legacy results without completion markers.
    if [[ -s "$out" && ! -e "$done_marker" ]]; then
        echo "  Legacy nonempty output present; skipping (check old log for completion)"
        ((skipped+=1)); continue
    fi


    # STEP 8: Run DIAMOND BLASTx
    # Remove previous temporary and timing files.
    rm -f "$temporary" "$timefile"

    # Record runtime and search translated nucleotide contigs
    # against the arbovirus protein database.
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

        # Main DIAMOND options:
        # --sensitive: Improve detection of distant protein matches.
        # --evalue 1e-5: Set maximum reported E-value.
        # --max-target-seqs 5: Report up to 5 target sequences per query.
        # --header simple: Include a simple column header.
        # --outfmt 6: Tab-separated results with 17 fields.
        #
        # Fields include sequence IDs, identity, coverage,
        # alignment positions, E-value, taxonomy and protein title.


        # STEP 9: Save successful results

        # Ensure an output file exists, even if no hits were found.
        [[ -f "$temporary" ]] || : > "$temporary"

        # Rename temporary output to final result.
        mv -f "$temporary" "$out"

        # Create a success marker.
        printf 'success\n' > "$done_marker"

        # Count successful searches with or without matches.
        if [[ -s "$out" ]]; then
            echo "  SUCCESS: matches found"
            ((success+=1))
        else
            echo "  SUCCESS: zero hits"
            ((nohits+=1))
        fi

    else

        # STEP 10: Handle failed searches
        # Display warning and remove incomplete output.
        echo "WARNING: DIAMOND failed for $acc; see $logfile" >&2
        rm -f "$temporary" "$done_marker"
        ((failed+=1))
    fi

# Read sequencing run IDs from the run list.
done < "$RUNLIST"


# STEP 11: Print summary
# Report the number of successful, skipped, missing and failed runs.
echo "Summary: hits=$success, zero_hits=$nohits, skipped=$skipped, missing_queries=$missing, failed=$failed"

# Return an error status if any run failed or had missing input.
if (( failed > 0 || missing > 0 )); then
    exit 1
fi

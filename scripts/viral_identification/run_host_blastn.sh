
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


# ============================================================
# STEP 1 — Enable strict Bash error handling
# ============================================================

# -e: Exit immediately if an unhandled command fails.
# -u: Treat references to undefined variables as errors.
# -o pipefail: Make a pipeline fail if any command in it fails,
#              rather than checking only the final command.
#
# These settings reduce the chance that the pipeline continues
# silently after encountering an error.

set -euo pipefail


# ============================================================
# STEP 2 — Define the project root and input/output directories
# ============================================================

# BASH_SOURCE[0] refers to the path of this script.
# dirname obtains the directory containing the script.
# /../.. moves two directories upward to the project root.
# cd changes to that directory, and pwd returns its full path.
#
# This allows the script to locate files regardless of the
# directory from which the user executes it.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

# Define the prefix of the Culicidae nucleotide BLAST database.
#
# Culicidae is the mosquito family.
# The database contains mosquito-associated reference sequences
# against which candidate contigs will be aligned.
#
# BLAST uses a database prefix, so an individual extension
# such as .nin or .nsq is not specified here.

DB="$ROOT_DIR/data/reference/culicidae_blast/culicidae"

# Directory containing the candidate FASTA files.
# These files were generated from the three alternative
# DIAMOND filtering branches: strict, loose and evalue.

HOST_DIR="$ROOT_DIR/results/viral_identification/host_filter"

# Directory for storing BLASTN alignment results.

OUTDIR="$HOST_DIR/blastn"

# Directory for storing BLASTN execution logs.

LOGDIR="$HOST_DIR/logs"

# Use eight CPU threads to accelerate the BLASTN searches.

THREADS=8


# ============================================================
# STEP 3 — Check whether the required software is installed
# ============================================================

# Iterate through the two required NCBI BLAST+ commands:
#
# blastn:
#   Aligns nucleotide query sequences against a nucleotide
#   reference database.
#
# blastdbcmd:
#   Inspects an existing BLAST database and retrieves
#   information about its contents.
#
# command -v checks whether each command is available.
#
# >/dev/null suppresses standard output.
# 2>&1 also suppresses standard error.
#
# If a command is unavailable, print an error message to
# standard error and terminate the script with exit code 1.

for cmd in blastn blastdbcmd; do
    command -v "$cmd" >/dev/null 2>&1 || {
        echo "ERROR: $cmd not installed" >&2
        exit 1
    }
done


# ============================================================
# STEP 4 — Verify the Culicidae BLAST database
# ============================================================

# blastdbcmd -db specifies the database prefix.
# -dbtype nucl specifies that this is a nucleotide database.
# -info requests information about the database.
#
# The command succeeds only if the database can be opened
# and inspected.
#
# Its normal output and error messages are suppressed.
#
# If the database is missing or unreadable, print an error
# and stop before attempting any sequence alignments.

blastdbcmd -db "$DB" -dbtype nucl -info >/dev/null 2>&1 || {
    echo "ERROR: Culicidae nucleotide BLAST database missing or unreadable: $DB" >&2
    exit 1
}


# ============================================================
# STEP 5 — Create the required output directories
# ============================================================

# mkdir creates directories.
# -p also creates missing parent directories and does not
# produce an error if the directories already exist.
#
# One directory stores the alignment results, while the
# other stores execution logs for troubleshooting.

mkdir -p "$OUTDIR" "$LOGDIR"


# ============================================================
# STEP 6 — Process the three candidate-filtering branches
# ============================================================

# The previous DIAMOND filtering stage generated three
# alternative candidate sets:
#
# strict:
#   Protein identity >= 70%, and query or subject coverage >= 80%.
#
# loose:
#   Protein identity >= 30%, and query or subject coverage >= 50%.
#
# evalue:
#   DIAMOND E-value < 1e-5.
#
# These branches are processed independently, not sequentially.
#
# The loop runs the same host-reference search for each branch.

for mode in strict loose evalue; do

    # --------------------------------------------------------
    # STEP 6.1 — Define mode-specific input and output files
    # --------------------------------------------------------

    # Input FASTA containing candidate contig nucleotide
    # sequences from the current filtering branch.
    #
    # Example for mode=strict:
    # strict_candidate_contigs.fa

    QUERY="$HOST_DIR/${mode}_candidate_contigs.fa"

    # Final tab-separated BLASTN alignment output.
    #
    # Example:
    # blastn/strict_vs_culicidae.tsv

    OUT="$OUTDIR/${mode}_vs_culicidae.tsv"

    # Temporary alignment file.
    #
    # BLASTN writes its results here first.
    # Only after successful completion is this file renamed
    # to the final output filename.
    #
    # This reduces the risk of treating an interrupted or
    # incomplete BLASTN search as a finished result.

    TMP="$OUT.tmp"

    # Log file containing BLASTN messages and possible errors.
    #
    # Example:
    # logs/strict_blastn.log

    LOG="$LOGDIR/${mode}_blastn.log"

    # Show which filtering branch is being processed.

    echo "Host BLASTN: $mode"


    # --------------------------------------------------------
    # STEP 6.2 — Check that the candidate FASTA file exists
    # --------------------------------------------------------

    # -f tests whether the path is an existing regular file.
    # ! negates the condition.
    #
    # If the expected FASTA is missing, the pipeline cannot
    # determine which contigs should be searched.
    #
    # Therefore, terminate rather than treating a missing
    # input as an empty candidate set.

    if [[ ! -f "$QUERY" ]]; then
        echo "ERROR: Missing extracted candidate FASTA: $QUERY" >&2
        exit 1
    fi


    # --------------------------------------------------------
    # STEP 6.3 — Remove any old temporary output
    # --------------------------------------------------------

    # -f forces removal without prompting and does not
    # report an error if the temporary file is absent.
    #
    # This prevents results from an earlier interrupted
    # attempt from being confused with the current search.

    rm -f "$TMP"


    # --------------------------------------------------------
    # STEP 6.4 — Handle an empty candidate FASTA file
    # --------------------------------------------------------

    # -s checks whether a file exists and has a size > 0 bytes.
    #
    # Because the file's existence was verified above,
    # ! -s means that the candidate FASTA is empty.
    #
    # An empty FASTA indicates that this filtering branch
    # contains no candidate nucleotide sequences to search.
    #
    # In this situation, BLASTN does not need to run.

    if [[ ! -s "$QUERY" ]]; then

        # Create or overwrite the result file as an empty file.
        #
        # The downstream removal script treats an existing
        # empty BLASTN file as indicating zero reported hits.

        : > "$OUT"

        # Record why no BLASTN search was performed.
        # %s is replaced by the current mode name.
        # \n inserts a newline.

        printf 'No candidate contigs for %s; created empty result.\n' "$mode" > "$LOG"

        # Skip the remaining commands in this iteration
        # and move to the next filtering branch.

        continue
    fi


    # ========================================================
    # STEP 7 — Search candidate contigs against Culicidae
    # ========================================================

    # BLASTN compares nucleotide sequences with nucleotide
    # sequences in the mosquito reference database.
    #
    # The purpose is to identify candidate viral contigs
    # that also show substantial nucleotide similarity
    # to mosquito-associated reference sequences.
    #
    # Such matches may indicate host-derived sequences or
    # other sequence similarities requiring further screening.
    #
    # BLASTN searches for local sequence alignments.
    # A reported match does not necessarily cover the full
    # length of the query contig.
    #
    # Therefore, the downstream Python script separately
    # evaluates query alignment coverage.


    if blastn \

        # -task blastn:
        # Use the standard BLASTN nucleotide alignment task.
        #
        # This is different from DIAMOND BLASTx, which
        # translates nucleotide queries into amino acids
        # before searching a protein database.

        -task blastn \

        # -query:
        # FASTA file containing candidate nucleotide contigs.

        -query "$QUERY" \

        # -db:
        # Culicidae nucleotide reference database prefix.

        -db "$DB" \

        # -num_threads:
        # Use eight CPU threads for the alignment search.

        -num_threads "$THREADS" \

        # -perc_identity 80:
        # Request alignments with at least approximately
        # 80% nucleotide sequence identity.
        #
        # Identity measures the proportion of identical
        # aligned nucleotide positions.
        #
        # This parameter does not require the alignment
        # to span 80% of the complete query contig.

        -perc_identity 80 \

        # -outfmt 6:
        # Produce a tab-delimited alignment table
        # without a header.
        #
        # The following twelve columns are requested:
        #
        # 1. qseqid   Query contig identifier.
        # 2. sseqid   Matching reference sequence identifier.
        # 3. pident   Percentage of identical aligned nucleotides.
        # 4. length   Alignment length, including gaps.
        # 5. qlen     Full length of the query contig.
        # 6. slen     Full length of the reference sequence.
        # 7. qstart   Start position of alignment on query.
        # 8. qend     End position of alignment on query.
        # 9. sstart   Start position of alignment on subject.
        # 10. send    End position of alignment on subject.
        # 11. evalue Statistical significance of the alignment.
        # 12. bitscore Normalised alignment score.
        #
        # These fields allow the downstream Python script
        # to calculate query alignment coverage:
        #
        # query_aligned_span = abs(qend - qstart) + 1
        #
        # query_coverage (%) =
        #     query_aligned_span / qlen * 100
        #
        # The downstream strong-match criterion is:
        #
        # nucleotide identity >= 80%
        # AND query alignment coverage >= 80%.

        -outfmt '6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore' \

        # -out:
        # Save the BLASTN alignment table to a temporary
        # file rather than directly to the final output.

        -out "$TMP" > "$LOG" 2>&1; then

        # ----------------------------------------------------
        # STEP 8 — Save the result after successful BLASTN
        # ----------------------------------------------------

        # mv renames the temporary output to the final
        # mode-specific output file.
        #
        # This happens only when BLASTN exits successfully.
        #
        # An empty output is still a valid result when
        # BLASTN succeeds but reports no alignments.

        mv "$TMP" "$OUT"

        # Print the location of the completed output.

        echo "SUCCESS: $mode -> $OUT"

    else

        # ----------------------------------------------------
        # STEP 9 — Handle a failed BLASTN search
        # ----------------------------------------------------

        # Remove the temporary and final output files.
        #
        # This prevents a previous or partially generated
        # result from being mistaken for a successful search.

        rm -f "$TMP" "$OUT"

        # Print the error and identify the log file
        # containing information for troubleshooting.

        echo "ERROR: Host BLASTN failed for $mode; see $LOG" >&2

        # Stop execution with a nonzero exit status.

        exit 1
    fi

# End of the current filtering branch.
# Continue until strict, loose and evalue have all been processed.

done


# ============================================================
# STEP 10 — Confirm that all three searches finished
# ============================================================

# The resulting TSV files can then be passed to
# remove_mosquito_contigs.py for host-reference filtering.

echo 'All three host BLASTN searches completed.'

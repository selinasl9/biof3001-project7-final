#!/usr/bin/env python3
"""
match_runs_to_logan.py

Find which mosquito RNA runs from our NCBI SRA list are also in Logan.

Inputs:
    SRR_Acc_List.csv                          our NCBI SRA run accessions
    logan_accessions_v1.2_SRA2025.csv.zst     Logan SRA run table, zstd-compressed

Outputs:
    logan_matches.tsv                         matched runs with Logan metadata
    logan_matches_by_bioproject.tsv           matched runs grouped by BioProject + link
    not_in_logan.txt                          NCBI runs Logan does not have
"""

import csv
import io
import os
import shutil
import subprocess
import sys
from collections import defaultdict

# Set file paths and the Logan columns we want to keep
# Paths are relative to the repo root; the Logan table sits outside the repo
NCBI_LIST = "data/runs/SRR_Acc_List.csv"
LOGAN_TABLE = os.path.expanduser("~/logan/logan_accessions_v1.2_SRA2025.csv.zst")
OUT_MATCHES = "results/run_selection/logan_matches.tsv"
OUT_BY_PROJECT = "results/run_selection/logan_matches_by_bioproject.tsv"
OUT_NOT_IN_LOGAN = "results/run_selection/not_in_logan.txt"
ACC_COLUMN = "acc"
WANTED_COLUMNS = [
    "acc",
    "bioproject",
    "biosample",
    "sample_acc",
    "experiment",
    "library_name",
    "organism",
    "assay_type",
    "librarysource",
    "libraryselection",
    "librarylayout",
    "platform",
    "instrument",
    "releasedate",
    "collection_date_sam",
    "geo_loc_name_country_calc",
    "geo_loc_name_country_continent_calc",
    "mbases",
    "mbytes",
    "avgspotlen",
    "jattr",
]

# Select columns for the summary in the BioProject table
BIOPROJECT_COLUMN = "bioproject"
SUMMARY_FIELDS = {
    "countries": "geo_loc_name_country_calc",
    "organisms": "organism",
    "librarysources": "librarysource",
    "assay_types": "assay_type",
    "platforms": "platform",
}

BIOPROJECT_URL = "https://www.ncbi.nlm.nih.gov/bioproject/{}"


def normalise(acc):
    """Clean accession by removing surrounding whitespace and stray quote marks."""
    return acc.strip().strip("\"'").strip()


def check_inputs():
    """Stop before doing any work if an input file or zstdcat is missing."""
    problems = []
    for path in (NCBI_LIST, LOGAN_TABLE):
        if not os.path.isfile(path):
            problems.append(f"input file not found: {path}")
    if shutil.which("zstdcat") is None:
        problems.append("zstdcat not found on PATH (install the zstd package)")
    if problems:
        for p in problems:
            print(f"ERROR: {p}", file=sys.stderr)
        sys.exit(1)


def load_ncbi_accessions():
    """Load the NCBI list into a set of cleaned accessions for fast look-ups.
    Duplicates become one entry."""
    wanted = set()
    with open(NCBI_LIST, encoding="utf-8-sig") as f:
        for line in f:
            acc = normalise(line)
            if acc:
                wanted.add(acc)
    return wanted


def remove_quietly(paths):
    """Delete leftover .tmp files after a failure."""
    for p in paths:
        try:
            os.remove(p)
        except FileNotFoundError:
            pass


def main():
    # Step 1: make sure all inputs are there
    check_inputs()

    # Step 2: load the NCBI accessions we are looking for
    wanted = load_ncbi_accessions()
    print(f"NCBI runs (unique): {len(wanted)}")

    # Step 3: stream the Logan table through zstdcat
    proc = subprocess.Popen(["zstdcat", LOGAN_TABLE], stdout=subprocess.PIPE)
    stream = io.TextIOWrapper(proc.stdout, encoding="utf-8", errors="replace", newline="")

    # Raise field size as some Logan fields may be long
    csv.field_size_limit(sys.maxsize)
    reader = csv.reader(stream)

    # Step 4: read the header and find where each column is
    header = next(reader, None)
    if header is None:
        proc.wait()
        print("ERROR: Logan table is empty or could not be read.", file=sys.stderr)
        sys.exit(1)
    header = [h.strip() for h in header]
    print("Logan header:")
    print("  " + ", ".join(header))

    if ACC_COLUMN not in header:
        proc.kill()
        proc.wait()
        print(f"ERROR: accession column '{ACC_COLUMN}' not in Logan header.", file=sys.stderr)
        sys.exit(1)

    missing = [c for c in WANTED_COLUMNS if c not in header]
    if missing:
        print("Wanted columns missing from Logan header: " + ", ".join(missing))
    else:
        print("All wanted columns present in Logan header.")

    # Keep the desired columns that exist, in order
    kept = [c for c in WANTED_COLUMNS if c in header]
    kept_idx = [header.index(c) for c in kept]
    acc_idx = header.index(ACC_COLUMN)

    # Step 5: scan every Logan row and keep only rows that match NCBI accession
    # Matched rows are written to a .tmp file + collect BioProject summaries as we go
    tmp_matches = OUT_MATCHES + ".tmp"
    tmp_by_project = OUT_BY_PROJECT + ".tmp"
    tmp_not_in = OUT_NOT_IN_LOGAN + ".tmp"
    tmp_files = [tmp_matches, tmp_by_project, tmp_not_in]

    found = set()
    project_runs = defaultdict(int)
    project_values = defaultdict(lambda: defaultdict(set))

    try:
        with open(tmp_matches, "w", encoding="utf-8", newline="") as out:
            writer = csv.writer(out, delimiter="\t", lineterminator="\n")
            writer.writerow(kept)

            for row in reader:
                if len(row) <= acc_idx:
                    continue  # short row, no accession to read
                acc = normalise(row[acc_idx])
                if acc not in wanted or acc in found:
                    continue  # not in our accession list or already recorded
                found.add(acc)

                values = [row[i] if i < len(row) else "" for i in kept_idx]
                writer.writerow(values)

                # Add this run to its BioProject's summary
                record = dict(zip(kept, values))
                project = record.get(BIOPROJECT_COLUMN, "").strip() or "(none)"
                project_runs[project] += 1
                for field, col in SUMMARY_FIELDS.items():
                    v = record.get(col, "").strip()
                    if v:
                        project_values[project][field].add(v)

        # Step 6: check zstdcat finished cleanly, delete .tmp files and stop
        # if stream cut short
        stream.close()
        rc = proc.wait()
        if rc != 0:
            remove_quietly(tmp_files)
            print(f"ERROR: zstdcat exited with code {rc}; output discarded.", file=sys.stderr)
            sys.exit(1)

        # Step 7: write BioProject summary with a link to each BioProject page
        with open(tmp_by_project, "w", encoding="utf-8", newline="") as out:
            writer = csv.writer(out, delimiter="\t", lineterminator="\n")
            writer.writerow(["bioproject", "n_runs", *SUMMARY_FIELDS, "link"])
            for project in sorted(project_runs, key=lambda p: (-project_runs[p], p)):
                link = BIOPROJECT_URL.format(project) if project != "(none)" else ""
                writer.writerow([
                    project,
                    project_runs[project],
                    *("; ".join(sorted(project_values[project][f])) for f in SUMMARY_FIELDS),
                    link,
                ])

        # Step 8: write the NCBI runs Logan does not have
        with open(tmp_not_in, "w", encoding="utf-8") as out:
            for acc in sorted(wanted - found):
                out.write(acc + "\n")

    except BaseException:
        # For any crash or Ctrl+C, don't leave half-written .tmp files behind
        if proc.poll() is None:
            proc.kill()
            proc.wait()
        remove_quietly(tmp_files)
        raise

    # Step 9: swap the .tmp files into place
    os.replace(tmp_matches, OUT_MATCHES)
    os.replace(tmp_by_project, OUT_BY_PROJECT)
    os.replace(tmp_not_in, OUT_NOT_IN_LOGAN)

    # Step 10: return summary in the terminal
    pct = (100 * len(found) / len(wanted)) if wanted else 0.0
    print(f"Total NCBI runs: {len(wanted)}")
    print(f"Matched in Logan: {len(found)}")
    print(f"Not in Logan:     {len(wanted) - len(found)}")
    print(f"% matched:        {pct:.1f}%")
    print(f"BioProjects with matches: {len(project_runs)}")
    print(f"Wrote {OUT_MATCHES}, {OUT_BY_PROJECT}, {OUT_NOT_IN_LOGAN}")


if __name__ == "__main__":
    main()
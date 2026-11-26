#!/usr/bin/env python3
"""
filter_runs_to_rna.py

Cut data/runs/logan_matches.tsv down to runs that are usable for virus screening,
and count how many runs each rule removes.

Rules (applied in this order):
  1. RNA or unclear   - the library must come from RNA (or be unclear, kept
                        for manual BioProject screening)
  2. Short-read       - platform must be a short-read platform, because
                        Logan's assembler is built for short reads
  3. Country recorded - geo_loc_name_country_calc must hold a real country,
                        because both research approaches need country info

Outputs
  data/runs/rna_runs.tsv               kept runs with rna_status column
  data/runs/bioprojects_to_screen.tsv  kept runs grouped by BioProject, with
                                       empty keep/reason columns for manual
                                       screening

"""

import pandas as pd

# File paths
IN_MATCHES = "data/runs/logan_matches.tsv"
OUT_RUNS = "data/runs/rna_runs.tsv"
OUT_SCREEN = "data/runs/bioprojects_to_screen.tsv"

BIOPROJECT_URL = "https://www.ncbi.nlm.nih.gov/bioproject/{}"

# Rule 1: RNA or unclear
#   - librarysource METATRANSCRIPTOMIC                        -> "RNA"
#   - librarysource METAGENOMIC and assay_type RNA-Seq        -> "RNA"
#   - librarysource METAGENOMIC and an RNA libraryselection   -> "RNA"
#   - librarysource METAGENOMIC and assay_type OTHER          -> "unclear"
#   - anything else                                           -> removed

RNA_LIBRARY_SOURCE = "METATRANSCRIPTOMIC"
METAGENOMIC_SOURCE = "METAGENOMIC"
RNA_ASSAY_TYPE = "RNA-Seq"
UNCLEAR_ASSAY_TYPE = "OTHER"
RNA_SELECTIONS = {"cDNA", "cDNA_randomPriming", "cDNA_oligo_dT", "RT-PCR",
                  "PolyA", "Inverse rRNA"}

# Rule 2: short-read platforms only
SHORT_READ_PLATFORMS = {"ILLUMINA", "BGISEQ", "DNBSEQ"}

# Rule 3: country must be recorded
COUNTRY_COLUMN = "geo_loc_name_country_calc"
COUNTRY_PLACEHOLDERS = {"uncalculated", "missing", "not collected", "not applicable"}

def rna_status(row):
    """Return "RNA", "unclear" or "" (= not RNA, remove) for one run."""
    source = row["librarysource"]
    if source == RNA_LIBRARY_SOURCE:
        return "RNA"
    if source == METAGENOMIC_SOURCE:
        if row["assay_type"] == RNA_ASSAY_TYPE:
            return "RNA"
        if row["libraryselection"] in RNA_SELECTIONS:
            return "RNA"
        if row["assay_type"] == UNCLEAR_ASSAY_TYPE:
            return "unclear"
    return ""

def join_unique(values):
    """Turn a column of values into one sorted "a; b; c" string, skipping blanks."""
    return "; ".join(sorted({v for v in values if v}))

def main():
    # Step 1: read the Logan matches table, keep accessions/IDs as text
    # Every empty cell is a blank string
    df = pd.read_csv(IN_MATCHES, sep="\t", dtype=str).fillna("")

    # Strip stray whitespace
    df = df.apply(lambda col: col.str.strip())

    print(f"Start (runs in {IN_MATCHES}): {len(df)}")

    # Step 2: Rule 1 - label each run RNA / unclear / not RNA, keep only RNA and unclear
    df["rna_status"] = df.apply(rna_status, axis=1)
    before = len(df)
    df = df[df["rna_status"] != ""]
    n_unclear = (df["rna_status"] == "unclear").sum()
    print(f"After rule 1 (RNA or unclear):      {len(df)}  (removed {before - len(df)}; "
          f"{len(df) - n_unclear} RNA, {n_unclear} unclear)")

    # Step 3: Rule 2 - keep only short-read platforms
    before = len(df)
    df = df[df["platform"].isin(SHORT_READ_PLATFORMS)]
    print(f"After rule 2 (short-read platform): {len(df)}  (removed {before - len(df)})")

    # Step 4: Rule 3 - turn placeholder country values into blanks, keep
    # only runs with a non-blank country
    df = df.copy()
    is_placeholder = df[COUNTRY_COLUMN].str.lower().isin(COUNTRY_PLACEHOLDERS)
    df.loc[is_placeholder, COUNTRY_COLUMN] = ""
    before = len(df)
    df = df[df[COUNTRY_COLUMN] != ""]
    print(f"After rule 3 (country recorded):    {len(df)}  (removed {before - len(df)})")

    # Step 5: write the kept runs (all original columns + rna_status)
    df.to_csv(OUT_RUNS, sep="\t", index=False)

    # Step 6: group kept runs by BioProject, one row per project
    screen = (
        df.groupby("bioproject")
        .agg(
            n_runs=("acc", "size"),
            rna_status=("rna_status", join_unique),
            countries=(COUNTRY_COLUMN, join_unique),
            organisms=("organism", join_unique),
        )
        .reset_index()
        .sort_values(["n_runs", "bioproject"], ascending=[False, True])
    )
    # Add a link to each project page and empty keep/reason column
    screen["link"] = screen["bioproject"].map(
        lambda p: BIOPROJECT_URL.format(p) if p else ""
    )
    screen["keep"] = ""
    screen["reason"] = ""
    screen.to_csv(OUT_SCREEN, sep="\t", index=False)

    # Step 7: return summary in the terminal
    print()
    print(f"Runs kept:        {len(df)}")
    print(f"BioProjects kept: {df['bioproject'].nunique()}")
    print(f"Countries kept:   {df[COUNTRY_COLUMN].nunique()}")
    print()
    print("Runs per country:")
    for country, n in df[COUNTRY_COLUMN].value_counts().items():
        print(f"  {country}\t{n}")
    print()
    print(f"Wrote {OUT_RUNS} and {OUT_SCREEN}")


if __name__ == "__main__":
    main()
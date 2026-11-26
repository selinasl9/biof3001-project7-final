#!/usr/bin/env python3
"""BIOF3001 Project 7 — Mosquito RNA Virome.

Script: make_full_summary.py

Objective:
    Aggregate per-run DIAMOND BLASTx results, select one best-scoring
    reference alignment per contig, attach NCBI taxonomy and run metadata,
    and generate run-level taxonomic count/detection summary tables.

Inputs (relative to the repository root):
    data/runs/all_project_runs_metadata.tsv
    results/viral_identification/diamond_available_runs.txt
    results/viral_identification/diamond/<accession>.blastx.tsv
    data/reference/taxonomy/nodes.dmp and data/reference/taxonomy/names.dmp

Outputs (in results/summary/):
    all_best_hits_with_metadata.tsv
    sample_metadata.tsv
    sample_virus_contig_counts.tsv
    sample_virus_presence_absence.tsv
    sample_family_presence_absence.tsv
    sample_genus_presence_absence.tsv
    zero_hit_samples.tsv
    missing_diamond_samples.tsv

Usage: python scripts/viral_identification/make_full_summary.py
Dependencies: Python 3, pandas; NCBI taxonomy dump.

Interpretation:
    - Scientific-name columns come from DIAMOND reference annotations,
      NOT CD-HIT/VSEARCH sequence clusters or confirmed virus species.
    - A detection value of 1 means a contig received that annotation;
      0 means no contig did, not confirmed biological absence.
    - Empty output files are interpreted as zero-hit searches only when
      successful completion is documented; new DIAMOND runs create a
      .blastx.done marker. Check legacy empty files against old logs.
    - Missing DIAMOND result files are excluded from detection matrices.
"""

import pandas as pd
from pathlib import Path

# ============================================================
# STEP 1: FILES
# Set paths and check required files.
# ============================================================

# Resolve paths from the script location, not the current working directory.
ROOT = Path(__file__).resolve().parents[2]
DIAMOND_DIR = ROOT / "results/viral_identification/diamond"
OUTDIR = ROOT / "results/summary"

METADATA_FILE = ROOT / "data/runs/all_project_runs_metadata.tsv"
RUNLIST_FILE = ROOT / "results/viral_identification/diamond_available_runs.txt"
NODES_FILE = ROOT / "data/reference/taxonomy/nodes.dmp"
NAMES_FILE = ROOT / "data/reference/taxonomy/names.dmp"

# Validate inputs before creating results (avoids misleading partial outputs).
for required in (METADATA_FILE, RUNLIST_FILE, NODES_FILE, NAMES_FILE):
    if not required.is_file():
        raise FileNotFoundError(f"Required input not found: {required}")

OUTDIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# STEP 2: DIAMOND COLUMNS
# Define the 17 DIAMOND output fields.
# ============================================================

diamond_cols = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "qlen",
    "slen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qcovhsp",
    "scovhsp",
    "staxids",
    "sscinames",
    "stitle"
]

# ============================================================
# STEP 3: READ SAMPLE METADATA
# Load metadata and select listed runs.
# ============================================================

metadata = pd.read_csv(
    METADATA_FILE,
    sep="\t"
)

metadata = metadata.rename(
    columns={
        "acc": "sample",
        "geo_loc_name_country_calc": "country",
        "organism": "host_organism"
    }
)

# Restrict to runs in the DIAMOND run list; the number is not hard-coded.
with open(RUNLIST_FILE) as f:
    analysed_samples = list(dict.fromkeys(
        line.strip() for line in f
        if line.strip() and not line.lstrip().startswith("#")
    ))

metadata = metadata[
    metadata["sample"].isin(analysed_samples)
].copy()

metadata = metadata.drop_duplicates("sample")

print("Analysed samples in metadata:", len(metadata))

# ============================================================
# STEP 4: READ NCBI TAXONOMY
# Load taxon names, parent IDs and ranks.
# ============================================================

print("Loading NCBI taxonomy...")

taxid_to_name = {}

with open(NAMES_FILE) as f:
    for line in f:
        parts = [x.strip() for x in line.split("|")]

        if len(parts) < 4:
            continue

        taxid = parts[0]
        name = parts[1]
        name_class = parts[3]

        if name_class == "scientific name":
            taxid_to_name[taxid] = name


parent = {}
rank = {}

with open(NODES_FILE) as f:
    for line in f:
        parts = [x.strip() for x in line.split("|")]

        if len(parts) < 3:
            continue

        taxid = parts[0]
        parent_taxid = parts[1]
        tax_rank = parts[2]

        parent[taxid] = parent_taxid
        rank[taxid] = tax_rank


# ============================================================
# STEP 5: TAXONOMY FUNCTION
# Trace each taxon to family, genus and species.
# ============================================================

def get_lineage(taxid_value):

    if pd.isna(taxid_value):
        return "", "", "", ""

    # DIAMOND may occasionally provide multiple taxids
    candidate_taxids = str(taxid_value).split(";")

    taxid = None

    for candidate in candidate_taxids:
        candidate = candidate.strip()

        if candidate in parent:
            taxid = candidate
            break

    if taxid is None:
        return "", "", "", ""

    lineage = []
    ranked_names = {}

    seen = set()

    while taxid in parent and taxid not in seen:

        seen.add(taxid)

        name = taxid_to_name.get(taxid, taxid)
        tax_rank = rank.get(taxid, "no rank")

        lineage.append(name)

        if tax_rank != "no rank":
            ranked_names[tax_rank] = name

        if taxid == "1":
            break

        next_taxid = parent[taxid]

        if next_taxid == taxid:
            break

        taxid = next_taxid

    lineage.reverse()

    family = ranked_names.get("family", "")
    genus = ranked_names.get("genus", "")
    species = ranked_names.get("species", "")
    full_lineage = "; ".join(lineage)

    return family, genus, species, full_lineage


# ============================================================
# STEP 6: READ ALL DIAMOND OUTPUTS
# Read each run and select its best hits.
# ============================================================

all_best_hits = []

zero_hit_samples = []
missing_samples = []

print("Reading DIAMOND outputs...")

for i, sample in enumerate(analysed_samples, start=1):

    file = DIAMOND_DIR / f"{sample}.blastx.tsv"

    # Missing result: track separately, not as a negative detection.
    if not file.exists():
        print("WARNING - DIAMOND output missing:", sample)
        missing_samples.append(sample)
        continue

    # Zero-byte outputs can be valid searches with no hits. For legacy
    # runs, verify successful completion from the upstream DIAMOND logs.
    # Empty result: record zero hits; warn if no success marker.
    if file.stat().st_size == 0:
        marker = DIAMOND_DIR / f"{sample}.blastx.done"
        if not marker.exists():
            print("WARNING - zero-byte result without success marker; verify log:", sample)
        zero_hit_samples.append(sample)
        continue

    # Read the tab-separated table; handle alternate header formats.
    try:
        df = pd.read_csv(
            file,
            sep="\t"
        )

        df.columns = [
            str(x).lstrip("#")
            for x in df.columns
        ]

        # Fallback if header not correctly detected
        if "qseqid" not in df.columns:
            df = pd.read_csv(
                file,
                sep="\t",
                comment="#",
                header=None,
                names=diamond_cols
            )

    except pd.errors.EmptyDataError:
        zero_hit_samples.append(sample)
        continue

    if len(df) == 0:
        zero_hit_samples.append(sample)
        continue

    # Fail clearly if DIAMOND outputs use unexpected column names.
    required_cols = {"qseqid", "bitscore", "evalue", "qcovhsp", "staxids", "sscinames"}
    absent = required_cols.difference(df.columns)
    if absent:
        raise ValueError(f"{file}: missing expected DIAMOND columns: {sorted(absent)}")

    # Numeric columns
    numeric_cols = [
        "pident",
        "length",
        "qlen",
        "slen",
        "evalue",
        "bitscore",
        "qcovhsp",
        "scovhsp"
    ]

    for col in numeric_cols:
        df[col] = pd.to_numeric(
            df[col],
            errors="coerce"
        )

    # --------------------------------------------------------
    # Keep ONE best hit per contig
    #
    # priority:
    # 1. highest bitscore
    # 2. lowest E-value
    # 3. highest query coverage
    # --------------------------------------------------------

    # Rank hits by bitscore, then E-value, then query coverage.
    df = df.sort_values(
        [
            "qseqid",
            "bitscore",
            "evalue",
            "qcovhsp"
        ],
        ascending=[
            True,
            False,
            True,
            False
        ]
    )

    # Retain only the first (best) reference hit per contig.
    best = (
        df
        .drop_duplicates("qseqid")
        .copy()
    )

    best["sample"] = sample

    all_best_hits.append(best)

    if i % 50 == 0:
        print(f"Processed {i} / {len(analysed_samples)} samples")


# ============================================================
# STEP 7: COMBINE BEST HITS
# Combine the runs into one table.
# ============================================================

if not all_best_hits:
    raise RuntimeError("No DIAMOND hits found in listed runs; cannot build best-hit tables.")

best_all = pd.concat(
    all_best_hits,
    ignore_index=True
)

# Only files actually present (including valid zero-hit files) can be
# represented as rows in the output matrices. Missing is not absence.
missing_set = set(missing_samples)
valid_samples = [s for s in analysed_samples if s not in missing_set]

print()
print("Samples with zero hits:", len(zero_hit_samples))
print("Best-hit rows:", len(best_all))


# ============================================================
# STEP 8: ADD TAXONOMY
# Attach taxonomy to each best hit.
# ============================================================

print("Adding taxonomy...")

# Look up family, genus, species and full lineage for each taxid.
taxonomy_results = best_all["staxids"].apply(
    get_lineage
)

best_all[
    [
        "family",
        "genus",
        "species",
        "full_lineage"
    ]
] = pd.DataFrame(
    taxonomy_results.tolist(),
    index=best_all.index
)


# ============================================================
# STEP 9: ADD COUNTRY / PROJECT / HOST / ASSEMBLY METADATA
# Attach run metadata to best hits.
# ============================================================

metadata_cols = [
    "sample",
    "bioproject",
    "country",
    "host_organism",
    "biosample",
    "collection_date_sam",
    "seqstats_contigs_n50",
    "seqstats_contigs_sumlen",
    "seqstats_contigs_nbseq",
    "seqstats_contigs_maxlen"
]

metadata_small = metadata[
    [
        c
        for c in metadata_cols
        if c in metadata.columns
    ]
].copy()

best_all = best_all.merge(
    metadata_small,
    on="sample",
    how="left"
)


# ============================================================
# STEP 10: REORDER COLUMNS
# Organise the output columns.
# ============================================================

preferred_order = [
    "sample",
    "bioproject",
    "country",
    "host_organism",
    "biosample",
    "collection_date_sam",

    "qseqid",
    "sseqid",

    "sscinames",
    "family",
    "genus",
    "species",
    "full_lineage",

    "pident",
    "length",
    "qlen",
    "slen",
    "qcovhsp",
    "scovhsp",
    "evalue",
    "bitscore",

    "staxids",
    "stitle",

    "seqstats_contigs_n50",
    "seqstats_contigs_sumlen",
    "seqstats_contigs_nbseq",
    "seqstats_contigs_maxlen",

    "qstart",
    "qend",
    "sstart",
    "send"
]

preferred_order = [
    c for c in preferred_order
    if c in best_all.columns
]

remaining = [
    c for c in best_all.columns
    if c not in preferred_order
]

best_all = best_all[
    preferred_order + remaining
]


# ============================================================
# STEP 11: SAVE DETAILED HIT TABLE
# Export the annotated best-hit table.
# ============================================================

best_all.to_csv(
    OUTDIR / "all_best_hits_with_metadata.tsv",
    sep="\t",
    index=False
)


# ============================================================
# STEP 12: SAVE SAMPLE METADATA TABLE
# Export metadata for runs with results.
# ============================================================

sample_metadata = metadata_small.copy()

sample_metadata = sample_metadata[
    sample_metadata["sample"].isin(valid_samples)
]

sample_metadata = sample_metadata.sort_values("sample")

sample_metadata.to_csv(
    OUTDIR / "sample_metadata.tsv",
    sep="\t",
    index=False
)


# ============================================================
# STEP 13: SAMPLE × VIRUS CONTIG COUNT MATRIX
# Count contigs per run and reference name.
# ============================================================

# Cross-tabulate contigs by run and DIAMOND reference scientific name.
counts = pd.crosstab(
    best_all["sample"],
    best_all["sscinames"]
)

# Critical:
# add samples with zero viral hits back into matrix
counts = counts.reindex(
    valid_samples,
    fill_value=0
)

counts.index.name = "sample"

counts.to_csv(
    OUTDIR / "sample_virus_contig_counts.tsv",
    sep="\t"
)


# ============================================================
# STEP 14: SAMPLE × VIRUS PRESENCE / ABSENCE
# Convert positive counts to 1 and zero counts to 0.
# ============================================================

# 1 = at least one annotated contig; 0 = none reported.
presence = (counts > 0).astype(int)

presence.to_csv(
    OUTDIR / "sample_virus_presence_absence.tsv",
    sep="\t"
)


# ============================================================
# STEP 15: FAMILY-LEVEL PRESENCE / ABSENCE
# Create a family detection matrix.
# ============================================================

family_hits = best_all[
    best_all["family"].notna()
    & (best_all["family"] != "")
].copy()

family_counts = pd.crosstab(
    family_hits["sample"],
    family_hits["family"]
)

family_counts = family_counts.reindex(
    valid_samples,
    fill_value=0
)

family_presence = (
    family_counts > 0
).astype(int)

family_presence.to_csv(
    OUTDIR / "sample_family_presence_absence.tsv",
    sep="\t"
)


# ============================================================
# STEP 16: GENUS-LEVEL PRESENCE / ABSENCE
# Create a genus detection matrix.
# ============================================================

genus_hits = best_all[
    best_all["genus"].notna()
    & (best_all["genus"] != "")
].copy()

genus_counts = pd.crosstab(
    genus_hits["sample"],
    genus_hits["genus"]
)

genus_counts = genus_counts.reindex(
    valid_samples,
    fill_value=0
)

genus_presence = (
    genus_counts > 0
).astype(int)

genus_presence.to_csv(
    OUTDIR / "sample_genus_presence_absence.tsv",
    sep="\t"
)


# ============================================================
# STEP 17: SAVE ZERO-HIT SAMPLES
# Save runs without reported hits.
# ============================================================

pd.DataFrame(
    {"sample": zero_hit_samples}
).to_csv(
    OUTDIR / "zero_hit_samples.tsv",
    sep="\t",
    index=False
)


# Track missing results separately so they are not mistaken for non-detections.
pd.DataFrame({"sample": missing_samples}).to_csv(
    OUTDIR / "missing_diamond_samples.tsv", sep="\t", index=False
)

# ============================================================
# STEP 18: FINAL SUMMARY
# Print counts and output file names.
# ============================================================

print()
print("=" * 60)
print("SUMMARY COMPLETE")
print("=" * 60)

print("DIAMOND runs listed:", len(analysed_samples))
print("DIAMOND files present:", len(valid_samples))
print("Missing DIAMOND files:", len(missing_samples))
print("Samples with zero hits:", len(zero_hit_samples))
print("Best-hit contigs:", len(best_all))
print("Unique scientific names:", best_all["sscinames"].nunique())
print("Unique families:", best_all["family"].replace("", pd.NA).nunique())
print("Unique genera:", best_all["genus"].replace("", pd.NA).nunique())

print()
print("Files created:")
print("results/summary/all_best_hits_with_metadata.tsv")
print("results/summary/sample_metadata.tsv")
print("results/summary/sample_virus_contig_counts.tsv")
print("results/summary/sample_virus_presence_absence.tsv")
print("results/summary/sample_family_presence_absence.tsv")
print("results/summary/sample_genus_presence_absence.tsv")
print("results/summary/zero_hit_samples.tsv")
print("results/summary/missing_diamond_samples.tsv")

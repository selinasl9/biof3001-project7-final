#!/usr/bin/env python3
"""BIOF3001 Project 7 — Putative reference annotation of viral clusters.

Objective: annotate 90%/70% sequence-based clusters with existing DIAMOND
reference species names, without performing additional alignments.
Inputs: original/China strict host-cleaned hit tables, cluster-membership
TSVs from build_matrices.py, and filtered matrix headers.
Outputs: CSVs for all/no-singleton/min5/min6 clusters and methods note.
Run: python scripts/viral_clustering/annotate_viral_clusters.py
Dependencies: Python 3, pandas.
Important: these are putative reference matches, NOT confirmed viral species.
Responsible contributor: [confirm before submission]
"""
"""BIOF3001 Project 7: annotate sequence clusters with EXISTING DIAMOND hits.

Outputs all/no-singletons/Min5/Min6 CSVs for 90% and 70% clustering.
No new alignments are performed. Virus labels are PUTATIVE, not confirmed species.
Run from repository root: python scripts/12_annotate_viral_clusters.py
"""
from pathlib import Path
import pandas as pd

ROOT = Path("results")
MAT = ROOT / "viral_clustering/matrices"
OUT = ROOT / "statistical_analysis/run_based_analysis/cluster_reference_virus_annotation"
HIT_FILES = [
    ROOT / "viral_identification/host_filter/strict_host_cleaned_hits.tsv",
    ROOT / "viral_identification/PRJNA1003643/host_filter/CN_PRJNA1003643_strict_host_cleaned_hits.tsv",
]
MIN_ALIGNED_AA = 50      # preliminary screen, NOT species identification
MAX_EVALUE = 1e-5
MIN_IDENT = 30.0
MIN_ANNOTATED_FRAC = 0.50


def load_hits():
    tables = []
    required = {"qseqid", "species", "length", "pident", "evalue", "bitscore"}
    for path in HIT_FILES:
        if not path.exists():
            raise FileNotFoundError(path)
        x = pd.read_csv(path, sep="\t", low_memory=False)
        if not required <= set(x.columns):
            raise ValueError(f"Missing annotation columns in {path}: {required - set(x.columns)}")
        tables.append(x[list(sorted(required))].copy())
    x = pd.concat(tables, ignore_index=True)
    for col in ("qseqid", "species"):
        x[col] = x[col].fillna("").astype(str).str.strip()
    for col in ("length", "pident", "evalue", "bitscore"):
        x[col] = pd.to_numeric(x[col], errors="coerce")
    x["qualified"] = (
        x["length"].ge(MIN_ALIGNED_AA)
        & x["pident"].ge(MIN_IDENT)
        & x["evalue"].le(MAX_EVALUE)
        & ~x["species"].str.casefold().isin(
            {"", "na", "nan", "none", "unknown", "unclassified"}
        )
    )
    if x["qseqid"].eq("").any():
        raise ValueError("Blank qseqid in annotation table")
    # Ambiguous contigs cannot be used for accepted cluster annotation.
    conflicts = set(x.loc[x.qualified].groupby("qseqid")["species"]
                    .nunique().loc[lambda z: z > 1].index)
    x = x.sort_values(["qseqid", "qualified", "bitscore", "evalue"],
                      ascending=[True, False, False, True], kind="stable")
    best = x.drop_duplicates("qseqid").copy()
    best.loc[best.qseqid.isin(conflicts), "qualified"] = False
    print(f"Annotation rows: {len(x)}; distinct contigs: {len(best)}; "
          f"contigs with conflicting labels: {len(conflicts)}")
    return best[["qseqid", "species", "qualified"]]


def filtered_features(identity, condition):
    if condition == "NoSingletons":
        path = (ROOT / f"statistical_analysis/run_based_analysis/{identity}_aS60_aL60/"
                "informative_runs_only.tsv")
    else:
        name = condition.lower()
        path = (ROOT / f"statistical_analysis/run_based_analysis/{name}_observations/"
                f"{identity}_aS60_aL60_informative_runs_{name}.tsv")
    if not path.exists():
        raise FileNotFoundError(path)
    header = pd.read_csv(path, sep="\t", nrows=0)
    cols = [c for c in header.columns if c.startswith("Cluster_")]
    if not cols:
        raise ValueError(f"No cluster columns: {path}")
    return cols


def annotate(identity, best):
    member_path = MAT / f"{identity}_aS60_aL60_cluster_membership.tsv"
    m = pd.read_csv(member_path, sep="\t", dtype=str)
    required = {"cluster", "metadata_qseqid", "sample", "country"}
    if not required <= set(m.columns):
        raise ValueError(f"Missing membership fields: {required - set(m.columns)}")
    m = m.rename(columns={"cluster": "cluster_id", "metadata_qseqid": "qseqid"})
    m = m.drop_duplicates(["cluster_id", "qseqid", "sample"])
    if m.qseqid.isna().any() or m.cluster_id.isna().any():
        raise ValueError("Missing membership identifiers")
    if m.groupby("qseqid").cluster_id.nunique().gt(1).any():
        raise ValueError("One contig occurs in multiple clusters")
    joined = m.merge(best, on="qseqid", how="left", validate="many_to_one")
    joined["qualified"] = joined.qualified.fillna(False).astype(bool)
    rows = []
    for cid, g in joined.groupby("cluster_id", sort=True):
        names = g.loc[g.qualified, "species"].value_counts()
        n, support = len(g), int(names.sum())
        frac = support / n
        top = str(names.index[0]) if len(names) else ""
        if support == 0:
            status = "unannotated_or_low_alignment_support"
        elif len(names) > 1:
            status = "conflicting_reference_species"
        elif frac < MIN_ANNOTATED_FRAC:
            status = "insufficient_annotated_members"
        elif n == 1:
            status = "single_contig_putative_match"
        else:
            status = "consistent_putative_reference_match"
        accepted = top if status in {
            "single_contig_putative_match", "consistent_putative_reference_match"
        } else ""
        # Geography is cluster occurrence by distinct SRRs (not confirmed virus presence).
        geography = (g[["sample", "country"]].drop_duplicates()
                     .groupby("country")["sample"].nunique().sort_index())
        rows.append(dict(
            cluster_id=cid, accepted_reference_species=accepted,
            top_reference_species=top, annotation_status=status,
            countries_detected="; ".join(str(c) for c in geography.index),
            country_run_counts="; ".join(f"{c}: {int(k)}" for c, k in geography.items()),
            n_member_contigs=n, n_member_runs=g["sample"].nunique(),
            n_supported_hits=support, annotated_member_fraction=round(frac, 6),
            supporting_species_counts="; ".join(f"{v}:{c}" for v, c in names.items())
        ))
    result = pd.DataFrame(rows).set_index("cluster_id", drop=False)
    fullpath = MAT / f"{identity}_aS60_aL60_run_with_singletons.tsv"
    fullcols = pd.read_csv(fullpath, sep="\t", nrows=0).columns
    fullfeatures = [c for c in fullcols if c.startswith("Cluster_")]
    if set(fullfeatures) != set(result.index):
        raise ValueError(f"Cluster membership IDs disagree with {fullpath}")
    conditions = {"All": fullfeatures}
    for name in ("NoSingletons", "Min5", "Min6"):
        conditions[name] = filtered_features(identity, name)
    for name, ids in conditions.items():
        missing = set(ids) - set(result.index)
        if missing:
            raise ValueError(f"Unknown {identity}% {name} cluster IDs: {sorted(missing)}")
        subset = result.loc[ids].reset_index(drop=True)
        path = OUT / f"{identity}_Cluster_Virus_Annotation_{name}.csv"
        subset.to_csv(path, index=False)
        n_accepted = subset.accepted_reference_species.ne("").sum()
        print(f"{identity}% {name}: {len(subset)} clusters, "
              f"{n_accepted} accepted putative reference labels → {path}")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    best = load_hits()
    for identity in (90, 70):
        annotate(identity, best)
    (OUT / "ANNOTATION_METHODS.txt").write_text(
        "Existing DIAMOND reference hits; aligned AA >=50, identity >=30%, "
        "E-value <=1e-5. Cluster assignment accepted only when all qualifying "
        "member hits agree on a reference label and >=50% of contigs qualify. "
        "Single-member clusters with qualifying hits are tentatively labelled. "
        "All other labels remain unresolved. These are PUTATIVE reference "
        "matches, not confirmed species. Min5 and Min6 use exact feature "
        "sets in the filtered matrices (prevalence across distinct SRRs). "
        "countries_detected and country_run_counts are based on cluster membership "
        "in the original runs, NOT confirmed presence of the labelled virus.\n"
    )


if __name__ == "__main__":
    main()

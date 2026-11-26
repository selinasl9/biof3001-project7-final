#!/usr/bin/env python3
"""BIOF3001 Project 7 — Filter run-level viral-cluster matrices.

Objective: derive no-singleton (>=2 SRRs), min5 (>=5 SRRs) and min6
(>=6 SRRs) detection matrices from the eight unfiltered cluster conditions.
Input: results/viral_clustering/matrices/*_run_with_singletons.tsv.
Outputs: results/statistical_analysis/run_based_analysis/ including both
all-SRR matrices and informative-only SRR subsets for ordination and tests.
Presence is binary (0/1), NOT viral abundance. The prevalence denominator
for analyses using informative-only runs differs from all-strict-positive runs.
Requirements: Python 3, pandas.
Usage: python scripts/viral_clustering/filter_matrices.py [--help]
Responsible contributor: [confirm before submission]
"""
import argparse
import os
from pathlib import Path
import pandas as pd

# STEP 1 — Define the 90%/70% and coverage combinations processed independently.
CONDITIONS=[f'{i}_{c}' for i in (90,70) for c in ('identity','aS60','aL60','aS60_aL60')]

# STEP 2 — Centralize TSV output writing and create directories when needed.
def save(d,path):
    path.parent.mkdir(parents=True,exist_ok=True)
    d.to_csv(path,sep='\t',index=False)

def main():
    # All documented relative paths are resolved from the repository root.
    os.chdir(Path(__file__).resolve().parents[2])
    ap=argparse.ArgumentParser()
    ap.add_argument('--input-dir',default='results/viral_clustering/matrices')
    ap.add_argument('--output-dir',default='results/statistical_analysis/run_based_analysis')
    args=ap.parse_args()
    input_dir=Path(args.input_dir); out=Path(args.output_dir)
    summary=[]
# STEP 3 — Read each unfiltered SRR × viral-cluster detection matrix.
    for condition in CONDITIONS:
        path=input_dir/f'{condition}_run_with_singletons.tsv'
        if not path.exists(): raise FileNotFoundError(path)
        df=pd.read_csv(path,sep='\t')
        if df['sample'].duplicated().any():raise ValueError(f'{path}: duplicate SRRs')
        clusters=[c for c in df if c.startswith('Cluster_')]
        if not clusters:raise ValueError(f'{path}: no clusters')
# STEP 4 — Enforce genuinely binary (0/1), nonmissing detection data.
        X=df[clusters].apply(pd.to_numeric,errors='raise')
        if X.isna().any().any() or not X.isin([0,1]).all().all():
            raise ValueError(f'{path}: expected clean binary cluster data')
# STEP 5 — Keep clusters detected in at least 2, 5 or 6 SRRs.
# The threshold counts RUNS, not the number of contigs or species.
        for threshold,label in [(2,'no_singletons'),(5,'min5'),(6,'min6')]:
            cols=X.columns[X.sum(axis=0)>=threshold].tolist()
            if not cols:
                print(f'WARNING: {condition} {label}: no retained clusters');continue
            filtered=df[['sample','country']+cols]
# STEP 6 — Distinguish the full retained-feature matrix from informative
# runs with >=1 detection; informative-only analyses have different denominators.
            nonzero=filtered[cols].sum(axis=1)>0
            info=filtered.loc[nonzero].copy()
# STEP 7 — Preserve output filenames expected by the statistical R scripts.
            if label=='no_singletons':
                # Legacy path compatible with original PERMANOVA / KW files.
                directory=out/condition
                save(filtered,directory/'run_without_singletons.tsv')
                save(info,directory/'informative_runs_only.tsv')
            else:
                directory=out/f'{label}_observations'
                save(filtered,directory/f'{condition}_run_matrix_{label}.tsv')
                save(info,directory/f'{condition}_informative_runs_{label}.tsv')
# STEP 8 — Audit the number of informative runs from each country.
            counts=info.groupby('country').size().reset_index(name='n_informative_SRRs')
            save(counts,directory/f'{condition}_country_counts_{label}.tsv')
            summary.append(dict(condition=condition,filter=label,clusters=len(cols),
                                all_SRRs=len(df),informative_SRRs=len(info),
                                informative_countries=info.country.nunique()))
            print(f'{condition:20s} {label:15s}: {len(cols):4} clusters, {len(info):4} informative SRRs')
# STEP 9 — Write one summary across all identity/coverage/threshold combinations.
    save(pd.DataFrame(summary),out/'cluster_filter_summary.tsv')

if __name__=='__main__':main()

#!/usr/bin/env python3
"""BIOF3001 Project 7 — Descriptive country-level viral-cluster overlap.

Objective: generate country shared-cluster barplots, pairwise overlaps and
UpSet-like intersection plots using country-level prevalence matrices.
Inputs: results/viral_clustering/matrices/*_country_with_singletons.tsv
Outputs: country overlap TSVs, barplot, heatmap and intersection PNGs.
Run: python scripts/statistical_analysis/plot_country_overlap.py
Dependencies: Python 3, pandas, numpy, matplotlib.
Shared means detected in >=2 countries before minimum-SRR country filter.
These are descriptive plots, not statistical hypothesis tests.
Responsible contributor: [confirm before submission]
"""
"""Country prevalence/overlap plots for all countries and both identities.

Sources are *COUNTRY-LEVEL prevalence matrices*: one row per country.
For shared-cluster plots, keep clusters observed in >=2 countries across the
FULL original set of countries before any group sample-size filtering.
This is not the SRR-level PERMANOVA dataset.
"""
from pathlib import Path
import argparse
from collections import Counter
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# STEP 1 — Analyse one clustering condition from the country-level
# prevalence table, not the SRR-level PERMANOVA matrix.
def plot_country(condition,indir,out,min_srr):
    path=indir/f'{condition}_country_with_singletons.tsv'
    if not path.exists():raise FileNotFoundError(path)
    df=pd.read_csv(path,sep='\t')
    cols=[col for col in df if col.startswith('Cluster_')]
    if not cols or not {'country','n_samples'}.issubset(df.columns):raise ValueError(f'Invalid country matrix: {path}')
    X=df[cols].apply(pd.to_numeric,errors='raise')
    # STEP 2 — Define globally shared clusters before filtering out
# countries with insufficient SRRs; detection means prevalence > 0.
    shared=X.columns[(X>0).sum(axis=0)>=2].tolist()
    if not shared:raise ValueError(f'{condition}: no clusters shared across at least 2 countries')
    d=df[df.n_samples>=min_srr].copy().sort_values('country')
    if d.empty:raise ValueError(f'{condition}: no eligible countries')
    # STEP 3 — Make a country × cluster detection matrix (0/1)
# for the countries that pass the minimum-run threshold.
    M=(d[shared]>0).astype(int)
    country_names=d.country.tolist()
    counts=M.sum(axis=1)
    dat=pd.DataFrame(dict(country=country_names,n_SRRs=d.n_samples,
                          n_shared_clusters_detected=counts.values))
    dat=dat.sort_values('n_shared_clusters_detected',ascending=False)
    dat.to_csv(out/f'{condition}_country_shared_cluster_counts.tsv',sep='\t',index=False)
    # STEP 4 — Plot the number of shared clusters found in each country.
    fig,ax=plt.subplots(figsize=(max(8,len(dat)*.75),5),dpi=140)
    ax.bar(dat.country,dat.n_shared_clusters_detected)
    ax.set(ylabel='Shared viral clusters detected',title=f'{condition}: shared clusters across original countries')
    ax.tick_params(axis='x',rotation=65);fig.tight_layout()
    fig.savefig(out/f'{condition}_shared_cluster_country_bar.png',dpi=250);plt.close(fig)
    # Pairwise overlap: #clusters present in both eligible countries.
    # STEP 5 — Multiply the binary matrix by its transpose to count
# clusters observed in BOTH countries for every pair of countries.
    overlap=M.to_numpy() @ M.to_numpy().T
    pd.DataFrame(overlap,index=country_names,columns=country_names).to_csv(
      out/f'{condition}_country_overlap.tsv',sep='\t')
    fig,ax=plt.subplots(figsize=(max(7,len(country_names)*.5),max(6,len(country_names)*.5)),dpi=140)
    im=ax.imshow(overlap,cmap='Blues')
    ax.set_xticks(range(len(country_names)),country_names,rotation=90)
    ax.set_yticks(range(len(country_names)),country_names)
    fig.colorbar(im,ax=ax,label='Shared clusters in both countries')
    ax.set_title(f'{condition}: pairwise viral-cluster overlap')
    fig.tight_layout();fig.savefig(out/f'{condition}_overlap_heatmap.png',dpi=200);plt.close(fig)
    # UpSet-like all-country intersection chart (exclusive combinations).
    # STEP 6 — Calculate exclusive presence patterns: for each cluster,
# record the precise set of eligible countries that detects it.
    patterns=Counter()
    for c in shared:
        key=tuple(country_names[i] for i,v in enumerate(M[c]) if v)
        # Retain intersections also if >=2 countries globally but 1 in displayed subset.
        if key:patterns[key]+=1
    ordered=sorted(patterns,key=lambda key:(-patterns[key],-len(key),key))
    top=ordered[:min(len(ordered),30)]
    fig,(ax,mat)=plt.subplots(2,1,figsize=(max(11,len(top)*.52),7),
              gridspec_kw=dict(height_ratios=[2,1.4]),sharex=True)
    ax.bar(range(len(top)),[patterns[k] for k in top]);ax.set(ylabel='Exclusive intersection size',
          title=f'{condition}: country intersections (top {len(top)} of {len(ordered)})')
    for j,k in enumerate(top):
        idxs=[country_names.index(g) for g in k]
        mat.scatter([j]*len(idxs),idxs,s=28,color='black')
        if len(idxs)>1:mat.plot([j,j],[min(idxs),max(idxs)],color='black',lw=1)
    mat.set_yticks(range(len(country_names)),country_names)
    mat.set_ylim(-.6,len(country_names)-.4);mat.invert_yaxis()
    mat.set_xlabel('Exclusive country combinations (columns)')
    mat.set_xticks(range(len(top)));mat.tick_params(axis='x',length=0,labelbottom=False)
    fig.tight_layout();fig.savefig(out/f'{condition}_all_country_upset.png',dpi=220);plt.close(fig)
    pd.DataFrame([dict(countries=';'.join(k),n_clusters=patterns[k]) for k in ordered]).to_csv(
        out/f'{condition}_all_country_intersections.tsv',sep='\t',index=False)
    print(condition,'all original countries:',len(df),'eligible:',len(d),'clusters shared across >=2 original countries:',len(shared))


# STEP 7 — Parse paths/threshold from command line and generate plots
# for both 90% CD-HIT and 70% VSEARCH conditions.
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--input-dir',type=Path,default=Path('results/viral_clustering/matrices'))
    ap.add_argument('--output-dir',type=Path,default=Path('results/statistical_analysis/country_overlap_final'))
    ap.add_argument('--min-srr-per-country',type=int,default=5)
    args=ap.parse_args();args.output_dir.mkdir(parents=True,exist_ok=True)
    for ident in (90,70):
        plot_country(f'{ident}_aS60_aL60',args.input_dir,args.output_dir,args.min_srr_per_country)

if __name__=='__main__':main()

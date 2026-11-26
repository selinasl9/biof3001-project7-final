#!/usr/bin/env python3
"""BIOF3001 Project 7 — Plot Bray–Curtis ordinations of viral-cluster detections.

Objective: visualise run-level binary cluster profiles by country or by
Chinese province using NMDS and/or t-SNE. Ordinations are exploratory and
must not be interpreted as inferential evidence of group differences.
Input: informative-run binary matrices and China province mapping.
Output: results/statistical_analysis/run_based_analysis/ordination_replots/
       (PNG plots, coordinate tables, summary table).
Requirements: Python 3, pandas, numpy, scipy, matplotlib, scikit-learn.
Usage: python scripts/statistical_analysis/plot_ordination.py --method both
A fixed seed (default 42) supports reproducible stochastic ordinations.
Responsible contributor: [confirm before submission]
"""
import argparse
import os
import inspect
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist,squareform
from sklearn.manifold import MDS,TSNE

CONDITIONS=[(i,lab) for i in (90,70) for lab in ('no_singletons','min5','min6')]

def source_path(base,ident,filter_name):
    cond=f'{ident}_aS60_aL60'
    if filter_name=='no_singletons': return base/cond/'informative_runs_only.tsv'
    return base/f'{filter_name}_observations'/f'{cond}_informative_runs_{filter_name}.tsv'

def fit(D,method,seed):
    if method=='nmds':
        p=inspect.signature(MDS).parameters
        kw=dict(n_components=2,n_init=4,max_iter=500,eps=1e-5,random_state=seed,n_jobs=1)
        if 'metric_mds' in p:kw.update(metric_mds=False,metric='precomputed')
        else:kw.update(metric=False,dissimilarity='precomputed')
        if 'normalized_stress' in p:kw['normalized_stress']=True
        m=MDS(**kw); coords=m.fit_transform(D)
        return coords,float(m.stress_)
    p=inspect.signature(TSNE).parameters
    kw=dict(n_components=2,metric='precomputed',init='random',
            perplexity=min(30.,max(2.,(len(D)-1)/3)),random_state=seed,learning_rate=200)
    kw['max_iter' if 'max_iter' in p else 'n_iter']=1000
    coords=TSNE(**kw).fit_transform(D)
    return coords,float('nan')

def analyse(source,analysis,prov,ident,fil,method,seed,out):
    if not source.is_file():raise FileNotFoundError(source)
    d=pd.read_csv(source,sep='\t')
    cols=[x for x in d if x.startswith('Cluster_')]
    if not cols or not {'sample','country'}.issubset(d.columns):raise ValueError(f'Invalid input matrix: {source}')
    if analysis=='china_province':
        d=d[d.country.astype(str).str.lower().eq('china')].copy()
        d['group']=d['sample'].astype(str).map(prov)
        d=d[d.group.notna() & d.group.ne('')].copy()
    else:d['group']=d.country.astype(str)
    X=d[cols].apply(pd.to_numeric,errors='raise').to_numpy(dtype=float)
    good=X.sum(axis=1)>0
    d=d.loc[good].reset_index(drop=True);X=X[good]
    if len(d)<5 or len(set(d.group))<2:
        print(f'SKIP {ident} {fil} {analysis}: insufficient samples/groups');return None
    if np.isnan(X).any() or (X<0).any():raise ValueError('Missing or negative cluster values')
    X=(X>0).astype(float)
    D=squareform(pdist(X,metric='braycurtis'))
    if not np.isfinite(D).all():raise ValueError('Invalid Bray-Curtis matrix')
    # Preserve zero dissimilarity between identical samples: sklearn's older
    # SMACOF variants interpreted these as missing, so small epsilon off-diagonal.
    z=(D==0)&~np.eye(len(D),dtype=bool)
    D[z]=1e-9
    coords,stress=fit(D,method,seed)
    stem=f'{ident}_{fil}_{analysis}_{method}_braycurtis'
    points=d[['sample','country','group']].copy()
    points['axis1']=coords[:,0];points['axis2']=coords[:,1]
    points.to_csv(out/f'{stem}_coordinates.tsv',sep='\t',index=False)
    fig,ax=plt.subplots(figsize=(9.5,7),dpi=125)
    cats=sorted(points.group.unique());palette=plt.get_cmap('tab20',max(20,len(cats)))
    for idx,c in enumerate(cats):
        mask=points.group.eq(c)
        ax.scatter(points.loc[mask,'axis1'],points.loc[mask,'axis2'],s=23,
                   alpha=.75,color=palette(idx),label=c,linewidths=.2,edgecolors='white')
    lab='NMDS' if method=='nmds' else 't-SNE'
    ax.set(xlabel=lab+'1',ylabel=lab+'2',title=f'{ident}% {fil} | {analysis} | {len(d)} SRRs, {len(cols)} clusters'+
           (f' | normalized stress {stress:.3f}' if np.isfinite(stress) else ''))
    ax.legend(loc='upper left',bbox_to_anchor=(1.02,1),title='Geography',fontsize=8,frameon=False)
    fig.tight_layout();fig.savefig(out/f'{stem}.png',dpi=250,bbox_inches='tight');plt.close(fig)
    print(stem,':',len(d),'SRRs,',len(cats),'groups')
    return dict(identity=ident,filter=fil,analysis=analysis,method=method,SRRs=len(d),
                groups=len(cats),viral_clusters=len(cols),normalized_stress=stress)

def main():
    # All documented relative paths are resolved from the repository root.
    os.chdir(Path(__file__).resolve().parents[2])
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--base',type=Path,default=Path('results/statistical_analysis/run_based_analysis'))
    ap.add_argument('--province-map',type=Path,default=Path('results/statistical_analysis/run_based_analysis/china_province_audit/china_srr_provinces.tsv'))
    ap.add_argument('--method',choices=['nmds','tsne','both'],default='nmds')
    ap.add_argument('--level',choices=['country','china_province','both'],default='both')
    ap.add_argument('--seed',type=int,default=42)
    args=ap.parse_args()
    if args.level in ('both', 'china_province') and not args.province_map.is_file():
        raise FileNotFoundError(f'Missing province mapping: {args.province_map}')
    prov_df = (pd.read_csv(args.province_map, sep='\t', dtype=str)
               if args.level in ('both', 'china_province') else pd.DataFrame(columns=['sample', 'province']))
    if prov_df['sample'].duplicated().any():raise ValueError('Duplicated SRRs in province map')
    prov=dict(zip(prov_df['sample'],prov_df['province']))
    out=args.base/'ordination_replots';out.mkdir(parents=True,exist_ok=True)
    methods=['nmds','tsne'] if args.method=='both' else [args.method]
    analyses=['country','china_province'] if args.level=='both' else [args.level]
    rows=[]
    for ident,fil in CONDITIONS:
        path=source_path(args.base,ident,fil)
        for analysis in analyses:
            for method in methods:
                r=analyse(path,analysis,prov,ident,fil,method,args.seed,out)
                if r:rows.append(r)
    pd.DataFrame(rows).to_csv(out/'ordination_summary.tsv',sep='\t',index=False)

if __name__=='__main__':main()

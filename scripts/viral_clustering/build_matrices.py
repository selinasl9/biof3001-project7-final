#!/usr/bin/env python3
"""BIOF3001 Project 7 — Build run-level viral-cluster detection matrices.

Objective: parse CD-HIT .clstr and VSEARCH .uc membership files, link contigs
with strict host-cleaned SRR metadata, and generate binary detection and
country-level cluster detection-prevalence matrices for eight conditions.

Inputs: strict host-cleaned hit tables for original and China supplement;
        results/viral_clustering/{cdhit,vsearch}/ cluster membership files.
Outputs: results/viral_clustering/matrices/{condition}_*.tsv.
Denominator: all SRRs with >=1 strict host-cleaned contig; zero rows are kept
for clusters not detected in a run. Representatives >=1000 nt are retained.
Requirements: Python 3, pandas.
Usage: python scripts/viral_clustering/build_matrices.py [--help]
Responsible contributor: [confirm before submission]
"""
import argparse
import os
import re
from pathlib import Path
import pandas as pd

META_DEFAULTS = [
    'results/viral_identification/host_filter/strict_host_cleaned_hits.tsv',
    'results/viral_identification/PRJNA1003643/host_filter/CN_PRJNA1003643_strict_host_cleaned_hits.tsv',
]
CONDITIONS = [f'{i}_{s}' for i in (90, 70)
              for s in ('identity', 'aS60', 'aL60', 'aS60_aL60')]


def load_metadata(paths):
    frames=[]
    for p in paths:
        if not Path(p).is_file():
            raise FileNotFoundError(f'Metadata file is missing: {p}')
        frames.append(pd.read_csv(p, sep='\t', dtype=str, usecols=['qseqid','sample','country']))
    meta=pd.concat(frames,ignore_index=True).drop_duplicates()
    for col in ['qseqid','sample','country']:
        meta[col]=meta[col].str.strip()
        if meta[col].isna().any() or meta[col].eq('').any():
            raise ValueError(f'Empty values in {col}')
    if meta.groupby('qseqid').agg(n_sample=('sample','nunique'),n_country=('country','nunique')).max().max()>1:
        raise ValueError('qseqid maps to multiple samples/countries')
    if meta.groupby('sample')['country'].nunique().max()>1:
        raise ValueError('A sample maps to multiple countries')
    return meta.drop_duplicates('qseqid')


def parse_clstr(path):
    rows=[]; cid=None
    with open(path,encoding='utf-8') as handle:
        for ln,line in enumerate(handle,1):
            line=line.strip()
            if not line: continue
            if line.startswith('>Cluster '):
                cid=f'Cluster_{int(line.split()[1]):06d}'; continue
            if cid is None: raise ValueError(f'Cluster member before cluster header at {path}:{ln}')
            m=re.search(r'\b(\d+)nt,\s*>(.+?)\.\.\.',line)
            if not m: raise ValueError(f'Invalid .clstr member at {path}:{ln}: {line}')
            rows.append((cid,m.group(2),int(m.group(1)),line.endswith('*')))
    if not rows: raise ValueError(f'No clusters found in {path}')
    df=pd.DataFrame(rows,columns=['cluster','qseqid','length_nt','representative'])
    if df.qseqid.duplicated().any(): raise ValueError('Duplicate contig IDs in .clstr')
    if not df.groupby('cluster').representative.sum().eq(1).all():
        raise ValueError('Every cluster must contain exactly one representative')
    return df


def condition_file(args,condition):
    id_str, setting=condition.split('_',1)
    suffix='' if setting=='identity' else '_'+setting
    if id_str=='90':
        return Path(args.cdhit_dir)/('strict_cdhit90'+suffix+'.fa.clstr')
    if id_str=='70':
        return Path(args.vsearch_dir)/('strict_vsearch70'+suffix+'.uc')
    raise ValueError(f'Unknown identity: {id_str}')


def parse_uc(path):
    # VSEARCH UC: S=seed/centroid, H=hit, C=summary (omit).
    rows=[]
    with open(path, encoding='utf-8') as handle:
        for ln,line in enumerate(handle,1):
            if not line.strip() or line.startswith('#'):
                continue
            fields=line.rstrip('\n').split('\t')
            if len(fields)<10: raise ValueError(f'Invalid VSEARCH UC row at {path}:{ln}')
            if fields[0] not in ('S','H'):continue
            cluster=f'Cluster_{int(fields[1]):06d}'
            query=fields[8].split()[0]
            rows.append((cluster,query,int(fields[2]),fields[0]=='S'))
    if not rows:raise ValueError(f'No VSEARCH S/H records in {path}')
    d=pd.DataFrame(rows,columns=['cluster','qseqid','length_nt','representative'])
    if d.qseqid.duplicated().any():raise ValueError(f'Duplicate VSEARCH member IDs: {path}')
    if not d.groupby('cluster').representative.sum().eq(1).all():
        raise ValueError(f'Each VSEARCH cluster needs exactly one S/seed row: {path}')
    return d


def save(frame, path):
    path.parent.mkdir(parents=True,exist_ok=True)
    frame.to_csv(path,sep='\t',index=False,float_format='%.8f')


def run_one(condition,meta,args):
    source=condition_file(args,condition)
    if not source.is_file(): raise FileNotFoundError(source)
    members=parse_clstr(source) if source.suffix=='.clstr' else parse_uc(source)
    # Reverse-complemented contigs in original project use _rev suffix.
    members['metadata_qseqid']=members.qseqid.str.replace(r'_rev$','',regex=True)
    if members.metadata_qseqid.duplicated().any():
        raise ValueError(f'{condition}: metadata IDs duplicated after removal of _rev')
    members=members.merge(meta,left_on='metadata_qseqid',right_on='qseqid',how='left',validate='one_to_one',suffixes=('','_meta'))
    if members['sample'].isna().any():
        missing=members.loc[members['sample'].isna(),'metadata_qseqid'].head(10).tolist()
        raise ValueError(f'{condition}: missing metadata for {missing}')
    if set(members.metadata_qseqid)!=set(meta.qseqid):
        raise ValueError(f'{condition}: .clstr member IDs do not exactly match metadata IDs')
    reps=members.loc[members.representative,['cluster','length_nt']].set_index('cluster').length_nt
    eligible=reps[reps>=args.min_rep_nt].index
    members=members[members.cluster.isin(eligible)].copy()
    if members.empty: raise ValueError(f'{condition}: no clusters with qualifying representative length')
    # Denominator: every strict-positive SRR, including runs with no retained clusters.
    runs=meta[['sample','country']].drop_duplicates().sort_values(['country','sample'])
    incid=members[['sample','cluster']].drop_duplicates().assign(present=1)
    wide=incid.pivot(index='sample',columns='cluster',values='present').fillna(0).astype('uint8')
    allcols=sorted(wide.columns,key=lambda x:int(x.split('_')[1]))
    mat=runs.merge(wide.reset_index(),on='sample',how='left',validate='one_to_one')
    mat[allcols]=mat[allcols].fillna(0).astype('uint8')
    counts=mat[allcols].sum(axis=0)
    no_singletons=counts.index[counts>=2].tolist()
    prefix=Path(args.out_dir)/condition
    save(mat,prefix.with_name(prefix.name+'_run_with_singletons.tsv'))
    save(mat[['sample','country']+no_singletons],prefix.with_name(prefix.name+'_run_without_singletons.tsv'))
    def country_matrix(cols):
        return mat.groupby('country',sort=True).agg(n_samples=('sample','nunique')).join(
            mat.groupby('country',sort=True)[cols].mean()).reset_index()
    save(country_matrix(allcols),prefix.with_name(prefix.name+'_country_with_singletons.tsv'))
    save(country_matrix(no_singletons),prefix.with_name(prefix.name+'_country_without_singletons.tsv'))
    summary=(members.groupby('cluster').agg(n_contigs=('metadata_qseqid','size'),
            n_SRRs=('sample','nunique'),n_countries=('country','nunique')).join(reps.rename('representative_length_nt')).reset_index())
    save(summary,prefix.with_name(prefix.name+'_cluster_summary.tsv'))
    save(members[['cluster','metadata_qseqid','sample','country','length_nt','representative']],
         prefix.with_name(prefix.name+'_cluster_membership.tsv'))
    print(f'{condition}: {len(mat)} SRRs, {len(allcols)} clusters with rep>={args.min_rep_nt}nt, {len(no_singletons)} observed in >=2 SRRs')


def main():
    # All documented relative paths are resolved from the repository root.
    os.chdir(Path(__file__).resolve().parents[2])
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--metadata',nargs='+',default=META_DEFAULTS)
    ap.add_argument('--cdhit-dir',default='results/viral_clustering/cdhit')
    ap.add_argument('--vsearch-dir',default='results/viral_clustering/vsearch')
    ap.add_argument('--out-dir',default='results/viral_clustering/matrices')
    ap.add_argument('--min-rep-nt',type=int,default=1000)
    args=ap.parse_args()
    meta=load_metadata(args.metadata)
    for c in CONDITIONS:run_one(c,meta,args)

if __name__=='__main__': main()

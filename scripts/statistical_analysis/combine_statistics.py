#!/usr/bin/env python3
"""BIOF3001 Project 7 — Combine pairwise PERMANOVA and PERMDISP.

Objective: merge pairwise multivariate composition and dispersion tests by
exact comparison keys, and summarise significance at BH-adjusted p < 0.05.
Inputs: pairwise_permanova_results.tsv and pairwise_permdisp_lingoes.tsv.
Outputs: combined pairwise table and significance-count summary.
Run: python scripts/statistical_analysis/combine_statistics.py
Dependencies: Python 3, pandas.
Responsible contributor: [confirm before submission]
"""
"""Combine pairwise PERMANOVA and Lingoes-corrected PERMDISP without row loss."""
from pathlib import Path
import pandas as pd

# STEP 1 — Locate the two pairwise result tables produced by the R tests.
out=Path('results/statistical_analysis/run_based_analysis/pairwise_permanova_braycurtis')
# STEP 2 — Read both results and define a unique geographic comparison
# key, including sample counts to detect inconsistent filtering.
a=pd.read_csv(out/'pairwise_permanova_results.tsv',sep='\t')
b=pd.read_csv(out/'pairwise_permdisp_lingoes.tsv',sep='\t')
keys=['identity','filter','analysis','group1','group2','n_group1','n_group2']
assert not a.duplicated(keys).any() and not b.duplicated(keys).any(),'Duplicate pairwise key'
# STEP 3 — Merge by exact keys. The outer join plus indicator ensures
# no pairwise comparison is silently dropped from either test.
merged=a.merge(b,on=keys,how='outer',validate='one_to_one',indicator=True)
if not merged['_merge'].eq('both').all():
    raise ValueError('PERMANOVA and PERMDISP pair sets do not match')
merged=merged.drop(columns='_merge')
# STEP 4 — Flag FDR-adjusted significance separately for community
# composition and dispersion. A PERMANOVA-only contrast is not proof
# that dispersion is identical, but avoids one common confound.
merged['PERMANOVA_significant_BH']=merged.p_BH<0.05
merged['PERMDISP_significant_BH']=merged.PERMDISP_p_BH<0.05
merged['PERMANOVA_sig_PERMDISP_not']=merged['PERMANOVA_significant_BH']&~merged['PERMDISP_significant_BH']
# STEP 5 — Export the complete comparison table and aggregate the
# number of significant pairs per clustering/filter/geography analysis.
merged.to_csv(out/'pairwise_permanova_permdisp_combined_lingoes.tsv',sep='\t',index=False)
summary=(merged.groupby(['identity','filter','analysis'],as_index=False)
         .agg(pairs_tested=('group1','size'),
              PERMANOVA_significant=('PERMANOVA_significant_BH','sum'),
              PERMDISP_significant=('PERMDISP_significant_BH','sum'),
              PERMANOVA_sig_PERMDISP_not=('PERMANOVA_sig_PERMDISP_not','sum')))
summary.to_csv(out/'pairwise_significance_summary_lingoes.tsv',sep='\t',index=False)
print(summary.to_string(index=False));print('Wrote combined pairwise table')

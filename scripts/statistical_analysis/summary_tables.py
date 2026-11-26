#!/usr/bin/env python3
"""BIOF3001 Project 7 — Final multivariate and binary-detection summaries.

Objective: merge pairwise PERMANOVA/PERMDISP results with the binary
any-retained-cluster detection Dunn/Fisher table, and summarise how many
geographical pairs show statistical significance across methods.
Inputs: combined PERMANOVA/PERMDISP table; ALL_dunn_fisher_pairs.tsv.
Outputs: pairwise_all_methods_comparison.tsv and condition_summary.tsv.
Run: python scripts/statistical_analysis/summary_tables.py
Dependencies: Python 3, pandas.
Caution: Dunn/Fisher here test ANY retained-cluster detection (0/1), NOT
individual-cluster enrichment. A significant PERMANOVA is different.
Responsible contributor: [confirm before submission]
"""
from pathlib import Path
import pandas as pd

BASE = Path('results/statistical_analysis/run_based_analysis')
COMBINED = BASE / 'pairwise_permanova_braycurtis/pairwise_permanova_permdisp_combined_lingoes.tsv'
DETECTION = BASE / 'kw_dunn_binary_detection/all_runs/ALL_dunn_fisher_pairs.tsv'
OUT = BASE / 'final_summary'


# STEP 1 — Read multivariate pairwise results and the separate binary
# any-retained-cluster detection comparisons.
def main():
    for path in (COMBINED, DETECTION):
        if not path.is_file():
            raise FileNotFoundError(f'Missing required result: {path}')
    multi = pd.read_csv(COMBINED, sep='\t')
    det = pd.read_csv(DETECTION, sep='\t')
    # The detection script uses "scope"; PERMANOVA uses "analysis".
    # STEP 2 — Harmonise key column names across the R scripts.
    # Distinguish binary-test sample counts from multivariate counts.
    det = det.rename(columns={'scope': 'analysis', 'n1': 'detection_n_group1',
                              'n2': 'detection_n_group2'})
    keys = ['identity','filter','analysis','group1','group2']
    select = keys + ['Dunn_p_BH','Fisher_p_BH',
                     'prevalence_pct1','prevalence_pct2',
                     'detection_n_group1','detection_n_group2']
    # STEP 3 — Verify geographic pair keys are complete and unique, so
    # the merge cannot duplicate or silently conflate contrasts.
    for label, frame in [('PERMANOVA/PERMDISP', multi), ('KW/Dunn/Fisher', det)]:
        missing = set(keys) - set(frame.columns)
        if missing:
            raise ValueError(f'{label} missing join keys: {sorted(missing)}')
        if frame.duplicated(keys).any():
            raise ValueError(f'{label} has duplicate geographical pairs')
    missing_cols = set(select) - set(det.columns)
    if missing_cols:
        raise ValueError(f'Detection table missing columns: {sorted(missing_cols)}')
    # STEP 4 — Retain every PERMANOVA/PERMDISP pair and attach
    # binary detection statistics only when the same pair was tested.
    merged = multi.merge(det[select], on=keys, how='left', validate='one_to_one',
                         indicator=True)
    # The methods use different inclusion rules; unmatched pairs are retained and flagged.
    merged['detection_pair_available'] = merged['_merge'].eq('both')
    merged = merged.drop(columns=['_merge'])
    # STEP 5 — Produce separate FDR-significance flags for each method.
    # These methods answer different scientific questions.
    merged['PERMANOVA_BH_sig'] = merged['p_BH'] < 0.05
    merged['PERMDISP_BH_sig'] = merged['PERMDISP_p_BH'] < 0.05
    merged['Dunn_BH_sig'] = merged['Dunn_p_BH'] < 0.05
    merged['Fisher_BH_sig'] = merged['Fisher_p_BH'] < 0.05
    OUT.mkdir(parents=True, exist_ok=True)
    merged.to_csv(OUT / 'pairwise_all_methods_comparison.tsv', sep='\t', index=False)
    # STEP 6 — Count significant pairs and successful matches by
    # identity, frequency filter and geographic analysis.
    summary = merged.groupby(['identity','filter','analysis'],as_index=False).agg(
        n_pairs=('group1','size'),
        PERMANOVA_sig=('PERMANOVA_BH_sig','sum'),
        PERMDISP_sig=('PERMDISP_BH_sig','sum'),
        detection_pairs_matched=('detection_pair_available','sum'),
        Dunn_sig=('Dunn_BH_sig','sum'),
        Fisher_sig=('Fisher_BH_sig','sum')
    )
    summary.to_csv(OUT / 'condition_summary.tsv', sep='\t', index=False)
    print(summary.to_string(index=False))
    print('NOTE: KW/Dunn/Fisher test detection of ANY retained cluster, not each cluster.')
    if not merged['detection_pair_available'].all():
        print('WARNING: some multivariate pairs have no matching binary-detection comparison.')


if __name__ == '__main__':
    main()

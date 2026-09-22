#!/usr/bin/env python3
"""Calculate Figure 1 CAP and DFS statistics for the 21-patient trial cohort."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from scipy.stats import rankdata, spearmanr


def cox(frame):
    fitter = CoxPHFitter()
    fitter.fit(frame[['dfs_months', 'dfs_event', 'ccdc3']],
               duration_col='dfs_months', event_col='dfs_event')
    row = fitter.summary.loc['ccdc3']
    return {key: float(row[key]) for key in ['exp(coef)', 'exp(coef) lower 95%',
                                             'exp(coef) upper 95%', 'p']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--trial', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--permutations', type=int, default=100_000)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.trial)
    # Accept either the dedicated Figure 1 table or the wider Figure 4 table.
    if {'ccdc3_z', 'time', 'event', 'response_score'}.issubset(df.columns):
        df['ccdc3'] = df.ccdc3_z
        df['dfs_months'] = df.time
        df['dfs_event'] = df.event
        df['cap_score'] = df.response_score
    if len(df) != 21 or df.patient.nunique() != 21 or df.dfs_event.sum() != 15:
        raise ValueError('Expected 21 unique trial patients and 15 DFS events')
    # Standardise explicitly using the sample SD, matching the reported unit.
    df['ccdc3'] = (df.ccdc3 - df.ccdc3.mean()) / df.ccdc3.std(ddof=1)
    cap = spearmanr(df.cap_score, df.ccdc3)
    cap_rank = rankdata(df.cap_score.to_numpy()).astype(float)
    ccdc3_rank = rankdata(df.ccdc3.to_numpy()).astype(float)
    cap_rank -= cap_rank.mean()
    ccdc3_rank -= ccdc3_rank.mean()
    norm = np.linalg.norm(cap_rank) * np.linalg.norm(ccdc3_rank)
    observed = float((cap_rank @ ccdc3_rank) / norm)
    rng = np.random.default_rng(20260916)
    extreme = 0
    remaining = args.permutations
    while remaining:
        size = min(remaining, 10_000)
        order = np.argsort(rng.random((size, len(df))), axis=1)
        null = (ccdc3_rank[order] @ cap_rank) / norm
        extreme += int(np.count_nonzero(np.abs(null) >= abs(observed) - 1e-12))
        remaining -= size
    full = cox(df)
    definite = cox(df[df.cap_score.isin([2, 3])])
    leave_one_out = [cox(df.drop(index=i))['exp(coef)'] for i in df.index]
    result = pd.DataFrame([{
        'n_patients': len(df), 'dfs_events': int(df.dfs_event.sum()),
        'cap_spearman_rho': cap.statistic, 'cap_asymptotic_p': cap.pvalue,
        'cap_permutation_p': (extreme + 1) / (args.permutations + 1),
        'permutation_seed': 20260916, 'permutations': args.permutations,
        'dfs_hr_per_sd': full['exp(coef)'],
        'dfs_hr_ci_lower': full['exp(coef) lower 95%'],
        'dfs_hr_ci_upper': full['exp(coef) upper 95%'], 'dfs_p': full['p'],
        'definite_residual_n': int(df.cap_score.isin([2, 3]).sum()),
        'definite_residual_hr': definite['exp(coef)'],
        'definite_residual_p': definite['p'],
        'loo_hr_min': min(leave_one_out), 'loo_hr_max': max(leave_one_out),
    }])
    result.to_csv(args.output / 'figure1_trial_results.csv', index=False)
    print(result.to_string(index=False))


if __name__ == '__main__':
    main()

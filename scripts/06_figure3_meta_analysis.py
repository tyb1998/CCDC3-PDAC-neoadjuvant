#!/usr/bin/env python3
"""Perform the Figure 3 inverse-variance meta-analysis."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm


def fixed_effect(rows):
    beta = np.log(rows.hr.to_numpy(float))
    # Recover the standard error from the reported confidence limits.
    se = (np.log(rows.ci_high.to_numpy(float)) -
          np.log(rows.ci_low.to_numpy(float))) / (2 * 1.959963984540054)
    weight = 1 / se**2
    estimate = np.sum(weight * beta) / np.sum(weight)
    pooled_se = np.sqrt(1 / np.sum(weight))
    z = estimate / pooled_se
    return {'hr': np.exp(estimate),
            'ci_low': np.exp(estimate - 1.959963984540054 * pooled_se),
            'ci_high': np.exp(estimate + 1.959963984540054 * pooled_se),
            'p': 2 * norm.sf(abs(z))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cohort-cox', required=True, type=Path)
    parser.add_argument('--adjusted-meta', required=True, type=Path)
    parser.add_argument('--vascular-meta', required=True, type=Path)
    parser.add_argument('--cin-meta', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    cohort = pd.read_csv(args.cohort_cox)
    adjusted = pd.read_csv(args.adjusted_meta)
    vascular = pd.read_csv(args.vascular_meta)
    cin = pd.read_csv(args.cin_meta)
    unadjusted = cohort.loc[cohort.model.str.lower() == 'unadjusted'].copy()
    if unadjusted.cohort.nunique() != 6:
        raise ValueError('Expected six unadjusted external cohorts')
    recalculated = fixed_effect(unadjusted)
    reported = adjusted.loc[adjusted.model == 'unadjusted'].iloc[0]
    results = pd.DataFrame([{
        'cohorts': unadjusted.cohort.nunique(),
        'patients_sum': int(unadjusted.n.sum()),
        'events_sum': int(unadjusted.events.sum()),
        **{f'recalculated_{key}': value for key, value in recalculated.items()},
        'reported_hr': reported.hr, 'reported_ci_low': reported.ci_low,
        'reported_ci_high': reported.ci_high, 'reported_p': reported.p,
        'absolute_hr_difference': abs(recalculated['hr'] - reported.hr)
    }])
    results.to_csv(args.output / 'figure3_fixed_effect_results.csv', index=False)
    adjusted.to_csv(args.output / 'figure3_adjusted_meta_copy.csv', index=False)
    vascular.to_csv(args.output / 'figure3_vascular_meta_copy.csv', index=False)
    cin.to_csv(args.output / 'figure3_cin_meta_copy.csv', index=False)
    print(results.to_string(index=False))


if __name__ == '__main__':
    main()

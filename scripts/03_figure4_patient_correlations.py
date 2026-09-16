#!/usr/bin/env python3
"""Recheck patient-level Figure 4 correlations against retained analysis tables.

These are intermediate, processed patient-level tables, not raw GEO reads. The
script validates reported unadjusted correlations and independent spatial
replication. Partial Spearman uses rank residuals and a covariate-adjusted t
test (n-3 degrees of freedom); its definition is made explicit to expose any
differences from the unpublished original workflow.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata, t


def partial_spearman(a, b, control):
    a, b, control = (rankdata(np.asarray(v, dtype=float)) for v in (a, b, control))
    if any(np.std(v) == 0 for v in (a, b, control)):
        raise ValueError('A variable has zero rank variance')
    design = np.column_stack([np.ones(len(a)), control])
    ra = a - design @ np.linalg.lstsq(design, a, rcond=None)[0]
    rb = b - design @ np.linalg.lstsq(design, b, rcond=None)[0]
    rho = float(np.corrcoef(ra, rb)[0, 1])
    degrees = len(a) - 3
    p = float(2 * t.sf(abs(rho) * np.sqrt(degrees / (1 - rho**2)), degrees))
    return rho, p


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--trial', required=True, type=Path)
    parser.add_argument('--spatial', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    trial = pd.read_csv(args.trial)
    spatial = pd.read_parquet(args.spatial)
    if len(trial) != 21 or trial.patient.nunique() != 21 or len(spatial) != 39:
        raise ValueError('Unexpected patient counts or duplicate trial patients')
    rows = []
    for gene, field in [('CDT1', 'cdt1'), ('CIN/replication', 'cin_replication'),
                        ('PKMYT1', 'pkmyt1'), ('PLK1', 'plk1'), ('CCNE1', 'ccne1')]:
        raw = spearmanr(trial.ccdc3, trial[field])
        adj_r, adj_p = partial_spearman(trial.ccdc3, trial[field], trial.generic_vascular)
        rows.append({'dataset': 'trial', 'target': gene, 'n_patients': len(trial),
                     'spearman_rho': raw.statistic, 'spearman_p': raw.pvalue,
                     'partial_spearman_rho': adj_r, 'partial_p_t_n_minus_3': adj_p})
    for gene in ['CCNE1', 'PKMYT1', 'PLK1']:
        raw = spearmanr(spatial.CCDC3, spatial[gene])
        rows.append({'dataset': 'spatial', 'target': gene, 'n_patients': len(spatial),
                     'spearman_rho': raw.statistic, 'spearman_p': raw.pvalue})
    pd.DataFrame(rows).to_csv(args.output / 'figure4_patient_correlations.csv', index=False)
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == '__main__':
    main()

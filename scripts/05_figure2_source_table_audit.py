#!/usr/bin/env python3
"""Audit Figure 2 values available in retained source/summary tables.

The Hwang single-nucleus h5ad is required to reproduce panels A-C. The spatial
tables accepted here are processed summaries and cannot regenerate spot-level
partial correlations.
"""
import argparse
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compartments', required=True, type=Path)
    parser.add_argument('--specificity', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    comp = pd.read_csv(args.compartments)
    spec = pd.read_csv(args.specificity)
    required = {'patient', 'compartment', 'n_spots',
                'mean_ccdc3_logcpm', 'ccdc3_detection'}
    if not required.issubset(comp.columns):
        raise ValueError(f'Missing compartment columns: {required - set(comp.columns)}')
    if comp.patient.nunique() != 39:
        raise ValueError(f'Expected 39 spatial patients; got {comp.patient.nunique()}')
    # The supplied CSV contains a blank separator followed by five compact
    # manuscript-display rows in its first three columns.
    full = spec.loc[pd.to_numeric(spec.n_patients, errors='coerce').eq(39) &
                    pd.to_numeric(spec.n_rois, errors='coerce').eq(108)].copy()
    display = spec.loc[pd.to_numeric(spec.n_patients, errors='coerce').notna() &
                       ~pd.to_numeric(spec.n_patients, errors='coerce').eq(39) &
                       spec.signature.notna(),
                       ['signature', 'n_patients', 'n_rois']].copy()
    display.columns = ['signature', 'displayed_partial_rho',
                       'displayed_fdr']
    if set(full.n_patients.astype(int)) != {39} or set(full.n_rois.astype(int)) != {108}:
        raise ValueError('Expected every specificity row to summarize 39 patients/108 ROIs')
    comp_summary = (comp.groupby('compartment', as_index=False)
                    .agg(n_patients=('patient', 'nunique'),
                         total_spots=('n_spots', 'sum'),
                         median_patient_mean_ccdc3=('mean_ccdc3_logcpm', 'median'),
                         median_patient_detection=('ccdc3_detection', 'median')))
    comp_summary.to_csv(args.output / 'figure2_compartment_audit.csv', index=False)
    columns = ['signature', 'n_patients', 'n_rois',
               'median_partial_rho_library_vascular',
               'wilcoxon_p_partial_rho_library_vascular',
               'wilcoxon_fdr_partial_rho_library_vascular']
    full[columns].to_csv(args.output / 'figure2_spatial_specificity_audit.csv',
                         index=False)
    display.to_csv(args.output / 'figure2_manuscript_display_values.csv',
                   index=False)
    print(comp_summary.to_string(index=False))
    print('\nFull specificity rows:', len(full))
    print('\nManuscript display rows:\n', display.to_string(index=False))


if __name__ == '__main__':
    main()

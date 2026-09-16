#!/usr/bin/env python3
"""Recalculate descriptive Figure 5 measures from counts and final annotation.

This script uses only the retained barcodes in the supplied annotation file.
It does not reconstruct the original filtering, clustering or annotation.
"""
import argparse
import csv
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread


def read_lines(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as stream:
        return [line.strip().split('\t')[0] for line in stream]


def load_sample(directory, genes):
    path = Path(directory)
    features = pd.read_csv(path / 'features.tsv.gz', sep='\t', header=None, compression='gzip')
    symbols = features.iloc[:, 1].astype(str).to_numpy()
    barcodes = read_lines(path / 'barcodes.tsv.gz')
    counts = mmread(str(path / 'matrix.mtx.gz')).tocsr()
    if counts.shape != (len(symbols), len(barcodes)):
        raise ValueError('Matrix and feature/barcode dimensions disagree')
    total = np.asarray(counts.sum(axis=0)).ravel()
    expr = {}
    for gene in genes:
        idx = np.flatnonzero(symbols == gene)
        if len(idx) != 1:
            raise ValueError(f'Expected one feature for {gene}, got {len(idx)}')
        expr[gene] = np.asarray(counts[idx[0], :].toarray()).ravel()
    return barcodes, total, expr


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--annotation', required=True, type=Path,
                        help='Figure5_PDAC_only_cell_metadata.csv.gz')
    parser.add_argument('--p1302539', required=True, type=Path)
    parser.add_argument('--p1368090', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(args.annotation, dtype={'patient_id': str})
    if 'Unnamed: 0' not in meta:
        raise ValueError('The annotation must contain the prefixed nucleus ID in Unnamed: 0')
    genes = ['CCDC3', 'CCNE1', 'PLK1']
    all_rows = []
    prefilter = []
    for sid, directory in [('P1302539', args.p1302539), ('P1368090', args.p1368090)]:
        barcodes, totals, expr = load_sample(directory, genes)
        patient = sid[1:]
        prefilter.append({'patient_id': patient, 'n_nuclei': len(barcodes),
                          'CCDC3_pseudobulk_cpm': 1e6 * expr['CCDC3'].sum() / totals.sum()})
        sample_meta = meta.loc[meta.patient_id == patient].copy()
        sample_meta['barcode'] = sample_meta['Unnamed: 0'].str.replace(
            f'^{patient}_', '', regex=True)
        if sample_meta['barcode'].duplicated().any():
            raise ValueError(f'Duplicate annotation barcodes for {sid}')
        indices = pd.Index(barcodes).get_indexer(sample_meta.barcode)
        if (indices < 0).any():
            raise ValueError(f'{sid}: annotated barcode absent from count matrix')
        sample_meta['total_umi'] = totals[indices]
        for gene in genes:
            sample_meta[f'{gene}_umi'] = expr[gene][indices]
            sample_meta[f'{gene}_cpm'] = 1e6 * expr[gene][indices] / np.maximum(totals[indices], 1)
        all_rows.append(sample_meta)
    cells = pd.concat(all_rows, ignore_index=True)
    if len(cells) != 18021:
        raise ValueError(f'Expected 18,021 retained nuclei; got {len(cells):,}')
    counts = (cells.groupby(['patient_id', 'response', 'manual_broad_celltype'], dropna=False)
              .size().rename('n_nuclei').reset_index())
    sample_totals = cells.groupby('patient_id').size()
    counts['fraction_of_retained_nuclei'] = counts['n_nuclei'] / counts['patient_id'].map(sample_totals)
    counts.to_csv(args.output / 'celltype_composition.csv', index=False)
    metrics = []
    for (patient, response, celltype), group in cells.groupby(
            ['patient_id', 'response', 'manual_broad_celltype'], dropna=False):
        row = {'patient_id': patient, 'response': response,
               'cell_type': celltype, 'n_nuclei': len(group)}
        for gene in genes:
            row[f'{gene}_mean_cpm'] = group[f'{gene}_cpm'].mean()
            row[f'{gene}_detected_fraction'] = (group[f'{gene}_umi'] > 0).mean()
        metrics.append(row)
    pd.DataFrame(metrics).to_csv(args.output / 'gene_expression_by_celltype.csv', index=False)
    pd.DataFrame(prefilter).to_csv(args.output / 'ccdc3_prefilter_sample_cpm.csv', index=False)
    pseudobulk = []
    for (patient, response, celltype), group in cells.groupby(
            ['patient_id', 'response', 'manual_broad_celltype'], dropna=False):
        row = {'patient_id': patient, 'response': response, 'cell_type': celltype,
               'n_nuclei': len(group)}
        for gene in genes:
            # Manuscript Fig. 5 lineage and malignant-gene CPM uses this definition.
            row[f'{gene}_pseudobulk_cpm'] = (
                1e6 * group[f'{gene}_umi'].sum() / group['total_umi'].sum())
            row[f'{gene}_detected_fraction'] = (group[f'{gene}_umi'] > 0).mean()
        pseudobulk.append(row)
    pd.DataFrame(pseudobulk).to_csv(args.output / 'gene_pseudobulk_by_celltype.csv', index=False)
    malignant = cells.loc[cells.manual_broad_celltype == 'Epithelial (malignant)']
    cycling = (malignant.groupby(['patient_id', 'response'])
               .agg(malignant_nuclei=('manual_cell_state', 'size'),
                    cycling_nuclei=('manual_cell_state',
                                    lambda s: int((s == 'Cycling malignant').sum())))
               .reset_index())
    cycling['cycling_fraction_of_malignant'] = (
        cycling.cycling_nuclei / cycling.malignant_nuclei)
    cycling.to_csv(args.output / 'malignant_cycling_fraction.csv', index=False)
    print('Retained nuclei:', sample_totals.to_dict(), 'total', len(cells))


if __name__ == '__main__':
    main()

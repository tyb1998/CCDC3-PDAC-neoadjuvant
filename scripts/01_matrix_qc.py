#!/usr/bin/env python3
"""Audit filtered CeleScope gene-by-nucleus matrices without applying secondary QC.

Example:
  python scripts/01_matrix_qc.py --sample P1302539=data/P1302539 \
      --sample P1368090=data/P1368090 --output qc_audit
"""
import argparse
import csv
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.io import mmread


def open_text(path):
    return gzip.open(path, 'rt') if path.suffix == '.gz' else path.open('r')


def find_file(directory, stem):
    for name in (stem + '.gz', stem):
        path = directory / name
        if path.is_file():
            return path
    raise FileNotFoundError(f'Missing {stem}(.gz) in {directory}')


def analyze(sample_id, directory, output_dir):
    matrix_path = find_file(directory, 'matrix.mtx')
    feature_path = find_file(directory, 'features.tsv')
    barcode_path = find_file(directory, 'barcodes.tsv')
    with open_text(feature_path) as handle:
        symbols = [line.rstrip('\n').split('\t')[1] for line in handle]
    with open_text(barcode_path) as handle:
        barcodes = [line.strip().split('\t')[0] for line in handle]
    matrix = mmread(str(matrix_path)).tocsr()
    if matrix.shape != (len(symbols), len(barcodes)):
        raise ValueError(f'{sample_id}: matrix shape does not match features/barcodes')
    nuclei = matrix.T.tocsr()
    umi = np.asarray(nuclei.sum(axis=1)).ravel().astype(np.int64)
    genes = np.diff(nuclei.indptr)
    mt_mask = np.array([name.startswith('MT-') for name in symbols])
    mt_umi = np.asarray(nuclei[:, mt_mask].sum(axis=1)).ravel()
    mt_percent = 100 * mt_umi / np.maximum(umi, 1)
    with (output_dir / f'{sample_id}_nucleus_metrics.csv').open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['sample_id', 'barcode', 'total_umi', 'detected_genes', 'mitochondrial_percent'])
        writer.writerows((sample_id, code, int(n), int(g), f'{p:.6f}')
                         for code, n, g, p in zip(barcodes, umi, genes, mt_percent))
    return dict(sample_id=sample_id, features=int(matrix.shape[0]),
                called_nuclei=int(matrix.shape[1]), nonzero_entries=int(matrix.nnz),
                total_umi=int(umi.sum()), median_umi=float(np.median(umi)),
                median_genes=float(np.median(genes)),
                median_mitochondrial_percent=float(np.median(mt_percent)),
                nuclei_with_fewer_than_200_genes=int((genes < 200).sum()),
                nuclei_above_20_percent_mitochondrial=int((mt_percent > 20).sum()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sample', action='append', required=True, metavar='ID=PATH',
                        help='Repeat for each sample; PATH contains matrix.mtx, features.tsv and barcodes.tsv (.gz accepted)')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    summary = []
    for spec in args.sample:
        if '=' not in spec:
            parser.error('--sample must have the form ID=PATH')
        sample_id, path = spec.split('=', 1)
        if not sample_id or not path:
            parser.error('--sample must have the form ID=PATH')
        summary.append(analyze(sample_id, Path(path), args.output))
    (args.output / 'matrix_qc_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()

# CCDC3 and neoadjuvant treated pancreatic cancer

This repository is being assembled for the CCDC3 PDAC manuscript. The input
libraries were generated from frozen tissue using single nucleus RNA sequencing.
The two manuscript samples are P1302539 and P1368090.

## Verified matrix audit

`scripts/01_matrix_qc.py` reads the filtered gene by nucleus matrices exported
by CeleScope 2.7.3. It checks matrix dimensions against the feature and barcode
files and calculates UMI counts, detected genes and mitochondrial transcript
percentages per nucleus. It does not repeat FASTQ processing, call doublets,
filter nuclei, integrate samples, annotate cell types or regenerate manuscript
figures. These steps require separate, validated code before the repository can
be described as the complete source code for the paper.

Place the three files for each sample in `data/P1302539` and `data/P1368090`:
`matrix.mtx.gz`, `features.tsv.gz`, `barcodes.tsv.gz`. Run from the repository root:

```bash
python -m pip install -r requirements.txt
python scripts/01_matrix_qc.py --sample P1302539=data/P1302539 \
  --sample P1368090=data/P1368090 --output qc_audit
```

The verified input dimensions are 38,606 features by 9,996 called nuclei for
P1302539 and 38,606 by 8,161 for P1368090 (18,157 nuclei total). These are
the matrices produced by the facility before the manuscript's downstream QC.
The 18,021 nucleus final atlas is not reproduced by this script.

No patient records, raw FASTQ files or restricted data are stored here.

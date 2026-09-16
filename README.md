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

## Descriptive Figure 5 audit

`scripts/02_figure5_descriptive.py` joins the original count matrices to the
previously saved 18,021 nucleus annotation file. It verifies that P1302539
contains 9,934 retained nuclei (62 fewer than its input matrix) and P1368090
contains 8,087 (74 fewer). It calculates cell type proportions and the mean
CPM and detection fractions of CCDC3, CCNE1 and PLK1 for each cell type. Its
output reproduces the reported malignant nucleus fractions of about 4.47% and
65.25%. Gene CPM estimates can differ slightly from values in an earlier figure
if its exact normalization and grouping rules differed.

The annotation file is an input to this script and is not stored in this
repository. After obtaining the data from the controlled access archive, run:

```bash
python scripts/02_figure5_descriptive.py \
  --annotation data/Figure5_PDAC_only_cell_metadata.csv.gz \
  --p1302539 data/P1302539 --p1368090 data/P1368090 \
  --output figure5_audit
```

No patient records, raw FASTQ files or restricted data are stored here.

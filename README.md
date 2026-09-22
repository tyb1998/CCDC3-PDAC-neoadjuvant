# CCDC3 in treated pancreatic ductal adenocarcinoma

Analysis code accompanying the manuscript **“A stromovascular CCDC3 state
marks tumour regression and low-CIN residual pancreatic cancer.”**

This repository provides documented workflows for the principal bulk
transcriptomic, spatial transcriptomic and in-house single-nucleus analyses
reported in the study. The scripts calculate sequencing-matrix quality
statistics, patient-level survival and correlation estimates, spatial summary
statistics, meta-analysis estimates, and the principal descriptive results of
the in-house cohort.

## Repository structure

- `scripts/01_matrix_qc.py`: per-sample single-nucleus matrix QC summaries.
- `scripts/02_figure5_descriptive.py`: in-house cell composition and
  pseudobulk CCDC3, CCNE1 and PLK1 summaries.
- `scripts/03_figure4_patient_correlations.py`: trial and spatial
  patient-level correlations.
- `scripts/04_figure1_trial_survival.py`: CAP association, permutation
  analysis, DFS Cox models and leave-one-out analysis.
- `scripts/05_figure2_spatial_summary.py`: spatial compartment and
  specificity summaries.
- `scripts/06_figure3_meta_analysis.py`: inverse-variance
  meta-analysis of cohort-specific Cox estimates.

## Environment

Python dependencies are listed in `requirements.txt`:

```bash
python -m pip install -r requirements.txt
```

## Input data

Public datasets used in the manuscript are available under GSE313101,
GSE202051 and GSE282302. In-house human sequencing data and associated
processed files are available through the controlled-access repository
specified in the manuscript Data Availability Statement.

For the in-house matrices, each sample directory contains
`matrix.mtx.gz`, `features.tsv.gz` and `barcodes.tsv.gz`. Figure source
tables use the filenames shown in the commands below. Large expression
objects, controlled-access data and patient-level source files are not
duplicated in this GitHub repository.

## Example commands

```bash
python scripts/01_matrix_qc.py \
  --sample P1302539=data/P1302539 \
  --sample P1368090=data/P1368090 \
  --output results/qc

python scripts/02_figure5_descriptive.py \
  --annotation data/Figure5_PDAC_only_cell_metadata.csv.gz \
  --p1302539 data/P1302539 \
  --p1368090 data/P1368090 \
  --output results/figure5

python scripts/03_figure4_patient_correlations.py \
  --trial data/Figure4_trial_residual_disease_scores.csv \
  --spatial data/Figure4_spatial_patient_expression.parquet \
  --output results/figure4

python scripts/04_figure1_trial_survival.py \
  --trial data/Figure1_trial_patient_response.csv \
  --output results/figure1

python scripts/05_figure2_spatial_summary.py \
  --compartments data/Figure2_spatial_compartment_source.csv \
  --specificity data/Figure2_spatial_specificity.csv \
  --output results/figure2

python scripts/06_figure3_meta_analysis.py \
  --cohort-cox data/Figure3_cohort_specific_cox.csv \
  --adjusted-meta data/Figure3_adjusted_meta.csv \
  --vascular-meta data/Figure3_vascular_meta.csv \
  --cin-meta data/Figure3_CIN_meta.csv \
  --output results/figure3
```

## Key analysis definitions

Gene-level pseudobulk CPM is calculated as the sum of gene UMIs divided by
the sum of all-gene UMIs within the specified sample or cell population,
multiplied by one million. Patient-level clinical associations use continuous
standardised CCDC3 expression without an outcome-derived cut point.
Cohort-specific Cox estimates are combined using inverse-variance
fixed-effect meta-analysis. Statistical analyses are conducted at the patient
level where multiple nuclei, tissue regions or sequencing blocks originate
from the same patient.

## Data protection

No human FASTQ files or directly identifiable clinical information are stored
in this repository. Users must obtain controlled-access data through the
repository and accession specified in the manuscript.

# CCDC3 in treated pancreatic ductal adenocarcinoma: reproducibility audit

This repository contains reconstructed analysis scripts for the manuscript. Scripts were checked against the original two in-house **single-nucleus** count matrices and retained Figure 1–5 source/summary tables. Original exploratory scripts, genomic-bin coordinates, doublet scores, the Hwang h5ad and six external cohort analysis objects have not been recovered. This is a reproducibility audit of specified results, not a claim that every figure can be regenerated from this repository.

## Inputs and commands

Create a local environment with `python -m pip install -r requirements.txt`. Obtain input files separately; do not put patient-level records or human FASTQ in this repository. The two input matrix directories contain `matrix.mtx.gz`, `features.tsv.gz` and `barcodes.tsv.gz` from CeleScope. The annotation CSV is a saved output from the original analysis, **not** independently recreated here. The trial CSV is a 21-patient previously processed expression/clinical table; the spatial parquet is previously processed expression for 39 patients.

```bash
python scripts/01_matrix_qc.py --sample P1302539=data/P1302539 --sample P1368090=data/P1368090 --output audit/qc
python scripts/02_figure5_descriptive.py --annotation data/Figure5_PDAC_only_cell_metadata.csv.gz --p1302539 data/P1302539 --p1368090 data/P1368090 --output audit/figure5
python scripts/03_figure4_patient_correlations.py --trial data/Figure4_trial_residual_disease_scores.csv --spatial data/Figure4_spatial_patient_expression.parquet --output audit/figure4
python scripts/04_figure1_trial_survival.py --trial data/Figure1_trial_patient_response.csv --output audit/figure1
python scripts/05_figure2_source_table_audit.py --compartments data/Figure2_spatial_compartment_source.csv --specificity data/Figure2_spatial_specificity.csv --output audit/figure2
python scripts/06_figure3_meta_summary_audit.py --cohort-cox data/Figure3_cohort_specific_cox.csv --adjusted-meta data/Figure3_adjusted_meta.csv --vascular-meta data/Figure3_vascular_meta.csv --cin-meta data/Figure3_CIN_meta.csv --output audit/figure3
```

## Checks performed

- `01`: input dimensions 38,606 genes × 9,996 nuclei (P1302539) and 38,606 × 8,161 (P1368090); per-nucleus QC from the *facility-filtered* matrices, not raw FASTQ.
- `02`: retained 9,934 and 8,087 nuclei using the saved annotation; malignant fractions 4.47% and 65.25%. **Pre-downstream-filtering** pseudobulk CCDC3 is 110.4248 vs 18.2627 CPM; **post-filtering, within-malignant** CCNE1 is 3.67896 vs 19.06836 CPM and PLK1 is 7.35792 vs 17.33195 CPM. Cell type pseudobulk CPM is summed gene UMIs divided by summed all-gene UMIs × 10^6, not the arithmetic average of per-cell CPM. The script outputs both measures explicitly; the original figure used pseudobulk. Cycling fractions derive from existing `manual_cell_state` annotation.
- `03`: trial unadjusted Spearman relationships and spatial CCDC3 correlations with CCNE1 (rho −0.45891, P 0.00330) and PKMYT1 (rho −0.54271, P 0.000358) match the manuscript. Our specified rank-residual partial Spearman correlations for some trial genes differ from the manuscript's reported partial coefficients; the original adjustment code was unavailable. Review these numbers before using them in the text.
- `04`: using the dedicated Figure 1 table, 21 trial patients and 15 DFS events; CAP rho −0.50301, asymptotic P 0.02011; DFS HR 0.568592 (95% CI 0.329123–0.982298; P 0.042971). A seeded, two-sided 100,000-permutation estimate is printed separately.
- `05`: validates the retained 39-patient/108-ROI spatial summary and manuscript display values. It does not recalculate spot-level partial correlations or Hwang atlas panels A–C.
- `06`: recombines the six unadjusted cohort estimates by inverse-variance fixed effects and exactly recovers 927 patients, 536 events and HR 0.834450. Adjusted, vascular and CIN tables are retained summaries and cannot substitute for refitting the cohort models.

## Current limits

These scripts do not regenerate the original 18-cluster embedding, Scrublet doublet calls, 177-gene classifier, inferCNV-like score, CIN score, Figure 2 cell atlas, Figure 2 spot-level analysis, Figure 3 cohort-specific Cox fits, or genome-wide candidate selection. The Figure 3 *fixed-effect pooling arithmetic* is reproduced, but the cohort fits feeding it are not. We cannot verify original filtering thresholds or classifier accuracy from the saved annotation alone. A patient-level or summary table does not independently reproduce its upstream processing from GEO.

Clinical conclusions drawn from the two in-house specimens are descriptive, since each response class contains one patient. The original sequencing data require deposition to an appropriate controlled-access archive; GitHub code alone does not satisfy the data-deposition request.

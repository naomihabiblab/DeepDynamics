# Synthetic tutorial cohort

**This directory contains a synthetic dataset.**

Donor identifiers and all numeric values are simulated. The files are provided
so the DeepDynamics tutorial notebook and method samples can be run without
study-restricted RNA-seq data. They are not a biological surrogate of the study
cohort and should not be interpreted as scientific results.

## Contents

| File | Description |
|------|-------------|
| `500.h5ad` | Synthetic AnnData object matching the tutorial schema |
| `shared_bulk_data_0.005.csv` | Filtered bulk features (`371 × 58`) |
| `shared_bulk_target_0.005.csv` | Targets `prAD`, `ABA`, `psuedotime` |
| `shared_bulk_data_mask.csv` | Same features under preprocess naming |
| `y.csv` | Same targets under preprocess naming |
| `features_names.csv` | Names of the 58 retained cell states |
| `X.npy` | Synthetic single-cell matrix after NaN filtering |

Shapes match the study-derived tutorial checks: 437 donors × 91 cell states in
the AnnData object, 419 shared bulk donors, 673 validation bulk donors, and 58
cell states retained by the correlation filter.

## Usage

The guided notebook and `examples/run_synthetic_pipeline.py` read these files
from `prediction/data/synthetic/`. Study-level data under `prediction/data/500.h5ad`
are not distributed here.

These files are committed for the tutorial and method samples; no regenerate
step is required to run the analyses.

## Analysis method samples

Additional CSVs (`donor_meta.csv`, `scrna_*.csv`, `proteomics_*.csv`,
`pathway_genes.csv`) support the R Markdown samples under
`../../analyses/notebooks/`. See `ANALYSIS_EXTRAS_README.md`.
These files are simulated and are not study results.

# DeepDynamics

DeepDynamics is a deep-learning framework that combines a small reference
single-cell or single-nucleus RNA-seq atlas with large bulk RNA-seq cohorts to
infer cell-type-specific and clinicopathological dynamics along annotated
trajectories.

The framework was developed to study processes such as disease progression and
aging, where bulk RNA-seq provides cohort-level statistical power but does not
directly retain cell-type resolution. DeepDynamics transfers trajectory
information learned from the single-cell reference to bulk cohorts and supports
downstream benchmarking, power analysis, scaling analysis, and interpretation.

In the accompanying study, DeepDynamics was applied to 1,092 cortical bulk
aging-brain profiles to investigate cellular, pathological, and molecular
mediators of Alzheimer's disease.

## Repository structure

```text
DeepDynamics/
├── DeepDynamics/
│   ├── DeepDynamics_example.ipynb   # End-to-end guided tutorial
│   ├── deepdynamics.yml             # Conda environment export
│   ├── prediction/                  # Model, loss, data structures, and training code
│   ├── explainability/              # Model-explainability notebook
│   ├── benchmarking/                # Baseline comparisons and DTW analyses
│   ├── power_analysis/              # Cohort and incremental-sampling power analyses
│   └── scaling_laws/                # Training-cohort-size scaling analysis
├── DeepDynamics_analyses/
│   ├── scripts/                     # Figure and validation analysis scripts
│   └── src/
│       ├── utils/                   # Shared R analysis utilities
│       └── visualization/           # Figure-generation functions
└── README.md
```

The main components are:

- `prediction/`: PyTorch implementation of DeepDynamics, including the model,
  loss function, dataset wrapper, preprocessing, and training utilities.
- `benchmarking/`: comparisons with linear regression, Elastic Net, XGBoost,
  and TabPFN, together with correlation, loss, visualization, and dynamic time
  warping utilities. Text files in this directory describe the associated
  methods.
- `power_analysis/`: analytic and simulation-based power estimates for group
  effects along the prAD and ABA trajectories.
- `scaling_laws/`: analysis of predictive performance as a function of bulk
  training-cohort size.
- `DeepDynamics_analyses/`: R scripts used for the paper-level downstream
  analyses and figure generation.

## Installation

### Prerequisites

- Git
- Conda or Mamba
- Linux is recommended because `deepdynamics.yml` was exported from a Linux
  environment.
- A CUDA-capable GPU is optional. CPU execution is supported but model fitting
  and repeated benchmarking will be slower.
- R is required only for the analyses under `DeepDynamics_analyses/` and for
  R-based trajectory-dynamics workflows.

### Clone the repository

```bash
git clone https://github.com/naomihabiblab/DeepDynamics.git
cd DeepDynamics
```

### Create the Python environment

Create an environment from the supplied Conda export. Passing `--name`
overrides the environment name stored in the exported YAML.

```bash
conda env create \
  --name deepdynamics \
  --file DeepDynamics/deepdynamics.yml
conda activate deepdynamics
```

The supplied environment includes Python, PyTorch, Scanpy, AnnData, NumPy,
Pandas, SciPy, scikit-learn, Matplotlib, Seaborn, Statsmodels, and the Jupyter
kernel components.

Install JupyterLab plus the additional benchmarking dependencies with:

```bash
python -m pip install jupyterlab optuna xgboost tabpfn
```

If GPU execution is required, install the PyTorch build matching the CUDA
version available on the system after creating the environment. See the
[PyTorch installation guide](https://pytorch.org/get-started/locally/) for the
appropriate command.

### Optional R environment

The R analysis utilities use CRAN and Bioconductor packages. A minimal setup can
be created from R as follows:

```r
install.packages(c(
  "anndata", "boot", "circlize", "colorspace", "cowplot", "dendsort",
  "dplyr", "ggnewscale", "ggplot2", "ggrepel", "gridExtra", "mgcv",
  "pheatmap", "progress", "purrr", "RColorBrewer", "reshape2",
  "reticulate", "stringr", "tidyr", "tidyverse"
))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c(
  "biomaRt", "ComplexHeatmap", "EnhancedVolcano", "SummarizedExperiment"
))
```

## Input data

Study-level input data are not distributed in this repository. Users should
provide their own authorized data or the synthetic tutorial cohort when it is
available. Do not commit participant-level or otherwise restricted data.

The guided notebook expects its inputs under `DeepDynamics/prediction/data/`.
The principal files used by the current workflows are:

```text
DeepDynamics/prediction/data/
├── 500.h5ad                       # Reference AnnData object
├── shared_bulk_data_mask.csv      # Bulk-derived input features
└── y.csv                          # Branch probabilities and pseudotime targets
```

The tutorial also refers to filtered synthetic tables named
`shared_bulk_data_0.005.csv` and `shared_bulk_target_0.005.csv`.

Power analysis uses a separate metadata table supplied with `--metadata`. The
required columns are documented in `DeepDynamics/power_analysis/config.py`.

## Getting started

Run commands from the repository root unless noted otherwise.

### Guided notebook

```bash
conda activate deepdynamics
jupyter lab DeepDynamics/DeepDynamics_example.ipynb
```

The notebook walks through loading the reference atlas and bulk-derived cell
state proportions, preparing trajectory targets, training DeepDynamics, and
inspecting predictions. The notebook is committed without saved outputs.

### Explainability notebook

```bash
jupyter lab DeepDynamics/explainability/Model_Explainability_analysis.ipynb
```

Update the input paths in the notebook for the local, authorized dataset before
running it.

### Power analysis

```bash
cd DeepDynamics
python -m power_analysis.power_analysis \
  --metadata /path/to/metadata.csv
```

Use `python -m power_analysis.power_analysis --help` to list the simulation,
bootstrap, significance-level, and sampling-fraction options.

### Scaling analysis

After placing the feature and target tables under `prediction/data/`:

```bash
cd DeepDynamics
python scaling_laws/scaling_laws.py
```

### Dynamic time warping

The DTW utility compares already-smoothed ground-truth and model dynamics:

```bash
cd DeepDynamics
python benchmarking/dtw_dynamics.py \
  --gt-pred-vals /path/to/ground_truth_dynamics_pred_vals.csv \
  --outputs-dir /path/to/model_dynamics_outputs \
  --out-csv /path/to/dtw_dynamics.csv
```

## Generated files

Analysis scripts may create `data/`, `outputs/`, `figures/`, model-weight, and
notebook-output artifacts. These files should remain local unless they are
explicitly approved for public release. In particular, notebooks should be
cleared before they are committed.

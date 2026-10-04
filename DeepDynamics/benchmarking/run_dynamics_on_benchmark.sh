#!/usr/bin/env bash
#
# Run the dynamics pipeline on every external_dynamics_input_*.csv and
# full_dynamics_input_*.csv in a given outputs directory, then generate dynamics plots.
# External runs: one dir per model (dynamics_$MODEL). Full runs: same dir, filenames with _full suffix.
#
# Usage (from project root DeepDynamics):
#   bash benchmarking/run_dynamics_on_benchmark.sh
#   bash benchmarking/run_dynamics_on_benchmark.sh [outputs_dir]
#
# Example:
#   bash benchmarking/run_dynamics_on_benchmark.sh benchmarking/outputs
#
# Requires: R with mgcv, dplyr, tidyr, progress; Python with pandas, matplotlib.
# If Rscript/python are not in PATH, the script will try to activate conda env
#   CONDA_ENV (default: dd_r). Set CONDA_ENV= to skip activation.
# For plotting, if your R env (dd_r) lacks pandas/matplotlib, set CONDA_ENV_PY
#   to an env that has them (e.g. dd): CONDA_ENV_PY=dd bash run_dynamics_on_benchmark.sh
#
# Dynamics are run only for these features (each model's prAD, ABA, pseudotime are used):
#   FEATURES (default: sqrt.amyloid_mf,sqrt.tangles_mf,cogng_demog_slope).
#   Missing features are merged by ID from the reference RDS (yifat_stuff/for_yuval.rds).
# In addition, a per-model 3x2 pathology grid PDF is written as
#   outputs/dynamics_$model/pathologies_grid_${model}.pdf (external)
#   outputs/dynamics_$model/pathologies_grid_${model}_full.pdf (full)

set -euo pipefail

# Pathology features to fit/plot: comma-separated; merged from reference RDS if missing in CSV
FEATURES="${FEATURES:-sqrt.amyloid_mf,sqrt.tangles_mf,cogng_demog_slope}"

# Cell cluster features to ALSO fit/plot for the clusters grid (SIG_CLUSTERS order).
# Set PLOT_CLUSTERS_GRID=0 to disable.
PLOT_CLUSTERS_GRID="${PLOT_CLUSTERS_GRID:-1}"
# CLUSTERS_FEATURES="${CLUSTERS_FEATURES:-Ast.1,Ast.10,Ast.2,Ast.4,Ast.5,Ast.6,Ast.7,End.1,End.3,Exc.1,Exc.12,Exc.3,Exc.8,Inh.12,Inh.15,Inh.16,Inh.5,Inh.6,Inh.7,Mic.1,Mic.7,Mic.12,Mic.13,Mic.14,Oli.3,Oli.4,Oli.5,Oli.7,Oli.9,Oli.11,OPC.1,OPC.2}"
CLUSTERS_FEATURES="${CLUSTERS_FEATURES:-Ast.1,Ast.10,Ast.5}"

# Activate conda env if Rscript or python not found (e.g. on cluster login node)
CONDA_ENV="${CONDA_ENV:-dd_r}"
if ! command -v Rscript &>/dev/null || ! command -v python &>/dev/null; then
  for conda_sh in "$HOME/miniconda3/etc/profile.d/conda.sh" "$HOME/anaconda3/etc/profile.d/conda.sh" "$CONDA_PREFIX/../etc/profile.d/conda.sh"; do
    if [[ -f "${conda_sh:-}" ]]; then
      source "$conda_sh"
      conda activate "$CONDA_ENV" 2>/dev/null || true
      break
    fi
  done
fi
if ! command -v Rscript &>/dev/null; then
  echo "Error: Rscript not found. Activate your R env first, e.g. conda activate dd_r" >&2
  exit 1
fi
if ! command -v python &>/dev/null; then
  echo "Error: python not found. Activate your env first (e.g. conda activate dd_r)" >&2
  exit 1
fi

# Project root: directory containing yifat_stuff and benchmarking
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
OUTPUTS_DIR="${1:-$SCRIPT_DIR/outputs}"

# Resolve to absolute path so R/Python see consistent paths
OUTPUTS_DIR="$(cd "$OUTPUTS_DIR" && pwd)"
YIFAT_DIR="$PROJECT_ROOT/yifat_stuff"
RUN_R="$YIFAT_DIR/run_dynamics.R"
PLOT_PY="$YIFAT_DIR/plot_dynamics.py"
PLOT_GRID_R="$YIFAT_DIR/plot_pathologies_grid.R"
PLOT_CLUSTERS_R="$YIFAT_DIR/plot_clusters_grid.R"

# Reference object for merging missing features by ID (pathologies + clusters).
# Default points at the CSV exported from bulk.prev.meta.RData (bulk.df).
REFERENCE_OBJ="${REFERENCE_OBJ:-$YIFAT_DIR/bulk.prev.meta.bulk.df.csv}"

if [[ ! -f "$RUN_R" ]]; then
  echo "Error: run_dynamics.R not found at $RUN_R" >&2
  exit 1
fi
if [[ ! -f "$PLOT_PY" ]]; then
  echo "Error: plot_dynamics.py not found at $PLOT_PY" >&2
  exit 1
fi
if [[ ! -f "$PLOT_GRID_R" ]]; then
  echo "Error: plot_pathologies_grid.R not found at $PLOT_GRID_R" >&2
  exit 1
fi
if [[ "$PLOT_CLUSTERS_GRID" != "0" ]] && [[ ! -f "$PLOT_CLUSTERS_R" ]]; then
  echo "Error: plot_clusters_grid.R not found at $PLOT_CLUSTERS_R" >&2
  exit 1
fi
if [[ ! -d "$OUTPUTS_DIR" ]]; then
  echo "Error: outputs directory not found: $OUTPUTS_DIR" >&2
  exit 1
fi
if [[ ! -f "$REFERENCE_OBJ" ]]; then
  echo "Error: reference object not found: $REFERENCE_OBJ" >&2
  exit 1
fi

# Python for plotting: use CONDA_ENV_PY if set (env with pandas/matplotlib), else current python
if [[ -n "${CONDA_ENV_PY:-}" ]] && command -v conda &>/dev/null; then
  PYTHON_CMD="conda run -n $CONDA_ENV_PY --no-capture-output python"
else
  PYTHON_CMD="python"
fi

echo "Outputs directory: $OUTPUTS_DIR"
echo "Running dynamics pipeline for external_dynamics_input_*.csv and full_dynamics_input_*.csv"
if [[ "$PLOT_CLUSTERS_GRID" != "0" ]]; then
  echo "Clusters grid: enabled (will also fit CLUSTERS_FEATURES; this can take longer)"
else
  echo "Clusters grid: disabled"
fi
echo "---"

shopt -s nullglob

# 1) External dynamics (validation-only samples)
for csv in "$OUTPUTS_DIR"/external_dynamics_input_*.csv; do
  base="$(basename "$csv" .csv)"
  model="${base#external_dynamics_input_}"
  out_subdir="$OUTPUTS_DIR/dynamics_$model"
  echo "Input: $csv"
  echo "Model: $model -> output dir: $out_subdir"

  mkdir -p "$out_subdir"
  if [[ "$PLOT_CLUSTERS_GRID" != "0" ]]; then
    ALL_FEATURES="${FEATURES},${CLUSTERS_FEATURES}"
  else
    ALL_FEATURES="${FEATURES}"
  fi
  Rscript "$RUN_R" "$csv" "$out_subdir" "$model" "$ALL_FEATURES" "$REFERENCE_OBJ" \
    || { echo "R failed for $model" >&2; exit 1; }
  $PYTHON_CMD "$PLOT_PY" "$out_subdir" "$model" || { echo "Python plot failed for $model" >&2; exit 1; }
  pred_csv="$out_subdir/dynamics_pred_vals_${model}.csv"
  grid_pdf="$out_subdir/pathologies_grid_${model}.pdf"
  clusters_pdf="$out_subdir/clusters_grid_${model}.pdf"
  if [[ -f "$pred_csv" ]]; then
    Rscript "$PLOT_GRID_R" "$pred_csv" "$grid_pdf" "$model (external)" \
      || { echo "Pathology grid plot failed for $model (external)" >&2; exit 1; }
    if [[ "$PLOT_CLUSTERS_GRID" != "0" ]]; then
      Rscript "$PLOT_CLUSTERS_R" "$pred_csv" "$clusters_pdf" "$model (external)" \
        || { echo "Clusters grid plot failed for $model (external)" >&2; exit 1; }
    fi
  else
    echo "Warning: expected $pred_csv not found; skipping pathology grid plot for $model (external)" >&2
  fi
  echo "Done: $model (external)"
  echo "---"
done

# 2) Full dynamics (train + test + external); same output dir per model, with _full suffix in filenames
for csv in "$OUTPUTS_DIR"/full_dynamics_input_*.csv; do
  base="$(basename "$csv" .csv)"
  model="${base#full_dynamics_input_}"
  out_subdir="$OUTPUTS_DIR/dynamics_$model"
  model_suffix="${model}_full"
  echo "Input: $csv"
  echo "Model: $model (full) -> output dir: $out_subdir (suffix: $model_suffix)"

  mkdir -p "$out_subdir"
  if [[ "$PLOT_CLUSTERS_GRID" != "0" ]]; then
    ALL_FEATURES="${FEATURES},${CLUSTERS_FEATURES}"
  else
    ALL_FEATURES="${FEATURES}"
  fi
  Rscript "$RUN_R" "$csv" "$out_subdir" "$model_suffix" "$ALL_FEATURES" "$REFERENCE_OBJ" \
    || { echo "R failed for $model (full)" >&2; exit 1; }
  $PYTHON_CMD "$PLOT_PY" "$out_subdir" "$model_suffix" || { echo "Python plot failed for $model (full)" >&2; exit 1; }
  pred_csv="$out_subdir/dynamics_pred_vals_${model_suffix}.csv"
  grid_pdf="$out_subdir/pathologies_grid_${model_suffix}.pdf"
  clusters_pdf="$out_subdir/clusters_grid_${model_suffix}.pdf"
  if [[ -f "$pred_csv" ]]; then
    Rscript "$PLOT_GRID_R" "$pred_csv" "$grid_pdf" "$model (full)" \
      || { echo "Pathology grid plot failed for $model (full)" >&2; exit 1; }
    if [[ "$PLOT_CLUSTERS_GRID" != "0" ]]; then
      Rscript "$PLOT_CLUSTERS_R" "$pred_csv" "$clusters_pdf" "$model (full)" \
        || { echo "Clusters grid plot failed for $model (full)" >&2; exit 1; }
    fi
  else
    echo "Warning: expected $pred_csv not found; skipping pathology grid plot for $model (full)" >&2
  fi
  echo "Done: $model (full)"
  echo "---"
done

echo "All dynamics runs and plots finished."

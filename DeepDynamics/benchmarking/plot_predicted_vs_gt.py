#!/usr/bin/env python
"""
Scatter plots of predicted pseudotime vs ground-truth pseudotime for each model,
colored by sqrt.tangles_mf, sqrt.amyloid_mf, braaksc, and cogng_demog_slope from the metadata.

Ground-truth pseudotime comes from prediction/data/y.csv (the target file used
during benchmarking). Pathology coloring comes from metadata.csv.

metadata.csv's own prAD, ABA, pseudotime are included as a "metadata" model
entry to serve as a reference (metadata pseudotime vs y.csv pseudotime).

Three outputs:
  1) *_our_labels.pdf + correlations_our_labels.csv
     Trajectory split uses metadata/GT prAD (consistent across all models).
  2) *.pdf + correlations.csv
     Trajectory split uses each model's own predicted prAD.
  3) correlations_gt_subset.csv
     GT pseudotime vs pathology + model pred pseudotime vs pathology,
     restricted to the ~371 GT subjects only.
"""

import argparse
import os
import glob

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import spearmanr, pearsonr

from viz import save_figure_for_illustrator

PROBABILITY_THRESHOLD = 0.5
PSEUDOTIME_THRESHOLD = 0.1

# ── paths ────────────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(SCRIPT_DIR)
OUTPUTS_DIR = os.path.join(SCRIPT_DIR, "outputs")
FIGURES_DIR = os.path.join(SCRIPT_DIR, "figures")
METADATA_PATH = os.path.join(BASE_DIR, "metadata.csv")
GT_TARGET_PATH = os.path.join(BASE_DIR, "prediction", "data", "y.csv")

os.makedirs(FIGURES_DIR, exist_ok=True)

# ── load metadata (for pathology coloring + as a reference "model") ───────────
meta_full = pd.read_csv(METADATA_PATH)
meta = meta_full[["ID", "sqrt.tangles_mf", "sqrt.amyloid_mf", "braaksc", "cogng_demog_slope"]].copy()
meta_with_traj = meta_full[["ID", "pseudotime", "prAD", "ABA"]].copy()
meta_with_traj = meta_with_traj.dropna(subset=["pseudotime"])

# ── load ground-truth targets (y.csv used in benchmarking) ───────────────────
gt_y = pd.read_csv(GT_TARGET_PATH, index_col=0)
# y.csv columns: prAD, ABA, psuedotime (note the typo in source)
gt_y = gt_y.rename(columns={"psuedotime": "pseudotime"})
gt_y["ID"] = gt_y.index.astype(int)
gt_y = gt_y.reset_index(drop=True)

# ── discover models + add metadata as a reference model ──────────────────────
pattern = os.path.join(OUTPUTS_DIR, "full_dynamics_input_*.csv")
model_files = sorted(glob.glob(pattern))

# Store each model as a DataFrame with columns: ID, pseudotime, prAD, ABA
model_data = {}
for fpath in model_files:
    fname = os.path.basename(fpath)
    model_name = fname.replace("full_dynamics_input_", "").replace(".csv", "")
    model_data[model_name] = pd.read_csv(fpath)

# Add metadata as a reference "model" using its own prAD, ABA, pseudotime
model_data["metadata"] = meta_with_traj[["ID", "pseudotime", "prAD", "ABA"]].copy()

print(f"Models (including metadata): {list(model_data.keys())}")

# ── Prepare merge table: GT pseudotime + GT prAD from y.csv + pathology ──────
meta_for_merge = gt_y[["ID", "pseudotime", "prAD"]].merge(meta, on="ID", how="left")
meta_for_merge = meta_for_merge.rename(columns={
    "pseudotime": "gt_pseudotime",
    "prAD": "gt_prAD",
})
meta_for_merge["gt_trajectory"] = np.where(
    meta_for_merge["gt_prAD"] > PROBABILITY_THRESHOLD, "prAD", "ABA"
)
meta_for_merge = meta_for_merge.merge(
    meta_full[["ID", "sex", "apoe_genotype"]], on="ID", how="left"
)

# ── Trajectory map from metadata prAD (pure metadata labels, no GT override) ──
meta_traj_all = meta_full[["ID", "prAD"]].copy().rename(columns={"prAD": "meta_prAD"})
meta_traj_all["trajectory"] = np.where(
    meta_traj_all["meta_prAD"] > PROBABILITY_THRESHOLD, "prAD", "ABA"
)
traj_map_all = meta_traj_all.set_index("ID")["trajectory"]

# ── constants ─────────────────────────────────────────────────────────────────
CMAP = "viridis"
PATHOLOGY_VARS = ["sqrt.amyloid_mf", "sqrt.tangles_mf", "braaksc", "cogng_demog_slope"]
TRAJECTORIES = ["prAD", "ABA"]


def _parse_args():
    parser = argparse.ArgumentParser(
        description="Predicted vs GT and pathology plots; optional pred-vs-GT and pred-vs-residuals with color-by and regression line.",
    )
    parser.add_argument(
        "--pred-vs-gt",
        action="store_true",
        help="Generate predicted vs ground-truth pseudotime plots (one per model).",
    )
    parser.add_argument(
        "--pred-vs-residuals",
        action="store_true",
        help="Generate predicted vs residuals (GT − predicted) plots (one per model).",
    )
    parser.add_argument(
        "--color-by",
        choices=["sex", "apoe"],
        default=None,
        help="Color points by this categorical variable (e.g. sex, apoe) in pred-vs-GT and pred-vs-residuals plots.",
    )
    parser.add_argument(
        "--regression-line",
        action="store_true",
        help="Add a regression line to pred-vs-GT and pred-vs-residuals plots.",
    )
    return parser.parse_args()


# ── helpers ───────────────────────────────────────────────────────────────────
def compute_corr(x, y):
    """Return (pearson_r, spearman_rho, n) or (NaN, NaN, 0) if insufficient data."""
    mask = np.isfinite(x) & np.isfinite(y)
    n = int(mask.sum())
    if n < 3:
        return np.nan, np.nan, n
    r_p, _ = pearsonr(x[mask], y[mask])
    r_s, _ = spearmanr(x[mask], y[mask])
    return r_p, r_s, n


def annotate_corr(ax, x, y):
    """Add Pearson & Spearman r to top-left corner."""
    r_p, r_s, n = compute_corr(x, y)
    if n < 3:
        return
    ax.text(
        0.05, 0.95,
        f"Pearson r = {r_p:.3f}\nSpearman rho = {r_s:.3f}",
        transform=ax.transAxes,
        fontsize=8,
        verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
    )


def add_regression_line(ax, x, y, color="red", linestyle="-", linewidth=1.5):
    """Fit y ~ x and plot the regression line over the x range."""
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 2:
        return
    x_vals = x[mask]
    y_vals = y[mask]
    slope, intercept, r, p, se = stats.linregress(x_vals, y_vals)
    x_line = np.array([x_vals.min(), x_vals.max()])
    y_line = slope * x_line + intercept
    ax.plot(x_line, y_line, color=color, linestyle=linestyle, linewidth=linewidth, label=f"Fit (r={r:.3f})")


def load_and_merge(model_name):
    """Load a model's predictions and merge with GT pseudotime + metadata."""
    pred = model_data[model_name].copy()
    pred = pred.rename(columns={"pseudotime": "pred_pseudotime"})
    merged = pred.merge(meta_for_merge, on="ID", how="inner")
    return merged


def load_and_merge_full(model_name):
    """Load a model's predictions and merge with metadata pathology (all subjects)."""
    pred = model_data[model_name].copy()
    pred = pred.rename(columns={"pseudotime": "pred_pseudotime"})
    merged = pred.merge(meta, on="ID", how="inner")
    return merged


# Categorical column name for color-by: CLI "apoe" -> apoe_genotype
COLOR_BY_COLUMN = {"sex": "sex", "apoe": "apoe_genotype"}


def _scatter_pred_gt_or_residuals(
    ax, merged, plot_type, color_by=None, add_regression_line_=False
):
    """
    Draw scatter on ax: either pred vs GT or pred vs residuals.
    plot_type in ("pred_vs_gt", "pred_vs_residuals").
    """
    x = merged["pred_pseudotime"].values
    if plot_type == "pred_vs_gt":
        y = merged["gt_pseudotime"].values
        ax.set_xlabel("Predicted pseudotime", fontsize=11)
        ax.set_ylabel("Ground-truth pseudotime", fontsize=11)
    else:
        y = merged["gt_pseudotime"].values - merged["pred_pseudotime"].values
        ax.set_xlabel("Predicted pseudotime", fontsize=11)
        ax.set_ylabel("Residual (GT − predicted)", fontsize=11)

    valid = np.isfinite(x) & np.isfinite(y)
    if not valid.any():
        return

    x_v, y_v = x[valid], y[valid]

    if color_by and color_by in merged.columns:
        col_vals = merged.loc[valid, color_by].astype(str)
        cats = col_vals.unique()
        colors = plt.cm.tab10(np.linspace(0, 1, max(len(cats), 1)))
        for i, cat in enumerate(cats):
            mask = col_vals == cat
            ax.scatter(
                x_v[mask], y_v[mask],
                label=cat, alpha=0.7, s=16, c=[colors[i % len(colors)]], edgecolors="none",
            )
        ax.legend(loc="best", fontsize=8)
    else:
        ax.scatter(x_v, y_v, alpha=0.6, s=16, c="steelblue", edgecolors="none")

    if plot_type == "pred_vs_gt":
        lims = [np.nanmin(np.r_[x_v, y_v]), np.nanmax(np.r_[x_v, y_v])]
        ax.plot(lims, lims, "--", color="grey", linewidth=0.8, alpha=0.6, label="y=x")
    else:
        ax.axhline(0, color="grey", linestyle="--", linewidth=0.8, alpha=0.6)

    if add_regression_line_:
        add_regression_line(ax, x_v, y_v)

    r_p, r_s, n = compute_corr(x_v, y_v)
    if n >= 3:
        ax.text(
            0.05, 0.95,
            f"Pearson r = {r_p:.3f}\nSpearman rho = {r_s:.3f}\nn = {n}",
            transform=ax.transAxes,
            fontsize=8,
            verticalalignment="top",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
        )


def generate_pred_vs_gt_plots(color_by=None, add_regression_line_=False):
    """One panel per model: predicted vs ground-truth pseudotime."""
    color_col = (COLOR_BY_COLUMN.get(color_by, color_by) if color_by else None)
    tag = []
    if color_col:
        tag.append(f"color_{color_by}")
    if add_regression_line_:
        tag.append("regression")
    suffix = "_" + "_".join(tag) if tag else ""

    for model_name in model_data:
        merged = load_and_merge(model_name)
        fig, ax = plt.subplots(1, 1, figsize=(6, 5), constrained_layout=True)
        _scatter_pred_gt_or_residuals(
            ax, merged, "pred_vs_gt",
            color_by=color_col,
            add_regression_line_=add_regression_line_,
        )
        ax.set_title(f"{model_name}: Predicted vs GT pseudotime", fontsize=12)
        safe_name = model_name.replace(" ", "_")
        out_path = os.path.join(
            FIGURES_DIR, f"pred_vs_gt_{safe_name}{suffix}.pdf",
        )
        save_figure_for_illustrator(fig, out_path)
        plt.close(fig)
        print(f"  {out_path}")


def generate_pred_vs_residuals_plots(color_by=None, add_regression_line_=False):
    """One panel per model: predicted pseudotime vs residuals (GT − predicted)."""
    color_col = (COLOR_BY_COLUMN.get(color_by, color_by) if color_by else None)
    tag = []
    if color_col:
        tag.append(f"color_{color_by}")
    if add_regression_line_:
        tag.append("regression")
    suffix = "_" + "_".join(tag) if tag else ""

    for model_name in model_data:
        merged = load_and_merge(model_name)
        fig, ax = plt.subplots(1, 1, figsize=(6, 5), constrained_layout=True)
        _scatter_pred_gt_or_residuals(
            ax, merged, "pred_vs_residuals",
            color_by=color_col,
            add_regression_line_=add_regression_line_,
        )
        ax.set_title(f"{model_name}: Predicted vs residuals", fontsize=12)
        safe_name = model_name.replace(" ", "_")
        out_path = os.path.join(
            FIGURES_DIR, f"pred_vs_residuals_{safe_name}{suffix}.pdf",
        )
        save_figure_for_illustrator(fig, out_path)
        plt.close(fig)
        print(f"  {out_path}")


def generate_plots_and_csv(suffix, use_our_labels):
    """
    Generate trajectory-split pathology plots and a correlations CSV.

    Parameters
    ----------
    suffix : str
        File name suffix, e.g. "" or "_our_labels".
    use_our_labels : bool
        If True, use metadata/GT-based trajectory assignment for all models.
        If False, use each model's own predicted prAD for trajectory assignment.
    """
    corr_rows = []

    for model_name in model_data:
        merged = load_and_merge(model_name)
        merged_full = load_and_merge_full(model_name)

        # ── Assign trajectory to the full set ──
        if use_our_labels:
            # Metadata / GT based (same for all models)
            merged_full["trajectory"] = merged_full["ID"].map(traj_map_all)
        else:
            # Model's own predicted prAD
            merged_full["trajectory"] = np.where(
                merged_full["prAD"] > PROBABILITY_THRESHOLD, "prAD", "ABA"
            )

        n_rows = len(PATHOLOGY_VARS)
        n_cols = len(TRAJECTORIES)

        fig, axes = plt.subplots(
            n_rows, n_cols,
            figsize=(6 * n_cols, 5 * n_rows),
            constrained_layout=True,
            squeeze=False,
        )

        for row_idx, cvar in enumerate(PATHOLOGY_VARS):
            all_c = merged[cvar].values.astype(float)
            valid_all = np.isfinite(all_c)
            if valid_all.any():
                vmin, vmax = np.nanmin(all_c), np.nanmax(all_c)
            else:
                vmin, vmax = 0, 1

            for col_idx, traj in enumerate(TRAJECTORIES):
                ax = axes[row_idx, col_idx]

                # Plot subset: only GT subjects, trajectory from GT
                subset = merged[merged["gt_trajectory"] == traj].copy()

                c_vals = subset[cvar].values.astype(float)
                valid = np.isfinite(c_vals)

                if (~valid).any():
                    ax.scatter(
                        subset.loc[~valid, "gt_pseudotime"],
                        subset.loc[~valid, "pred_pseudotime"],
                        s=14, alpha=0.35, color="lightgrey", edgecolors="none",
                        label="NA",
                    )

                if valid.any():
                    sc = ax.scatter(
                        subset.loc[valid, "gt_pseudotime"],
                        subset.loc[valid, "pred_pseudotime"],
                        c=c_vals[valid],
                        cmap=CMAP,
                        vmin=vmin, vmax=vmax,
                        s=16, alpha=0.7, edgecolors="none",
                    )
                    plt.colorbar(sc, ax=ax, label=cvar, shrink=0.85)

                all_vals = np.concatenate([
                    subset["gt_pseudotime"].values,
                    subset["pred_pseudotime"].values,
                ])
                if len(all_vals) > 0:
                    lims = [np.nanmin(all_vals), np.nanmax(all_vals)]
                    ax.plot(lims, lims, "--", color="grey", linewidth=0.8, alpha=0.6)

                # Correlation on FULL set of subjects
                # Require both trajectory label and minimum predicted pseudotime
                full_subset = merged_full[
                    (merged_full["trajectory"] == traj)
                    & (merged_full["pred_pseudotime"] > PSEUDOTIME_THRESHOLD)
                ]
                full_c = full_subset[cvar].values.astype(float)
                full_valid = np.isfinite(full_c)
                if full_valid.any():
                    x_corr = full_subset.loc[full_valid, "pred_pseudotime"].values
                    y_corr = full_c[full_valid]
                    annotate_corr(ax, x_corr, y_corr)
                    r_p, r_s, n = compute_corr(x_corr, y_corr)
                    corr_rows.append({
                        "model": model_name,
                        "trajectory": traj,
                        "pathology": cvar,
                        "pearson_r": r_p,
                        "spearman_rho": r_s,
                        "n": n,
                    })

                n_plot = len(subset)
                n_corr = int(full_valid.sum())
                ax.set_xlabel("Ground-truth pseudotime", fontsize=11)
                ax.set_ylabel("Predicted pseudotime", fontsize=11)
                ax.set_title(
                    f"{traj} trajectory (plot n={n_plot}, corr n={n_corr}) - {cvar}",
                    fontsize=10,
                )

        label_tag = " [our labels]" if use_our_labels else " [model labels]"
        fig.suptitle(
            f"{model_name}: Predicted vs GT Pseudotime by Trajectory{label_tag}",
            fontsize=13, fontweight="bold",
        )

        safe_name = model_name.replace(" ", "_")
        out_path = os.path.join(
            FIGURES_DIR, f"pathology_by_trajectory_{safe_name}{suffix}.pdf",
        )
        save_figure_for_illustrator(fig, out_path)
        plt.close(fig)
        print(f"  {out_path}")

    # Save correlations CSV
    corr_df = pd.DataFrame(corr_rows)
    corr_csv_path = os.path.join(FIGURES_DIR, f"correlations{suffix}.csv")
    corr_df.to_csv(corr_csv_path, index=False, float_format="%.6f")
    print(f"  CSV -> {corr_csv_path}")
    return corr_df


# ══════════════════════════════════════════════════════════════════════════════
# 1) Our labels (metadata / GT prAD trajectory assignment)
# ══════════════════════════════════════════════════════════════════════════════
# Commented out per request: disable generation of `*_our_labels` figures/CSVs.
# print("\n=== Generating plots with OUR labels (metadata/GT prAD) ===")
# corr_our = generate_plots_and_csv(suffix="_our_labels", use_our_labels=True)

# ══════════════════════════════════════════════════════════════════════════════
# 2) Model labels (each model's own predicted prAD)
# ══════════════════════════════════════════════════════════════════════════════
print("\n=== Generating plots with MODEL labels (model predicted prAD) ===")
corr_model = generate_plots_and_csv(suffix="", use_our_labels=False)

# ══════════════════════════════════════════════════════════════════════════════
# 3) GT-subset CSV: ground-truth + model correlations on GT subjects only (~371)
# ══════════════════════════════════════════════════════════════════════════════
print("\n=== Computing GT-subset correlations (only ~371 GT subjects) ===")
gt_subset_rows = []

# Ground-truth pseudotime vs pathology (GT subjects only, GT trajectory)
gt_merged = gt_y[["ID", "pseudotime"]].merge(meta, on="ID", how="inner")
gt_merged["gt_trajectory"] = gt_merged["ID"].map(traj_map_all)

for cvar in PATHOLOGY_VARS:
    for traj in TRAJECTORIES:
        sub = gt_merged[gt_merged["gt_trajectory"] == traj]
        c_vals = sub[cvar].values.astype(float)
        valid = np.isfinite(c_vals)
        if valid.any():
            r_p, r_s, n = compute_corr(
                sub.loc[valid, "pseudotime"].values, c_vals[valid],
            )
            gt_subset_rows.append({
                "model": "ground_truth",
                "trajectory": traj,
                "pathology": cvar,
                "pearson_r": r_p,
                "spearman_rho": r_s,
                "n": n,
            })
        else:
            print(f"No valid data for {cvar} and {traj}")

# Each model's predicted pseudotime vs pathology (GT subjects only, GT trajectory)
for model_name in model_data:
    merged = load_and_merge(model_name)  # already restricted to GT subjects
    for cvar in PATHOLOGY_VARS:
        for traj in TRAJECTORIES:
            # Restrict to this trajectory and minimum predicted pseudotime
            sub = merged[
                (merged["gt_trajectory"] == traj)
                & (merged["pred_pseudotime"] > PSEUDOTIME_THRESHOLD)
            ]
            c_vals = sub[cvar].values.astype(float)
            valid = np.isfinite(c_vals)
            if valid.any():
                r_p, r_s, n = compute_corr(
                    sub.loc[valid, "pred_pseudotime"].values, c_vals[valid],
                )
                gt_subset_rows.append({
                    "model": model_name,
                    "trajectory": traj,
                    "pathology": cvar,
                    "pearson_r": r_p,
                    "spearman_rho": r_s,
                    "n": n,
                })

gt_subset_df = pd.DataFrame(gt_subset_rows)
gt_csv_path = os.path.join(FIGURES_DIR, "correlations_gt_subset.csv")
gt_subset_df.to_csv(gt_csv_path, index=False, float_format="%.6f")
print(f"  CSV -> {gt_csv_path}")
print(gt_subset_df.to_string(index=False))

# ══════════════════════════════════════════════════════════════════════════════
# 4) Optional: Predicted vs GT and Predicted vs residuals (with --pred-vs-gt / --pred-vs-residuals)
# ══════════════════════════════════════════════════════════════════════════════
args = _parse_args()
if args.pred_vs_gt:
    print("\n=== Predicted vs GT pseudotime (optional color-by and regression line) ===")
    generate_pred_vs_gt_plots(color_by=args.color_by, add_regression_line_=args.regression_line)
if args.pred_vs_residuals:
    print("\n=== Predicted vs residuals (optional color-by and regression line) ===")
    generate_pred_vs_residuals_plots(color_by=args.color_by, add_regression_line_=args.regression_line)

print("\nDone.")

#!/usr/bin/env python
"""
Violin plots of pathology values across trajectory groups from predicted dynamics.

For probmodel, tabpfn, and xgboost predictions, plus the ground-truth targets
used by plot_predicted_vs_gt.py and the labels stored in metadata.csv, subjects
are split into: healthy, prAD, and ABA. One grouped violin plot is saved per
pathology.
"""

import glob
import os
from typing import Dict, Iterable, List

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from viz import save_figure_for_illustrator


PROBABILITY_THRESHOLD = 0.5
PSEUDOTIME_THRESHOLD = 0.1

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(SCRIPT_DIR)
OUTPUTS_DIR = os.path.join(SCRIPT_DIR, "outputs")
FIGURES_DIR = os.path.join(SCRIPT_DIR, "figures", "violin_plots")
METADATA_PATH = os.path.join(BASE_DIR, "metadata.csv")
GT_TARGET_PATH = os.path.join(BASE_DIR, "prediction", "data", "y.csv")

REQUIRED_DYNAMICS_COLUMNS = ["ID", "pseudotime", "prAD", "ABA"]
MODEL_ORDER = ["probmodel", "tabpfn", "xgboost"]
PATHOLOGY_COLUMNS = ["sqrt.amyloid", "sqrt.tangles", "sqrt.amyloid_mf", "sqrt.tangles_mf", "cogng_demog_slope"]
GROUP_ORDER = ["healthy", "prAD", "ABA"]
GROUP_COLORS = {
    "healthy": "#2ecc71",
    "prAD": "#e74c3c",
    "ABA": "#3498db",
}


def _require_columns(df: pd.DataFrame, columns: Iterable[str], source: str) -> None:
    missing = [col for col in columns if col not in df.columns]
    if missing:
        raise ValueError(f"{source} is missing required columns: {missing}")


def _coerce_id_column(df: pd.DataFrame, source: str) -> pd.DataFrame:
    df = df.copy()
    df["ID"] = pd.to_numeric(df["ID"], errors="coerce")
    df = df.dropna(subset=["ID"])
    df["ID"] = df["ID"].astype(int)
    if df.empty:
        raise ValueError(f"{source} has no valid IDs after coercion.")
    return df


def load_metadata() -> pd.DataFrame:
    metadata = pd.read_csv(METADATA_PATH)
    _require_columns(metadata, ["ID"] + PATHOLOGY_COLUMNS, METADATA_PATH)
    metadata = _coerce_id_column(metadata, METADATA_PATH)
    return metadata[["ID"] + PATHOLOGY_COLUMNS].copy()


def load_ground_truth() -> pd.DataFrame:
    gt = pd.read_csv(GT_TARGET_PATH, index_col=0)
    gt = gt.rename(columns={"psuedotime": "pseudotime"})
    _require_columns(gt, ["pseudotime", "prAD", "ABA"], GT_TARGET_PATH)
    gt["ID"] = gt.index
    gt = gt.reset_index(drop=True)
    gt = _coerce_id_column(gt, GT_TARGET_PATH)
    return gt[REQUIRED_DYNAMICS_COLUMNS].copy()


def load_metadata_labels() -> pd.DataFrame:
    metadata = pd.read_csv(METADATA_PATH)
    _require_columns(metadata, REQUIRED_DYNAMICS_COLUMNS, METADATA_PATH)
    metadata = _coerce_id_column(metadata, METADATA_PATH)
    for column in ["pseudotime", "prAD", "ABA"]:
        metadata[column] = pd.to_numeric(metadata[column], errors="coerce")
    return metadata[REQUIRED_DYNAMICS_COLUMNS].copy()


def load_model_predictions() -> Dict[str, pd.DataFrame]:
    model_files = sorted(glob.glob(os.path.join(OUTPUTS_DIR, "full_dynamics_input_*.csv")))
    if not model_files:
        raise FileNotFoundError(
            f"No full_dynamics_input_*.csv files found in {OUTPUTS_DIR}"
        )

    model_data = {}
    for path in model_files:
        model_name = (
            os.path.basename(path)
            .replace("full_dynamics_input_", "")
            .replace(".csv", "")
        )
        if model_name not in MODEL_ORDER:
            continue
        pred = pd.read_csv(path)
        _require_columns(pred, REQUIRED_DYNAMICS_COLUMNS, path)
        pred = _coerce_id_column(pred, path)
        for column in ["pseudotime", "prAD", "ABA"]:
            pred[column] = pd.to_numeric(pred[column], errors="coerce")
        model_data[model_name] = pred[REQUIRED_DYNAMICS_COLUMNS].copy()

    missing_models = [model for model in MODEL_ORDER if model not in model_data]
    if missing_models:
        raise FileNotFoundError(
            f"Missing requested model dynamics files for: {missing_models}"
        )

    return {model: model_data[model] for model in MODEL_ORDER}


def assign_groups(df: pd.DataFrame) -> pd.Series:
    pseudotime = df["pseudotime"]
    groups = pd.Series(pd.NA, index=df.index, dtype="object")
    groups.loc[pseudotime < PSEUDOTIME_THRESHOLD] = "healthy"
    groups.loc[
        (df["prAD"] > PROBABILITY_THRESHOLD) & (pseudotime >= PSEUDOTIME_THRESHOLD)
    ] = "prAD"
    groups.loc[
        (df["ABA"] > PROBABILITY_THRESHOLD) & (pseudotime >= PSEUDOTIME_THRESHOLD)
    ] = "ABA"
    return groups


def build_plot_data(
    source_data: Dict[str, pd.DataFrame],
    metadata: pd.DataFrame,
    pathology: str,
) -> pd.DataFrame:
    rows = []
    for source, dynamics in source_data.items():
        merged = dynamics.merge(metadata[["ID", pathology]], on="ID", how="inner")
        merged["group"] = assign_groups(merged)
        merged["pathology_value"] = pd.to_numeric(merged[pathology], errors="coerce")
        merged["source"] = source
        merged = merged.dropna(subset=["group", "pathology_value"])
        rows.append(merged[["source", "group", "pathology_value"]])

    if not rows:
        return pd.DataFrame(columns=["source", "group", "pathology_value"])
    return pd.concat(rows, ignore_index=True)


def _safe_filename(name: str) -> str:
    return name.replace(".", "_").replace("/", "_").replace("\\", "_")


def _source_label(source: str, plot_data: pd.DataFrame) -> str:
    counts = (
        plot_data.loc[plot_data["source"] == source, "group"]
        .value_counts()
        .reindex(GROUP_ORDER, fill_value=0)
    )
    return (
        f"{source}\n"
        f"H={counts['healthy']} P={counts['prAD']} A={counts['ABA']}"
    )


def plot_pathology_violin(
    plot_data: pd.DataFrame,
    pathology: str,
    source_order: List[str],
) -> None:
    if plot_data.empty:
        print(f"Skipping {pathology}: no data after filtering.")
        return

    fig_width = max(10, 1.8 * len(source_order))
    fig, ax = plt.subplots(figsize=(fig_width, 6))

    sns.violinplot(
        data=plot_data,
        x="source",
        y="pathology_value",
        hue="group",
        order=source_order,
        hue_order=GROUP_ORDER,
        palette=GROUP_COLORS,
        cut=0,
        inner="quartile",
        linewidth=0.8,
        density_norm="width",
        ax=ax,
    )

    for collection in ax.collections:
        collection.set_alpha(0.65)

    ax.set_title(f"{pathology} by predicted trajectory group", fontsize=14, fontweight="bold")
    ax.set_xlabel("Source", fontsize=12)
    ax.set_ylabel(pathology, fontsize=12)
    ax.set_xticks(range(len(source_order)))
    ax.set_xticklabels(
        [_source_label(source, plot_data) for source in source_order],
        rotation=25,
        ha="right",
        fontsize=9,
    )
    ax.grid(axis="y", alpha=0.3)
    ax.legend(title="Group", loc="best")

    os.makedirs(FIGURES_DIR, exist_ok=True)
    output_stem = os.path.join(FIGURES_DIR, f"{_safe_filename(pathology)}_violin")
    fig.tight_layout()
    save_figure_for_illustrator(fig, f"{output_stem}.pdf")
    plt.close(fig)
    print(f"Saved {output_stem}.pdf and {output_stem}.svg")


def main() -> None:
    metadata = load_metadata()
    model_data = load_model_predictions()
    source_data = {
        "Gilad's": load_ground_truth(),
        "metadata": load_metadata_labels(),
        **model_data,
    }
    source_order = list(source_data.keys())

    print(f"Sources: {source_order}")
    for pathology in PATHOLOGY_COLUMNS:
        plot_data = build_plot_data(source_data, metadata, pathology)
        plot_pathology_violin(plot_data, pathology, source_order)


if __name__ == "__main__":
    main()

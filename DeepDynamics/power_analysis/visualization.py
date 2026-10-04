"""Save power-analysis tables, figures, and terminal summaries."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import pandas as pd

from power_analysis.config import (
    COHORT_SCOPES,
    FIGURE_ROOT,
    FIGURE_VECTOR_FORMATS,
    OUTPUT_ROOT,
    PATHOLOGIES,
    TRAJECTORIES,
)


def configure_matplotlib_for_illustrator() -> None:
    """Use vector output with editable text in Adobe Illustrator."""
    mpl.rcParams.update(
        {
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "axes.unicode_minus": False,
        }
    )


def disable_figure_rasterization(fig: Figure) -> None:
    """Keep artists as vectors instead of embedded bitmaps."""
    for ax in fig.axes:
        ax.set_rasterization_zorder(None)
        for artist in ax.get_children():
            if hasattr(artist, "set_rasterized"):
                artist.set_rasterized(False)


def save_figure_for_illustrator(fig: Figure, output_path: Path, formats: Iterable[str] = FIGURE_VECTOR_FORMATS) -> None:
    """Save fully vector figures with editable text for Illustrator."""
    disable_figure_rasterization(fig)
    base_path = output_path.with_suffix("")
    for fmt in formats:
        target = base_path.with_suffix(f".{fmt}")
        fig.savefig(
            target,
            format=fmt,
            bbox_inches="tight",
            facecolor="white",
            edgecolor="none",
        )
        print(f"Saved {target}")


configure_matplotlib_for_illustrator()


def readable_label(value: str) -> str:
    return value.replace("_", " ")


def ensure_output_dirs() -> Dict[str, Path]:
    paths = {
        "analytic_outputs": OUTPUT_ROOT / "analytic",
        "simulation_outputs": OUTPUT_ROOT / "simulation",
        "analytic_figures": FIGURE_ROOT / "analytic",
        "simulation_figures": FIGURE_ROOT / "simulation",
    }
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)
    return paths


def save_tables(
    cohort_analytic: pd.DataFrame,
    cohort_simulation: pd.DataFrame,
    incremental_analytic: pd.DataFrame,
    incremental_simulation: pd.DataFrame,
    cohort_analytic_ci: pd.DataFrame,
    cohort_simulation_ci: pd.DataFrame,
    incremental_analytic_ci: pd.DataFrame,
    incremental_simulation_ci: pd.DataFrame,
    paths: Dict[str, Path],
) -> None:
    cohort_analytic.to_csv(paths["analytic_outputs"] / "cohort_power_summary.csv", index=False)
    cohort_simulation.to_csv(paths["simulation_outputs"] / "cohort_power_summary.csv", index=False)
    incremental_analytic.to_csv(
        paths["analytic_outputs"] / "incremental_power_summary.csv",
        index=False,
    )
    incremental_simulation.to_csv(
        paths["simulation_outputs"] / "incremental_power_summary.csv",
        index=False,
    )
    cohort_analytic_ci.to_csv(
        paths["analytic_outputs"] / "cohort_power_bootstrap_ci.csv",
        index=False,
    )
    cohort_simulation_ci.to_csv(
        paths["simulation_outputs"] / "cohort_power_bootstrap_ci.csv",
        index=False,
    )
    incremental_analytic_ci.to_csv(
        paths["analytic_outputs"] / "incremental_power_bootstrap_ci.csv",
        index=False,
    )
    incremental_simulation_ci.to_csv(
        paths["simulation_outputs"] / "incremental_power_bootstrap_ci.csv",
        index=False,
    )


def merge_bootstrap_ci(summary: pd.DataFrame, ci: pd.DataFrame) -> pd.DataFrame:
    if ci.empty:
        return summary.copy()
    summary_keyed = summary.copy()
    ci_keyed = ci.copy()
    summary_keyed["sample_fraction_unshared"] = pd.to_numeric(
        summary_keyed["sample_fraction_unshared"],
        errors="coerce",
    ).fillna(-1.0)
    ci_keyed["sample_fraction_unshared"] = pd.to_numeric(
        ci_keyed["sample_fraction_unshared"],
        errors="coerce",
    ).fillna(-1.0)

    return summary_keyed.merge(
        ci_keyed,
        on=[
            "comparison",
            "pathology",
            "trajectory",
            "effect",
            "cohort_scope",
            "sample_fraction_unshared",
        ],
        how="left",
    )


def plot_grid(
    data: pd.DataFrame,
    method: str,
    effect: str,
    output_path: Path,
    x_column: str,
    x_label: str,
    x_order: tuple[str, ...] | None = None,
) -> None:
    fig, axes = plt.subplots(len(PATHOLOGIES), len(TRAJECTORIES), figsize=(12, 10), sharey=True)

    for row, pathology in enumerate(PATHOLOGIES):
        for col, trajectory in enumerate(TRAJECTORIES):
            ax = axes[row, col]
            plot_data = data[
                (data["effect"] == effect)
                & (data["pathology"] == pathology)
                & (data["trajectory"] == trajectory)
            ].copy()
            if x_order is not None:
                plot_data = plot_data.set_index(x_column).reindex(x_order).reset_index()
            else:
                plot_data = plot_data.sort_values(x_column)
            x_positions = list(range(len(plot_data)))
            (line,) = ax.plot(x_positions, plot_data["power"], marker="o", linewidth=2)
            line.set_rasterized(False)
            if {"power_ci_lower", "power_ci_upper"}.issubset(plot_data.columns):
                lower = plot_data["power_ci_lower"].to_numpy(dtype=float)
                upper = plot_data["power_ci_upper"].to_numpy(dtype=float)
                if pd.notna(lower).any() and pd.notna(upper).any():
                    ci_band = ax.fill_between(x_positions, lower, upper, alpha=0.2)
                    ci_band.set_rasterized(False)
            ax.axhline(0.8, color="gray", linestyle="--", linewidth=0.8)
            ax.set_ylim(0, 1.02)
            ax.set_title(f"{pathology} | {trajectory}")
            ax.set_xlabel(x_label)
            if col == 0:
                ax.set_ylabel("Power")
            labels = [f"{value:.2f}" if isinstance(value, float) else str(value) for value in plot_data[x_column]]
            ax.set_xticks(x_positions)
            ax.set_xticklabels(labels)
            ax.tick_params(axis="x", rotation=25)
            ax.grid(axis="y", alpha=0.3)

    fig.suptitle(f"{method.title()} power: {readable_label(effect)}")
    fig.tight_layout()
    save_figure_for_illustrator(fig, output_path)
    plt.close(fig)


def plot_power_figures(cohort: pd.DataFrame, incremental: pd.DataFrame, method: str, output_dir: Path) -> None:
    ok_cohort = cohort[cohort["status"] == "ok"]
    ok_incremental = incremental[incremental["status"] == "ok"]
    if ok_cohort.empty or ok_incremental.empty:
        print(f"Skipping {method} figures: no successful power estimates.")
        return

    for comparison in ok_cohort["comparison"].unique():
        comparison_cohort = ok_cohort[ok_cohort["comparison"] == comparison]
        comparison_incremental = ok_incremental[ok_incremental["comparison"] == comparison]
        for effect in comparison_cohort["effect"].unique():
            plot_grid(
                comparison_cohort,
                method,
                effect,
                output_dir / f"cohort_power_{comparison}_{effect}.pdf",
                "cohort_scope",
                "Cohort",
                COHORT_SCOPES,
            )
            plot_grid(
                comparison_incremental,
                method,
                effect,
                output_dir / f"incremental_power_{comparison}_{effect}.pdf",
                "sample_fraction_unshared",
                "Fraction of unshared samples added",
            )


def plot_effect_diagnostics(summary: pd.DataFrame, method: str, output_dir: Path) -> None:
    ok = summary[(summary["status"] == "ok") & (summary["cohort_scope"].isin(COHORT_SCOPES))]
    if ok.empty:
        return

    fig, axes = plt.subplots(len(PATHOLOGIES), len(TRAJECTORIES), figsize=(12, 10))
    for row, pathology in enumerate(PATHOLOGIES):
        for col, trajectory in enumerate(TRAJECTORIES):
            ax = axes[row, col]
            plot_data = ok[(ok["pathology"] == pathology) & (ok["trajectory"] == trajectory)]
            labels = [f"{r.cohort_scope}\n{readable_label(r.effect)}" for r in plot_data.itertuples()]
            errorbars = ax.errorbar(
                range(len(plot_data)),
                plot_data["estimate"],
                yerr=1.96 * plot_data["standard_error"],
                fmt="o",
                capsize=3,
            )
            errorbars[0].set_rasterized(False)
            for cap in errorbars[1]:
                cap.set_rasterized(False)
            errorbars[2].set_rasterized(False)
            ax.axhline(0, color="gray", linestyle="--", linewidth=0.8)
            ax.set_xticks(range(len(plot_data)))
            ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
            ax.set_title(f"{pathology} | {trajectory}")
            if col == 0:
                ax.set_ylabel("Estimate +/- 1.96 SE")
            ax.grid(axis="y", alpha=0.3)

    fig.suptitle(f"{method.title()} fitted effect diagnostics")
    fig.tight_layout()
    output_path = output_dir / "effect_size_diagnostics.pdf"
    save_figure_for_illustrator(fig, output_path)
    plt.close(fig)


def save_figures(
    cohort_analytic: pd.DataFrame,
    cohort_simulation: pd.DataFrame,
    incremental_analytic: pd.DataFrame,
    incremental_simulation: pd.DataFrame,
    cohort_analytic_ci: pd.DataFrame,
    cohort_simulation_ci: pd.DataFrame,
    incremental_analytic_ci: pd.DataFrame,
    incremental_simulation_ci: pd.DataFrame,
    paths: Dict[str, Path],
) -> None:
    cohort_analytic_with_ci = merge_bootstrap_ci(cohort_analytic, cohort_analytic_ci)
    cohort_simulation_with_ci = merge_bootstrap_ci(cohort_simulation, cohort_simulation_ci)
    incremental_analytic_with_ci = merge_bootstrap_ci(incremental_analytic, incremental_analytic_ci)
    incremental_simulation_with_ci = merge_bootstrap_ci(incremental_simulation, incremental_simulation_ci)

    plot_power_figures(
        cohort_analytic_with_ci,
        incremental_analytic_with_ci,
        "analytic",
        paths["analytic_figures"],
    )
    plot_effect_diagnostics(cohort_analytic, "analytic", paths["analytic_figures"])
    plot_power_figures(
        cohort_simulation_with_ci,
        incremental_simulation_with_ci,
        "simulation",
        paths["simulation_figures"],
    )
    plot_effect_diagnostics(cohort_simulation, "simulation", paths["simulation_figures"])


def print_power_gain(summary: pd.DataFrame, method: str) -> None:
    ok = summary[summary["status"] == "ok"]
    if ok.empty:
        print(f"\n{method.title()} summary: no successful power estimates.")
        return

    pivot = ok.pivot_table(
        index=["comparison", "pathology", "trajectory", "effect"],
        columns="cohort_scope",
        values="power",
        aggfunc="first",
    )
    if {"shared", "all"}.issubset(pivot.columns):
        pivot["all_minus_shared"] = pivot["all"] - pivot["shared"]

    columns = [col for col in ("shared", "unshared", "all", "all_minus_shared") if col in pivot]
    print(f"\n{method.title()} power gain summary:")
    print(pivot[columns].round(3).to_string())


def print_summary(cohort_analytic: pd.DataFrame, cohort_simulation: pd.DataFrame) -> None:
    print_power_gain(cohort_analytic, "analytic")
    print_power_gain(cohort_simulation, "simulation")

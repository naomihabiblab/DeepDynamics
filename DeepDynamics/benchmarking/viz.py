"""
Visualization utilities for benchmarking results.
"""

import os
from typing import Dict, List

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

mpl.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"})


def save_figure_for_illustrator(fig, path: str) -> None:
    """Save vector PDF/SVG with editable text for Adobe Illustrator."""
    base, _ = os.path.splitext(path)
    for ext in (".pdf", ".svg"):
        fig.savefig(f"{base}{ext}", bbox_inches="tight", facecolor="white")


def plot_benchmark_results(results: Dict[str, List[dict]], figures_path: str) -> None:
    """
    Create box plots for each metric showing distribution across bootstrap iterations.

    Args:
        results: dict mapping model_name -> list of metric dicts from each bootstrap iteration
        figures_path: directory to save figures into
    """
    os.makedirs(figures_path, exist_ok=True)

    models = list(results.keys())
    metric_names = list(results[models[0]][0].keys())

    metric_groups = {
        "Loss Metrics": ["total_loss", "kl_loss", "ce_loss", "mse_loss"],
        "Spearman Correlations": [
            "spearman_ABA",
            "spearman_prAD",
            "spearman_time",
        ],
        "Pearson Correlations": [
            "pearson_ABA",
            "pearson_prAD",
            "pearson_time",
        ],
    }

    colors = [
        "#2ecc71",
        "#3498db",
        "#e74c3c",
        "#9b59b6",
        "#f39c12",
    ]

    # Individual metric plots
    for metric in metric_names:
        fig, ax = plt.subplots(figsize=(8, 5))

        data = [[r[metric] for r in results[model]] for model in models]
        bp = ax.boxplot(data, patch_artist=True, tick_labels=models)

        for patch, color in zip(bp["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)

        ax.set_xlabel("Model", fontsize=12)
        ax.set_ylabel(metric, fontsize=12)
        ax.set_title(f"{metric} by Model", fontsize=14, fontweight="bold")
        ax.set_xticklabels(models, rotation=15, ha="right")
        ax.grid(axis="y", alpha=0.3)

        plt.tight_layout()
        save_figure_for_illustrator(
            fig,
            os.path.join(figures_path, f"benchmark_{metric}.pdf"),
        )
        plt.close()

    # Grouped plots
    for group_name, group_metrics in metric_groups.items():
        fig, axes = plt.subplots(
            1, len(group_metrics), figsize=(5 * len(group_metrics), 5)
        )
        if len(group_metrics) == 1:
            axes = [axes]

        for ax, metric in zip(axes, group_metrics):
            data = [[r[metric] for r in results[model]] for model in models]
            bp = ax.boxplot(data, patch_artist=True, tick_labels=models)

            for patch, color in zip(bp["boxes"], colors):
                patch.set_facecolor(color)
                patch.set_alpha(0.7)

            ax.set_xlabel("Model", fontsize=11)
            ax.set_ylabel(metric, fontsize=11)
            ax.set_title(metric, fontsize=12, fontweight="bold")
            ax.set_xticklabels(models, rotation=15, ha="right", fontsize=9)
            ax.grid(axis="y", alpha=0.3)

        fig.suptitle(group_name, fontsize=14, fontweight="bold", y=1.02)
        plt.tight_layout()
        safe_group_name = group_name.replace(" ", "_").lower()
        save_figure_for_illustrator(
            fig,
            os.path.join(figures_path, f"benchmark_{safe_group_name}.pdf"),
        )
        plt.close()

    print(f"Plots saved to {figures_path}")


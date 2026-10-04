"""
Scaling Laws Analysis for Bulk Sample Size.

This script systematically varies the number of training bulk samples and measures
how prediction performance scales. Results are visualized as scatter plots with
number of samples on the x-axis and various performance metrics on the y-axis.
"""

import sys
import os

# Add prediction folder to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'prediction'))

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import DataLoader
from sklearn.model_selection import train_test_split
from scipy.stats import spearmanr, pearsonr
from tqdm import tqdm
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"})


def save_figure_for_illustrator(fig, path: str) -> None:
    """Save vector PDF/SVG with editable text for Adobe Illustrator."""
    base, _ = os.path.splitext(path)
    for ext in (".pdf", ".svg"):
        fig.savefig(f"{base}{ext}", bbox_inches="tight", facecolor="white")


from model import ProbModel
from data_structures import SCData
from loss import ProbLoss


# Data paths (relative to scaling_laws folder)
DATA_PATH = '../prediction/data/shared_bulk_data_mask.csv'
TARGET_PATH = '../prediction/data/y.csv'
OUTPUT_PATH = 'outputs/scaling_laws_results.csv'
FIGURES_PATH = 'figures/'

# Scaling laws settings
SAMPLE_SIZES = [20, 35, 60, 100, 175, 278]  # Logarithmically spaced sample sizes (max ~278 with 25% test split)
N_BOOTSTRAP = 10  # Number of repetitions per sample size
TEST_SIZE = 0.25  # Fixed test set size (consistent with benchmark.py)

# ProbModel training settings (from benchmark.py)
PROBMODEL_EPOCHS = 2000
PROBMODEL_LR = 1e-5
PROBMODEL_BATCH_SIZE = 10
PROBMODEL_PATIENCE = 10
PROBMODEL_CONFIG = {"l1_lambda": 5.541763126909059e-05, "l2_lambda": 0}


def _get_device():
    return torch.device("cuda" if torch.cuda.is_available() else "cpu")


class ProbModelWrapper:
    """Wrapper to give ProbModel a sklearn-like interface (aligned with benchmarking/wrappers.py)."""

    def __init__(self):
        self.model = None
        self.device = _get_device()

    def fit(self, X, y):
        self.model = ProbModel(input_size=X.shape[1]).to(self.device)

        # Split into train/val for validation-based early stopping (matches benchmark)
        X_train, X_val, y_train, y_val = train_test_split(
            X, y, test_size=0.25, random_state=42,
        )

        train_data = SCData(
            torch.tensor(X_train, dtype=torch.float32),
            torch.tensor(y_train, dtype=torch.float32),
        )
        val_data = SCData(
            torch.tensor(X_val, dtype=torch.float32),
            torch.tensor(y_val, dtype=torch.float32),
        )
        train_loader = DataLoader(
            train_data,
            batch_size=min(PROBMODEL_BATCH_SIZE, len(X_train)),
            shuffle=True,
        )
        val_loader = DataLoader(
            val_data,
            batch_size=min(PROBMODEL_BATCH_SIZE, len(X_val)),
            shuffle=False,
        )

        optimizer = optim.Adam(self.model.parameters(), lr=PROBMODEL_LR)
        loss_func = ProbLoss()

        best_loss = np.inf
        counter = 0

        for _ in range(PROBMODEL_EPOCHS):
            # --- train ---
            self.model.train()
            for batch_X, batch_y in train_loader:
                batch_X = batch_X.to(self.device)
                batch_y = batch_y.to(self.device)
                output = self.model(batch_X)
                loss, _, _, _ = loss_func(output, batch_y)

                l1_reg = sum(param.abs().sum() for param in self.model.parameters())
                l2_reg = sum(param.pow(2).sum() for param in self.model.parameters())
                loss += (
                    PROBMODEL_CONFIG["l1_lambda"] * l1_reg
                    + PROBMODEL_CONFIG["l2_lambda"] * l2_reg
                )

                optimizer.zero_grad()
                loss.backward()
                optimizer.step()

            # --- validation-based early stopping ---
            self.model.eval()
            val_loss = 0.0
            with torch.no_grad():
                for batch_X, batch_y in val_loader:
                    batch_X = batch_X.to(self.device)
                    batch_y = batch_y.to(self.device)
                    output = self.model(batch_X)
                    loss, _, _, _ = loss_func(output, batch_y)
                    val_loss += loss.item()
            val_loss /= len(val_loader)

            if val_loss < best_loss:
                best_loss = val_loss
                counter = 0
            else:
                counter += 1

            if counter >= PROBMODEL_PATIENCE:
                break

        return self

    def predict(self, X):
        self.model.eval()
        with torch.no_grad():
            X_tensor = torch.tensor(X, dtype=torch.float32, device=self.device)
            output = self.model(X_tensor)
        return output.cpu().numpy()


def compute_metrics(predictions, targets):
    """
    Compute all evaluation metrics (logic aligned with benchmarking/metrics.py).
    Uses ABA/prAD naming: col 0 = ABA, col 1 = prAD.

    Args:
        predictions: numpy array (N, 3) - [ABA_prob, prAD_prob, pseudotime]
        targets: numpy array (N, 3) - ground truth values

    Returns:
        dict with all metrics
    """
    output = torch.tensor(predictions, dtype=torch.float32)
    target = torch.tensor(targets, dtype=torch.float32)

    # Loss components (same as ProbLoss)
    log_output = torch.log_softmax(output[:, :2], dim=1)
    kl_loss = F.kl_div(log_output, target[:, :2], reduction="batchmean").item()
    ce_loss = F.cross_entropy(output[:, :2], target[:, :2]).item()
    mse_loss = F.mse_loss(output[:, 2], target[:, 2]).item()

    total_loss = ce_loss + kl_loss + mse_loss

    pred_probs = F.softmax(output[:, :2], dim=1).numpy()
    pred_time = predictions[:, 2]

    # Global correlations (ABA = col 0, prAD = col 1)
    spearman_ABA, _ = spearmanr(pred_probs[:, 0], targets[:, 0])
    spearman_prAD, _ = spearmanr(pred_probs[:, 1], targets[:, 1])
    spearman_time, _ = spearmanr(pred_time, targets[:, 2])

    pearson_ABA, _ = pearsonr(pred_probs[:, 0], targets[:, 0])
    pearson_prAD, _ = pearsonr(pred_probs[:, 1], targets[:, 1])
    pearson_time, _ = pearsonr(pred_time, targets[:, 2])

    # Branch-specific pseudotime correlations
    ABA_branch_mask = pred_probs[:, 0] > 0.5
    prAD_branch_mask = pred_probs[:, 1] > 0.5

    if np.sum(ABA_branch_mask) >= 2:
        spearman_time_ABA, _ = spearmanr(
            pred_time[ABA_branch_mask], targets[ABA_branch_mask, 2]
        )
        pearson_time_ABA, _ = pearsonr(
            pred_time[ABA_branch_mask], targets[ABA_branch_mask, 2]
        )
    else:
        spearman_time_ABA = np.nan
        pearson_time_ABA = np.nan

    if np.sum(prAD_branch_mask) >= 2:
        spearman_time_prAD, _ = spearmanr(
            pred_time[prAD_branch_mask], targets[prAD_branch_mask, 2]
        )
        pearson_time_prAD, _ = pearsonr(
            pred_time[prAD_branch_mask], targets[prAD_branch_mask, 2]
        )
    else:
        spearman_time_prAD = np.nan
        pearson_time_prAD = np.nan

    mean_spearman_time_branch = np.nanmean(
        [spearman_time_ABA, spearman_time_prAD]
    )
    mean_pearson_time_branch = np.nanmean(
        [pearson_time_ABA, pearson_time_prAD]
    )

    return {
        "total_loss": total_loss,
        "kl_loss": kl_loss,
        "ce_loss": ce_loss,
        "mse_loss": mse_loss,
        "spearman_ABA": spearman_ABA,
        "spearman_prAD": spearman_prAD,
        "spearman_time": spearman_time,
        "pearson_ABA": pearson_ABA,
        "pearson_prAD": pearson_prAD,
        "pearson_time": pearson_time,
        "spearman_time_ABA": spearman_time_ABA,
        "spearman_time_prAD": spearman_time_prAD,
        "pearson_time_ABA": pearson_time_ABA,
        "pearson_time_prAD": pearson_time_prAD,
        "mean_spearman_time_branch": mean_spearman_time_branch,
        "mean_pearson_time_branch": mean_pearson_time_branch,
    }


# All metrics returned by compute_metrics (single source of truth for plotting)
ALL_METRICS = [
    "total_loss",
    "kl_loss",
    "ce_loss",
    "mse_loss",
    "spearman_ABA",
    "spearman_prAD",
    "spearman_time",
    "pearson_ABA",
    "pearson_prAD",
    "pearson_time",
    "spearman_time_ABA",
    "spearman_time_prAD",
    "pearson_time_ABA",
    "pearson_time_prAD",
    "mean_spearman_time_branch",
    "mean_pearson_time_branch",
]


def _safe_errorbar(ax, results_df, metric, color, **kwargs):
    """Plot one metric with NaN-safe error bars. No-op if metric missing or empty."""
    metric_data = results_df[results_df["metric"] == metric]
    if metric_data.empty:
        ax.text(0.5, 0.5, f"No data: {metric}", ha="center", va="center", transform=ax.transAxes)
        return
    x = metric_data["n_samples"].values
    y = metric_data["mean"].values
    # Replace NaN std with 0 so errorbar does not fail (e.g. branch metrics with no samples)
    yerr = np.nan_to_num(metric_data["std"].values, nan=0.0, posinf=0.0, neginf=0.0)
    ax.errorbar(x, y, yerr=yerr, color=color, **kwargs)
    ax.set_xscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(x)
    ax.grid(True, alpha=0.3)


def plot_scaling_laws(results_df):
    """
    Create scatter plots with error bars showing scaling laws.
    Plots all metrics from compute_metrics; handles NaN in std (e.g. branch-specific metrics).
    """
    os.makedirs(FIGURES_PATH, exist_ok=True)

    # Define metric groups (ABA/prAD naming, aligned with benchmarking)
    metric_groups = {
        "Loss Metrics": ["total_loss", "kl_loss", "ce_loss", "mse_loss"],
        "Spearman Correlations": ["spearman_ABA", "spearman_prAD", "spearman_time"],
        "Pearson Correlations": ["pearson_ABA", "pearson_prAD", "pearson_time"],
        "Branch-Specific Pseudotime Correlations": [
            "spearman_time_ABA",
            "spearman_time_prAD",
            "pearson_time_ABA",
            "pearson_time_prAD",
            "mean_spearman_time_branch",
            "mean_pearson_time_branch",
        ],
    }

    colors = {
        "total_loss": "#e74c3c",
        "kl_loss": "#3498db",
        "ce_loss": "#2ecc71",
        "mse_loss": "#9b59b6",
        "spearman_ABA": "#e74c3c",
        "spearman_prAD": "#3498db",
        "spearman_time": "#2ecc71",
        "pearson_ABA": "#e74c3c",
        "pearson_prAD": "#3498db",
        "pearson_time": "#2ecc71",
        "spearman_time_ABA": "#e67e22",
        "spearman_time_prAD": "#16a085",
        "pearson_time_ABA": "#e67e22",
        "pearson_time_prAD": "#16a085",
        "mean_spearman_time_branch": "#f39c12",
        "mean_pearson_time_branch": "#1abc9c",
    }
    # Default color for any metric not in dict
    default_color = "#333333"

    # Plot each metric group (only metrics that exist in results_df)
    for group_name, metrics in metric_groups.items():
        present = [m for m in metrics if (results_df["metric"] == m).any()]
        if not present:
            continue
        fig, axes = plt.subplots(1, len(present), figsize=(5 * len(present), 4.5))
        if len(present) == 1:
            axes = [axes]
        for ax, metric in zip(axes, present):
            c = colors.get(metric, default_color)
            _safe_errorbar(
                ax, results_df, metric, c,
                fmt="o-", capsize=4, capthick=1.5, markersize=8, linewidth=2, alpha=0.8,
            )
            ax.set_xlabel("Number of Training Samples", fontsize=11)
            ax.set_ylabel(metric, fontsize=11)
            ax.set_title(metric.replace("_", " ").title(), fontsize=12, fontweight="bold")
        fig.suptitle(f"Scaling Laws: {group_name}", fontsize=14, fontweight="bold", y=1.02)
        plt.tight_layout()
        safe_group_name = group_name.replace(" ", "_").lower()
        outpath = f"{FIGURES_PATH}scaling_laws_{safe_group_name}.pdf"
        try:
            save_figure_for_illustrator(fig, outpath)
        except Exception as e:
            print(f"Warning: could not save {outpath}: {e}")
        plt.close()

    # One figure with every metric (guarantees all metrics are plotted)
    n_metrics = len(ALL_METRICS)
    n_cols = 4
    n_rows = (n_metrics + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 3.5 * n_rows))
    axes = np.atleast_2d(axes)
    for idx, metric in enumerate(ALL_METRICS):
        r, c = idx // n_cols, idx % n_cols
        ax = axes[r, c]
        c = colors.get(metric, default_color)
        _safe_errorbar(
            ax, results_df, metric, c,
            fmt="o-", capsize=3, markersize=5, linewidth=1.5, alpha=0.8,
        )
        ax.set_xlabel("N samples", fontsize=9)
        ax.set_ylabel(metric, fontsize=9)
        ax.set_title(metric.replace("_", " ").title(), fontsize=10)
    # Hide unused subplots
    for idx in range(n_metrics, n_rows * n_cols):
        r, c = idx // n_cols, idx % n_cols
        axes[r, c].set_visible(False)
    fig.suptitle("Scaling Laws: All Metrics", fontsize=14, fontweight="bold", y=1.001)
    plt.tight_layout()
    try:
        save_figure_for_illustrator(fig, f"{FIGURES_PATH}scaling_laws_all_metrics.pdf")
    except Exception as e:
        print(f"Warning: could not save scaling_laws_all_metrics.pdf: {e}")
    plt.close()

    # Combined overview plot (4 panels)
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    overview = [
        ("total_loss", "#e74c3c", "Total Loss"),
        ("mse_loss", "#9b59b6", "MSE Loss"),
        ("pearson_time", "#2ecc71", "Pearson (Pseudotime)"),
        ("spearman_time", "#3498db", "Spearman (Pseudotime)"),
    ]
    for ax, (metric, color, title) in zip(axes.flat, overview):
        _safe_errorbar(
            ax, results_df, metric, color,
            fmt="o-", capsize=4, markersize=8, linewidth=2,
        )
        ax.set_xlabel("Number of Training Samples", fontsize=11)
        ax.set_ylabel(metric if "time" not in metric else "Correlation", fontsize=11)
        ax.set_title(title, fontsize=12, fontweight="bold")
    fig.suptitle("Scaling Laws Overview", fontsize=14, fontweight="bold")
    plt.tight_layout()
    try:
        save_figure_for_illustrator(fig, f"{FIGURES_PATH}scaling_laws_overview.pdf")
    except Exception as e:
        print(f"Warning: could not save scaling_laws_overview.pdf: {e}")
    plt.close()

    print(f"Plots saved to {FIGURES_PATH}")


def run_scaling_laws():
    """Run the scaling laws analysis."""
    print("Loading data...")
    X_df = pd.read_csv(DATA_PATH, index_col=0)
    y_df = pd.read_csv(TARGET_PATH, index_col=0)

    X = X_df.values.astype(np.float32)
    y = y_df.values.astype(np.float32)

    print(f"Data loaded: X shape = {X.shape}, y shape = {y.shape}")

    # Initial split: fixed holdout test set (same as benchmark.py, no leakage)
    print("\n" + "=" * 60)
    print("INITIAL DATA SPLIT")
    print("=" * 60)
    print("Splitting data into train/test (random_state=42)...")
    X_train_full, X_test_heldout, y_train_full, y_test_heldout = train_test_split(
        X, y, test_size=TEST_SIZE, random_state=42,
    )
    print(f"Training pool: {X_train_full.shape[0]} samples")
    print(f"Test set (held out): {X_test_heldout.shape[0]} samples")

    # Filter sample sizes that are feasible given training pool size
    max_train_size = len(X_train_full)
    valid_sample_sizes = [s for s in SAMPLE_SIZES if s <= max_train_size]

    if len(valid_sample_sizes) < len(SAMPLE_SIZES):
        print(f"Note: Some sample sizes were too large. Using: {valid_sample_sizes}")

    print("\n" + "=" * 60)
    print("SCALING LAWS ANALYSIS")
    print("=" * 60)
    print(f"Sample sizes to test: {valid_sample_sizes}")
    print(f"Bootstrap iterations per size: {N_BOOTSTRAP}")
    print(f"Test set size: {TEST_SIZE * 100:.0f}% (fixed holdout)")

    # Store results
    all_results = []

    # For each sample size
    for n_samples in valid_sample_sizes:
        print(f"\n--- Training with {n_samples} samples ---")

        sample_results = []

        for i in tqdm(range(N_BOOTSTRAP), desc=f"n={n_samples}"):
            # Subsample n_samples from training pool (different subsample each bootstrap i)
            indices = np.random.RandomState(i).choice(
                len(X_train_full), size=n_samples, replace=False
            )
            X_train = X_train_full[indices]
            y_train = y_train_full[indices]

            # Train (with validation-based early stopping inside fit) and evaluate on fixed test set
            model = ProbModelWrapper()
            model.fit(X_train, y_train)
            predictions = model.predict(X_test_heldout)
            metrics = compute_metrics(predictions, y_test_heldout)

            sample_results.append(metrics)
        
        # Aggregate results for this sample size
        metric_names = list(sample_results[0].keys())
        for metric_name in metric_names:
            values = [r[metric_name] for r in sample_results]
            all_results.append({
                'n_samples': n_samples,
                'metric': metric_name,
                'mean': np.mean(values),
                'std': np.std(values),
            })
    
    # Save to CSV
    results_df = pd.DataFrame(all_results)
    results_df.to_csv(OUTPUT_PATH, index=False)
    print(f"\nResults saved to {OUTPUT_PATH}")
    
    # Print summary
    print("\n" + "=" * 60)
    print("SCALING LAWS SUMMARY")
    print("=" * 60)
    
    for n_samples in valid_sample_sizes:
        print(f"\nn_samples = {n_samples}:")
        sample_data = results_df[results_df['n_samples'] == n_samples]
        for _, row in sample_data.iterrows():
            print(f"  {row['metric']}: {row['mean']:.4f} ± {row['std']:.4f}")
    
    # Generate plots
    plot_scaling_laws(results_df)
    
    return results_df


if __name__ == "__main__":
    # Change to script directory for relative paths
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    run_scaling_laws()

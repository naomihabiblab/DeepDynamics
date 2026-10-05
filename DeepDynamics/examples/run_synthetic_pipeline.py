"""Run a compact, end-to-end DeepDynamics example on synthetic data.

The script keeps training data and model weights in memory. Generated files are
written only to an explicit output directory, or to a temporary directory when
``--output-dir`` is omitted.
"""

from __future__ import annotations

import argparse
import json
import math
import random
import re
import sys
import tempfile
from copy import deepcopy
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from scipy import stats
from sklearn.model_selection import train_test_split
from torch.utils.data import DataLoader, TensorDataset

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from prediction.loss import ProbLoss
from prediction.model import ProbModel


PROVIDED_SYNTHETIC_CONTRACT = {
    "observations": 437,
    "all_features": 91,
    "shared_samples": 419,
    "validation_samples": 673,
    "retained_features": 58,
    "labelled_samples": 371,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Train DeepDynamics and plot synthetic cell-state dynamics."
    )
    parser.add_argument(
        "--data",
        type=Path,
        default=PROJECT_ROOT / "prediction" / "data" / "synthetic" / "500.h5ad",
        help="Synthetic AnnData input.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Local output directory. Defaults to a new system temporary directory.",
    )
    parser.add_argument(
        "--adj-pval-threshold",
        type=float,
        default=0.005,
        help="Maximum adjusted p-value for retaining a cell state (default: 0.005).",
    )
    parser.add_argument(
        "--min-correlation",
        type=float,
        default=0.0,
        help="Minimum stored correlation for retaining a cell state (default: 0).",
    )
    parser.add_argument(
        "--trajectory-columns",
        nargs=2,
        metavar=("TRAJECTORY_1", "TRAJECTORY_2"),
        default=["prAD", "ABA"],
        help="Two branch-probability columns in the AnnData trajectory table.",
    )
    parser.add_argument(
        "--pseudotime-column",
        default="psuedotime",
        help="Name assigned to the pseudotime target and prediction column.",
    )
    parser.add_argument(
        "--plot-cell-states",
        nargs="*",
        help="Exactly three retained cell states to plot. Defaults to the top three correlations.",
    )
    parser.add_argument("--epochs", type=int, default=200)
    parser.add_argument("--patience", type=int, default=10)
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--learning-rate", type=float, default=1e-4)
    parser.add_argument("--test-size", type=float, default=0.25)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--device",
        choices=["auto", "cpu", "cuda"],
        default="auto",
        help="Training device (default: auto).",
    )
    parser.add_argument(
        "--smoothing-bandwidth",
        type=float,
        default=0.10,
        help="Gaussian-kernel bandwidth as a fraction of the plotted pseudotime range.",
    )
    return parser.parse_args()


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    if hasattr(torch.backends, "cudnn"):
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False


def resolve_device(requested: str) -> torch.device:
    if requested == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("--device cuda was requested, but CUDA is unavailable")
    if requested == "auto":
        requested = "cuda" if torch.cuda.is_available() else "cpu"
    return torch.device(requested)


def prepare_output_dir(path: Path | None) -> Path:
    if path is None:
        return Path(tempfile.mkdtemp(prefix="deepdynamics-toy-"))
    output_dir = path.expanduser().resolve()
    if output_dir == PROJECT_ROOT or PROJECT_ROOT in output_dir.parents:
        raise ValueError(
            "Toy outputs must be written outside the repository; choose a temporary "
            "or private analysis directory"
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


def load_pipeline_data(
    path: Path,
    adj_pval_threshold: float,
    min_correlation: float,
    trajectory_columns: list[str],
    pseudotime_column: str,
) -> dict:
    if not 0 <= adj_pval_threshold <= 1:
        raise ValueError("--adj-pval-threshold must be between 0 and 1")

    data = ad.read_h5ad(path)
    if "synthetic" not in data.uns:
        raise ValueError(
            "Input is not marked as a synthetic cohort; refusing to create toy-run "
            "artifacts from potentially restricted data"
        )
    celmod = data.uns["celmod"]
    correlations = celmod["test.corrs"]
    required_corr_columns = {"adj.pval", "corr"}
    missing_corr_columns = required_corr_columns.difference(correlations.columns)
    if missing_corr_columns:
        raise ValueError(
            f"Correlation table is missing columns: {sorted(missing_corr_columns)}"
        )

    retained = correlations.index[
        (correlations["adj.pval"] < adj_pval_threshold)
        & (correlations["corr"] > min_correlation)
    ].astype(str)
    if len(retained) < 3:
        raise ValueError(
            f"Only {len(retained)} cell states passed the filter; at least 3 are required"
        )

    bulk_tables = celmod["avg.predicted.prop"]
    shared = bulk_tables["train"].copy()
    validation = bulk_tables["validation"].copy()
    shared.columns = shared.columns.astype(str)
    validation.columns = validation.columns.astype(str)
    retained_list = retained.tolist()
    shared = shared.loc[:, retained_list]
    validation = validation.loc[:, retained_list]

    trajectory = data.uns["trajectories"]["palantir"]
    branch_probabilities = trajectory["branch.probs"].copy()
    missing_targets = set(trajectory_columns).difference(branch_probabilities.columns)
    if missing_targets:
        raise ValueError(
            f"Trajectory table is missing requested columns: {sorted(missing_targets)}"
        )
    pseudotime = np.asarray(trajectory["pseudotime"], dtype=float)
    if len(pseudotime) != len(branch_probabilities):
        raise ValueError("Pseudotime and branch-probability row counts differ")

    target = branch_probabilities.loc[:, trajectory_columns].copy()
    target[pseudotime_column] = pseudotime
    shared.index = shared.index.astype(str)
    validation.index = validation.index.astype(str)
    target.index = target.index.astype(str)
    labelled_ids = shared.index.intersection(target.index[target[pseudotime_column].notna()])
    labelled_features = shared.loc[labelled_ids].astype(float)
    labelled_target = target.loc[labelled_ids].astype(float)

    if labelled_features.isna().any().any() or labelled_target.isna().any().any():
        raise ValueError("Filtered training data contain missing values")
    if not np.isfinite(labelled_features.to_numpy()).all():
        raise ValueError("Filtered training features contain non-finite values")

    all_bulk = pd.concat([shared, validation], axis=0)
    if all_bulk.index.has_duplicates:
        raise ValueError("Shared and validation bulk tables contain duplicate sample IDs")

    observed_contract = {
        "observations": int(data.n_obs),
        "all_features": int(data.n_vars),
        "shared_samples": int(len(shared)),
        "validation_samples": int(len(validation)),
        "retained_features": int(len(retained_list)),
        "labelled_samples": int(len(labelled_features)),
    }
    matches_provided_contract = observed_contract == PROVIDED_SYNTHETIC_CONTRACT

    return {
        "data": data,
        "correlations": correlations,
        "retained_features": retained_list,
        "labelled_features": labelled_features,
        "labelled_target": labelled_target,
        "all_bulk": all_bulk.astype(float),
        "observed_contract": observed_contract,
        "matches_provided_contract": matches_provided_contract,
    }


def make_loader(
    features: pd.DataFrame,
    target: pd.DataFrame,
    batch_size: int,
    shuffle: bool,
    seed: int,
) -> DataLoader:
    dataset = TensorDataset(
        torch.tensor(features.to_numpy(), dtype=torch.float32),
        torch.tensor(target.to_numpy(), dtype=torch.float32),
    )
    generator = torch.Generator().manual_seed(seed)
    return DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=shuffle,
        generator=generator,
    )


def evaluate_loss(
    model: ProbModel,
    loader: DataLoader,
    loss_function: ProbLoss,
    device: torch.device,
) -> float:
    model.eval()
    total = 0.0
    sample_count = 0
    with torch.no_grad():
        for features, target in loader:
            features = features.to(device)
            target = target.to(device)
            loss, _, _, _ = loss_function(model(features), target)
            batch_count = len(features)
            total += float(loss.item()) * batch_count
            sample_count += batch_count
    return total / sample_count


def train_model(
    model: ProbModel,
    train_loader: DataLoader,
    test_loader: DataLoader,
    device: torch.device,
    epochs: int,
    patience: int,
    learning_rate: float,
) -> tuple[list[float], list[float]]:
    loss_function = ProbLoss().to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
    train_history: list[float] = []
    test_history: list[float] = []
    best_loss = math.inf
    best_state = deepcopy(model.state_dict())
    stale_epochs = 0

    for epoch in range(epochs):
        model.train()
        total = 0.0
        sample_count = 0
        for features, target in train_loader:
            features = features.to(device)
            target = target.to(device)
            optimizer.zero_grad()
            loss, _, _, _ = loss_function(model(features), target)
            loss.backward()
            optimizer.step()
            batch_count = len(features)
            total += float(loss.item()) * batch_count
            sample_count += batch_count

        train_loss = total / sample_count
        test_loss = evaluate_loss(model, test_loader, loss_function, device)
        train_history.append(train_loss)
        test_history.append(test_loss)

        if test_loss < best_loss:
            best_loss = test_loss
            best_state = deepcopy(model.state_dict())
            stale_epochs = 0
        else:
            stale_epochs += 1

        if (epoch + 1) % 10 == 0 or epoch == 0:
            print(
                f"Epoch {epoch + 1:>3}: train_loss={train_loss:.5f}, "
                f"test_loss={test_loss:.5f}"
            )
        if stale_epochs >= patience:
            print(f"Early stopping after epoch {epoch + 1}")
            break

    model.load_state_dict(best_state)
    return train_history, test_history


def predict(model: ProbModel, features: pd.DataFrame, device: torch.device) -> np.ndarray:
    model.eval()
    tensor = torch.tensor(features.to_numpy(), dtype=torch.float32, device=device)
    with torch.no_grad():
        raw = model(tensor)
        branch_probabilities = torch.softmax(raw[:, :2], dim=1)
        output = torch.cat([branch_probabilities, raw[:, 2:3]], dim=1)
    return output.cpu().numpy()


def prediction_metrics(prediction: np.ndarray, target: pd.DataFrame) -> dict:
    truth = target.to_numpy(dtype=float)
    pseudotime_true = truth[:, 2]
    pseudotime_pred = prediction[:, 2]
    pearson = stats.pearsonr(pseudotime_true, pseudotime_pred)
    return {
        "pseudotime_mse": float(np.mean((pseudotime_pred - pseudotime_true) ** 2)),
        "pseudotime_pearson": float(pearson.statistic),
        "dominant_branch_accuracy": float(
            np.mean(np.argmax(prediction[:, :2], axis=1) == np.argmax(truth[:, :2], axis=1))
        ),
        "trajectory_1_mae": float(np.mean(np.abs(prediction[:, 0] - truth[:, 0]))),
        "trajectory_2_mae": float(np.mean(np.abs(prediction[:, 1] - truth[:, 1]))),
        "max_branch_sum_error": float(
            np.max(np.abs(prediction[:, :2].sum(axis=1) - 1.0))
        ),
    }


def choose_plot_features(
    requested: list[str] | None,
    retained_features: list[str],
    correlations: pd.DataFrame,
) -> list[str]:
    if requested:
        if len(requested) != 3:
            raise ValueError("--plot-cell-states requires exactly three names")
        missing = set(requested).difference(retained_features)
        if missing:
            raise ValueError(
                f"Requested plot cell states were not retained: {sorted(missing)}"
            )
        return requested
    ranked = correlations.loc[retained_features, "corr"].sort_values(ascending=False)
    return ranked.index[:3].astype(str).tolist()


def weighted_kernel_curve(
    pseudotime: np.ndarray,
    values: np.ndarray,
    trajectory_probability: np.ndarray,
    grid: np.ndarray,
    bandwidth: float,
) -> tuple[np.ndarray, np.ndarray]:
    means = np.full_like(grid, np.nan, dtype=float)
    confidence = np.full_like(grid, np.nan, dtype=float)
    for index, point in enumerate(grid):
        local_weight = np.exp(-0.5 * ((pseudotime - point) / bandwidth) ** 2)
        weight = trajectory_probability * local_weight
        total_weight = weight.sum()
        if total_weight <= 1e-10:
            continue
        mean = np.sum(weight * values) / total_weight
        effective_n = total_weight**2 / np.sum(weight**2)
        variance = np.sum(weight * (values - mean) ** 2) / total_weight
        means[index] = mean
        confidence[index] = 1.96 * np.sqrt(variance / max(effective_n, 1.0))
    return means, confidence


def safe_filename(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._")
    return cleaned or "cell_state"


def plot_dynamics(
    dynamics_input: pd.DataFrame,
    cell_states: list[str],
    trajectory_columns: list[str],
    pseudotime_column: str,
    bandwidth_fraction: float,
    output_dir: Path,
) -> list[Path]:
    if bandwidth_fraction <= 0:
        raise ValueError("--smoothing-bandwidth must be positive")
    plot_dir = output_dir / "dynamics"
    plot_dir.mkdir(parents=True, exist_ok=True)

    pseudotime = dynamics_input[pseudotime_column].to_numpy(dtype=float)
    finite = np.isfinite(pseudotime)
    low, high = np.quantile(pseudotime[finite], [0.01, 0.99])
    if high <= low:
        raise ValueError("Predicted pseudotime has no usable range for dynamics plots")
    grid = np.linspace(low, high, 200)
    bandwidth = max((high - low) * bandwidth_fraction, 1e-6)
    colors = ["#1f77b4", "#d62728"]
    written: list[Path] = []

    for cell_state in cell_states:
        values = dynamics_input[cell_state].to_numpy(dtype=float)
        fig, axis = plt.subplots(figsize=(7, 4.5))
        for trajectory, color in zip(trajectory_columns, colors, strict=True):
            probability = dynamics_input[trajectory].to_numpy(dtype=float)
            keep = finite & np.isfinite(values) & np.isfinite(probability)
            mean, confidence = weighted_kernel_curve(
                pseudotime[keep],
                values[keep],
                probability[keep],
                grid,
                bandwidth,
            )
            axis.plot(grid, mean, color=color, linewidth=2, label=trajectory)
            axis.fill_between(
                grid,
                mean - confidence,
                mean + confidence,
                color=color,
                alpha=0.18,
                linewidth=0,
            )
        axis.set_title(f"Synthetic {cell_state} dynamics")
        axis.set_xlabel(f"Predicted {pseudotime_column}")
        axis.set_ylabel("Synthetic cell-state value")
        axis.legend(title="Trajectory")
        fig.tight_layout()
        output_path = plot_dir / f"{safe_filename(cell_state)}_dynamics.png"
        fig.savefig(output_path, dpi=160)
        plt.close(fig)
        written.append(output_path)
    return written


def main() -> None:
    args = parse_args()
    set_seed(args.seed)
    device = resolve_device(args.device)
    trajectory_columns = list(args.trajectory_columns)

    bundle = load_pipeline_data(
        args.data,
        args.adj_pval_threshold,
        args.min_correlation,
        trajectory_columns,
        args.pseudotime_column,
    )
    output_dir = prepare_output_dir(args.output_dir)
    features = bundle["labelled_features"]
    target = bundle["labelled_target"]
    train_ids, test_ids = train_test_split(
        features.index,
        test_size=args.test_size,
        random_state=args.seed,
    )
    x_train, x_test = features.loc[train_ids], features.loc[test_ids]
    y_train, y_test = target.loc[train_ids], target.loc[test_ids]

    train_loader = make_loader(
        x_train, y_train, args.batch_size, shuffle=True, seed=args.seed
    )
    test_loader = make_loader(
        x_test, y_test, args.batch_size, shuffle=False, seed=args.seed
    )
    model = ProbModel(input_size=features.shape[1]).to(device)
    train_history, test_history = train_model(
        model,
        train_loader,
        test_loader,
        device,
        args.epochs,
        args.patience,
        args.learning_rate,
    )

    test_prediction = predict(model, x_test, device)
    metrics = prediction_metrics(test_prediction, y_test)
    all_prediction = predict(model, bundle["all_bulk"], device)
    prediction_columns = trajectory_columns + [args.pseudotime_column]
    predictions = pd.DataFrame(
        all_prediction,
        index=bundle["all_bulk"].index,
        columns=prediction_columns,
    )
    predictions.index.name = "synthetic_sample_id"

    plot_cell_states = choose_plot_features(
        args.plot_cell_states,
        bundle["retained_features"],
        bundle["correlations"],
    )
    dynamics_input = pd.concat(
        [bundle["all_bulk"].loc[:, plot_cell_states], predictions], axis=1
    )
    plot_paths = plot_dynamics(
        dynamics_input,
        plot_cell_states,
        trajectory_columns,
        args.pseudotime_column,
        args.smoothing_bandwidth,
        output_dir,
    )

    predictions.to_csv(output_dir / "predictions.csv")
    dynamics_input.to_csv(output_dir / "dynamics_input.csv")
    summary = {
        "input": str(args.data.resolve()),
        "synthetic_only": True,
        "device": str(device),
        "seed": args.seed,
        "filter": {
            "adj_pval_threshold": args.adj_pval_threshold,
            "min_correlation": args.min_correlation,
        },
        "target_columns": prediction_columns,
        "observed_contract": bundle["observed_contract"],
        "matches_provided_synthetic_contract": bundle["matches_provided_contract"],
        "train_samples": len(x_train),
        "test_samples": len(x_test),
        "epochs_completed": len(train_history),
        "best_test_loss": float(min(test_history)),
        "test_metrics": metrics,
        "plot_cell_states": plot_cell_states,
        "plot_files": [str(path) for path in plot_paths],
    }
    with (output_dir / "run_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    print("\nSynthetic DeepDynamics toy run complete")
    print(f"  input contract: {bundle['observed_contract']}")
    print(
        "  provided-data contract matched: "
        f"{bundle['matches_provided_contract']}"
    )
    print(f"  train/test: {len(x_train)} / {len(x_test)}")
    print(f"  retained features: {features.shape[1]}")
    print(f"  pseudotime MSE: {metrics['pseudotime_mse']:.5f}")
    print(f"  pseudotime Pearson: {metrics['pseudotime_pearson']:.5f}")
    print(f"  branch accuracy: {metrics['dominant_branch_accuracy']:.5f}")
    print(f"  dynamics cell states: {', '.join(plot_cell_states)}")
    print(f"  local outputs: {output_dir}")


if __name__ == "__main__":
    main()

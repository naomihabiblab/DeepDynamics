#!/usr/bin/env python
"""
Compute DTW distance between *ground-truth dynamics* and model-predicted dynamics.

This script expects dynamics predictions produced by `yifat_stuff/run_dynamics.R`,
specifically the `dynamics_pred_vals_*.csv` files containing smoothed trajectories
with columns like:
  - x (pseudotime grid)
  - fit (smoothed mean)
  - se.fit (optional)
  - feature (pathology/cluster feature name)
  - trajectory ("prAD" / "ABA")

It computes DTW distances between the GT curve and each model curve per
(trajectory, feature), and writes a single CSV summary.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import tempfile
from dataclasses import dataclass
from shutil import which
from typing import Dict, Iterable, List, Optional

import numpy as np
import pandas as pd


REQUIRED_PRED_COLS = {"x", "fit", "feature", "trajectory"}


def _dtw_distance(y_a: np.ndarray, y_b: np.ndarray) -> float:
    """
    Classic DTW distance (L1 local cost) with O(N*M) DP.
    For our use-case N,M are small (~50), so this is fine and dependency-free.
    """
    y_a = np.asarray(y_a, dtype=float)
    y_b = np.asarray(y_b, dtype=float)
    if y_a.size == 0 or y_b.size == 0:
        return float("nan")

    n, m = int(y_a.size), int(y_b.size)
    dp = np.full((n + 1, m + 1), np.inf, dtype=float)
    dp[0, 0] = 0.0

    for i in range(1, n + 1):
        a = y_a[i - 1]
        for j in range(1, m + 1):
            cost = abs(a - y_b[j - 1])
            dp[i, j] = cost + min(dp[i - 1, j],
                                  dp[i, j - 1],
                                  dp[i - 1, j - 1])

    return float(dp[n, m])


def _interp_to_grid(x: np.ndarray, y: np.ndarray, grid: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    grid = np.asarray(grid, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    if x.size < 2:
        return np.full_like(grid, np.nan, dtype=float)
    return np.interp(grid, x, y)


def _read_pred_vals(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    missing = REQUIRED_PRED_COLS.difference(df.columns)
    if missing:
        raise ValueError(f"Missing required columns {sorted(missing)} in {path}")
    return df


@dataclass(frozen=True)
class DtwRow:
    model: str
    trajectory: str
    feature: str
    dtw_l1: float
    dtw_l1_per_step: float
    n_grid: int
    x_min: float
    x_max: float
    gt_path: str
    model_path: str


def _assert_monotonic_x(
    x: np.ndarray,
    *,
    model: str,
    trajectory: str,
    feature: str,
    source_path: str,
) -> None:
    """
    Sanity check that x is monotonic non-decreasing.

    We do NOT sort, because that could hide upstream issues and change the implied
    mapping between points and pseudotime order. `np.interp` requires sorted x.
    """
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < 2:
        return
    d = np.diff(x)
    if np.any(d < 0):
        first = int(np.where(d < 0)[0][0])
        raise SystemExit(
            "Dynamics CSV is not sorted by pseudotime `x` within (trajectory, feature).\n"
            f"- model: {model}\n"
            f"- trajectory: {trajectory}\n"
            f"- feature: {feature}\n"
            f"- file: {os.path.abspath(source_path)}\n"
            f"- first decrease at indices {first}->{first+1}: x={x[first]} then x={x[first+1]}\n"
            "Fix: ensure `dynamics_pred_vals_*.csv` rows are ordered by increasing `x` within each group."
        )


def compute_dtw_table(
    *,
    gt_pred_vals_csv: str,
    model_pred_vals: Dict[str, str],
    features: Optional[Iterable[str]] = None,
    trajectories: Iterable[str] = ("prAD", "ABA"),
    n_grid: int = 50,
) -> pd.DataFrame:
    """
    Compute DTW distances between GT dynamics and each model dynamics.

    We interpolate both GT and model curves onto a *shared, evenly spaced pseudotime grid*
    over the overlap of their x ranges for that (trajectory, feature).
    """
    gt_df = _read_pred_vals(gt_pred_vals_csv)
    if features is None:
        features = sorted(gt_df["feature"].dropna().unique().tolist())
    else:
        features = list(features)

    rows: List[DtwRow] = []
    trajectories = list(trajectories)

    for model_name, model_path in model_pred_vals.items():
        m_df = _read_pred_vals(model_path)

        for traj in trajectories:
            for feat in features:
                gt_sub = gt_df[(gt_df["trajectory"] == traj) & (gt_df["feature"] == feat)]
                m_sub = m_df[(m_df["trajectory"] == traj) & (m_df["feature"] == feat)]
                if gt_sub.empty or m_sub.empty:
                    continue

                x_gt = gt_sub["x"].to_numpy(dtype=float)
                y_gt = gt_sub["fit"].to_numpy(dtype=float)
                x_m = m_sub["x"].to_numpy(dtype=float)
                y_m = m_sub["fit"].to_numpy(dtype=float)

                _assert_monotonic_x(
                    x_gt,
                    model="ground_truth",
                    trajectory=traj,
                    feature=feat,
                    source_path=gt_pred_vals_csv,
                )
                _assert_monotonic_x(
                    x_m,
                    model=model_name,
                    trajectory=traj,
                    feature=feat,
                    source_path=model_path,
                )

                min_x = max(np.nanmin(x_gt), np.nanmin(x_m))
                max_x = min(np.nanmax(x_gt), np.nanmax(x_m))
                if not np.isfinite(min_x) or not np.isfinite(max_x) or max_x <= min_x:
                    continue

                grid = np.linspace(min_x, max_x, int(n_grid))
                y_gt_g = _interp_to_grid(x_gt, y_gt, grid)
                y_m_g = _interp_to_grid(x_m, y_m, grid)
                valid = np.isfinite(y_gt_g) & np.isfinite(y_m_g)
                if int(valid.sum()) < 3:
                    continue

                y_gt_g = y_gt_g[valid]
                y_m_g = y_m_g[valid]

                dtw = _dtw_distance(y_gt_g, y_m_g)
                rows.append(
                    DtwRow(
                        model=model_name,
                        trajectory=traj,
                        feature=feat,
                        dtw_l1=dtw,
                        dtw_l1_per_step=float(dtw / max(len(y_gt_g), 1)),
                        n_grid=int(len(y_gt_g)),
                        x_min=float(min_x),
                        x_max=float(max_x),
                        gt_path=os.path.abspath(gt_pred_vals_csv),
                        model_path=os.path.abspath(model_path),
                    )
                )

    out = pd.DataFrame([r.__dict__ for r in rows])
    if not out.empty:
        out = out.sort_values(["model", "trajectory", "feature"]).reset_index(drop=True)
    return out


def _discover_model_pred_vals(outputs_dir: str, full_only: bool = True) -> Dict[str, str]:
    """
    Discover per-model dynamics_pred_vals CSVs under benchmarking/outputs.

    Expected layout:
      outputs/dynamics_<model>/dynamics_pred_vals_<model>_full.csv
    """
    outputs_dir = os.path.abspath(outputs_dir)
    model_map: Dict[str, str] = {}
    if not os.path.isdir(outputs_dir):
        return model_map

    for name in os.listdir(outputs_dir):
        if not name.startswith("dynamics_"):
            continue
        subdir = os.path.join(outputs_dir, name)
        if not os.path.isdir(subdir):
            continue

        for fname in os.listdir(subdir):
            if not (fname.startswith("dynamics_pred_vals_") and fname.endswith(".csv")):
                continue
            if full_only and not fname.endswith("_full.csv"):
                continue
            path = os.path.join(subdir, fname)

            # Infer model name from filename: dynamics_pred_vals_<model>.csv
            model = fname.replace("dynamics_pred_vals_", "").replace(".csv", "")
            model_map[model] = path

    return dict(sorted(model_map.items(), key=lambda kv: kv[0].lower()))


def generate_gt_dynamics_pred_vals(
    *,
    project_root: str,
    out_dir: str,
    features_csv: str,
    reference_obj: str,
    rscript: str = "Rscript",
    conda_env_r: Optional[str] = None,
) -> str:
    """
    Generate GT dynamics using y.csv (prAD, ABA, pseudotime) by calling `yifat_stuff/run_dynamics.R`.

    `y.csv` lacks pathology columns, so we rely on run_dynamics.R's merge-by-ID feature:
    it will pull the requested features from `reference_obj` (CSV/RDS/RData) by ID.
    """
    project_root = os.path.abspath(project_root)
    out_dir = os.path.abspath(out_dir)
    os.makedirs(out_dir, exist_ok=True)

    y_path = os.path.join(project_root, "prediction", "data", "y.csv")
    run_r = os.path.join(project_root, "yifat_stuff", "run_dynamics.R")
    if not os.path.isfile(y_path):
        raise FileNotFoundError(f"GT y.csv not found: {y_path}")
    if not os.path.isfile(run_r):
        raise FileNotFoundError(f"run_dynamics.R not found: {run_r}")
    if not os.path.isfile(reference_obj):
        raise FileNotFoundError(f"Reference object not found: {reference_obj}")

    gt_y = pd.read_csv(y_path, index_col=0)
    if "psuedotime" in gt_y.columns and "pseudotime" not in gt_y.columns:
        gt_y = gt_y.rename(columns={"psuedotime": "pseudotime"}) # awful typo in source data need to fix
    need = {"pseudotime", "prAD", "ABA"}
    missing = need.difference(gt_y.columns)
    if missing:
        raise ValueError(f"Missing columns {sorted(missing)} in {y_path}")

    gt_df = gt_y[["prAD", "ABA", "pseudotime"]].copy()
    gt_df.insert(0, "ID", gt_y.index.astype(str))

    with tempfile.TemporaryDirectory() as tmp:
        gt_in = os.path.join(tmp, "gt_dynamics_input.csv")
        gt_df.to_csv(gt_in, index=False)

        run_r_args = [
            run_r,
            gt_in,
            out_dir,
            "ground_truth",
            features_csv,
            os.path.abspath(reference_obj),
        ]

        # Prefer direct Rscript if available; otherwise optionally run via conda.
        cmd: List[str]
        if which(rscript) is not None:
            cmd = [rscript, *run_r_args]
        elif conda_env_r and which("conda") is not None:
            cmd = ["conda", "run", "-n", conda_env_r, "--no-capture-output", "Rscript", *run_r_args]
        else:
            raise SystemExit(
                "Could not find `Rscript` on PATH (needed for --generate-gt).\n"
                "- Option A: activate an environment that has R (so `Rscript` works)\n"
                "- Option B: pass `--conda-env-r dd_r` (or your R env name) so we run `conda run -n <env> Rscript ...`\n"
                "- Option C: skip generation and pass `--gt-pred-vals path/to/dynamics_pred_vals_*.csv`"
            )

        try:
            subprocess.run(cmd, check=True)
        except FileNotFoundError as e:
            raise SystemExit(
                f"Failed to execute R runner: {cmd[0]!r} not found.\n"
                "Make sure R is installed/activated (Rscript on PATH), or use --conda-env-r."
            ) from e

    return os.path.join(out_dir, "dynamics_pred_vals_ground_truth.csv")


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute DTW between GT and model dynamics curves.")
    parser.add_argument(
        "--project-root",
        default=os.path.abspath(os.path.join(os.path.dirname(__file__), "..")),
        help="Project root (contains prediction/ and yifat_stuff/).",
    )
    parser.add_argument(
        "--outputs-dir",
        default=os.path.join(os.path.dirname(__file__), "outputs"),
        help="benchmarking/outputs directory (contains dynamics_<model>/).",
    )
    parser.add_argument(
        "--gt-pred-vals",
        default=None,
        help="Path to GT dynamics_pred_vals CSV. If omitted, pass --generate-gt to create it.",
    )
    parser.add_argument(
        "--generate-gt",
        action="store_true",
        help="Generate GT dynamics_pred_vals by running yifat_stuff/run_dynamics.R on prediction/data/y.csv.",
    )
    parser.add_argument(
        "--features",
        default="sqrt.amyloid_mf,sqrt.tangles_mf,braaksc,cogng_demog_slope",
        help="Comma-separated feature list to DTW and (if generating GT) to request from reference object.",
    )
    parser.add_argument(
        "--reference-obj",
        default=os.path.join(
            os.path.abspath(os.path.join(os.path.dirname(__file__), "..")),
            "yifat_stuff",
            "bulk.prev.meta.bulk.df.csv",
        ),
        help="Reference object (CSV/RDS/RData) used by run_dynamics.R to merge missing features by ID.",
    )
    parser.add_argument(
        "--rscript",
        default="Rscript",
        help="Rscript executable (only used with --generate-gt).",
    )
    parser.add_argument(
        "--conda-env-r",
        default="dd_r",
        help="If Rscript is not on PATH, run via `conda run -n <env> Rscript ...` (only used with --generate-gt).",
    )
    parser.add_argument(
        "--full-only",
        action="store_true",
        default=True,
        help="Use only *_full dynamics_pred_vals files when discovering models (default: true).",
    )
    parser.add_argument(
        "--n-grid",
        type=int,
        default=50,
        help="Interpolation grid size per (feature, trajectory).",
    )
    parser.add_argument(
        "--out-csv",
        default=os.path.join(os.path.dirname(__file__), "figures", "dtw_dynamics.csv"),
        help="Output CSV path.",
    )
    args = parser.parse_args()

    features = [f.strip() for f in args.features.split(",") if f.strip()]
    os.makedirs(os.path.dirname(os.path.abspath(args.out_csv)), exist_ok=True)

    gt_path = args.gt_pred_vals
    if args.generate_gt:
        gt_out_dir = os.path.join(os.path.abspath(args.outputs_dir), "dynamics_ground_truth")
        gt_path = generate_gt_dynamics_pred_vals(
            project_root=args.project_root,
            out_dir=gt_out_dir,
            features_csv=",".join(features),
            reference_obj=args.reference_obj,
            rscript=args.rscript,
            conda_env_r=args.conda_env_r,
        )
    if not gt_path:
        raise SystemExit("Provide --gt-pred-vals or pass --generate-gt.")
    if not os.path.isfile(gt_path):
        raise SystemExit(f"GT pred vals CSV not found: {gt_path}")

    model_map = _discover_model_pred_vals(args.outputs_dir, full_only=bool(args.full_only))
    if not model_map:
        raise SystemExit(f"No model dynamics_pred_vals CSVs found under {args.outputs_dir}")

    df = compute_dtw_table(
        gt_pred_vals_csv=gt_path,
        model_pred_vals=model_map,
        features=features,
        n_grid=int(args.n_grid),
    )
    df.to_csv(args.out_csv, index=False, float_format="%.6f")
    print(f"Wrote {args.out_csv} ({len(df)} rows)")


if __name__ == "__main__":
    main()


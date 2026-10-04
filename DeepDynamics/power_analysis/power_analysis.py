#!/usr/bin/env python
"""Run APOE4-vs-APOE3/APOE2 power analyses along prAD and ABA trajectories."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from power_analysis.analysis import (
    analyze_cohort_scopes,
    analyze_incremental_samples,
    bootstrap_power_intervals,
)
from power_analysis.config import (
    DEFAULT_ALPHA,
    DEFAULT_BOOTSTRAP_REPETITIONS,
    DEFAULT_INCREMENT_FRACTIONS,
    DEFAULT_N_SIMULATIONS,
    DEFAULT_SEED,
    METADATA_PATH,
)
from power_analysis.data import load_metadata
from power_analysis.visualization import ensure_output_dirs, print_summary, save_figures, save_tables


def parse_increment_fractions(value: str) -> List[float]:
    """Parse fractions of unshared samples to add to the shared cohort."""
    fractions = sorted({float(item.strip()) for item in value.split(",") if item.strip()})
    invalid = [fraction for fraction in fractions if fraction < 0 or fraction > 1]
    if invalid:
        raise ValueError(f"Increment fractions must be between 0 and 1: {invalid}")
    return fractions


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Estimate one-sided APOE4 power for pathology means and "
            "pseudotime dynamics along prAD and ABA."
        )
    )
    parser.add_argument("--metadata", type=Path, default=METADATA_PATH, help="Path to metadata.csv.")
    parser.add_argument("--alpha", type=float, default=DEFAULT_ALPHA, help="One-sided alpha.")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED, help="Random seed.")
    parser.add_argument(
        "--n-simulations",
        type=int,
        default=DEFAULT_N_SIMULATIONS,
        help="Number of simulations per simulation-based power estimate.",
    )
    parser.add_argument(
        "--increment-fractions",
        default=",".join(str(value) for value in DEFAULT_INCREMENT_FRACTIONS),
        help="Comma-separated fractions of unshared samples to add to shared samples.",
    )
    parser.add_argument(
        "--bootstrap-repetitions",
        type=int,
        default=DEFAULT_BOOTSTRAP_REPETITIONS,
        help="Number of stratified bootstrap repetitions for power confidence intervals.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.n_simulations < 0:
        raise ValueError("--n-simulations must be non-negative.")
    if args.bootstrap_repetitions < 0:
        raise ValueError("--bootstrap-repetitions must be non-negative.")

    rng = np.random.default_rng(args.seed)
    increment_fractions = parse_increment_fractions(args.increment_fractions)
    metadata = load_metadata(args.metadata)
    print(f"Loaded {len(metadata)} APOE3/APOE4 trajectory-assigned rows from {args.metadata}")
    print(metadata.groupby(["cohort", "trajectory", "apoe_group"]).size().to_string())

    cohort_analytic, cohort_simulation = analyze_cohort_scopes(
        metadata,
        args.alpha,
        args.n_simulations,
        rng,
    )
    incremental_analytic, incremental_simulation = analyze_incremental_samples(
        metadata,
        increment_fractions,
        args.alpha,
        args.n_simulations,
        rng,
    )
    (
        cohort_analytic_ci,
        cohort_simulation_ci,
        incremental_analytic_ci,
        incremental_simulation_ci,
    ) = bootstrap_power_intervals(
        metadata,
        increment_fractions,
        args.alpha,
        args.n_simulations,
        args.bootstrap_repetitions,
        rng,
    )

    paths = ensure_output_dirs()
    save_tables(
        cohort_analytic,
        cohort_simulation,
        incremental_analytic,
        incremental_simulation,
        cohort_analytic_ci,
        cohort_simulation_ci,
        incremental_analytic_ci,
        incremental_simulation_ci,
        paths,
    )
    save_figures(
        cohort_analytic,
        cohort_simulation,
        incremental_analytic,
        incremental_simulation,
        cohort_analytic_ci,
        cohort_simulation_ci,
        incremental_analytic_ci,
        incremental_simulation_ci,
        paths,
    )
    print_summary(cohort_analytic, cohort_simulation)


if __name__ == "__main__":
    main()

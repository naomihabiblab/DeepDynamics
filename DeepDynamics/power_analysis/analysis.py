"""Analysis loops for cohort and incremental APOE power estimates."""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from power_analysis.config import (
    COMPARISONS,
    COHORT_SCOPES,
    EFFECTS,
    PATHOLOGIES,
    PATHOLOGY_ALTERNATIVES,
    TRAJECTORIES,
)
from power_analysis.data import cohort_subset, stratified_bootstrap_sample, stratified_sample_unshared
from power_analysis.models import (
    analytic_power,
    build_design,
    fit_ols,
    p_value,
    simulation_power,
)


def group_counts(df: pd.DataFrame, group_column: str, labels: tuple[str, str]) -> Dict[str, int]:
    ref_label, ind_label = labels
    counts = df[group_column].value_counts()
    return {
        f"n_{ref_label}": int(counts.get(ref_label, 0)),
        f"n_{ind_label}": int(counts.get(ind_label, 0)),
    }


def skipped_row(base: Dict[str, object], status: str) -> Dict[str, object]:
    return {
        **base,
        "estimate": np.nan,
        "standard_error": np.nan,
        "p_value": np.nan,
        "residual_sigma": np.nan,
        "residual_df": np.nan,
        "power": np.nan,
        "status": status,
    }


def analyze_subset(
    df: pd.DataFrame,
    pathology: str,
    trajectory: str,
    cohort_scope: str,
    sample_fraction_unshared: Optional[float],
    alpha: float,
    n_simulations: int,
    rng: np.random.Generator,
) -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    """Analyze one pathology/trajectory subset for both power methods."""
    subset = df[df["trajectory"] == trajectory].copy()
    subset[pathology] = pd.to_numeric(subset[pathology], errors="coerce")
    subset = subset.dropna(subset=[pathology, "pseudotime"])

    analytic_rows: List[Dict[str, object]] = []
    simulation_rows: List[Dict[str, object]] = []
    if subset.empty:
        return analytic_rows, simulation_rows

    for comparison_name, spec in COMPARISONS.items():
        comparison_subset = subset
        if spec["subset_query"]:
            comparison_subset = comparison_subset.query(spec["subset_query"]).copy()

        indicator_column = spec["indicator"]
        group_column = spec["group_column"]
        labels = spec["group_labels"]
        sidedness = spec["test_sidedness"]

        comparison_subset = comparison_subset.dropna(subset=[indicator_column, group_column])
        counts = group_counts(comparison_subset, group_column, labels)

        for effect_name in EFFECTS:
            # For APOE comparisons we keep the pathology-specific one-sided direction;
            # sex is analyzed two-sided by default.
            if sidedness == "two-sided":
                alternative = "two-sided"
            else:
                alternative = PATHOLOGY_ALTERNATIVES[pathology]

            include_interaction = effect_name == "dynamic_difference"
            coefficient_name = (
                f"{indicator_column}_x_pseudotime" if include_interaction else indicator_column
            )

            base = {
                "comparison": comparison_name,
                "indicator": indicator_column,
                "group_column": group_column,
                "pathology": pathology,
                "trajectory": trajectory,
                "cohort_scope": cohort_scope,
                "sample_fraction_unshared": sample_fraction_unshared,
                "effect": effect_name,
                "coefficient": coefficient_name,
                "test_alternative": alternative,
                "n": int(len(comparison_subset)),
                **counts,
            }

            ref_label, ind_label = labels
            if counts[f"n_{ref_label}"] < 2 or counts[f"n_{ind_label}"] < 2:
                row = skipped_row(base, "skipped_too_few_group_samples")
                analytic_rows.append(row.copy())
                simulation_rows.append(row.copy())
                continue

            x, names = build_design(
                comparison_subset,
                indicator_column=indicator_column,
                include_interaction=include_interaction,
            )
            result = fit_ols(comparison_subset[pathology].to_numpy(dtype=float), x)
            if result is None or coefficient_name not in names:
                row = skipped_row(base, "skipped_singular_or_underpowered_design")
                analytic_rows.append(row.copy())
                simulation_rows.append(row.copy())
                continue

            idx = names.index(coefficient_name)
            estimate = result.coefficients[idx]
            standard_error = result.standard_errors[idx]
            t_value = estimate / standard_error
            common = {
                **base,
                "estimate": float(estimate),
                "standard_error": float(standard_error),
                "p_value": p_value(t_value, result.residual_df, alternative),
                "residual_sigma": float(result.residual_sigma),
                "residual_df": int(result.residual_df),
                "status": "ok",
            }

            analytic_rows.append(
                {
                    **common,
                    "power": analytic_power(
                        estimate,
                        result.residual_sigma,
                        result.xtx_inv,
                        idx,
                        result.residual_df,
                        alpha,
                        alternative,
                    ),
                }
            )
            simulation_rows.append(
                {
                    **common,
                    "power": simulation_power(
                        x,
                        result.coefficients,
                        result.residual_sigma,
                        idx,
                        alpha,
                        alternative,
                        n_simulations,
                        rng,
                    ),
                }
            )

    return analytic_rows, simulation_rows


def analyze_dataset(
    df: pd.DataFrame,
    cohort_scope: str,
    sample_fraction_unshared: Optional[float],
    alpha: float,
    n_simulations: int,
    rng: np.random.Generator,
) -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    analytic_rows: List[Dict[str, object]] = []
    simulation_rows: List[Dict[str, object]] = []

    for trajectory in TRAJECTORIES:
        for pathology in PATHOLOGIES:
            analytic, simulation = analyze_subset(
                df,
                pathology,
                trajectory,
                cohort_scope,
                sample_fraction_unshared,
                alpha,
                n_simulations,
                rng,
            )
            analytic_rows.extend(analytic)
            simulation_rows.extend(simulation)

    return analytic_rows, simulation_rows


def analyze_cohort_scopes(
    metadata: pd.DataFrame,
    alpha: float,
    n_simulations: int,
    rng: np.random.Generator,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Compare shared, unshared, and all samples."""
    analytic_rows: List[Dict[str, object]] = []
    simulation_rows: List[Dict[str, object]] = []

    for scope in COHORT_SCOPES:
        analytic, simulation = analyze_dataset(
            cohort_subset(metadata, scope),
            scope,
            None,
            alpha,
            n_simulations,
            rng,
        )
        analytic_rows.extend(analytic)
        simulation_rows.extend(simulation)

    return pd.DataFrame(analytic_rows), pd.DataFrame(simulation_rows)


def analyze_incremental_samples(
    metadata: pd.DataFrame,
    fractions: Sequence[float],
    alpha: float,
    n_simulations: int,
    rng: np.random.Generator,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Start with shared rows and add increasing fractions of unshared rows."""
    analytic_rows: List[Dict[str, object]] = []
    simulation_rows: List[Dict[str, object]] = []
    shared = cohort_subset(metadata, "shared")
    unshared = cohort_subset(metadata, "unshared")

    for fraction in fractions:
        combined = pd.concat(
            [
                shared,
                stratified_sample_unshared(
                    unshared,
                    fraction,
                    rng,
                    stratify_columns=["trajectory", "apoe_group", "sex_group"],
                ),
            ],
            ignore_index=True,
        )
        analytic, simulation = analyze_dataset(
            combined,
            f"shared_plus_{fraction:.2f}_unshared",
            fraction,
            alpha,
            n_simulations,
            rng,
        )
        analytic_rows.extend(analytic)
        simulation_rows.extend(simulation)

    return pd.DataFrame(analytic_rows), pd.DataFrame(simulation_rows)


BOOTSTRAP_KEY_COLUMNS = [
    "comparison",
    "pathology",
    "trajectory",
    "effect",
    "cohort_scope",
    "sample_fraction_unshared",
]


def summarize_bootstrap_ci(bootstrap_rows: List[Dict[str, object]]) -> pd.DataFrame:
    """Aggregate bootstrap runs into percentile confidence intervals for power."""
    if not bootstrap_rows:
        return pd.DataFrame(columns=[*BOOTSTRAP_KEY_COLUMNS, "power_ci_lower", "power_ci_upper"])

    bootstrap_df = pd.DataFrame(bootstrap_rows)
    ok = bootstrap_df[bootstrap_df["status"] == "ok"].copy()
    if ok.empty:
        return pd.DataFrame(columns=[*BOOTSTRAP_KEY_COLUMNS, "power_ci_lower", "power_ci_upper"])

    grouped = (
        ok.groupby(BOOTSTRAP_KEY_COLUMNS, dropna=False)["power"]
        .agg(
            power_ci_lower=lambda values: float(np.quantile(values, 0.025)),
            power_ci_upper=lambda values: float(np.quantile(values, 0.975)),
        )
        .reset_index()
    )
    return grouped


def bootstrap_power_intervals(
    metadata: pd.DataFrame,
    fractions: Sequence[float],
    alpha: float,
    n_simulations: int,
    n_bootstrap_repetitions: int,
    rng: np.random.Generator,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Estimate power confidence intervals from stratified bootstrap resampling."""
    if n_bootstrap_repetitions <= 0:
        empty = pd.DataFrame(columns=[*BOOTSTRAP_KEY_COLUMNS, "power_ci_lower", "power_ci_upper"])
        return empty, empty.copy(), empty.copy(), empty.copy()

    cohort_analytic_bootstrap_rows: List[Dict[str, object]] = []
    cohort_simulation_bootstrap_rows: List[Dict[str, object]] = []
    incremental_analytic_bootstrap_rows: List[Dict[str, object]] = []
    incremental_simulation_bootstrap_rows: List[Dict[str, object]] = []

    shared = cohort_subset(metadata, "shared")
    unshared = cohort_subset(metadata, "unshared")

    for _ in range(n_bootstrap_repetitions):
        for scope in COHORT_SCOPES:
            scope_data = cohort_subset(metadata, scope)
            scope_bootstrap = stratified_bootstrap_sample(scope_data, rng)
            analytic, simulation = analyze_dataset(
                scope_bootstrap,
                scope,
                None,
                alpha,
                n_simulations,
                rng,
            )
            cohort_analytic_bootstrap_rows.extend(analytic)
            cohort_simulation_bootstrap_rows.extend(simulation)

        for fraction in fractions:
            combined = pd.concat(
                [
                    shared,
                    stratified_sample_unshared(
                        unshared,
                        fraction,
                        rng,
                        stratify_columns=["trajectory", "apoe_group", "sex_group"],
                    ),
                ],
                ignore_index=True,
            )
            combined_bootstrap = stratified_bootstrap_sample(combined, rng)
            analytic, simulation = analyze_dataset(
                combined_bootstrap,
                f"shared_plus_{fraction:.2f}_unshared",
                fraction,
                alpha,
                n_simulations,
                rng,
            )
            incremental_analytic_bootstrap_rows.extend(analytic)
            incremental_simulation_bootstrap_rows.extend(simulation)

    return (
        summarize_bootstrap_ci(cohort_analytic_bootstrap_rows),
        summarize_bootstrap_ci(cohort_simulation_bootstrap_rows),
        summarize_bootstrap_ci(incremental_analytic_bootstrap_rows),
        summarize_bootstrap_ci(incremental_simulation_bootstrap_rows),
    )

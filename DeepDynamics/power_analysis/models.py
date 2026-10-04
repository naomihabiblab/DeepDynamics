"""OLS fitting and power calculations."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats

from power_analysis.config import PSEUDOTIME_BIN_ORIGIN, PSEUDOTIME_BIN_WIDTH


@dataclass(frozen=True)
class OLSResult:
    coefficients: np.ndarray
    standard_errors: np.ndarray
    p_values: np.ndarray
    residual_sigma: float
    residual_df: int
    xtx_inv: np.ndarray


def build_design(
    df: pd.DataFrame,
    indicator_column: str,
    include_interaction: bool,
) -> Tuple[np.ndarray, List[str]]:
    """Build the APOE4 + pseudotime design matrix.

    Notes:
        - Pseudotime is first quantized to bins of width `PSEUDOTIME_BIN_WIDTH`
          (numeric, not categorical), then mean-centered.
        - The "dynamic difference" effect corresponds to the interaction term
          `{indicator_column}_x_pseudotime = indicator * centered(pseudotime_binned)`.
    """
    pseudotime = df["pseudotime"].to_numpy(dtype=float)
    if PSEUDOTIME_BIN_WIDTH and PSEUDOTIME_BIN_WIDTH > 0:
        # Quantize to a fixed 0.1 grid (numeric), anchored at PSEUDOTIME_BIN_ORIGIN.
        pseudotime = (
            np.round((pseudotime - PSEUDOTIME_BIN_ORIGIN) / PSEUDOTIME_BIN_WIDTH) * PSEUDOTIME_BIN_WIDTH
            + PSEUDOTIME_BIN_ORIGIN
        )
    pseudotime = pseudotime - np.mean(pseudotime)
    indicator = df[indicator_column].to_numpy(dtype=float)

    columns = [np.ones(len(df)), indicator, pseudotime]
    names = ["intercept", indicator_column, "pseudotime"]
    if include_interaction:
        columns.append(indicator * pseudotime)
        names.append(f"{indicator_column}_x_pseudotime")

    return np.column_stack(columns), names


def fit_ols(y: np.ndarray, x: np.ndarray) -> Optional[OLSResult]:
    """Fit OLS with numpy and return coefficient uncertainty."""
    n_samples, n_parameters = x.shape
    residual_df = n_samples - n_parameters
    if residual_df <= 0 or n_samples == 0:
        return None

    try:
        xtx_inv = np.linalg.inv(x.T @ x)
    except np.linalg.LinAlgError:
        return None

    coefficients = xtx_inv @ x.T @ y
    residuals = y - x @ coefficients
    residual_variance = float(residuals.T @ residuals) / residual_df
    if not np.isfinite(residual_variance) or residual_variance <= 0:
        return None

    residual_sigma = math.sqrt(residual_variance)
    standard_errors = np.sqrt(np.diag(xtx_inv) * residual_variance)
    with np.errstate(divide="ignore", invalid="ignore"):
        t_values = coefficients / standard_errors
    p_values = 2 * stats.t.sf(np.abs(t_values), residual_df)

    return OLSResult(
        coefficients=coefficients,
        standard_errors=standard_errors,
        p_values=p_values,
        residual_sigma=residual_sigma,
        residual_df=residual_df,
        xtx_inv=xtx_inv,
    )


def p_value(t_value: float, residual_df: int, alternative: str) -> float:
    """Return a p-value under the requested alternative.

    Supported alternatives:
        - "greater": one-sided, \(H_1: \\beta > 0\)
        - "less": one-sided, \(H_1: \\beta < 0\)
        - "two-sided": two-sided, \(H_1: \\beta \\neq 0\)
    """
    if alternative == "greater":
        return float(stats.t.sf(t_value, residual_df))
    if alternative == "less":
        return float(stats.t.cdf(t_value, residual_df))
    if alternative == "two-sided":
        return float(2 * stats.t.sf(abs(t_value), residual_df))
    raise ValueError(f"Unknown alternative: {alternative}")


def analytic_power(
    coefficient: float,
    residual_sigma: float,
    xtx_inv: np.ndarray,
    coefficient_index: int,
    residual_df: int,
    alpha: float,
    alternative: str,
) -> float:
    """Power from the noncentral t distribution for the fitted design."""
    variance_multiplier = xtx_inv[coefficient_index, coefficient_index]
    if residual_df <= 0 or residual_sigma <= 0 or variance_multiplier <= 0:
        return np.nan

    standard_error = residual_sigma * math.sqrt(variance_multiplier)
    noncentrality = coefficient / standard_error
    critical = stats.t.ppf(1 - alpha, residual_df)

    if alternative == "greater":
        return float(stats.nct.sf(critical, residual_df, noncentrality))
    if alternative == "less":
        return float(stats.nct.cdf(-critical, residual_df, noncentrality))
    if alternative == "two-sided":
        critical_two_sided = stats.t.ppf(1 - alpha / 2, residual_df)
        return float(
            stats.nct.sf(critical_two_sided, residual_df, noncentrality)
            + stats.nct.cdf(-critical_two_sided, residual_df, noncentrality)
        )
    raise ValueError(f"Unknown alternative: {alternative}")


def simulation_power(
    x: np.ndarray,
    coefficients: np.ndarray,
    residual_sigma: float,
    coefficient_index: int,
    alpha: float,
    alternative: str,
    n_simulations: int,
    rng: np.random.Generator,
) -> float:
    """Estimate power by simulating outcomes from the fitted OLS model."""
    n_samples, n_parameters = x.shape
    residual_df = n_samples - n_parameters
    if n_simulations <= 0 or residual_sigma <= 0 or residual_df <= 0:
        return np.nan

    try:
        xtx_inv = np.linalg.inv(x.T @ x)
    except np.linalg.LinAlgError:
        return np.nan

    expected = x @ coefficients
    simulated_y = expected[:, None] + rng.normal(
        0,
        residual_sigma,
        size=(n_samples, n_simulations),
    )
    beta_hat = xtx_inv @ x.T @ simulated_y
    residuals = simulated_y - x @ beta_hat
    residual_variance = np.sum(residuals**2, axis=0) / residual_df
    se = np.sqrt(np.maximum(residual_variance, 0) * xtx_inv[coefficient_index, coefficient_index])

    with np.errstate(divide="ignore", invalid="ignore"):
        t_values = beta_hat[coefficient_index] / se

    if alternative == "greater":
        p_values = stats.t.sf(t_values, residual_df)
    elif alternative == "less":
        p_values = stats.t.cdf(t_values, residual_df)
    elif alternative == "two-sided":
        p_values = 2 * stats.t.sf(np.abs(t_values), residual_df)
    else:
        raise ValueError(f"Unknown alternative: {alternative}")
    return float(np.mean(p_values < alpha))

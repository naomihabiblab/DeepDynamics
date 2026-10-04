"""Data loading and cohort construction for the APOE power analysis."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from power_analysis.config import (
    PATHOLOGIES,
    PROBABILITY_THRESHOLD,
    PSEUDOTIME_THRESHOLD,
    REQUIRED_COLUMNS,
)


def require_columns(df: pd.DataFrame, columns: Iterable[str], source: Path) -> None:
    """Validate that the metadata file matches the expected schema."""
    missing = [column for column in columns if column not in df.columns]
    if missing:
        raise ValueError(f"{source} is missing required columns: {missing}")


def assign_trajectory(df: pd.DataFrame) -> pd.Series:
    """Assign subjects to healthy, prAD, or ABA using the project thresholds."""
    trajectory = pd.Series(pd.NA, index=df.index, dtype="object")
    healthy = df["pseudotime"] < PSEUDOTIME_THRESHOLD
    trajectory.loc[healthy] = "healthy"
    trajectory.loc[(~healthy) & (df["prAD"] > PROBABILITY_THRESHOLD)] = "prAD"
    trajectory.loc[(~healthy) & (df["ABA"] > PROBABILITY_THRESHOLD)] = "ABA"
    return trajectory


def assign_apoe_group(df: pd.DataFrame) -> pd.Series:
    """Create the APOE contrast used by the models.

    APOE4 excludes APOE2 carriers. The comparison group combines APOE3 and
    APOE2 subjects as long as they are not APOE4 carriers.
    """
    group = pd.Series(pd.NA, index=df.index, dtype="object")
    apoe4 = (df["apoe_4"] == 1) & (df["apoe_2"] != 1)
    apoe3_or_2 = (df["apoe_4"] == 0) & ((df["apoe_3"] == 1) | (df["apoe_2"] == 1))
    group.loc[apoe3_or_2] = "APOE3"
    group.loc[apoe4] = "APOE4"
    return group


def load_metadata(path: Path) -> pd.DataFrame:
    """Load metadata and keep rows usable for APOE trajectory contrasts."""
    metadata = pd.read_csv(path)
    require_columns(metadata, REQUIRED_COLUMNS, path)

    numeric_columns = [
        "ID",
        "apoe_2",
        "apoe_3",
        "apoe_4",
        "msex",
        "prAD",
        "ABA",
        "pseudotime",
        *PATHOLOGIES,
    ]
    metadata[numeric_columns] = metadata[numeric_columns].apply(
        pd.to_numeric,
        errors="coerce",
    )

    metadata["trajectory"] = assign_trajectory(metadata)
    metadata["apoe_group"] = assign_apoe_group(metadata)
    metadata["apoe4"] = (metadata["apoe_group"] == "APOE4").astype(int)
    metadata["msex"] = pd.to_numeric(metadata["msex"], errors="coerce")
    metadata["sex_group"] = metadata["msex"].map({0: "female", 1: "male"})
    metadata["male"] = (metadata["msex"] == 1).astype(int)

    needed = ["ID", "trajectory", "apoe_group", "pseudotime", "sex_group", "male"]
    return (
        metadata[metadata["cohort"].isin(["shared", "unshared"])]
        .dropna(subset=needed)
        .copy()
    )


def cohort_subset(df: pd.DataFrame, scope: str) -> pd.DataFrame:
    """Return shared, unshared, or all rows."""
    if scope == "all":
        return df.copy()
    if scope in {"shared", "unshared"}:
        return df[df["cohort"] == scope].copy()
    raise ValueError(f"Unknown cohort scope: {scope}")


def stratified_sample_unshared(
    unshared: pd.DataFrame,
    fraction: float,
    rng: np.random.Generator,
    stratify_columns: list[str] | None = None,
) -> pd.DataFrame:
    """Sample a fraction of unshared rows within trajectory/APOE strata."""
    if fraction <= 0 or unshared.empty:
        return unshared.iloc[0:0].copy()
    if fraction >= 1:
        return unshared.copy()

    sampled = []
    stratify_columns = stratify_columns or ["trajectory", "apoe_group", "sex_group"]
    for _, stratum in unshared.groupby(stratify_columns, dropna=False):
        n = int(round(len(stratum) * fraction))
        if n:
            sampled.append(stratum.loc[rng.choice(stratum.index, size=n, replace=False)])

    return pd.concat(sampled, ignore_index=False) if sampled else unshared.iloc[0:0].copy()


def stratified_bootstrap_sample(
    df: pd.DataFrame,
    rng: np.random.Generator,
    stratify_columns: list[str] | None = None,
) -> pd.DataFrame:
    """Draw a within-stratum bootstrap sample with replacement."""
    if df.empty:
        return df.copy()

    sampled = []
    stratify_columns = stratify_columns or ["trajectory", "apoe_group", "sex_group"]
    for _, stratum in df.groupby(stratify_columns, dropna=False):
        sampled_idx = rng.choice(stratum.index.to_numpy(), size=len(stratum), replace=True)
        sampled.append(stratum.loc[sampled_idx].copy())

    return pd.concat(sampled, ignore_index=True)

"""Shared configuration for APOE power analyses."""

from __future__ import annotations

from pathlib import Path


PROBABILITY_THRESHOLD = 0.5
PSEUDOTIME_THRESHOLD = 0.1
PSEUDOTIME_BIN_WIDTH = 0.1
PSEUDOTIME_BIN_ORIGIN = 0.1
DEFAULT_ALPHA = 0.05
DEFAULT_N_SIMULATIONS = 1000
DEFAULT_BOOTSTRAP_REPETITIONS = 0 #1000
DEFAULT_SEED = 0
DEFAULT_INCREMENT_FRACTIONS = (0.0, 0.1, 0.25, 0.5, 0.75, 1.0)

SCRIPT_DIR = Path(__file__).resolve().parent
BASE_DIR = SCRIPT_DIR.parent
METADATA_PATH = BASE_DIR / "metadata.csv"
OUTPUT_ROOT = SCRIPT_DIR / "outputs"
FIGURE_ROOT = SCRIPT_DIR / "figures"
# Vector formats for Illustrator: TrueType text in PDF, editable text in SVG.
FIGURE_VECTOR_FORMATS = ("pdf", "svg")

PATHOLOGIES = ("amyloid_mf", "tangles_mf", "cogng_demog_slope")
TRAJECTORIES = ("prAD", "ABA")
COHORT_SCOPES = ("shared", "unshared", "all")

EFFECTS = ("mean_difference", "dynamic_difference")

# Which comparisons to run. Each comparison is defined by:
# - indicator: binary 0/1 column used in the regression (1 = group of interest)
# - group_column: categorical label column used for counts/stratification
# - group_labels: (reference_label, indicator_label) for reporting
# - subset_query: optional pandas query string applied before analysis
COMPARISONS = {
    # Baseline: APOE4 carriers vs APOE3/APOE2 (non-APOE4) across all sexes.
    "apoe4_vs_noncarrier": {
        "indicator": "apoe4",
        "group_column": "apoe_group",
        "group_labels": ("APOE3", "APOE4"),
        "subset_query": None,
        "test_sidedness": "one-sided",
    },
    # Sex effect: males vs females.
    "male_vs_female": {
        "indicator": "male",
        "group_column": "sex_group",
        "group_labels": ("female", "male"),
        "subset_query": None,
        "test_sidedness": "one-sided",
    },
    # Female-only APOE4 vs noncarrier.
    "female_apoe4_vs_noncarrier": {
        "indicator": "apoe4",
        "group_column": "apoe_group",
        "group_labels": ("APOE3", "APOE4"),
        "subset_query": "sex_group == 'female'",
        "test_sidedness": "one-sided",
    },
}

# One-sided alternatives for the APOE4 coefficient.
PATHOLOGY_ALTERNATIVES = {
    "amyloid_mf": "greater",
    "tangles_mf": "greater",
    "cogng_demog_slope": "less",
}

REQUIRED_COLUMNS = (
    "ID",
    "apoe_2",
    "apoe_3",
    "apoe_4",
    "msex",
    "prAD",
    "ABA",
    "pseudotime",
    "cohort",
    *PATHOLOGIES,
)

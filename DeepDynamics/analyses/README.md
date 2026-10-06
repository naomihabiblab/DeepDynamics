# DeepDynamics analyses (synthetic method samples)

Short runnable demos of the paper's downstream methods on the synthetic cohort
under `../prediction/data/synthetic/`. Outputs are **not** scientific results.

**Disclosure:** inputs are simulated (`donor_meta` IDs `80000001+`, toy gene
counts, toy proteomics). Extras are derived only from the synthetic tutorial
`500.h5ad`. No study donor matrices or `data/sync/` paths are read. Public
gene/cluster names may appear as labels; abundances are RNG draws.

There is no Harmony or Seurat-anchor integration. When two cohorts are present,
`batch` / `dataset` is a covariate in the merged models (levels `ROSMAP` /
`cuimc2`).

## Setup

Synthetic cohort and analysis extras are already committed under
`../prediction/data/synthetic/`. From `DeepDynamics/analyses/`:

```bash
Rscript notebooks/01_pathway_gam_dataset_batch.R
Rscript notebooks/02_deg_and_fig5_heatmap.R
Rscript notebooks/03_trait_proteomics.R
```

To knit HTML with figures inline (from `analyses/`):

```bash
Rscript -e 'rmarkdown::render("notebooks/01_pathway_gam_dataset_batch.Rmd")'
Rscript -e 'rmarkdown::render("notebooks/02_deg_and_fig5_heatmap.Rmd")'
Rscript -e 'rmarkdown::render("notebooks/03_trait_proteomics.Rmd")'
```

Open the resulting `.html` files to review plots without re-running. Synthetic
extras (beyond the tutorial `500.h5ad`): `donor_meta.csv`, `scrna_*.csv`,
`proteomics_*.csv`, `pathway_genes.csv`. See
`../prediction/data/synthetic/ANALYSIS_EXTRAS_README.md`.

## Notebooks

Grouped by analysis family (three samples). See
[notebooks/README.md](notebooks/README.md) for Method alignment and how to
run / knit.

| Sample | Method |
|--------|--------|
| `01_pathway_gam_dataset_batch` | **Dynamics.** Part A: merge → `AddModuleScore` → donor mean → GAM `mean_scaled ~ s(pseudotime, k=7) + batch` (levels `ROSMAP`/`cuimc2`, ribbon at `ROSMAP`) + ANOVA with `k=7`. Part B: original-cohort `Ast.10` and `sqrt.amyloid_mf` GAMs (`s(pseudotime)` / `apoe_4 + s(pseudotime, by=apoe_4)`; no batch, no `k`) |
| `02_deg_and_fig5_heatmap` | **Merged DEGs.** Map `batchA`/`batchB` → `dataset` `ROSMAP`/`cuimc2`; healthy-window edgeR `~ apoe_4 + dataset` with `glmQLFTest(..., coef=2)`; Poisson `FindMarkers` (`mean.fxn`, `logfc=0.25`, `min.pct=0.1`, latents `dataset.num`+`projid.num`+`sex_num`); Fig 5-style volcano |
| `03_trait_proteomics` | **Trait association + proteomics.** `trait ~ covariate + age_death + pmi + RIN` (pathologies vs sex/APOE, then vs paper `SIG.CLUSTERS` states) with Fig 3-style signed-FDR heatmaps; proteomics density Wilcoxon (`23`/`33` vs `34`/`44`). No ANOVA |

Each sample has a `.R` CLI script, an `.Rmd` source, and a knitted `.html`
(with figures). Prefer the HTML for viewing; use `.R` / `.Rmd` to re-run.

## Packages

Beyond `dd_r`: Seurat, edgeR, SummarizedExperiment, ggrepel, ComplexHeatmap,
reshape2, gridExtra, rmarkdown, knitr.

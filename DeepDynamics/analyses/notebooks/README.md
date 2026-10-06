# Method sample notebooks

See the parent [analyses/README.md](../README.md) for setup. Synthetic inputs
only — outputs are **not** paper results.

| Sample | Family | View |
|--------|--------|------|
| `01_pathway_gam_dataset_batch` | Dynamics | `.html` / `.Rmd` / `.R` |
| `02_deg_and_fig5_heatmap` | Merged DEGs | `.html` / `.Rmd` / `.R` |
| `03_trait_proteomics` | Trait + proteomics | `.html` / `.Rmd` / `.R` |

Open the `.html` files to see figures. Re-run from `DeepDynamics/analyses/`:

```bash
Rscript notebooks/01_pathway_gam_dataset_batch.R
Rscript notebooks/02_deg_and_fig5_heatmap.R
Rscript notebooks/03_trait_proteomics.R

# Or knit HTML
Rscript -e 'rmarkdown::render("notebooks/01_pathway_gam_dataset_batch.Rmd")'
Rscript -e 'rmarkdown::render("notebooks/02_deg_and_fig5_heatmap.Rmd")'
Rscript -e 'rmarkdown::render("notebooks/03_trait_proteomics.Rmd")'
```

PDFs/CSVs under `outputs/` are local/gitignored demos and are not paper results.

## Method alignment

| Sample | What it fits | Cohort covariate |
|--------|--------------|------------------|
| 01 Part A | Pathway GAM `mean_scaled ~ s(pseudotime, k=7) + batch`; ANOVA with `k=7`; ribbon at batch `ROSMAP` | `batch` (`ROSMAP` / `cuimc2`) |
| 01 Part B | `Ast.10` / `sqrt.amyloid_mf` GAMs by APOE4 (no `k`, no batch) | none |
| 02 edgeR | Healthy-window `~ apoe_4 + dataset`; `glmQLFTest(..., coef = 2)` | `dataset` (`ROSMAP` / `cuimc2`) |
| 02 Poisson | `FindMarkers` (`mean.fxn`, `logfc=0.25`, `min.pct=0.1`) | `dataset.num` + `projid.num` + `sex_num` |
| 03 trait | `trait ~ covariate + age_death + pmi + RIN`; signed-FDR heatmaps | technical covariates only |
| 03 proteomics | Density Wilcoxon (`23`/`33` vs `34`/`44`) | genotype bins |

Synthetic `batchA` / `batchB` map to paper cohort levels `ROSMAP` / `cuimc2`.

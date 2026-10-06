# Synthetic analysis extras

Synthetic analysis extras. Donor covariates, gene counts, and proteomics values
are simulated for method samples and are not study results.

Donor IDs match the tutorial synthetic cohort (`80000001+`). Gene counts and
proteomics abundances are independent draws.

**Disclosure:** these files are simulated. They are not study donor matrices, do
not contain ROSMAP participant identifiers, and are derived only from the
synthetic tutorial `500.h5ad`. Public gene and cluster *names* may appear as
labels; abundances are RNG draws.

| File | Role |
|------|------|
| `donor_meta.csv` | Covariates, pathologies, trajectories, dataset, SIG.CLUSTERS |
| `scrna_counts.csv` / `scrna_cells.csv` | Toy gene-level cells |
| `pathway_genes.csv` | Gene set for AddModuleScore sample |
| `proteomics_*.csv` | Toy proteomics assay + metadata |

Method samples: `../../analyses/notebooks/`.

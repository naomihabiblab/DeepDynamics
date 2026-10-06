#!/usr/bin/env Rscript
# Sample notebook 02 — Merged DEGs
#   edgeR: ~ apoe_4 + dataset, glmQLFTest(..., coef = 2)
#   Poisson: latent.vars = dataset.num + projid.num + sex_num
# Synthetic data only — not paper results.
#
# Run from DeepDynamics/analyses/:
#   Rscript notebooks/02_deg_and_fig5_heatmap.R

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(edgeR)
  library(ggplot2)
  library(ggrepel)
})

# --- Paths and healthy-window subset ------------------------------------------
# Healthy window: pseudotime < 0.1 (same cut used for Fig 5 DEG tables).
# Prefer sibling of analyses/ (.. /prediction/...); fall back to old sibling layout.
.root_candidates <- c(
  file.path("..", "prediction", "data", "synthetic"),
  file.path("..", "DeepDynamics", "prediction", "data", "synthetic"),
  file.path("..", "..", "DeepDynamics", "prediction", "data", "synthetic")
)
root <- NULL
for (cand in .root_candidates) {
  if (dir.exists(cand)) {
    root <- normalizePath(cand)
    break
  }
}
if (is.null(root)) stop("Synthetic data directory not found. Tried: ", paste(.root_candidates, collapse = ", "))
message("Synthetic root: ", root)

E4POS_COL <- "#ef6c00"
DS.LEVELS <- c("ROSMAP", "cuimc2")

.map_to_paper_dataset <- function(x) {
  x <- as.character(x)
  dplyr::case_when(
    x %in% c("batchA", "500", "ROSMAP") ~ "ROSMAP",
    x %in% c("batchB", "cuimc2") ~ "cuimc2",
    TRUE ~ x
  )
}

counts <- as.matrix(read.csv(file.path(root, "scrna_counts.csv"), row.names = 1, check.names = FALSE))
cells <- read.csv(file.path(root, "scrna_cells.csv"), row.names = 1, check.names = FALSE) %>%
  filter(is.finite(pseudotime), pseudotime < 0.1) %>%
  mutate(dataset = factor(.map_to_paper_dataset(dataset), levels = DS.LEVELS))
counts <- counts[rownames(cells), , drop = FALSE]
message(
  "Cells in healthy window: ", nrow(cells),
  "; donors: ", dplyr::n_distinct(cells$projid)
)
stopifnot(nlevels(droplevels(cells$dataset)) >= 2L)

# --- Merge two cohort objects (paper: ROSMAP + cuimc2 → dataset covariate) -----
# Synthetic batchA/batchB map to dataset levels "ROSMAP" / "cuimc2".
# No Harmony/anchors; dataset is the covariate in edgeR and Poisson FindMarkers.
obj_list <- lapply(DS.LEVELS, function(ds) {
  idx <- rownames(cells)[as.character(cells$dataset) == ds]
  o <- CreateSeuratObject(
    counts = t(counts[idx, , drop = FALSE]),
    meta.data = cells[idx, , drop = FALSE],
    project = ds
  )
  o$dataset <- ds
  o
})
names(obj_list) <- DS.LEVELS
obj <- merge(
  obj_list[[1L]],
  y = obj_list[-1L],
  add.cell.ids = names(obj_list),
  project = "merged_synthetic"
)
obj$apoe_4 <- factor(as.character(obj$apoe_4), levels = c("0", "1"))
obj$dataset <- factor(as.character(obj$dataset), levels = DS.LEVELS)
# Latents for run.2group.poisson.merged:
#   dataset.num = 1 if cuimc2; projid.num = integer factor; sex_num = 1 if male
obj$dataset.num <- as.integer(as.character(obj$dataset) == "cuimc2")
obj$projid.num <- as.integer(factor(as.character(obj$projid)))
obj$sex_num <- as.integer(tolower(as.character(obj$sex)) == "male")
DefaultAssay(obj) <- "RNA"
# Seurat v5 keeps per-object layers after merge; join before Normalize/FindMarkers.
obj <- tryCatch(JoinLayers(obj), error = function(e) obj)
obj <- NormalizeData(obj, verbose = FALSE)
obj <- tryCatch(JoinLayers(obj), error = function(e) obj)
message(
  "Merged Seurat object: ", ncol(obj), " cells from ",
  paste(DS.LEVELS, collapse = " + "), " (no integration)"
)

# --- Pseudobulk edgeR (run_pseudobulk_2groups_edger_merged) --------------------
# Aggregate counts per donor, then edgeR QL:
#   design ~ apoe_4 + dataset
#   glmQLFTest(..., coef = 2)  # apoe_41; positive logFC = higher in APOE4
# No |logFC| drop before saving; volcano applies FDR and |logFC| for display only.
pb_list <- AggregateExpression(
  obj,
  assays = "RNA",
  group.by = "projid",
  return.seurat = FALSE
)
pb <- as.matrix(pb_list$RNA)
colnames(pb) <- sub("^g", "", colnames(pb))

donor_anno <- slot(obj, "meta.data") %>%
  group_by(projid) %>%
  summarise(
    apoe_4 = dplyr::first(apoe_4),
    dataset = dplyr::first(dataset),
    .groups = "drop"
  ) %>%
  mutate(
    projid = as.character(projid),
    apoe_4 = factor(as.character(apoe_4), levels = c("0", "1")),
    dataset = factor(as.character(dataset), levels = DS.LEVELS)
  )
donor_anno <- donor_anno[match(colnames(pb), donor_anno$projid), ]
stopifnot(!any(is.na(donor_anno$projid)))
stopifnot(all(donor_anno$projid == colnames(pb)))

keep_valid <- !is.na(donor_anno$apoe_4)
pb <- pb[, keep_valid, drop = FALSE]
donor_anno <- donor_anno[keep_valid, , drop = FALSE]

y <- DGEList(counts = pb, group = donor_anno$apoe_4)
design <- model.matrix(~ apoe_4 + dataset, data = donor_anno)
keep <- filterByExpr(y, design)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2)
tt <- topTags(qlf, n = Inf)$table
e4.10 <- tt %>%
  tibble::rownames_to_column("gene") %>%
  transmute(
    gene,
    avg_log2FC = logFC,
    p_val = PValue,
    p_val_adj = FDR
  )

# --- Poisson companion (run.2group.poisson.merged) -----------------------------
# Same Seurat FindMarkers settings as the paper 2-group Poisson:
#   mean.fxn = log2(rowMeans + 1e-6), logfc.threshold = 0.25, min.pct = 0.1,
#   recorrect_umi = FALSE, ident.1 = 1 vs ident.2 = 0,
#   latent.vars = c("dataset.num", "projid.num", "sex_num").
Idents(obj) <- "apoe_4"
mean.fxn <- function(x) log(rowMeans(x) + 1e-6, base = 2)
poisson_res <- FindMarkers(
  obj,
  ident.1 = 1,
  ident.2 = 0,
  test.use = "poisson",
  mean.fxn = mean.fxn,
  latent.vars = c("dataset.num", "projid.num", "sex_num"),
  logfc.threshold = 0.25,
  min.pct = 0.1,
  recorrect_umi = FALSE,
  verbose = FALSE
)
pois_df <- poisson_res %>%
  tibble::rownames_to_column("gene") %>%
  transmute(
    gene,
    avg_log2FC = avg_log2FC,
    p_val = p_val,
    p_val_adj = p_val_adj,
    method = "poisson"
  )

out_dir <- file.path("notebooks", "outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(
  bind_rows(e4.10 %>% mutate(method = "edgeR"), pois_df),
  file.path(out_dir, "02_deg_results.csv"),
  row.names = FALSE
)

# --- Fig 5 volcano (.volcano_panel) --------------------------------------------
# x = avg_log2FC, y = -log10(p_val_adj). Genes pass "sig" when FDR < 0.05 and
# |logFC| >= 0.25. Pathway genes are highlighted like the IFN panel in Fig 5.
highlight_genes <- read.csv(file.path(root, "pathway_genes.csv"))$gene
p_cutoff <- 0.05
fc_cutoff <- 0.25

df <- e4.10 %>%
  mutate(
    sig_fc = !is.na(p_val_adj) & p_val_adj < p_cutoff &
      !is.na(avg_log2FC) & abs(avg_log2FC) >= fc_cutoff,
    is_ifn = gene %in% highlight_genes,
    color_grp = case_when(
      is_ifn & sig_fc ~ "IFN",
      sig_fc ~ "Significant",
      TRUE ~ "NS"
    ),
    neg_log10_p = pmin(-log10(pmax(p_val_adj, .Machine$double.xmin)), 50)
  )

color_vals <- c(
  "IFN" = E4POS_COL,
  "Significant" = "black",
  "NS" = "grey70"
)
size_vals <- c("IFN" = 0.9, "Significant" = 0.55, "NS" = 0.25)

ifn_lab <- df %>% filter(is_ifn, sig_fc) %>% arrange(p_val_adj) %>% slice_head(n = 4L) %>% pull(gene)
top_up <- df %>% filter(sig_fc, !is_ifn, avg_log2FC > 0) %>% arrange(p_val_adj) %>% slice_head(n = 2L) %>% pull(gene)
top_down <- df %>% filter(sig_fc, !is_ifn, avg_log2FC < 0) %>% arrange(p_val_adj) %>% slice_head(n = 2L) %>% pull(gene)
label_set <- unique(c(ifn_lab, top_up, top_down))
df$label <- ifelse(df$gene %in% label_set, df$gene, "")

p <- ggplot(df, aes(x = avg_log2FC, y = neg_log10_p)) +
  geom_point(aes(color = color_grp, size = color_grp), alpha = 0.6) +
  geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", linewidth = 0.2, color = "grey40") +
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", linewidth = 0.2, color = "grey40") +
  ggrepel::geom_text_repel(
    data = dplyr::filter(df, nzchar(label)),
    aes(label = label, color = color_grp),
    size = 2.2,
    segment.size = 0.2,
    max.overlaps = 12,
    show.legend = FALSE,
    box.padding = 0.25,
    point.padding = 0.15,
    seed = 1
  ) +
  scale_color_manual(values = color_vals, name = NULL) +
  scale_size_manual(values = size_vals, guide = "none") +
  labs(
    title = "Toy healthy-window DEGs (synthetic)",
    subtitle = "Pseudobulk · healthy window · edgeR ~ apoe_4 + dataset (coef = 2)",
    x = "log2 FC (APOE4+ vs APOE4-)",
    y = "-log10 padj"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.25)
  )

ggsave(file.path(out_dir, "02_deg_volcano.pdf"), p, width = 5, height = 4.5)
old_hm <- file.path(out_dir, "02_fig5_style_heatmap.pdf")
if (file.exists(old_hm)) file.remove(old_hm)

message("Wrote ", file.path(out_dir, "02_deg_results.csv"))
message("Wrote ", file.path(out_dir, "02_deg_volcano.pdf"))
message("edgeR significant (FDR<0.05 & |logFC|>=0.25): ", sum(df$sig_fc, na.rm = TRUE))
message("poisson significant (FDR<0.05): ", sum(pois_df$p_val_adj < 0.05, na.rm = TRUE))
message("Done notebook 02")

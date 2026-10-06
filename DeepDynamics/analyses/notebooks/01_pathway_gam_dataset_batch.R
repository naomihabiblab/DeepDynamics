#!/usr/bin/env Rscript
# Sample notebook 01 — Dynamics (two parts)
#   Part A: merged pathway (k=7 + batch) — fit_split_dynamics_batch / ED Stage-3
#   Part B: original-cohort cell state + pathology (no k, no batch) — ANOVA.dyn.R
# Synthetic data only — not paper results.
#
# Run from DeepDynamics/analyses/:
#   Rscript notebooks/01_pathway_gam_dataset_batch.R

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(mgcv)
  library(ggplot2)
  library(patchwork)
})

# --- Paths and constants -------------------------------------------------------
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

PT.MAX <- 0.34
PT.HEALTHY <- 0.1
GAM.K <- 7L
E4POS_COL <- "#ef6c00"
E4NEG_COL <- "#4db6ac"
BATCH.LEVELS <- c("ROSMAP", "cuimc2")

out_dir <- file.path("notebooks", "outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

.map_to_paper_batch <- function(x) {
  x <- as.character(x)
  dplyr::case_when(
    x %in% c("batchA", "500", "ROSMAP") ~ "ROSMAP",
    x %in% c("batchB", "cuimc2") ~ "cuimc2",
    TRUE ~ x
  )
}

.fmt_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 0.001) "p < 0.001" else paste0("p = ", formatC(p, format = "f", digits = 3))
}

.dynamics_panel <- function(pv, feat_pos, feat_neg, anova_p, title, subtitle, ylab) {
  pv <- pv %>%
    filter(feature %in% c(feat_pos, feat_neg), is.finite(x), x >= 0, x <= PT.MAX) %>%
    mutate(feature = factor(feature, levels = c(feat_pos, feat_neg)))
  cols <- setNames(c(E4POS_COL, E4NEG_COL), c(feat_pos, feat_neg))
  labs_map <- setNames(c("APOE4+", "APOE4-"), c(feat_pos, feat_neg))
  ggplot(pv, aes(x = x, y = fit, color = feature, fill = feature)) +
    annotate(
      "rect", xmin = 0, xmax = PT.HEALTHY, ymin = -Inf, ymax = Inf,
      fill = "grey90", alpha = 0.5, colour = NA
    ) +
    geom_ribbon(
      aes(ymin = fit - se.fit, ymax = fit + se.fit),
      alpha = 0.18, linetype = "dashed", linewidth = 0.25, show.legend = FALSE
    ) +
    geom_line(linewidth = 0.45) +
    annotate(
      "text", x = PT.MAX, y = -Inf, label = .fmt_p(anova_p),
      hjust = 1.05, vjust = -0.6, size = 2.4, colour = "grey20"
    ) +
    scale_color_manual(values = cols, labels = labs_map, name = NULL) +
    scale_fill_manual(values = cols, labels = labs_map, name = NULL) +
    scale_x_continuous(limits = c(0, PT.MAX), expand = c(0.01, 0)) +
    labs(
      title = title, subtitle = subtitle,
      x = "Progression to AD (Pseudotime)", y = ylab
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.25)
    )
}

##############################################################################
# Part A — Merged pathway dynamics (k = 7 + batch)
##############################################################################
message("=== Part A: merged pathway dynamics (k=7 + batch) ===")

counts <- as.matrix(read.csv(file.path(root, "scrna_counts.csv"), row.names = 1, check.names = FALSE))
cells <- read.csv(file.path(root, "scrna_cells.csv"), row.names = 1, check.names = FALSE)
pathway <- read.csv(file.path(root, "pathway_genes.csv"))$gene
donor_meta <- read.csv(file.path(root, "donor_meta.csv"), check.names = FALSE)

# Merge two cohort objects (paper: ROSMAP + cuimc2). No Harmony/anchors.
cells$batch <- factor(.map_to_paper_batch(cells$dataset), levels = BATCH.LEVELS)
stopifnot(nlevels(droplevels(cells$batch)) >= 2L)

obj_list <- lapply(BATCH.LEVELS, function(b) {
  idx <- rownames(cells)[as.character(cells$batch) == b]
  o <- CreateSeuratObject(
    counts = t(counts[idx, , drop = FALSE]),
    meta.data = cells[idx, , drop = FALSE],
    project = b
  )
  o$batch <- b
  o
})
names(obj_list) <- BATCH.LEVELS
obj <- merge(
  obj_list[[1L]], y = obj_list[-1L],
  add.cell.ids = names(obj_list), project = "merged_synthetic"
)
obj$batch <- factor(as.character(obj$batch), levels = BATCH.LEVELS)
obj <- tryCatch(JoinLayers(obj), error = function(e) obj)
message("Merged Seurat object: ", ncol(obj), " cells from ROSMAP + cuimc2")

# AddModuleScore defaults only (features + name); donor-mean → mean_scaled.
obj <- NormalizeData(obj, verbose = FALSE)
obj <- AddModuleScore(obj, features = list(pathway = pathway), name = "toy_pathway")
meta <- slot(obj, "meta.data")
score_col <- grep("^toy_pathway", colnames(meta), value = TRUE)
stopifnot(length(score_col) == 1)

donor_scores <- meta %>%
  group_by(projid) %>%
  summarise(mean_scaled = mean(.data[[score_col]], na.rm = TRUE), .groups = "drop") %>%
  mutate(projid = as.character(projid))

bulk <- donor_meta %>%
  mutate(
    ID = as.character(ID),
    apoe_4 = factor(as.character(apoe_4), levels = c("0", "1")),
    batch = factor(.map_to_paper_batch(dataset), levels = BATCH.LEVELS)
  ) %>%
  inner_join(donor_scores, by = c("ID" = "projid")) %>%
  filter(
    is.finite(pseudotime), pseudotime <= PT.MAX,
    is.finite(prAD), is.finite(mean_scaled),
    !is.na(apoe_4), !is.na(batch)
  )
message(
  "Part A donors: ", nrow(bulk),
  " (E4-=", sum(bulk$apoe_4 == "0"), ", E4+=", sum(bulk$apoe_4 == "1"), ")"
)

# mean_scaled ~ s(pseudotime, k=7) + batch; predict at batch = "ROSMAP"
fit_one_group_batch <- function(df_sub, label) {
  df_sub$batch <- droplevels(factor(df_sub$batch, levels = BATCH.LEVELS))
  if (nlevels(df_sub$batch) < 2L) stop("Need both batch levels in group ", label)
  mgcv::gam(
    mean_scaled ~ s(pseudotime, k = GAM.K) + batch,
    weights = prAD, data = df_sub
  )
}
fit.e3 <- fit_one_group_batch(bulk %>% filter(apoe_4 == "0"), "APOE4-")
fit.e4 <- fit_one_group_batch(bulk %>% filter(apoe_4 == "1"), "APOE4+")

ref_lvl <- "ROSMAP"
xs_a <- seq(min(bulk$pseudotime), PT.MAX, length.out = 50L)
newdata_a <- data.frame(
  pseudotime = xs_a,
  batch = factor(ref_lvl, levels = BATCH.LEVELS)
)
pred_one_a <- function(fit, feature_label) {
  preds <- predict(fit, newdata = newdata_a, se.fit = TRUE)
  data.frame(
    x = xs_a, fit = as.numeric(preds$fit), se.fit = as.numeric(preds$se.fit),
    feature = feature_label, stringsAsFactors = FALSE
  )
}
pred.vals.a <- bind_rows(
  pred_one_a(fit.e3, "toy_pathway.E4-"),
  pred_one_a(fit.e4, "toy_pathway.E4+")
)

null_a <- mgcv::gam(
  mean_scaled ~ s(pseudotime, k = GAM.K) + batch,
  weights = prAD, data = bulk
)
full_a <- mgcv::gam(
  mean_scaled ~ apoe_4 + batch + s(pseudotime, by = apoe_4, k = GAM.K),
  weights = prAD, data = bulk
)
anova_a <- mgcv::anova.gam(null_a, full_a, test = "Chisq")
print(anova_a)
anova_p_a <- anova_a[["Pr(>Chi)"]][[2L]]
if (is.null(anova_p_a) || is.na(anova_p_a)) stop("Part A ANOVA p is NA")
message("Part A ANOVA p: ", signif(anova_p_a, 4))

write.csv(
  data.frame(
    part = "A_merged_pathway",
    feature = "toy_pathway",
    anova_p_value = anova_p_a,
    anova_deviance = anova_a[["Deviance"]][[2L]],
    n = nrow(bulk),
    ref_batch = ref_lvl,
    formula_note = "s(pseudotime,k=7)+batch; full +apoe_4 + by-smooth k=7"
  ),
  file.path(out_dir, "01_pathway_gam_anova.csv"),
  row.names = FALSE
)

p_a <- .dynamics_panel(
  pred.vals.a, "toy_pathway.E4+", "toy_pathway.E4-", anova_p_a,
  title = "Part A: toy pathway (merged)",
  subtitle = paste0(
    "mean_scaled ~ s(pseudotime, k=7) + batch; ribbon at batch = ", ref_lvl
  ),
  ylab = "mean_scaled"
)
ggsave(file.path(out_dir, "01_pathway_gam.pdf"), p_a, width = 5.5, height = 3.5)
message("Wrote Part A outputs")

##############################################################################
# Part B — Original-cohort cell state + pathology (no batch, no k)
# Matches ANOVA.dyn.R fit.split.GAM
##############################################################################
message("=== Part B: cohort cell-state + pathology dynamics (no batch, no k) ===")

cohort <- donor_meta %>%
  mutate(
    ID = as.character(ID),
    apoe_4 = factor(as.character(apoe_4), levels = c("0", "1"))
  ) %>%
  filter(
    is.finite(pseudotime), pseudotime <= PT.MAX,
    is.finite(prAD), !is.na(apoe_4),
    is.finite(Ast.10), is.finite(sqrt.amyloid_mf)
  )
message("Part B donors: ", nrow(cohort))

# Per-genotype curve (no batch): feature ~ s(pseudotime)
.fit_split_and_anova <- function(df, feature) {
  df <- df %>% mutate(y = .data[[feature]])
  fit_e3 <- mgcv::gam(y ~ s(pseudotime), weights = prAD, data = df %>% filter(apoe_4 == "0"))
  fit_e4 <- mgcv::gam(y ~ s(pseudotime), weights = prAD, data = df %>% filter(apoe_4 == "1"))
  xs <- seq(min(df$pseudotime), PT.MAX, length.out = 50L)
  newdata <- data.frame(pseudotime = xs)
  pred <- bind_rows(
    {
      pr <- predict(fit_e3, newdata = newdata, se.fit = TRUE)
      data.frame(
        x = xs, fit = as.numeric(pr$fit), se.fit = as.numeric(pr$se.fit),
        feature = paste0(feature, ".E4-"), stringsAsFactors = FALSE
      )
    },
    {
      pr <- predict(fit_e4, newdata = newdata, se.fit = TRUE)
      data.frame(
        x = xs, fit = as.numeric(pr$fit), se.fit = as.numeric(pr$se.fit),
        feature = paste0(feature, ".E4+"), stringsAsFactors = FALSE
      )
    }
  )
  # ANOVA.dyn.R fit.split.GAM: null s(pseudotime) vs apoe_4 + s(pseudotime, by=apoe_4)
  null_mod <- mgcv::gam(y ~ s(pseudotime), weights = prAD, data = df)
  full_mod <- mgcv::gam(
    y ~ apoe_4 + s(pseudotime, by = apoe_4),
    weights = prAD, data = df
  )
  ar <- mgcv::anova.gam(null_mod, full_mod, test = "Chisq")
  pval <- ar[["Pr(>Chi)"]][[2L]]
  if (is.null(pval) || is.na(pval)) stop("Part B ANOVA p NA for ", feature)
  list(pred = pred, anova_p = pval, anova = ar, n = nrow(df))
}

res_ast <- .fit_split_and_anova(cohort, "Ast.10")
res_amy <- .fit_split_and_anova(cohort, "sqrt.amyloid_mf")
print(res_ast$anova)
print(res_amy$anova)
message("Part B Ast.10 ANOVA p: ", signif(res_ast$anova_p, 4))
message("Part B amyloid ANOVA p: ", signif(res_amy$anova_p, 4))

write.csv(
  data.frame(
    part = c("B_cohort_state", "B_cohort_pathology"),
    feature = c("Ast.10", "sqrt.amyloid_mf"),
    anova_p_value = c(res_ast$anova_p, res_amy$anova_p),
    n = c(res_ast$n, res_amy$n),
    formula_note = "s(pseudotime); full apoe_4 + s(pseudotime, by=apoe_4); no batch, no k"
  ),
  file.path(out_dir, "01_cohort_dynamics_anova.csv"),
  row.names = FALSE
)

p_ast <- .dynamics_panel(
  res_ast$pred, "Ast.10.E4+", "Ast.10.E4-", res_ast$anova_p,
  title = "Part B: Ast.10 (cohort)",
  subtitle = "Ast.10 ~ s(pseudotime)",
  ylab = "Ast.10"
)
p_amy <- .dynamics_panel(
  res_amy$pred, "sqrt.amyloid_mf.E4+", "sqrt.amyloid_mf.E4-", res_amy$anova_p,
  title = "Part B: amyloid (cohort)",
  subtitle = "sqrt.amyloid_mf ~ s(pseudotime)",
  ylab = "sqrt.amyloid_mf"
)
ggsave(
  file.path(out_dir, "01_cohort_dynamics.pdf"),
  p_ast + p_amy,
  width = 10, height = 3.5
)
message("Wrote Part B outputs under ", out_dir)
message("Done notebook 01")

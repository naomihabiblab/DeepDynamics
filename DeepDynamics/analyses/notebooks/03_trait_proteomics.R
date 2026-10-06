#!/usr/bin/env Rscript
# Sample: trait association (associate.traits) + Fig 3-style TA heatmaps +
# proteomics density Wilcoxon. Matches utils.TA.R / fig.TA.funs.R
# (pathology.sex.apoe, states.pathologies) and utils.prot.R.
# Synthetic data only — not paper results.
#
# Run from DeepDynamics/analyses/:
#   Rscript notebooks/03_trait_proteomics.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(reshape2)
  library(grid)
  library(gridExtra)
  library(ComplexHeatmap)
  library(SummarizedExperiment)
})

# --- Paths and donor table -----------------------------------------------------
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

DEFAULT.PATHOLOGIES <- c("sqrt.amyloid_mf", "sqrt.tangles_mf", "cogng_demog_slope")
DEFAULT.CONTROLS <- c("age_death", "pmi", "RIN")
orange <- "#ef6c00"
turqouise <- "#4db6ac"

donor <- read.csv(file.path(root, "donor_meta.csv"), check.names = FALSE) %>%
  mutate(
    apoe_4 = as.numeric(as.character(apoe_4)),
    n.apoe_4 = ifelse(apoe_4 == 1, 0, 1),
    msex = as.numeric(sex == "Male"),
    fsex = as.numeric(sex == "Female")
  )

# --- associate.traits / create.sets.and.run.SA (utils.TA.R) --------------------
# For each trait × covariate: lm(trait ~ covariate + age_death + pmi + RIN),
# then BH FDR within trait. Stars mark FDR bins (***/**/*).
associate.traits <- function(traits, covariates, controls, p.adjust.method = "BH") {
  df <- data.frame(traits, covariates, controls, check.names = FALSE)
  out <- do.call(rbind, lapply(colnames(traits), function(trait) {
    do.call(rbind, lapply(colnames(covariates), function(covariate) {
      control <- colnames(controls)
      control <- if (trait == "cww") {
        paste(setdiff(control, "msex"), collapse = " + ")
      } else {
        paste(control, collapse = " + ")
      }
      formula <- stringr::str_interp("${trait} ~ ${covariate} + ${control}")
      m <- summary(lm(formula, df[!is.na(df[[trait]]), ]))
      data.frame(
        trait = trait,
        covariate = covariate,
        beta = m$coefficients[covariate, 1],
        se = m$coefficients[covariate, 2],
        tstat = m$coefficients[covariate, 3],
        pval = m$coefficients[covariate, 4],
        r.sq = m$adj.r.squared,
        n = sum(!is.na(df[[trait]])),
        formula = formula,
        stringsAsFactors = FALSE
      )
    }))
  }))
  out %>%
    group_by(trait) %>%
    mutate(
      adj.pval = p.adjust(pval, method = p.adjust.method),
      sig = cut(adj.pval, c(-0.1, 0.001, 0.01, 0.05, Inf), c("***", "**", "*", ""))
    ) %>%
    ungroup()
}

create.sets.and.run.SA <- function(data, set.definitions, trait.names, control.names) {
  sets <- lapply(set.definitions, function(def) {
    list(
      covariates = data[, def$covariate.names, drop = FALSE],
      traits = data[, trait.names, drop = FALSE],
      controls = data[, control.names, drop = FALSE]
    )
  })
  names(sets) <- sapply(set.definitions, function(def) def$name)
  sapply(
    names(sets),
    function(n) associate.traits(sets[[n]]$traits, sets[[n]]$covariates, sets[[n]]$controls),
    simplify = FALSE,
    USE.NAMES = TRUE
  )
}

# --- Set 1: pathology.sex.apoe — separate sex and APOE covariate sets ----------
set.sex.apoe <- list(
  list(name = "sex", covariate.names = c("msex", "fsex")),
  list(name = "apoe", covariate.names = c("apoe_4", "n.apoe_4"))
)
ta.sex.apoe <- create.sets.and.run.SA(
  donor, set.sex.apoe, DEFAULT.PATHOLOGIES, DEFAULT.CONTROLS
)

# --- Set 2: states.pathologies — SIG.CLUSTERS covariates -----------------------
# Same fallback cluster list as paper utils.TA.R; values are synthetic proportions.
SIG.CLUSTERS <- c(
  "Oli.11", "Mic.12", "Oli.7", "Ast.10", "Mic.13",
  "Ast.7", "Oli.3", "Oli.4", "Ast.4", "Inh.6",
  "Exc.1", "Exc.3", "Inh.15", "OPC.1", "Inh.16",
  "Inh.5", "Exc.8", "End.3", "Inh.12", "Oli.5",
  "Mic.1", "Exc.12", "Inh.7", "Ast.6", "End.1",
  "Ast.2", "OPC.2", "Mic.14", "Ast.1", "Oli.9"
)
state.cols <- intersect(SIG.CLUSTERS, colnames(donor))
if (length(state.cols) < 2L) {
  stop("Need >=2 SIG.CLUSTERS columns in donor_meta.csv")
}
message("Trait association states: ", length(state.cols), " clusters")
set.states <- list(
  list(name = "states", covariate.names = state.cols)
)
ta.states <- create.sets.and.run.SA(
  donor, set.states, DEFAULT.PATHOLOGIES, DEFAULT.CONTROLS
)

out_dir <- file.path("notebooks", "outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ta.all <- bind_rows(
  ta.sex.apoe$sex %>% mutate(set = "sex"),
  ta.sex.apoe$apoe %>% mutate(set = "apoe"),
  ta.states$states %>% mutate(set = "states")
)
write.csv(ta.all, file.path(out_dir, "03_trait_association.csv"), row.names = FALSE)
print(ta.all %>% arrange(adj.pval) %>% head(12))

# --- Trait-association heatmaps (ta.heatmap / pathology.sex.apoe style) --------
# Cell color = -log10(FDR) * sign(beta); stars from FDR bins; turquoise–white–
# orange scale. Sex and APOE panels side-by-side (paper Fig 3 A–B layout);
# states panel written separately.
color.pheatmap.neg.pos <- function(df, colors = c(turqouise, "white", orange)) {
  paletteLength <- 100
  colors_step <- colorRampPalette(colors)(n = paletteLength)
  breaks <- c(
    seq(-max(abs(df), na.rm = TRUE), 0, length.out = 50),
    seq(0.01, max(abs(df), na.rm = TRUE), length.out = 50)
  )
  list(colors = colors_step, breaks = breaks)
}

return.signed.pvalues <- function(ta) {
  beta <- as.matrix(
    dcast(as.data.frame(ta), trait ~ covariate, value.var = "beta") %>%
      tibble::column_to_rownames("trait")
  )
  adj.pval <- as.matrix(
    dcast(as.data.frame(ta), trait ~ covariate, value.var = "adj.pval") %>%
      tibble::column_to_rownames("trait")
  )
  -log10(adj.pval) * sign(beta)
}

return.stars.pv <- function(ta) {
  as.matrix(
    dcast(as.data.frame(ta), trait ~ covariate, value.var = "sig") %>%
      tibble::column_to_rownames("trait")
  )
}

ta.heatmap <- function(ta, title, color.legend,
                       cluster.rows = FALSE, cluster.cols = FALSE,
                       angle.col = "0") {
  p.signed <- return.signed.pvalues(ta)
  stars <- return.stars.pv(ta)
  cb.list <- color.pheatmap.neg.pos(p.signed)

  row_dict <- c(
    "sqrt.amyloid_mf" = "Amyloid",
    "sqrt.tangles_mf" = "Tangles",
    "cogng_demog_slope" = "Cognitive Decline"
  )
  col_dict <- c(
    "msex" = "Males",
    "fsex" = "Females",
    "apoe_4" = "ApoE4\nCarriers",
    "n.apoe_4" = "ApoE4\nNon-Carriers"
  )

  if (all(rownames(p.signed) %in% names(row_dict))) {
    rownames(p.signed) <- row_dict[rownames(p.signed)]
    rownames(stars) <- rownames(p.signed)
  }
  if (all(colnames(p.signed) %in% names(col_dict))) {
    colnames(p.signed) <- col_dict[colnames(p.signed)]
    colnames(stars) <- colnames(p.signed)
  }

  ComplexHeatmap::pheatmap(
    p.signed,
    cell_fun = function(j, i, x, y, w, h, fill) grid.text(stars[i, j], x, y),
    color = cb.list$colors,
    breaks = cb.list$breaks,
    show_colnames = TRUE,
    cluster_rows = cluster.rows,
    cluster_cols = cluster.cols,
    show_column_dend = FALSE,
    heatmap_legend_param = list(
      title_position = "topcenter",
      legend_direction = "horizontal",
      title = color.legend
    ),
    na_col = "grey",
    main = title,
    angle_col = angle.col
  )
}

color.label <- "-log(FDR) * sign(beta)"
p.sex <- ta.heatmap(
  ta.sex.apoe$sex,
  "Sex and Pathologies (synthetic)",
  color.label
)
p.apoe <- ta.heatmap(
  ta.sex.apoe$apoe,
  "APOE4 and Pathologies (synthetic)",
  color.label
)
p.states <- ta.heatmap(
  ta.states$states,
  "States and Pathologies (synthetic)",
  color.label,
  cluster.cols = TRUE
)

pdf(file.path(out_dir, "03_trait_association_heatmap.pdf"), width = 10, height = 4)
grid.arrange(
  grid.grabExpr(draw(p.sex, heatmap_legend_side = "bottom")),
  grid.grabExpr(draw(p.apoe, heatmap_legend_side = "bottom")),
  ncol = 2
)
dev.off()

pdf(file.path(out_dir, "03_trait_states_heatmap.pdf"), width = 14, height = 4)
draw(p.states, heatmap_legend_side = "bottom")
dev.off()
message("Wrote ", file.path(out_dir, "03_trait_association_heatmap.pdf"))
message("Wrote ", file.path(out_dir, "03_trait_states_heatmap.pdf"))

# --- Proteomics: prepare.prot.for.density / calc.genotype.pval -----------------
# Density of protein abundance by genotype; Wilcoxon of 23/33 vs 34/44.
assay <- as.matrix(read.csv(file.path(root, "proteomics_assay.csv"), row.names = 1, check.names = FALSE))
coldata <- read.csv(file.path(root, "proteomics_coldata.csv"), row.names = 1, check.names = FALSE)
rowdata <- read.csv(file.path(root, "proteomics_rowdata.csv"), check.names = FALSE)
rownames(rowdata) <- rowdata$Symbol

proteomics <- SummarizedExperiment(
  assays = list(abundance = assay),
  colData = S4Vectors::DataFrame(coldata),
  rowData = S4Vectors::DataFrame(rowdata)
)

calc.genotype.pval <- function(df) {
  apoe33 <- as.numeric(df$gene[df$apoe_genotype %in% c("33", "23")])
  apoe34 <- as.numeric(df$gene[df$apoe_genotype %in% c("34", "44")])
  wilcox.test(x = apoe33, y = apoe34)
}

prepare.prot.for.density <- function(proteomics, name) {
  gene.mapping <- rowData(proteomics)$Symbol == name
  if (!any(gene.mapping)) return(NULL)
  idx <- which(gene.mapping)[1]
  df <- cbind(
    gene = assay(proteomics)[idx, ],
    as.data.frame(colData(proteomics))
  )
  p.val <- calc.genotype.pval(df)$p.value
  list(df = df, p.val = p.val, uni.name = rowData(proteomics)$UniProt[idx])
}

marker <- "SPP1"
var.lst <- prepare.prot.for.density(proteomics, marker)
stopifnot(!is.null(var.lst))

plot_df <- var.lst$df %>%
  mutate(
    genotype_group = ifelse(apoe_genotype %in% c("34", "44"), "APOE4+", "APOE4-")
  )

p <- ggplot(plot_df, aes(x = gene, fill = genotype_group)) +
  geom_density(alpha = 0.45) +
  labs(
    title = paste0(marker, " abundance (synthetic)"),
    subtitle = paste0("Wilcoxon p = ", signif(var.lst$p.val, 3), " — not a paper result"),
    x = "Abundance",
    y = "Density"
  ) +
  theme_minimal()
ggsave(file.path(out_dir, "03_proteomics_density.pdf"), p, width = 5, height = 3.5)

write.csv(
  data.frame(gene = marker, wilcox_p = var.lst$p.val, uniprot = var.lst$uni.name),
  file.path(out_dir, "03_proteomics_wilcox.csv"),
  row.names = FALSE
)

message("Wrote trait association + proteomics outputs under ", out_dir)

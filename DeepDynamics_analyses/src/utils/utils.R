##############################################################################
#                              Utils for Paper                               #
##############################################################################
##############################################################################
#                               Initialization                               #
##############################################################################
# Constant:
packages.to.import <- c(
    'anndata', 'dplyr', 'mgcv', 'tidyverse', 'ggplot2',
    'boot', 'reticulate', "colorspace", "ggrepel", "ggnewscale", 
    "progress", "RColorBrewer", 'pheatmap', 'ComplexHeatmap', 'circlize',
    "gridExtra", "reshape2", "dendsort", "cowplot", "tidyr", 
    "purrr", "stringr", "SummarizedExperiment", "EnhancedVolcano", "ggrepel",
    "biomaRt"
)

# Import packages
for (package in packages.to.import) {
    suppressMessages(library(package, character.only = TRUE))
}

# Set seed
set.seed(2024)

# Initalization of working directory
# Get the current working directory
cur.dir <- getwd()

# Function to check if a path ends with a specific directory
ends.with.dir <- function(path, dir) {
    grepl(paste0(dir, "$"), path)
}

# Check the current directory and act accordingly
if (ends.with.dir(cur.dir, "paper")) {
    message("Already in the 'paper' folder. No change in working directory.")
} else if (ends.with.dir(cur.dir, "AD-Metamodelling")) {
    if (dir.exists("paper")) {
        setwd("paper")
        message("Working directory set to 'paper' folder.")
    } else {
        stop("'paper' folder not found in AD-Metamodelling directory.")
    }
} else {
    # Check if we're one level deeper in the paper folder
    parent.dir <- dirname(cur.dir)
    if (ends.with.dir(parent.dir, "AD-Metamodelling") && ends.with.dir(cur.dir, "paper")) {
        message("Already in a subdirectory of 'paper'. No change in working directory.")
    } else {
        stop("Not in AD-Metamodelling or paper folder. Aborting.")
    }
}


# Paths for data files
PATH.500had <- '../data/sync/500.updated.h5ad'
PATH.BULK.META <- '../data/sync/dataset_707_basic_04-25-2022.csv'
PATH.RIN <- '../data/sync/RIN_n1089_08222023.csv'
PATH.bulk.df <- '../data/sync/bulk.prev.meta.RData'
PATH.prot <- '../data/sync/proteomics.rds'
PATH.BULK.SEQ <- '../data/sync/geneTpmResidualsAgeGenderUnadj'


# Paths for saving results
PATH.FIGS <- 'results/figures'
PATH.SUPP <- 'results/supp.figures'
PATH.DRAFTS <- 'results/drafts'

# Paths for other util files
PATH.TA <- 'src/utils/utils.TA.R'
PATH.DYN <- 'src/utils/utils.dynamics.R'

# Used colors
orange <- "#ef6c00"
turquoise <- "#4db6ac"
border.col <- "grey25"
apoe4.pallet.col <- "Pastel1"


# Default width and height for figures:
DYN.SIZE <- list(w = 3.5, h = 2)  # For two trajectories
SINGLE.DYN.SIZE <- list(w = 1.75, h = 2)  # For one trajectory
TA.SIZE <- list(w = 11, h = 3.4)  # For the TA plots
PROT.SINGLE.SIZE <- list(w = 3.5, h = 2)  # For one protein density plot
VOLCANO.SIZE <- list(w = 7, h = 7.5) # For volcano plots
SPIDER.SIZE <- list(w = 2.5, h = 2.5) # For spider plots
HEATMAP.SIZE <- list(w = 10, h = 3.4) # For heatmaps
COMPARE.HEATMAP.SIZE <- list(w = 15, h = 10) # For heatmaps

# Default sizes for axis ticks and title:
AXIS.TICK.SIZE <- element_text(size = 5)
AXIS.TITLE.SIZE <- element_text(size = 7)
FACET.TITLE.SIZE <- element_text(size = 8)
PLOT.TITLE.SIZE <- element_text(size = 8)
LEGEND.TITLE.SIZE <- element_text(size = 6)
LEGEND.TEXT.SIZE <- element_text(size = 5)


DYN.THEME <- theme_minimal() +
    theme(
        axis.line = element_line(color = "black", linewidth = 0.25),
        axis.text = AXIS.TICK.SIZE,
        axis.text.y = element_text(angle=0),
        axis.title = AXIS.TITLE.SIZE,
        strip.text = FACET.TITLE.SIZE,
        aspect.ratio = 0.8,
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none"
    )

DENSITY.THEME <- theme_minimal() +
    theme(
        axis.line = element_line(color = "black", linewidth = 0.25),
        axis.text = AXIS.TICK.SIZE,
        axis.text.y = element_text(angle=0),
        axis.title = AXIS.TITLE.SIZE,
        strip.text = FACET.TITLE.SIZE,
        plot.title = PLOT.TITLE.SIZE, 
        legend.text = LEGEND.TEXT.SIZE,
        legend.title = LEGEND.TITLE.SIZE,
        legend.key.size = unit(0.3, "cm"),
        aspect.ratio = 0.8,
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom"
    )

SPIDER.THEME <- theme_minimal() +
    theme(
        axis.line = element_line(color = "black", linewidth = 0.25),
        axis.text = AXIS.TICK.SIZE,
        axis.text.y = element_text(angle=0),
        axis.title = AXIS.TITLE.SIZE,
        strip.text = FACET.TITLE.SIZE,
        plot.title = element_text(size = 5), 
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.key.size = unit(0.3, "cm"),
        aspect.ratio = 0.8,
        panel.background = element_blank(),
        plot.background = element_blank(),
        legend.position = "bottom"
    )

# TODO decide if it should be here or not
# Source other util files
# source(PATH.TA)
# source(PATH.DYN)


# General Functions for Loading

#' load.metadata
#' @description
#' Function to load the whole metadata for the many people there are.
load.metadata <- function(bulk.dir = PATH.BULK.META) {
    read.csv(bulk.dir, header = T) %>%
        dplyr::mutate(
            projid = as.character(projid),
            fsex = abs(1 - msex),
            sex = factor(plyr::mapvalues(msex, c(F, T), c(
                "Female", "Male"
            ))),
            # Cognition state
            cogdx = factor(cogdx, levels = 1:6),
            cogdx_ad = as.numeric(
                recode(
                    cogdx,
                    "1" = "1",
                    "2" = "2",
                    "3" = NA_character_,
                    "4" = "3",
                    "5" = NA_character_,
                    "6" = NA_character_
                )
            ),
            cogdx_grouped = as.numeric(
                recode(
                    cogdx,
                    "3" = "2",
                    "4" = "3",
                    "5" = "3",
                    "6" = NA_character_
                )
            ),
            # Semi-quantitative pathology measuresceradsc = factor(ceradsc, levels = 4:1, ordered=T),braaksc = factor(case_when(braaksc == 0 ~ NA_integer_, T~braaksc), levels = 1:6, ordered = T),niareagansc = factor(niareagansc, levels = 4:1, ordered=T),pAD = factor(plyr::mapvalues(as.numeric(as.character(niareagansc)), from = c(1,2,3,4), to=c(1,1,0,0)), labels = c(0,1)),
            dlbdx = factor(dlbdx, levels = 0:3, ordered = T),
            tdp_stage4 = factor(tdp_stage4, levels = 0:3, ordered = T),
            # Quantitative pathology measures
            across(matches("amyloid|tangles|nft|plaq"), sqrt, .names = "sqrt.{.col}"),
            apoe_genotype = factor(apoe_genotype),
            apoe_2 = as.numeric(apoe_genotype %in% c(22, 23)),
            apoe_3 = as.numeric(apoe_genotype %in% c(33)),
            apoe_4 = as.numeric(apoe_genotype %in% c(34, 44))
        ) %>%
        dplyr::mutate(
            cerad.txt = factor(
                ceradsc,
                levels = 4:1,
                labels = c("No AD", "Possible", "Probable", "Definite")
            ),
            braak.txt = factor(
                braaksc,
                levels = 1:6,
                labels = c("I", "II", "III", "IV", "V", "VI")
            ),
            braak.grouped.txt = factor(
                braaksc,
                levels = 0:6,
                labels = c("0", "I", "II", "III", "IV", "V+VI", "V+VI")
            ),
            cdx.txt = factor(
                cogdx,
                levels = c(1, 2, 3, 4, 5, 6),
                labels = c(
                    "No Cognitive\nImpairment",
                    "Mild\nCognitive Impairment",
                    "Mild Cognitive\nImpairment\nand Another\nCause of CI",
                    "AD",
                    "AD\nand Another\nCause of CI",
                    "Other Dementia"
                )
            ),
            niareagen.txt = factor(
                niareagansc,
                levels = 4:1,
                labels = c("No AD", "Low", "Intermediate", "High")
            )
        ) %>%
        tibble::column_to_rownames("projid")
}

# Load metadata of all bulk samples together with RIN
# NOTE: must be used after sourcing utils.R.
load.full.bulk <- function() {
    if (!exists("bulk.metadata", envir = .GlobalEnv)) {
        bulk.metadata <- merge(
            load.metadata(),  # the bulk metadata
            read.csv(PATH.RIN, row.names = 1),  # The RIN of the bulk
            by.x = "row.names",
            by.y = "row.names",
            all.x = TRUE
        )
        assign("bulk.metadata", bulk.metadata, envir = .GlobalEnv)
        message("bulk.metadata has been created and assigned to the global environment.")
    } else {
        message("bulk.metadata already exists in the global environment. No changes were made.")
    }
}

load.celmod.bulk <- function(){
  if (!exists("bulk.df", envir = .GlobalEnv)) {
    load(PATH.bulk.df)
    
    bulk.df <- merge(
      bulk.df,
      read.csv(PATH.RIN, row.names = 1),  # The RIN of the bulk
      by.x = "ID",
      by.y = "row.names",
      all.x = TRUE
    )
    
    assign("bulk.df", bulk.df, envir = .GlobalEnv)
    
    message("bulk.df has been created and assigned to the global environment.")
  } else {
    message("bulk.df already exists in the global environment. No changes were made.")
  }
}


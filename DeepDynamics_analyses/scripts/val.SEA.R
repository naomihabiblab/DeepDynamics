##############################################################################
# Validation - Seattle AD #
##############################################################################
# This file is a script that is supposed to except arguments from sbatch
# script and run the validation pipeline for the Seattle AD dataset.
#' The validation pipline concludes the following steps:
#' 1. Create a Seurat object from the Seattle AD dataset.
#' 2. Create desired metadata features
#' 3. Run FindMarkers with MAST test for the object, ApoE4 vs the others
#' 4. Save the results in a file
#' 5. Create a volcano plot for the results


# Import libraries:
invisible(lapply(
    c(
        "anndata",
        "dplyr",
        "Seurat",
        "cowplot",
        "ggplot2",
        "reshape2",
        "reticulate",
        "circlize",
        "RColorBrewer",
        "EnhancedVolcano",
        "Matrix",
        "tidyverse"
    ),
    library,
    character.only = T
))

# Source the val.utils.R file:
source('./src/r/val.utils.R')

# Set seed:
set.seed(0)

# Constants:
PLOT.DIR <- './results/SEA/graphs/'
DATA.MIC <- '../../Shared/SeaAD/sc_data_files/cell_type_files/sea_ad_microglia.h5ad'
DATA.AST <- '../../Shared/SeaAD/sc_data_files/cell_type_files/single_cell_data/Astrocyte.h5ad'#'../../Shared/SeaAD/sc_data_files/cell_type_files/sea_ad_Astrocyte.h5ad'
DATA.OLI <- '../../Shared/SeaAD/sc_data_files/cell_type_files/sea_ad_Oli.h5ad'
DATA.OPC <- '../../Shared/SeaAD/sc_data_files/cell_type_files/sea_ad_opc.h5ad'
DATA.SST <- '../../Shared/SeaAD/sc_data_files/cell_type_files/sea_ad_sst.h5ad'
DATA.L2L3 <- '../../Shared/SeaAD/sc_data_files/cell_type_files/single_cell_data/L2_3 IT_no_UMIs.h5ad'
DATA.OLI.NEW <- 'oligodendrocyte_subset.h5ad'

# Define dictionary for args to cell types:
cell.types.arg <- c(
    'mic' = DATA.MIC,
    'ast' = DATA.AST,
    'oli' = DATA.OLI,
    'oli.new' = DATA.OLI.NEW,
    'opc' = DATA.OPC,
    'sst' = DATA.SST,
    'l2l3' = DATA.L2L3
)

cell.types.str <- c(
    'mic' = 'Microglia',
    'ast' = 'Astrocytes',
    'oli' = 'Oligodendrocytes',
    'oli.new' = 'Oligodendrocytes (new subset)',
    'opc' = 'OPCs',
    'sst' = 'SST',
    'l2l3' = 'L2/3 IT Neurons'
)

# Retrieve arguments from sbatch script:
args <- commandArgs(trailingOnly = TRUE)
cell.type <- args[1]
data.path <- cell.types.arg[cell.type]

min.logfc <- args[2]
min.pct <- args[3]
test.use <- args[4]
run.type <- args[5]

# Load data:
obj <- create.seurat.from.h5ad(data.path)

# Add ApoE4 & Sex metadata:
obj <- add.meta.essentials(obj)

if (run.type == "degs") {
    # Run DEGs analysis:
    e4.10 <- run.apoe4.vs.non(obj, cell.type = cell.types.str[cell.type], 
                              min.logfc = min.logfc, min.pct = min.pct, 
                              test = test.use)
    
    # Run the Volcano plot funciton:
    create.volcano.plot(e4.10, cell.type = cell.types.str[cell.type])
    
    message('Validation pipeline for DEGs completed successfully!')
} else if (run.type == "pathways") {
    message("Running pathway module score analysis...")
    
    pathway.results <- run.pathway.module.analysis(
        obj = obj,
        cell.type = cell.type,
        min.pct = as.numeric(min.pct),
        min.logfc = as.numeric(min.logfc),
        test.use = test.use
    )
    
    message('Validation pipeline for pathway module analysis completed successfully!')
}











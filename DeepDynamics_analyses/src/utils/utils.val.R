##############################################################################
#                         Validation Util Functions                          #
##############################################################################
create.seurat.from.h5ad <- function(data.path) {
    # Read the AnnData file
    adata <- read_h5ad(data.path)
    message("Successfully read h5ad file: ", data.path)
    
    # Extract the normalized counts
    norm.counts <- Matrix::t(adata$X)
    message("Successfully extracted normalized counts from h5ad file")
    
    # Extract metadata
    metadata <- as.data.frame(adata$obs)
    message("Successfully extracted metadata from h5ad file")
    
    # Create a Seurat object with the normalized data in the correct slot
    obj <- CreateSeuratObject(counts = norm.counts, 
                              meta.data = metadata,
                              project = "SEA-AD")
    message("Successfully created Seurat object from normalized counts")
    # Move the normalized data to the correct slot
    obj[["RNA"]]$data <- obj[["RNA"]]$counts
    obj[["RNA"]]$counts <- NULL
    message("Successfully moved normalized data to the RNA slot in Seurat object")
    # Set the default assay
    DefaultAssay(obj) <- "RNA"
    message("Successfully created Seurat object from h5ad file")
    
    return(obj)
}

# import h5py
# 
# with h5py.File('../../Shared/SeaAD/sc_data_files/cell_type_files/single_cell_data/L2_3 IT.h5ad', "a") as f:
#     if "/layers/UMIs" in f:
#     del f["/layers/UMIs"]
# print("Deleted 'UMIs' directly from disk")

add.meta.essentials <- function(obj) {
    # Remove people with APoE genotype 2/4:
    Idents(obj) <- "APOE Genotype"
    obj <- subset(obj, idents = c("2/2", "2/3", "3/3", "3/4", "4/4"))
    
    # Remove people that are not males or females:
    Idents(obj) <- "Sex"
    obj <- subset(obj, idents = c("Female", "Male"))
    
    # Define genotype groups
    obj@meta.data <- obj@meta.data %>%
        mutate(APOE_group = if_else(`APOE Genotype` %in% c("3/4", "4/4"), 
                                    "Risk", "NonRisk"))
    
    # Create interaction in metadata of sex and APOE risk:
    obj <- Seurat::AddMetaData(object = obj,
                               metadata = interaction(obj@meta.data$Sex, obj@meta.data$APOE_group, sep = "."),
                               col.name = "sex.apoe4")
    
    message("Successfully added essential metadata to Seurat object, filtered NaNs in sex and removed 2/2 individuals. ")
    
    return(obj)
}


run.apoe4.vs.non <- function(obj, cell.type, test = "MAST", 
                             min.logfc = 0.25, min.pct = 0.1,
                             save.path = './results/SEA/DEGs/',
                             add.fname = '') {
    # Change the idents to ApoE group:
    Idents(obj) <- "APOE_group"
    
    # Print a message with number of individuals in each group:
    message('Number of cells in each group:')
    print(table(obj$APOE_group))
    
    
    degs <- FindMarkers(obj, ident.2 = "Risk", ident.1 = "NonRisk",
                        test.use = test, min.pct = min.pct, logfc.threshold = min.logfc)
    
    e4.10 <- degs %>% rownames_to_column("gene")
    
    # Save degs as RData:
    fname = paste0(save.path, cell.type,'.',
                   test, '.', add.fname, ".e4.10.RData")
    save(e4.10, file = fname)
    
    message('Successfully ran DEGs for APOE4 vs Non-Risk individuals,
            and also saved the results as', fname)
    
    return(e4.10)
    
}


create.volcano.plot <- function(e4.10, cell.type, 
                                save.path = "./results/SEA/graphs/") {
    # Use EnhancedVolcano:
    p.co <- 0.05
    f.co <- 0.24
    p.all <- EnhancedVolcano(e4.10,
                             lab = e4.10$gene,
                             x = 'avg_log2FC',
                             y = 'p_val_adj',
                             title = paste(cell.type, 'APOE4 Carriers vs. Non-Carriers'),
                             pCutoff = p.co,
                             FCcutoff = f.co,
                             pointSize = 3.0,
                             labSize = 4.0,
                             legendLabels = c("NS", expression(Log[2] ~ FC), 
                                              "p-value (adjusted)", 
                                              expression(p - value ~ and ~ log[2] ~ FC)),
                             drawConnectors = TRUE)
    p.all$labels$subtitle <- paste('Adjusted p-value cutoff:', p.co, '. Log2FC cutoff:', f.co)
    
    # save plot in pdf:
    fname <- paste0(save.path, cell.type, ".SEA.e4.10.volcano.pdf")
    ggsave(fname, p.all, width = 10, height = 10)
    
    # Create another one with only HSP and DNAJB1 genes labeled:
    genes.to.label <- e4.10$gene[grepl("^HSP|^DNAJ", e4.10$gene)]
    if (cell.type == "SST") {
        genes.to.label <- c(genes.to.label, "SST")
    }
    
    p.hsp <- EnhancedVolcano(e4.10,
                             lab = e4.10$gene,
                             x = 'avg_log2FC',
                             y = 'p_val_adj',
                             title = paste(cell.type, 'ApoE4 Carriers vs. Non-Carriers'),
                             pCutoff = p.co,
                             FCcutoff = f.co,
                             pointSize = 3.0,
                             labSize = 4.0,
                             legendLabels = c("NS", expression(Log[2] ~ FC), 
                                              "p-value (adjusted)", 
                                              expression(p - value ~ and ~ log[2] ~ FC)),
                             selectLab = genes.to.label,
                             drawConnectors = TRUE)
    p.hsp$labels$subtitle <- paste('Adjusted p-value cutoff:', p.co, '. Log2FC cutoff:', f.co)
    
    # save plot in pdf:
    fname <- paste0(save.path, cell.type, ".SEA.e4.10.volcano.HSP.pdf")
    ggsave(fname, p.hsp, width = 10, height = 10)
    
    # Save RData with these two plots:
    fname <- paste0(save.path, cell.type, ".SEA.e4.10.volcano.plots.RData")
    save(p.all, p.hsp, file = fname)
    
    message('Successfully created volcano plots for', cell.type, 'and saved them as', fname)
} 


# Genesets to check
ast.met.ion.trans.synaptic <- c("ATP1A2", "ATP1B3", "CACNA1D", "CACNB2", 
                                "KCNE4", "KCNJ10", "KCNN3", "SLC30A1", 
                                "SLC38A2", "SLC39A11", "SLC39A12", "SLC5A3", 
                                "TRPM3")

ast.react.oxy.metabolic <- c("CLU", "CRYAB", "DDIT4", "EGFR", 
                             "MT3", "PDGFRB", "PRCP", "RORA", 
                             "SESN1", "SLC5A3", "SLC7A2")

mic.lipid.catabolic <- c("APOE", "ASAH1", "LIPA", "LPL", 
                         "PDE3B", "PRKCE", "SPP1") 

oli.morphogenesis <- c("ADAM10", "DOCK1", "DOCK5", "FMNL2", 
                       "KANK1", "KIF13B", "LRP4", "MKLN1", 
                       "RDX", "RHOBTB1", "SEMA3C", "SPP1", 
                       "TANC2", "TRPC5", "ZSWIM6")

#' Compute module scores and differential analysis for gene sets by cell type
#' 
run.pathway.module.analysis <- function(
        obj, 
        cell.type, 
        min.pct = 0.1, 
        min.logfc = 0.25, 
        test.use = "wilcox", 
        output.dir = "./results/pathways.sea-ad"
) {
    message(paste0("Running pathway module analysis for cell type: ", cell.type))
    
    # Create output dir if missing
    if (!dir.exists(output.dir)) {
        dir.create(output.dir, recursive = TRUE)
    }
    
    # Get all gene sets in global env matching this cell type
    gene.set.vars <- ls(envir = .GlobalEnv, pattern = paste0("^", cell.type, "\\."))
    
    if (length(gene.set.vars) == 0) {
        stop(paste0("No gene sets found for cell type: ", cell.type))
    }
    
    results <- list()
    pseudocount <- 1e-6
    
    for (gene.set.var in gene.set.vars) {
        genes <- get(gene.set.var, envir = .GlobalEnv)
        module.name <- paste0("module_", gene.set.var)
        
        # Ensure genes exist in data
        genes.present <- intersect(genes, rownames(obj))
        if (length(genes.present) < 3) {
            message(paste0("Skipping ", gene.set.var, ": fewer than 3 genes found in dataset"))
            next
        }
        
        # Add module score
        obj <- AddModuleScore(obj, features = list(genes.present), name = module.name)
        score.col <- paste0(module.name, "1")
        
        # Extract scores for each group
        risk.scores <- obj@meta.data[obj$APOE_group == "Risk", score.col]
        nonrisk.scores <- obj@meta.data[obj$APOE_group == "NonRisk", score.col]
        
        # Check for sufficient cells
        if (length(risk.scores) < 3 || length(nonrisk.scores) < 3) {
            message(paste0("Skipping ", gene.set.var, ": insufficient cells in a group"))
            next
        }
        
        # Calculate log fold change (add pseudocount to avoid log(0))
        avg_logFC <- log2(mean(risk.scores) + pseudocount) - log2(mean(nonrisk.scores) + pseudocount)
        
        
        # Wilcoxon test (two-sided)
        wilcox.res <- wilcox.test(risk.scores, nonrisk.scores)
        
        # Wilcoxon test (directional: Risk > NonRisk)
        wilcox.dir <- wilcox.test(risk.scores, nonrisk.scores, alternative = "greater")
        
        wilcox.df <- data.frame(
            p_val = wilcox.res$p.value,
            avg_diff = mean(risk.scores) - mean(nonrisk.scores),
            avg_logFC = avg_logFC,
            median_diff = median(risk.scores) - median(nonrisk.scores),
            p_val_directional = wilcox.dir$p.value,  # <-- Directional test result
            test = "Wilcoxon",
            row.names = gene.set.var
        )
        
        # t-test
        ttest.res <- t.test(risk.scores, nonrisk.scores)
        ttest.df <- data.frame(
            p_val = ttest.res$p.value,
            avg_diff = mean(risk.scores) - mean(nonrisk.scores),
            avg_logFC = avg_logFC,
            median_diff = median(risk.scores) - median(nonrisk.scores),
            test = "t-test",
            row.names = gene.set.var
        )
        
        # Save both results as a list
        both.results <- list(wilcoxon = wilcox.df, ttest = ttest.df)
        results[[gene.set.var]] <- both.results
        
        # Save both results to file
        save(both.results, file = file.path(output.dir, paste0(gene.set.var, "_both_tests.RData")))
        
        message(paste0("Saved Wilcoxon and t-test results for ", gene.set.var))
    }
    
    return(results)
}


# Create table of SEA results of needed:
# Set your folder path
folder_path <- "../results/pathways.sea-ad"

# List all RData files with the relevant suffix
rdata_files <- list.files(folder_path, pattern = "both_tests( copy)?\\.RData$", full.names = TRUE)

# Initialize a list to store results
summary_list <- list()

for (file in rdata_files) {
    # Load the file (loads 'both_results' into the environment)
    load(file)
    
    # Extract the wilcox dataframe from both_results
    wilcox_df <- both.results$wilcox
    
    # If wilcox_df is a dataframe, extract the first row (or modify as needed)
    # Adjust the row/column names if needed based on your data structure
    avg_logFC <- wilcox_df$avg_logFC[1]
    p_val <- wilcox_df$p_val[1]
    p_val_directional <- wilcox_df$p_val_directional[1]
    if (length(wilcox_df$p_val_directional[1]) ==0) {
        p_val_directional <- NA  # Handle NA values if they exist
        message("Warning: p_val_directional is NA for file:", file)
    }
    
    # Clean up the filename
    file_base <- basename(file)
    file_base <- sub(" both_tests\\.RData$", "", file_base)
    
    # Store results in a named vector
    summary_list[[file_base]] <- c(
        avg_logFC = avg_logFC,
        wilcox_p_val = p_val,
        wilcox_p_val_directional = p_val_directional
    )
    message(paste("Processed file:", file_base))
}

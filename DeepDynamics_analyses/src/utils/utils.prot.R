#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#                      Utils for Proteomics Analyses                      ####
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# -------------------------------------------------------------------------- #
##                     Constants and Global Variables                     ####
# -------------------------------------------------------------------------- #
EXCLUDED.GENOTYPES <- c("22", "24", "44")
AST.MARKERS <- c('SLC38A2', 'QDPR')
MIC.MARKERS <- c('CPM', 'SPP1', 'GPNMB')
NEUR.MARKERS <- c('SST','NPTX2')

# Check if proteomics is loaded:
if (!exists("proteomics", envir = .GlobalEnv)) {
    proteomics <- readRDS(PATH.prot)
    message("Loaded proteomics data.")
} else {
    message("Proteomics data already loaded.")
}
# -------------------------------------------------------------------------- #
##             Functions for Proteomics Distribution Analysis             ####
# -------------------------------------------------------------------------- #
calc.genotype.pval <- function(df, opt = "two.sided") {
    # Calc p-value for each genotype group:
    apoe33 <- df$gene[df$apoe_genotype %in% c('33', '23')] %>% as.numeric()
    apoe34 <- df$gene[df$apoe_genotype %in% c('34', '44')] %>% as.numeric()
    wilcox.test(x = apoe33, y = apoe34, opt = opt)
}


prepare.prot.for.density <- function(proteomics, name) {
    gene.mapping <- rowData(proteomics)$Symbol == name
    if(!(any(gene.mapping))) return(NULL)
    
    # In case there is more then one corresponding uniprot, use of lapply:
    prot.ls <- lapply(which(gene.mapping), function(uniprot) {
        df <-  cbind(gene = (assay(proteomics)[uniprot, ]), 
                     colData(proteomics)) %>%
            as.data.frame()
        uniprot <- rowData(proteomics)$UniProt[uniprot]
        #' TODO Notice it is not working now for more than once at the current
        #' structure but it does the last option and returns it,
        #' which is OK for the figures.
        return(list(df = df, uniprot = uniprot))
    })
    
    p.val <- calc.genotype.pval(prot.ls[[1]]$df)
    
    return(list(df = prot.ls[[1]]$df, p.val = p.val$p.value, uni.name = prot.ls[[1]]$uniprot))
}


prepare.prot.group.for.density <- function(proteomics, group.genes, group.name) {
    # Get gene symbols from dataset
    available.genes <- rowData(proteomics)$Symbol
    
    # Match and find missing
    matched.idx <- which(available.genes %in% group.genes)
    matched.genes <- available.genes[matched.idx]
    missing.genes <- setdiff(group.genes, matched.genes)
    
    # Report missing genes
    if (length(missing.genes) > 0) {
        message("In group '", group.name, "', the following genes are missing or NA in proteomics data: ",
                paste(missing.genes, collapse = ", "))
    }
    
    if (length(matched.idx) == 0) {
        warning("None of the genes in the group '", group.name, "' were found in the proteomics dataset.")
        return(NULL)
    }
    
    # Extract expression matrix
    expr.mat <- assay(proteomics)[matched.idx, , drop = FALSE]
    
    # Normalize each gene (row) — z-score
    expr.mat.z <- t(apply(expr.mat, 1, function(x) {
        if (sd(x, na.rm = TRUE) == 0) {
            rep(0, length(x))  # Avoid division by zero
        } else {
            scale(x)
        }
    }))
    
    # Compute mean score per sample
    mean.abundance <- colMeans(expr.mat.z, na.rm = TRUE)
    
    # Combine with sample metadata
    df <- cbind(gene = mean.abundance, colData(proteomics)) %>%
        as.data.frame()
    
    # Calculate p-value
    p.val <- calc.genotype.pval(df)$p.value
    
    return(list(
        df = df,
        p.val = p.val,
        uni.name = "mean score"
    ))
}

prepare.prot.group.for.density <- function(proteomics, gene.set, group.name = NULL) {
    gene.mapping <- rowData(proteomics)$Symbol %in% gene.set
    if (!(any(gene.mapping))) {
        warning("No matches found for: ", paste(gene.set, collapse = ", "))
        return(NULL)
    }
    
    # Extract matching assays
    matched.symbols <- rowData(proteomics)$Symbol[gene.mapping]
    matched.uniprot <- rowData(proteomics)$UniProt[gene.mapping]
    
    expr <- assay(proteomics)[gene.mapping, , drop = FALSE]
    
    # Normalize each protein across samples
    expr.norm <- t(scale(t(expr)))  # row-wise z-score
    
    # RowMeans (mean z-score per sample for the group)
    group.score <- colMeans(expr.norm, na.rm = TRUE)
    
    df <- cbind(gene = group.score, colData(proteomics)) %>%
        as.data.frame()
    
    # P-value (using overall score)
    p.val <- calc.genotype.pval(df)
    
    list(
        df = df,
        uni.name = "mean score",
        p.val = p.val$p.value,
        missing = setdiff(gene.set, matched.symbols),
        present = matched.symbols
    )
}

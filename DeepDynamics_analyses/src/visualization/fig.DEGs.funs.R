#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#                       Functions to Create DEGs Plots                    ####
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# -------------------------------------------------------------------------- #
## Constants to be used: ####
# -------------------------------------------------------------------------- #
heatshock_groups <- list(
    HSP70 = c("HSPA1A", "HSPA1B", "HSPA6", "HSPA5", "HSPA8", "HSPA4", "HSPA4L", "HSPA12A"),
    HSP40 = c("DNAJA1", "DNAJA4", "DNAJB1", "DNAJB2", "DNAJB4", "DNAJB5", 
              "DNAJC2", "DNAJC3", "DNAJC5", "DNAJC6", "DNAJC9", "DNAJB11"),
    HSP90 = c("HSP90AA1", "HSP90AB1", "HSP90B1"),
    HSP60_10 = c("HSPD1", "HSPE1"),
    Small_HSPs = c("HSPB1", "HSPB3", "HSPBAP1"),
    HSP110 = c("HSPH1")
)

# -------------------------------------------------------------------------- #
##                            2-Group Plots                               ####
# -------------------------------------------------------------------------- #
# Volcano plots:
create.volcano.plot <- function(e4.10, cell.type, 
                                p.co = 0.05, f.co = 0.25,
                                degs.title = 'APOE4 Carriers vs. Non-Carriers',
                                xlim = NULL) {
    #' This function assumes to accepts in e4.10 the results of Seurat::FindMarkers
    #' Function to create volcano plots for DEGs.
    #' 
    if (!'gene' %in% colnames(e4.10)) {
        e4.10 <- rownames_to_column(e4.10, 'gene')
    }
    
    if (is.null(xlim)) {
        xlim <- c(min(e4.10$avg_log2FC, na.rm = TRUE) - 1.5, 
                  max(e4.10$avg_log2FC, na.rm = TRUE) + 1.5)
    }
    
    p.degs <- EnhancedVolcano(e4.10,
                             lab = e4.10$gene,
                             x = 'avg_log2FC',
                             y = 'p_val_adj',
                             xlim = xlim,
                             title = paste(cell.type, degs.title),
                             pCutoff = p.co,
                             FCcutoff = f.co,
                             pointSize = 3.0,
                             labSize = 4.0,
                             legendLabels = c("NS", expression(Log[2] ~ FC), 
                                              "p-value (adjusted)", 
                                              expression(p - value ~ and ~ log[2] ~ FC)),
                             drawConnectors = TRUE)
    p.degs$labels$subtitle <- paste('Adjusted p-value cutoff:', p.co, '. Log2FC cutoff:', f.co)
    
    return(p.degs)
}

# -------------------------------------------------------------------------- #
##                            4-Group Plots                               ####
# -------------------------------------------------------------------------- #
create.spider.graph <- function(fname.rd,
                                 top.degs.num = 10,
                                 min.abs.log2FC = 0,
                                 max.adj.pval = 0.05,
                                 add.title = '') {
    # Read RData file
    load(fname.rd)
    
    # Format df of the loaded results and filter bu inputs #TODO function
    f0f1 <- f0f1 %>%
        rownames_to_column('gene') %>%
        dplyr::select(gene, p_val_adj, avg_log2FC) %>%
        filter(p_val_adj < max.adj.pval) %>%
        `colnames<-`(c('gene', 'p.val.adj.f0f1', 'avg.log2FC.f0f1'))
    
    m0m1 <- m0m1 %>%
        rownames_to_column('gene') %>%
        dplyr::select(gene, p_val_adj, avg_log2FC) %>%
        filter(p_val_adj < max.adj.pval) %>%
        `colnames<-`(c('gene', 'p.val.adj.m0m1', 'avg.log2FC.m0m1'))
    
    m0f0 <- m0f0 %>%
        rownames_to_column('gene') %>%
        dplyr::select(gene, p_val_adj, avg_log2FC) %>%
        filter(p_val_adj < max.adj.pval) %>%
        `colnames<-`(c('gene', 'p.val.adj.m0f0', 'avg.log2FC.m0f0'))
    
    m1f1 <- m1f1 %>%
        rownames_to_column('gene') %>%
        dplyr::select(gene, p_val_adj, avg_log2FC) %>%
        filter(p_val_adj < max.adj.pval) %>%
        `colnames<-`(c('gene', 'p.val.adj.m1f1', 'avg.log2FC.m1f1'))
    
    print("m0m1 HSP")
    print(m0m1%>%filter(grepl('^HSP', gene)) %>% pull(gene))
    
    print("f0f1 HSP")
    print(f0f1%>%filter(grepl('^HSP', gene)) %>% pull(gene))
    
    # Merge into 1 df by gene
    APOE.effect.df <- merge(m0m1, f0f1, by = 'gene', all = TRUE)
    sex.effect.df <- merge(m0f0, m1f1, by = 'gene', all = TRUE)
    
    # APOE4 effect:
    top.f.up <- get.top.degs(APOE.effect.df,
                             "avg.log2FC.f0f1",
                             min.abs.log2FC,
                             top.degs.num,
                             "up")
    top.f.down <- get.top.degs(APOE.effect.df,
                               "avg.log2FC.f0f1",
                               min.abs.log2FC,
                               top.degs.num,
                               "down")
    top.m.up <- get.top.degs(APOE.effect.df,
                             "avg.log2FC.m0m1",
                             min.abs.log2FC,
                             top.degs.num,
                             "up")
    top.m.down <- get.top.degs(APOE.effect.df,
                               "avg.log2FC.m0m1",
                               min.abs.log2FC,
                               top.degs.num,
                               "down")
    
    # Put the names of the genes if they are in one of the top n genes.
    APOE.effect.df$label.flag <- ifelse(
        APOE.effect.df$gene %in%
            c(top.f.up, top.m.up, top.f.down, top.m.down) | grepl("^(HSP|DNAJ|NAMPT|FLT1|PTPRG)", APOE.effect.df$gene), # If needed extra conditions: | grepl("^HSP", APOE.effect.df$gene),
        APOE.effect.df$gene,
        ""
    )
    
    # APOE.effect.df$label.flag <- APOE.effect.df$gene
    
    
    # Identigy geens that are present in the top n of both comparisons
    APOE.effect.df$in.both.topn <- ifelse(APOE.effect.df$gene %in% intersect(c(top.f.up, top.f.down), c(top.m.up, top.m.down)), "Both", "One")
    
    APOE.p <- ggplot(APOE.effect.df, aes(
        ifelse((
            is.na(avg.log2FC.m0m1) |
                is.na(p.val.adj.m0m1)
        ), 0, avg.log2FC.m0m1),
        ifelse((
            is.na(avg.log2FC.f0f1) |
                is.na(p.val.adj.f0f1)
        ), 0, avg.log2FC.f0f1)
    )) +
        geom_point(aes(color = in.both.topn, size = in.both.topn)) +
        scale_color_manual(values = c("Both" = "darkred", "One" = "black")) +
        scale_size_manual(values = c("Both" = .2, "One" = .2)) +
        labs(
            x = paste('Males APOE3 vs. APOE4 ', "logFC"),
            y = paste('Females APOE3 vs. APOE4 ', "logFC"),
            color = paste("Top", top.degs.num, "Overlap"),
            size = paste("Top", top.degs.num, "Overlap")
        ) +
        geom_text_repel(aes(label = label.flag),
                        max.overlaps = 700000,
                        size = 1.7,
                        segment.size = 0.25) +
        geom_abline(intercept = 0, slope = 0) + geom_vline(xintercept = 0) +
        theme_minimal() +
        ggtitle(paste('APOE4 Effect on Males vs. Females', add.title)) +
        theme(plot.title = element_text(hjust = 0.5)) +
        SPIDER.THEME
    
    
    return(APOE.p)
    
}


# -------------------------------------------------------------------------- #
##                             Heatmap Plots                              ####
# -------------------------------------------------------------------------- #

# Process before plotting:
create.HM.order.annotation <- function(mat, split.by = 'celltype') {
    col.info <- colnames(mat) %>%
        str_split("_") %>%
        do.call(rbind, .) %>%
        as.data.frame() %>%
        `colnames<-`(c('celltype', 'description', 'datatype', 'sex'))
    
    col.info$order.factor <- with(
        col.info, interaction(sex, description, sep = "_")
    )
    
    ordered.cols <- order(col.info$order.factor)
    mat <- mat[, ordered.cols] # Reorder the matrix 
    
    col.split.vec <- col.info[[split.by]][ordered.cols] %>% # Split by celltype
        as.factor()
    
    col.anno.ordered <- HeatmapAnnotation(
        celltype = anno_text(col.info$celltype[ordered.cols], rot = 90, gp = gpar(fontsize = 8)),
        description = anno_text(col.info$description[ordered.cols], rot = 90, gp = gpar(fontsize = 8)),
        sex = anno_text(col.info$sex[ordered.cols], rot = 90, gp = gpar(fontsize = 8)),
        show_annotation_name = TRUE,
        annotation_name_side = "left"
    )
    return(list("anno" = col.anno.ordered, "split.vec" = col.split.vec, "mat.ordered" = mat))
}

# Create the heatmap of desired genes
create.genes.heatmap <- function(dir.path,
                                 split.by = 'celltype',
                                 title = "APOE4 Carriers vs Non-Carriers - HSP Genes",
                                 pattern.genes = "^(HSP|DNAJ)",
                                 filter.genes = NULL,
                                 color.sat = T) {

    # Create a matrix of all DEGs:
    mat <- create.degs.matrix(dir.path)

    # If filter.genes not null, filter the matrix:
    if (!is.null(filter.genes)) {
        # Find which genes from filter.genes are actually in the matrix
        existing.genes <- intersect(rownames(mat), filter.genes)

        # Check if any genes were found
        if (length(existing.genes) > 0) {
            mat <- mat[existing.genes, , drop = FALSE]
            cat("Selected", length(existing.genes), "genes out of", length(filter.genes), "requested.\n")
        } else {
            warning("None of the specified genes were found in the matrix.")
            # Optionally, you could return the original matrix or an empty matrix here
            # mat <- matrix(nrow = 0, ncol = ncol(mat))
        }
    } else {

        mat <- mat[grepl("^(HSP|DNAJ)", rownames(mat)), , drop = FALSE] #TODO make this a parameter

        print(grepl(pattern.genes, rownames(mat)))
        mat <- mat[grepl(pattern.genes, rownames(mat)), , drop = FALSE]

        cat("No genes specified for filtering. Using all genes.\n")
    }

    # # Create heatmap parameters:
    hm.params <- create.HM.order.annotation(mat, split.by)

    # Color map:
    col.map = colorRamp2(c(min(mat, na.rm = TRUE), 0, max(mat, na.rm = TRUE)),
                     c(turquoise, "white", orange))
    if (min(mat, na.rm = TRUE) >= 0) {
        col.map = colorRamp2(c(min(mat, na.rm = TRUE), 0, max(mat, na.rm = TRUE)),
                             c("white", "white", orange))
    }
    if (color.sat){
        col.map = colorRamp2(c(-2, 0, 2),
                             c(turquoise, "white", orange)) 
    }
    # Create heatmap:
    heatmap <- Heatmap(hm.params$mat.ordered,
                       name = "avg_log2FC",
                       top_annotation = hm.params$anno,
                       show_row_names = TRUE,
                       show_column_names = FALSE,
                       column_names_gp = gpar(fontsize = 4),
                       column_title = title,
                       col = col.map,
                       na_col = "grey",
                       cluster_rows = TRUE,
                       cluster_columns = FALSE,
                       column_split = hm.params$split.vec,
                       gap = unit(2, "mm"),
                       border = TRUE,
                       column_gap = unit(2, "mm"),
                       row_gap = unit(2, "mm"))

    return(heatmap)
}

# --------------------------------------------------------------------------------- #
#                           Validation Heatmap Functions                         ####
# --------------------------------------------------------------------------------- #
process.single.e4.10.file <- function(file.path) {
    load(file.path)  # loads object called e4.10
    if (!exists("e4.10")) stop(paste("e4.10 object not found in", file.path))
    
    celltype <- str_extract(basename(file.path), "^[^\\.]+")  # extract prefix before first dot
    
    tryCatch({
        e4.10 %>%
            filter(p_val_adj <= 0.05) %>%
            # rownames_to_column("gene") %>%
            dplyr::select(gene, avg_log2FC) %>%
            mutate(celltype = celltype)
    }, error = function(e) {
        e4.10 %>%
            filter(p_val_adj <= 0.05) %>%
            rownames_to_column("gene") %>%
            dplyr::select(gene, avg_log2FC) %>%
            mutate(celltype = celltype)
    })
}

# Create combined DEGs matrix
create.e4.10.degs.matrix <- function(dir.path, gene.pattern = "^(HSP|DNAJ)", file.pattern = "e4\\.10\\.RData$") {
    files <- list.files(dir.path, pattern = file.pattern, full.names = TRUE)
    
    degs <- map(files, safely(process.single.e4.10.file)) %>%
        map("result") %>%
        compact() %>%
        bind_rows()
    
    # Filter genes by regex
    degs <- degs %>%
        filter(str_detect(gene, gene.pattern))
    
    all.genes <- unique(degs$gene)
    
    # Complete matrix
    degs.mat <- degs %>%
        complete(gene = all.genes, celltype, fill = list(avg_log2FC = 0)) %>%
        pivot_wider(
            names_from = celltype,
            values_from = avg_log2FC,
            values_fill = 0
        )
    
    mat <- as.matrix(degs.mat[, -1])
    rownames(mat) <- degs.mat$gene
    
    return(mat)
}

create.simple.genes.heatmap.annotated <- function(
        dir.path,
        gene.pattern = "^(HSP|DNAJ)",
        file.pattern = "e4\\.10\\.RData$",
        title = "APOE4+ vs APOE4- (Simplified)",
        color.sat = TRUE,
        cluster.columns = TRUE,
        cluster.rows = TRUE
) {
    mat <- create.e4.10.degs.matrix(dir.path, gene.pattern, file.pattern = file.pattern)
    
    if (nrow(mat) == 0) {
        warning("No matching genes found.")
        return(NULL)
    }
    
    genes <- rownames(mat)
    
    # Flatten the groups to a data frame
    gene_to_complex <- stack(heatshock_groups)
    colnames(gene_to_complex) <- c("gene", "complex")
    
    # Map genes to complexes
    row_anno <- data.frame(gene = genes)
    row_anno$complex <- gene_to_complex$complex[match(row_anno$gene, gene_to_complex$gene)]
    row_anno$complex[is.na(row_anno$complex)] <- "Other"
    row_anno$complex <- factor(row_anno$complex, levels = unique(row_anno$complex))
    row_anno_vec <- setNames(row_anno$complex, row_anno$gene)
    
    # Define pastel discrete colors for the annotation
    anno_colors <- list(
        Complex = c(
            HSP70      = "#FFB3BA",  # Pastel red
            HSP40      = "#BAE1FF",  # Pastel blue
            HSP90      = "#BAFFB9",  # Pastel green
            HSP60_10   = "#DABFFF",  # Pastel purple
            Small_HSPs = "#FFDFBA",  # Pastel orange
            HSP110     = "#FFC8FF",  # Pastel pink
            Other      = "#E0E0E0"   # Light grey
        )
    )
    
    # Create the row annotation
    ha_row <- rowAnnotation(
        Complex = row_anno_vec,
        col = anno_colors,
        show_annotation_name = TRUE
    )
    
    if (color.sat) {
        col.map <- colorRamp2(c(-1, 0, 1), c(turquoise, "white", orange))
    } else {
        col.map <- colorRamp2(c(min(mat), 0, max(mat)), c(turquoise, "white", orange))
    }
    
    # Order the matrix by annotation group (optional, but helps with split)
    mat <- mat[order(row_anno$complex), , drop = FALSE]
    row_anno_vec <- row_anno_vec[rownames(mat)]
    
    # Define row_split
    row_split <- row_anno_vec
    
    heatmap <- Heatmap(
        mat,
        name = "avg_log2FC",
        column_title = title,
        col = col.map,
        show_row_names = TRUE,
        show_column_names = TRUE,
        cluster_columns = cluster.columns,
        cluster_rows = cluster.rows,
        row_split = row_split,  # Split by annotation group
        row_title_rot = 0,
        border = TRUE,
        na_col = "grey",
        column_names_gp = gpar(fontsize = 8),
        row_names_gp = gpar(fontsize = 8),
        right_annotation = ha_row,
        show_row_dend =F
    )
    
    return(heatmap)
}

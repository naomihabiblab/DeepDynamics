#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#                     Functions of TA & Heatmap Figures                   ####
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# NOTE - All the functions here assume utils.R was sourced before.

pathology.sex.apoe <- function(sex.title = "Bulk, All 1061, Pathologies Vs. Sex",
                               apoe.title = "Bulk, All 1061, Pathologies Vs. APOE4",
                               color.label = "-log(FDR) * sign(beta)",                               
                               control.names = DEFAULT.CONTROLS,
                               trait.names = DEFAULT.PATHOLOGIES) {
    
    # Load bulk (if needed)
    load.full.bulk()
    
    bulk.df <- bulk.df %>% mutate(
        fsex = ifelse(msex, 0, 1),
        n.apoe_4 = ifelse(apoe_4, 0, 1)
    )
    # Example usage
    set.definitions <- list(
        list(name = "sex", covariate.names = c("msex", "fsex")),
        list(name = "apoe", covariate.names = c("apoe_4", "n.apoe_4"))
        # Add more sets as needed
    )
    
    # Run the analysis
    # results <- create.sets.and.run.SA(bulk.metadata, set.definitions, trait.names, control.names)
    results <- create.sets.and.run.SA(bulk.df, set.definitions, trait.names, control.names) #TODO remove
    
    ta.sex <- results$sex
    ta.apoe <- results$apoe
    
    # Create Heatmaps and save:
    p.sex <- ta.heatmap(
        ta.sex,
        sex.title,
        color.label,
        cluster.rows = F,
        cluster.cols = F,
        leg.pos = "bottom",
        fontface = 'bold',
        angle.col = "0"
    )
    
    p.apoe <- ta.heatmap(
        ta.apoe,
        apoe.title,
        color.label,
        cluster.rows = F,
        cluster.cols = F,
        fontface = 'bold',
        leg.pos = "bottom",
        angle.col = "0"
    )
    
    # Change order  of p.sex
    current_order <- p.sex@row_order
    current_names <- rownames(p.sex@matrix)
    
    desired.names <- DEFAULT.PATHOLOGIES
    # After change of names:
    desired.names <- c("Amyloid", "Tau", "Cognitive decline rate")
    new_order <- match(desired.names, current_names)
    
    p.sex@row_order <- new_order
    p.sex@column_order <- c(2,1)
    
    # Change order of p.apoe
    current_order <- p.apoe@row_order
    current_names <- rownames(p.apoe@matrix)
    
    new_order <- match(desired.names, current_names)
    
    p.apoe@row_order <- new_order
    
    # Create the original heatmap plot
    p1 <- grid.grabExpr(draw(p.sex, heatmap_legend_side = "bottom"))
    
    # Create the modified heatmap plot with legend at the bottom
    p2 <- grid.grabExpr(draw(p.apoe, heatmap_legend_side = "bottom"))
    
    message('Full bulk pathology trait association plots have been created.')
    # Arrange both plots side by side
    grid.arrange(p1, p2, ncol = 2)
    
    grid.arrange(
        p1, p2, 
        ncol = 2,
        widths = c(6, 6),  # Adjust the width of each plot
        heights = 5,       # Adjust the height of the plots
        padding = unit(3, "cm")  # Add padding between plots
        #  top = textGrob("Your Overall Title", gp = gpar(fontsize = 16, font = 2))  # Add a title if needed
        # bottom = textGrob("Your footnote", gp = gpar(fontsize = 10, font = 3))  # Add a footnote if needed
    )
}

states.pathologies <- function(data,
                               states = SIG.CLUSTERS,
                               color.label = "-log(adj.pval) * sign(beta)",
                               control.names = DEFAULT.CONTROLS,
                               trait.names = DEFAULT.PATHOLOGIES,
                               title = "Blk, 1,067 Individuals, Pathologies Vs. States",
                               add2title = "",
                               cluster.rows = F,
                               cluster.cols = T,
                               cluster.cols.dist = 'euclidean',
                               return.flag =  F,
                               ...) {
    
    
    # Example usage
    set.definitions <- list(
        list(name = "states", covariate.names = states)
        # Add more sets as needed
    )
    
    # Run the analysis
    results <- create.sets.and.run.SA(data, set.definitions, trait.names, control.names)
    
    ta.states <- results$states
    
    # Create Heatmaps and save:
    p.states <- ta.heatmap(
        ta.states,
        paste(title, add2title),
        color.label,
        cluster.rows = cluster.rows,
        cluster.cols = cluster.cols,
        cluster.cols.dist = cluster.cols.dist,
        leg.pos = "bottom",
        ...
    )
    
    # Change order of p.states
    current_order <- p.states@row_order
    current_names <- rownames(p.states@matrix)
    
    if (all(trait.names == DEFAULT.PATHOLOGIES)) {
        desired.names <- DEFAULT.PATHOLOGIES
        new_order <- match(desired.names, current_names)
        
        p.states@row_order <- new_order
    } else {
        desired.names <- trait.names
        new_order <- match(desired.names, current_names)
        
        p.states@row_order <- new_order
    }
    
    # Change order of rows
    if (cluster.cols == F) {
        current_order <- p.states@column_order
        current_names <- colnames(p.states@matrix) 
        
        desired.names <- states
        new_order <- match(desired.names, current_names)
        
        p.states@column_order <- new_order
    }
    
    
    # Create the original heatmap plot
    p1 <- grid.grabExpr(draw(p.states, heatmap_legend_side = "bottom"))
    
    message('Bulk states - pathology trait association plots have been created.')
    # Arrange both plots side by side
    grid.arrange(p1)
    
    # if (return.flag) {
    #     return(p.states)
    # }
}



# -------------------------------------------------------------------------- #
##                Function for States with SHAP direction                 ####
# -------------------------------------------------------------------------- #
create.SHAP.states.HM <- function(cells, prAD.ls, ABA.ls,
                                  col.fun = colorRamp2(c(-1, 0, 1), c("cyan4", "grey", "darkorange3")),
                                  cluster.rows = F, cluster.cols = F,
                                  show.row.names = T, show.col.names = T) {
    hm.shap.df <- data.frame(
        Cell = cells
    )
    
    hm.shap.df <- hm.shap.df %>%
        mutate(
            prAD = case_when(
                Cell %in% prAD.ls$pos ~ 1, 
                Cell %in% prAD.ls$neg ~ -1, 
                TRUE ~ 0
            ),
            ABA = case_when(
                Cell %in% ABA.ls$pos ~ 1, 
                Cell %in% ABA.ls$neg ~ -1, 
                TRUE ~ 0
            )
        ) %>%
        column_to_rownames('Cell')
    
    # Create the heatmap
    heatmap <- Heatmap(
        as.matrix(hm.shap.df),   # Convert data frame to matrix
        name = "Values",
        col = col.fun,
        cluster_rows = cluster.rows,
        cluster_columns = cluster.cols,
        show_row_names = show.row.names,
        show_column_names = show.col.names,
        heatmap_legend_param = list(title = "Color Legend"),
        border = TRUE,             # Add borders between cells
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.rect(x, y, width, height, gp = gpar(col = "black", fill = NA)) # Draw black rectangles as borders
        }
    )
    
    # Draw the heatmap
    draw(heatmap)
}




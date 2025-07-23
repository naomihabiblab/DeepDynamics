#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#                    Functions to Create Dynamics Figures                 ####
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# _________________________________________________________________________ #
##                       Displaying Temporal KS Scores                   ####
# _________________________________________________________________________ #
# TODO take colors out of function/ as defaults
create.temporal.ks.heatmap <- function(res, 
                                       features,
                                       trajs = c('prAD', 'ABA'),
                                       title = NULL) {
    
    mat_prAD <- prepare.data(res, "prAD")
    mat_ABA <- prepare.data(res, "ABA")
     
     # Create color function
    col_fun <- circlize::colorRamp2(c(min(c(mat_prAD, mat_ABA)), 
                                      median(c(mat_prAD, mat_ABA)), 
                                      max(c(mat_prAD, mat_ABA))), 
                                    c("#4db6ac", "white", "#ef6c00"))
    
    col_fun_prAD <- circlize::colorRamp2(c(min(c(mat_prAD)), 
                                           median(c(mat_prAD)), 
                                           max(c(mat_prAD))), 
                                         c("#4db6ac", "white", "#ef6c00"))
    
    col_fun_ABA <- circlize::colorRamp2(c(min(c(mat_ABA)), 
                                          median(c(mat_ABA)), 
                                          max(c(mat_ABA))), 
                                        c("#4db6ac", "white", "#ef6c00"))
    
    # Step 3: Create heatmaps
    ht_prAD <- Heatmap(mat_prAD,
                       name = "KS statistics (prAD)",
                       col = col_fun_prAD,
                       cell_fun = function(j, i, x, y, width, height, fill) {
                           if (is.na(x)) return(NULL) # Skip NA values
                           #grid.rect(x, y, width, height, gp = gpar(fill = fill)) # Draw the rectangle
                           value = mat_prAD[i, j] # Get the value
                           threshold = median(c(mat_prAD)) # Set the threshold
                           
                           # Set the text properties
                           text_gp = gpar(fontsize = 10)
                           if (value > threshold) {
                               text_gp$fontface = "bold"
                           }
                           
                           # Draw the text
                           grid.text(sprintf("%.2f", value), x, y, 
                                     gp = text_gp)
                       },
                       cluster_rows = TRUE,
                       cluster_columns = FALSE,
                       column_title = "prAD",
                       row_names_side = "left",
                       column_names_rot = 45)
    
    ht_ABA <- Heatmap(mat_ABA,
                      name = "KS statistics (ABA)",
                      col = col_fun_ABA,
                      cluster_rows = TRUE,
                      cluster_columns = FALSE,
                      column_title = "ABA",
                      row_names_side = "left",
                      column_names_rot = 45)
    
    # Step 4: Combine heatmaps
    combined_ht <- ht_prAD + ht_ABA
    
    # Step 5: Draw the combined heatmap
    return(draw(combined_ht, ht_gap = unit(1, "cm")))
    

 }

# Function to prepare the data for the heatmap
prepare.data <- function(res, trajectory, features = c(DEFAULT.PATHOLOGIES, SIG.CLUSTERS)) {
    res %>%
        filter(trajectory == !!trajectory) %>%
        filter(feature %in% features) %>%
        dplyr::select(feature, split, result) %>%
        pivot_wider(names_from = split, values_from = result) %>%
        column_to_rownames("feature") %>%
        dplyr::select(start, mid, end) %>%
        as.matrix()
}

create.ks.score.barplot <- function(ks.res,
                                    features = c(DEFAULT.PATHOLOGIES, SIG.CLUSTERS),
                                    quantile = 0.5,
                                    title = "Scores by Feature",
                                    x_label = "States",
                                    y_label = "KS Temporal Score",
                                    annotation_text = "Quantile:") {
    # Calculate the desired quantile
    quant_val <- quantile(ks.res$result, probs = quantile)
    
    # Filter and order the data
    ks.res <- ks.res %>% 
        filter(feature %in% features) %>%
        mutate(split = factor(split, levels = c( "end", "mid","start")))
    
    # Define pastel colors for the splits
    pastel_colors <- c(
        "start" = start.col,
        "mid" = middle.col,
        "end" = end.col
    )
    
    # Create the bar plot
    # Desired sizes: h=5, w=7
    p.ks.score <- ggplot(ks.res, aes(x = feature, y = result, fill = split)) +
        geom_bar(stat = "identity",
                 color = border.col,
                 position = position_dodge()) +
        coord_flip() +  # Flip coordinates for horizontal bars
        geom_hline(
            yintercept = quant_val,
            linetype = "dotted",
            color = border.col,
            size = 1
        ) +
        labs(title = title, x = x_label, y = y_label) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(
                angle = 45,
                hjust = 1,
                size = 12
            ),
            axis.text.y = element_text(size = 12),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14)
        ) +
        scale_fill_manual(values = pastel_colors) +
        annotate(
            "text",
            x = quant_val,
            y = 1.22,
            label = paste(annotation_text, round(quant_val, 2)),
            vjust = -0.5,
            hjust = 1.1,
            color = border.col
        )
    
    return(p.ks.score)
}


# _________________________________________________________________________ #
##                             Dynamics Heatmap                          ####
# _________________________________________________________________________ #
create.single.heatmap <- function(dyn, traj, features = SIG.CLUSTERS,
                                  type.str = 'snuc', scale.flag = T,
                                  width = 9, cluster.rows = T,
                                  lfc.opt = 'min.max') {
    
    # Set str variables:
    name = paste(traj, type.str, "heatmap")
    col.title = paste("Dynamics", type.str, "in", traj, "Heatmap")
    
    # Select desired columns:
    dyn <- dyn %>%
        py_to_r() %>%
        filter(feature %in% features & trajectory == traj) %>%
        dplyr::select(x, fit, trajectory, feature)
    
    # Arrange for heatmap
    dyn.wide <- dyn %>%
        pivot_wider(names_from = feature, values_from = fit) %>%
        arrange(., x) %>% 
        mutate(type = type.str)
    
    # Set columns to scale
    scale.cols <- setdiff(names(dyn.wide), c("x", "type", "trajectory"))
    
    for.scale <- dyn.wide %>%
        mutate(across(all_of(scale.cols), scale))
    
    # Scale if needed:
    dyn.wide <- dyn.wide %>% 
        dplyr::select(-x, -type, -trajectory) %>%
        scale() %>%
        t()
    
    # Calculate a signed log fold change
    calculate_signed_log_fold_change <- function(start, end) {
        # Use log(abs()) to handle negative values
        log_start <- sign(start) * log(abs(start) + 1e-6)
        log_end <- sign(end) * log(abs(end) + 1e-6)
        
        # Calculate the difference
        return(log_end - log_start)
    }
    
    # Apply the function to your data
    lfc.se <- calculate_signed_log_fold_change(dyn.wide[, 1], dyn.wide[, ncol(dyn.wide)])
    
    calculate_signed_log_fold_change_min_max <- function(row) {
        min_val <- min(row)
        max_val <- max(row)
        
        log_min <- sign(min_val) * log(abs(min_val) + 1)
        log_max <- sign(max_val) * log(abs(max_val) + 1)
        
        return(log_max - log_min)
    }
    
    # Apply the function to each row of dyn.wide
    lfc.mm <- apply(dyn.wide, 1, calculate_signed_log_fold_change_min_max)
    
    row_ha <- rowAnnotation(
        LFC = anno_barplot(ifelse(lfc.opt == 'min.max', lfc.mm, lfc.se), 
                           gp = gpar(fill = ifelse(lfc > 0, "orange3", "steelblue")),
                           width = unit(2, "cm"))
    )
    
    # Create color function:
    col_fun <- circlize::colorRamp2(c(min(dyn.wide), 
                                      median(dyn.wide), 
                                      max(dyn.wide)), 
                                    c(low.col, mid.col, high.col))
    
    if (typeof(cluster.rows) != "logical") {
        flat_order = unlist(cluster.rows)
        
        # Create a dendrogram
        dend = as.dendrogram(hclust(dist(dyn.wide[flat_order, ])))
        
        # Optionally, sort the dendrogram for better visualization
        sorted_dend = dendsort(dend)
        cluster.rows <- sorted_dend
    }
    
    # Create heatmap
    h.dyn <- Heatmap(dyn.wide,
                     cluster_columns = F,
                     col = col_fun,
                     name = name,
                     column_title = col.title,
                     width = unit(width, "cm"),
                     cluster_rows = cluster.rows,
                     right_annotation = row_ha) # Add the annotation here
    
    return(h.dyn)
}


create.snuc.vs.bulk <- function(dyn.s, dyn.b, traj, features = SIG.CLUSTERS,
                                scale.flag = T, width = 9, lfc.opt = 'min.max') {
    
    # Calculate a signed log fold change
    calculate_signed_log_fold_change <- function(start, end) {
        # Use log(abs()) to handle negative values
        log_start <- sign(start) * log(abs(start) + 1e-6)
        log_end <- sign(end) * log(abs(end) + 1e-6)
        
        # Calculate the difference
        return(log_end - log_start)
    }
    
    calculate_signed_log_fold_change_min_max <- function(row) {
        min_val <- min(row)
        max_val <- max(row)
        
        log_min <- sign(min_val) * log(abs(min_val) + 1)
        log_max <- sign(max_val) * log(abs(max_val) + 1)
        
        return(log_max - log_min)
    }
    
    # Set str variables:
    name = paste(traj, "snuc.vs.bulk", "heatmap")
    col.title = paste(traj, "Scaled")
    
    # Select desired columns:
    dyn.s <- dyn.s %>%
        py_to_r() %>%
        filter(feature %in% features & trajectory == traj) %>%
        dplyr::select(x, fit, trajectory, feature)
    
    dyn.b <- dyn.b %>%
        py_to_r() %>%
        filter(feature %in% features & trajectory == traj) %>%
        dplyr::select(x, fit, trajectory, feature)
    
    # Arrange for heatmap
    dyn.s.wide <- dyn.s %>%
        pivot_wider(names_from = feature, values_from = fit) %>%
        arrange(., x) %>% 
        mutate(type = 'snuc')
    
    dyn.b.wide <- dyn.b %>%
        pivot_wider(names_from = feature, values_from = fit) %>%
        arrange(., x) %>% 
        mutate(type = 'bulk')
    
    # Set columns to scale
    state.cols <- setdiff(names(dyn.s.wide), c("x", "type", "trajectory"))
    
    lfc.se.s <- calculate_signed_log_fold_change(dyn.s.wide[1, state.cols], dyn.s.wide[ncol(dyn.s.wide), state.cols])
    lfc.mm.s <- apply(dyn.s.wide[, state.cols], 2, calculate_signed_log_fold_change_min_max)    
    
    # Bulk logfc calculations:
    lfc.se.b <- calculate_signed_log_fold_change(dyn.b.wide[1, state.cols], dyn.b.wide[ncol(dyn.b.wide), state.cols])
    lfc.mm.b <- apply(dyn.b.wide[, state.cols], 2, calculate_signed_log_fold_change_min_max)  
    
    for.scale <- dyn.s.wide %>%
        mutate(across(all_of(state.cols), scale))
    
    # Scale if needed:
    dyn.s.wide <- dyn.s.wide %>% 
        dplyr::select(-x, -type, -trajectory) %>%
        scale() %>%
        t()
    
    dyn.b.wide <- dyn.b.wide %>% 
        dplyr::select(-x, -type, -trajectory) %>%
        scale() %>%
        t()
   
    
    lfc.s <- if(lfc.opt == 'min.max') lfc.mm.s else lfc.se.s
    lfc.b <- if(lfc.opt == 'min.max') lfc.mm.b else lfc.se.b
    
    
    row_ha.s <- rowAnnotation(
        LFC = anno_barplot(lfc.s, 
                           gp = gpar(fill = ifelse(lfc.s > 0, "orange3", "steelblue")),
                           width = unit(2, "cm"))
    )
  
    row_ha.b <- rowAnnotation(
        LFC = anno_barplot(lfc.b, 
                           gp = gpar(fill = ifelse(lfc.b > 0, "orange3", "steelblue")),
                           width = unit(2, "cm"))
    )
    
    combined.dyn <- cbind(dyn.s.wide, dyn.b.wide)
    
    row.clust <- hclust(dist(combined.dyn))
    
    h.s <- Heatmap(dyn.s.wide,
                   cluster_columns = F,
                   col = circlize::colorRamp2(c(min(dyn.s.wide), 
                                                median(dyn.s.wide), 
                                                max(dyn.s.wide)), 
                                              c(low.col, mid.col, high.col)),
                   name = paste("sNuc", col.title),
                   column_title = paste("Dynamics sNuc in", traj, "Heatmap"),
                   width = unit(width, "cm"),
                   cluster_rows = row.clust,
                   left_annotation = row_ha.s,
                   use_raster = T)
    
    h.b <- Heatmap(dyn.b.wide,
                   cluster_columns = F,
                   col = circlize::colorRamp2(c(-4, 
                                                0, 
                                                4), 
                                              c(low.col, mid.col, high.col)),
                   name = paste("Bulk", col.title),
                   column_title = paste("Dynamics Bulk in", traj, "Heatmap"),
                   width = unit(width, "cm"),
                   cluster_rows = row.clust,
                   right_annotation = row_ha.b)
    
    h.b <- Heatmap(dyn.b.wide,
                   cluster_columns = F,
                   col = circlize::colorRamp2(c(-4,
                                                0,
                                                4),
                                              c(low.col, mid.col, high.col)),
                   name = paste("Bulk", col.title),
                   column_title = paste("Dynamics Bulk in", traj, "Heatmap"),
                   width = unit(width, "cm"),
                   cluster_rows = row.clust,
                   use_raster = T)
    
    return(list(h.s, h.b))
    
}
    
    
    
plot.genes.dynamcis.heatmap <- function(dyn, traj = 'prAD', row.clust = T, 
                                        low.col = "#023047", mid.col = "#219ebc", high.col = "#fb8500", 
                                        width = 9, col.title = "Gene Dynamics Heatmap") {
    # Define features as a string list that includes _scaled in them, without mean_scaled in them:
    
    # Select desired columns:
    dyn <- dyn %>%
        py_to_r() %>%
        filter(grepl("_scaled", feature) & !grepl("mean_scaled", feature) & trajectory == traj) %>%
        dplyr::select(x, fit, trajectory, feature, split)
    
    # Arrange for heatmap
    dyn.wide <- dyn %>%
        pivot_wider(names_from = feature, values_from = fit) %>%
        arrange(., x) 
    
    # Select columns for each APOE4 group
    gene.cols.0 <- grep("_scaled.*E4\\-$", colnames(dyn.wide), value = TRUE)
    gene.cols.1 <- grep("_scaled.*E4\\+$", colnames(dyn.wide), value = TRUE)
    
    # Create matrices for the heatmaps
    mat.0 <- as.matrix(dyn.wide[, gene.cols.0] %>% remove_missing() %>% t())
    mat.1 <- as.matrix(dyn.wide[, gene.cols.1] %>% remove_missing() %>% t())
    
    
    # Heatmap for APOE4 non-carriers
    h.0 <- Heatmap(mat.0,
                   cluster_columns = FALSE,
                   col = circlize::colorRamp2(
                       c(min(mat.0, mat.1, na.rm = TRUE), 
                         median(mat.0, mat.1, na.rm = TRUE), 
                         max(mat.0, mat.1, na.rm = TRUE)), 
                       c(low.col, mid.col, high.col)),
                   name = paste("Non-carriers", col.title),
                   column_title = paste("Gene dynamics", traj, "non-carriers"),
                   width = unit(width, "cm"),
                   cluster_rows = row.clust,
                   use_raster = TRUE)
    
    # Heatmap for APOE4 carriers
    h.1 <- Heatmap(mat.1,
                   cluster_columns = FALSE,
                   col = circlize::colorRamp2(
                       c(min(mat.1, mat.0, na.rm = TRUE), 
                         median(mat.1, mat.0, na.rm = TRUE), 
                         max(mat.1, mat.0, na.rm = TRUE)), 
                       c(low.col, mid.col, high.col)),
                   name = paste("Carriers", col.title),
                   column_title = paste("Gene dynamics", traj, "carriers"),
                   width = unit(width, "cm"),
                   cluster_rows = row.clust,
                   use_raster = TRUE)
    
    return(list(h.0, h.1))
}
    



# _________________________________________________________________________ #
##                Functions to plot2 or 4 groups of dynamics             ####
# _________________________________________________________________________ #

plot.4g.split <- function(state, pts.flag = FALSE, 
                          trajs = c('prAD', 'ABA'),
                          feature.name.flag = TRUE,
                          orientation = 'row') {
    
    orientation.lst <- list(
        'row' = list(nrow = 1),
        'col' = list(ncol = 1)
    )
    
    args <- orientation.lst[[orientation]]
    
    new.ylab <- if (state %in% c('sqrt.amyloid_mf', 'sqrt.tangles_mf', 'cogng_demog_slope')) 'Abundance' else 'Prevalence'
    new.xlab <- 'Pseudotime'
    
    plot.single.state <- function(f, risk.name = 'APOE4_genotype') {  # In old - ApoE4
        fts.e4 <- paste0(f, '.', gsub("_genotype", "", risk.name), c('M+', 'F+'))
        fts.e3 <- paste0(f, '.', gsub("_genotype", "", risk.name), c('M-', 'F-'))
        
        get.preds <- function(features) {
            data$uns$trajectories$celmod.e4.sex.new.splits$pred.vals %>%
                py_to_r() %>%
                filter((feature %in% features) & (trajectory %in% trajs))
        }
        
        preds.4 <- get.preds(fts.e4)
        preds.3 <- get.preds(fts.e3)
        
        ymin <- min(preds.4$fit - 2*preds.4$se.fit, 
                    preds.3$fit - 2*preds.3$se.fit, 
                    na.rm = TRUE)
        ymax <- max(preds.4$fit + 2*preds.4$se.fit, 
                    preds.3$fit + 2*preds.3$se.fit, 
                    na.rm = TRUE)
        
        # Calculate overall x-axis range
        xmin <- min(preds.4$x, preds.3$x, na.rm = TRUE)
        xmax <- max(preds.4$x, preds.3$x, na.rm = TRUE)
        
        
        get.minmax <- function(features) {
            data$uns$trajectories$celmod.e4.sex.new.splits$fitted.vals %>%
                py_to_r() %>%
                filter((feature %in% features) & (trajectory %in% trajs)) %>%
                summarise(min = min(fit), max = max(fit))
        }
        
        minmax.e4 <- get.minmax(fts.e4)
        minmax.e3 <- get.minmax(fts.e3)
        
        overall.min <- min(minmax.e4$min, minmax.e3$min)
        overall.max <- max(minmax.e4$max, minmax.e3$max)
        
        plot.dyn.couple <- function(features, colors) {
            plot.dynamics.wrapper(
                data$uns$trajectories$celmod.e4.sex.new.splits,
                features,
                scale = "free_x",
                include.points = pts.flag,
                overlap.pseudotime = .1,
                cols = setNames(colors, features),
                label = FALSE,
                ymin = ymin,
                ymax = ymax,
                xmin = xmin,
                xmax = xmax,
                trajectories = trajs
            ) +
                geom_hline(yintercept = c(overall.min, overall.max), 
                           linetype = "longdash", 
                           color = "darkgrey", size = 0.5) +
                labs(x = new.xlab, y = new.ylab) +
                DYN.THEME +
                coord_cartesian(xlim = c(xmin, xmax)) 
        }
        
        p.plus <- plot.dyn.couple(fts.e4, c('#fb8500','#ffb703'))
        p.minus <- plot.dyn.couple(fts.e3, c('#023047','#219ebc'))
        
        f <- switch(f,
                    'sqrt.tangles_mf' = 'Tangles',
                    'sqrt.amyloid_mf' = 'Amyloid',
                    'cogng_demog_slope' = 'Cognitive Decline',
                    f)
        
        if (feature.name.flag) {
            
            do.call(
                plot_grid, 
                c(list(p.plus, p.minus), args, labels = paste(f, "Sex and E4", c('+', '-')))
                )
            
        } else {
            do.call(
                plot_grid, 
                c(list(p.plus, p.minus), args)
                )
            
        }
    }
    
    lapply('ApoE4_genotype', function(risk.name) {
        plot_grid(plotlist = lapply(state, plot.single.state), nrow = 1)
    })
}




plot.2g.split <- function(features, split = 'E4', pts.flag = F, trajs = c('prAD', 'ABA'), 
                          fit.section.e4 = "celmod.all.e4.new.split", fit.section.sex = "celmod.all.sex.new.split",
                          feature.name.flag = T, nrow = 1) {
    new.xlab = 'Pseudotime'
    
    p <- plot_grid(
        plotlist = lapply(c(features), function(f) {
            if (split == 'E4') {
                flist <- list(f = paste0(f, c('.E4+', '.E4-'))) 
                data.split <- data$uns$trajectories[[fit.section.e4]]
            } else if (split == 'Sex') {
                flist <- list(f = paste0(f, c('.Female', '.Male')))
                data.split <- data$uns$trajectories[[fit.section.sex]]
            }
            if (f %in% c('sqrt.amyloid_mf', 'sqrt.tangles_mf', 'cogng_demog_slope')) {
                new.ylab = 'Abundance'
            } else {
                new.ylab = 'Prevalance'
            }
            
            p.celmod <-
                plot.dynamics.wrapper(
                    data.split,
                    flist$f,
                    scale = "free_x",
                    include.points = pts.flag,
                    overlap.pseudotime = .1,
                    trajectories = trajs,
                    label = F,
                    cols = setNames(c("#ef6c00","#4db6ac"), flist$f)
                ) +
                labs(x = new.xlab, y = new.ylab) +
                DYN.THEME +
                theme(axis.text.y = element_text(angle=90))  #To get the width the same for each sub figure
            
            if (f == 'sqrt.tangles_mf') {
                f <- 'Tangles'
            } else if (f == 'sqrt.amyloid_mf') {
                f <- 'Amyloid'
            } else if (f == 'cogng_demog_slope') {
                f <- 'Cognitive Decline'
            }
            
            if (feature.name.flag) {
                return(plot_grid(p.celmod,
                                 labels = c(paste(f, ifelse(split == 'E4', "ApoE4", "Sex")))))
            } else {
                return(p.celmod)
            }
            
        }), nrow = nrow
    )
    return(p)
}


# _________________________________________________________________________ #
##                   Dynamics without split (communities)               ####
# _________________________________________________________________________ #
plot.nosplit.bulk <- function(names, colors, 
                                  pts.flag = F, trajs = c('prAD', 'ABA')) {
    THEME <- DYN.THEME
    # Create the plot
    p.dyn <- plot.dynamics.wrapper(
        data$uns$trajectories$celmod.fits.nosplit,
        names,
        facet.by = 'trajectory',
        scale = "free_x",
        include.points = pts.flag,
        trajectories = trajs,
        overlap.pseudotime = 0.1,
        label = FALSE,
        cols = setNames(colors, names)
    ) +
        labs(x = 'Pseudotime', y = 'Prevalence') +
        THEME +
        theme(axis.text.y = element_text(angle=90))
    
    return(p.dyn)
}

plot.nosplit.snuc <- function(names, colors,
                              pts.flag = F, trajs = c('prAD', 'ABA'),
                              data.slot = "snuc.fits.nosplit") {
    if (length(trajs) == 2) {
        THEME <- DYN.THEME
    } else {
        THEME <- DYN.THEME
    }
    # Create the plot
    p.dyn <- plot.dynamics.wrapper(
        data$uns$trajectories[[data.slot]],
        names,
        facet.by = 'trajectory',
        scale = "free_x",
        include.points = pts.flag,
        trajectories = trajs,
        overlap.pseudotime = 0.1,
        label = FALSE,
        cols = setNames(colors, names)
    ) +
        labs(x = 'Pseudotime', y = 'Prevalence') +
        THEME +
        theme(axis.text.y = element_text(angle=90))
    
    return(p.dyn)
}




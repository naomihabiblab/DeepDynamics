##############################################################################
#                       Trait Association Functions                          #
##############################################################################
# Consts (used as default values):
orange <- "#ef6c00"
turqouise <- "#4db6ac"

DEFAULT.PATHOLOGIES <- c('sqrt.amyloid_mf', 'sqrt.tangles_mf', 'cogng_demog_slope')
DEFAULT.CONTROLS <- c("age_death", "pmi", "RIN")

LABEL.SIZE = 10

# To define SIG.CLUSTERS, we need to first load the 500 and then delete it.
data <- anndata::read_h5ad(PATH.500had)

tryCatch({
  SIG.CLUSTERS <- data$uns$celmod$test.corrs %>% #old.data$uns$celmod$cor
    py_to_r() %>% filter(adj.pval < .005 & corr > 0) %>%
    rownames()
}, error = function(e) {
  SIG.CLUSTERS <- c('Oli.11', 'Mic.12', 'Oli.7', 'Ast.10', 'Mic.13',
                    'Ast.7', 'Oli.3', 'Oli.4', 'Ast.4', 'Inh.6',
                    'Exc.1', 'Exc.3', 'Inh.15', 'OPC.1', 'Inh.16',
                    'Inh.5', 'Exc.8', 'End.3', 'Inh.12', 'Oli.5',
                    'Mic.1', 'Exc.12', 'Inh.7', 'Ast.6', 'End.1',
                    'Ast.2', 'OPC.2', 'Mic.14', 'Ast.1', 'Oli.9') # (only SIG for TA)
})

rm(data)

##############################################################################
#                       Statistical Assocaition                              #
##############################################################################
#' Perform trait associations
#'
#' Associate co-variates to traits will controlling for technical variables
#'
#' @param traits dataframe of shape (n_samples, n_traits)
#' @param covariates dataframe of shape (n_samples, n_covariates)
#' @param controls dataframe of shape (n_samples, n_controls)
#'
#' @return dataframe with trait association results for every n_traits*n_covariates
#' while correcting for multiple hypothesis testing within each of the given traits
#'
associate.traits <- function(traits, covariates, controls, p.adjust.method = "BH") {
    df <- data.frame(traits, covariates, controls)
    
    df <- do.call(rbind, lapply(colnames(traits), function(trait)
        do.call(rbind, lapply(colnames(covariates), function(covariate) {
            
            # Set controlled covariates
            print(paste(trait, covariate))
            control <- colnames(controls)
            control <- case_when(trait == "cww" ~ paste(setdiff(control, "msex"), collapse = " + "),
                                 T ~ paste(control, collapse = " + "))
            
            formula <- stringr::str_interp("${trait} ~ ${covariate} + ${control}")
            m <- summary(lm(formula, df[!is.na(df[, trait]),]))
            data.frame(trait = trait,
                       covariate = covariate,
                       beta = m$coefficients[covariate, 1],
                       se = m$coefficients[covariate, 2],
                       tstat = m$coefficients[covariate, 3],
                       pval = m$coefficients[covariate, 4],
                       r.sq = m$adj.r.squared,
                       n = sum(!is.na(df[, trait])),
                       formula = formula)
        }))))
    print(paste('number of rows:', nrow(df)))
    return(df %>%
               group_by(trait) %>%
               mutate(adj.pval = p.adjust(pval, method = p.adjust.method),
                      sig = cut(adj.pval, c(-.1, 0.001, 0.01, 0.05, Inf), c("***", "**", "*", ""))) %>%
               ungroup())
}


# Function for specific statistical association:
create.sets.and.run.SA <- function(data, set.definitions, trait.names, control.names) {
    # Create sets dynamically based on input
    sets <- lapply(set.definitions, function(def) {
        list(
            covariates = data[, def$covariate.names, drop = FALSE],
            traits = data[, trait.names, drop = FALSE],
            controls = data[, control.names, drop = FALSE]
        )
    })
    
    # Assign names to the sets
    names(sets) <- sapply(set.definitions, function(def) def$name)
    
    # Calculate trait association
    ta <- sapply(names(sets), function(n)
        associate.traits(sets[[n]]$traits, sets[[n]]$covariates, sets[[n]]$controls),
        simplify = FALSE, USE.NAMES = TRUE)
    
    message("Statistical association has been calculated successfully.")
    return(ta)
}


##############################################################################
#                            Functions for Plots                             #
##############################################################################
# Assign color scale
color.pheatmap.neg.pos <- function(df, colors = c(turqouise, "white", orange)) {
    paletteLength <- 100
    colors_step <- colorRampPalette(colors)(n = paletteLength)
    breaks <- c(
        seq(-max(abs(df), na.rm = T), 0, length.out = 50),
        seq(0.01, max(abs(df), na.rm = T), length.out = 50)
    )
    return(list(colors = colors_step, breaks = breaks))
}

# TWO OPTIONS FOR COLORS OF HEATMAP
# Return log of p-value with sign, where significant in the map:
return.signed.pvalues <- function(ta) {
    beta <- dcast(ta, trait ~ covariate, value.var = "beta") %>% column_to_rownames("trait")
    sig <- dcast(ta, trait ~ covariate, value.var = "sig") %>% column_to_rownames("trait")
    adj.pval <- dcast(ta, trait ~ covariate, value.var = "adj.pval") %>% column_to_rownames("trait")
    p.signed <- -log10(adj.pval) * sign(beta)
    return(p.signed)
}

# Return the coefs:
return.coefs <- function(ta) {
    beta <- dcast(ta, trait ~ covariate, value.var = "beta") %>% column_to_rownames("trait")
    return(beta)
}

# Return where to put stars and how many
return.stars.pv <- function(ta) {
    stars <- dcast(ta, trait ~ covariate, value.var = "sig") %>% column_to_rownames("trait")
    return(stars)
}

# Plot heatmap of trait associations
ta.heatmap <- function(ta, title, color.legend, 
                       sigonly = F, color.beta = F,
                       cluster.rows = T, cluster.cols = T,
                       leg.pos = "right", cluster.cols.dist = 'euclidean',
                       fontface = 'plain', hm.width = 5, hm.height = 7.5,
                       angle.col = "90", fontsize = LABEL.SIZE, ...){
    return.psigned <- ifelse(color.beta, 
                             return.coefs, 
                             return.signed.pvalues) # Choose which function of color to do
    if (sigonly) {
        p.signed <- return.psigned(ta %>% filter(sig != "")) %>% replace(is.na(.), 0)
        cb.list <- color.pheatmap.neg.pos(p.signed)
        stars <- return.stars.pv(ta %>% filter(sig != ""))
        stars <- stars %>% replace(is.na(.), "")
    } else {
        p.signed <- return.psigned(ta)
        cb.list <- color.pheatmap.neg.pos(p.signed)
        stars <- return.stars.pv(ta)
    }
    
    # Define dictionaries for renaming
    row_dict <- c(
        "sqrt.amyloid_mf" = "Amyloid",
        "sqrt.tangles_mf" = "Tangles",
        "cogng_demog_slope" = "Cognitive Decline"
    )
    
    col_dict <- c(
        "msex" = "Males",
        "fsex" = "Females",
        "apoe_4" = "ApoE4\nCarriers",
        "apoe_3" = "ApoE4\nNon-Carriers"
    )
    
    # Check if rownames match DEFAULT.PATHOLOGIES and colnames match col_dict names
    if (all(rownames(p.signed) %in% names(row_dict)) && 
        all(colnames(p.signed) %in% names(col_dict))) {
        
        # Rename rows
        rownames(p.signed) <- row_dict[rownames(p.signed)]
        
        # Rename columns
        colnames(p.signed) <- col_dict[colnames(p.signed)]
        
        # Also update the stars matrix
        rownames(stars) <- rownames(p.signed)
        colnames(stars) <- colnames(p.signed)
    }
    
    p <- ComplexHeatmap::pheatmap(
        p.signed,
        cell_fun = function(j, i, x, y, w, h, fill) grid.text(stars[i,j], x,y),
        color = cb.list$colors,
        breaks = cb.list$breaks,
        show_colnames = TRUE,
        cluster_rows = cluster.rows,
        cluster_cols = cluster.cols,
        show_column_dend = FALSE,
        clustering_distance_cols = cluster.cols.dist,
        heatmap_legend_param = list(
            title_position = "topcenter",
            legend_direction = "horizontal",
            legend_width = unit(0.5, "npc"),
            legend_height = unit(0.5, "npc"),
            title_gp = gpar(fontsize = 10),
            title = color.legend
        ),
        # cellwidth = unit(4, "cm"),
        # cellheight = unit(4, "cm"),
        na_col = 'grey',
        main = title,
        treeheight_col = unit(0.5, "cm"),
        width = hm.width,
        height = hm.height,
        fontface = fontface,
        angle_col = angle.col
    )
    
    p@matrix_color_mapping@name <- color.legend
    return(p)   
}



#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#           Functions to plot proteomics data in different ways           ####
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# -------------------------------------------------------------------------- #
##                             Plot Functions                             ####
# -------------------------------------------------------------------------- #
create.density.plot <- function(df, gene, p.val, data.type = "prot", uniprot = NULL) {
    # Check if uniprot is NULL and data.type is prot, stop and raise error:
    if (is.null(uniprot) && data.type == "prot") {
        stop("Error: uniprot cannot be NULL when data.type is 'prot'")
    } else if (data.type == "bulk" && !is.null(uniprot)) {
        warning("Warning: uniprot is not NULL when data.type is 'bulk'")
    }
    
    # calc apoe4 genotype conuts:
    genotype.counts <- df %>%
        filter(!(apoe_genotype %in% EXCLUDED.GENOTYPES)) %>%
        group_by(apoe_genotype) %>%
        summarise(count = n()) %>%
        mutate(label = paste0(apoe_genotype, " (n=", count, ")"))
    
    # Create title:
    if (data.type == "prot") {
        p.title <- paste0(gene, 
                          " (", uniprot, ")", 
                          "\np-value = ", 
                          format(p.val, scientific = TRUE, digits = 3))
        # Create density plot:
        ggplot(df %>%
                   filter(!(apoe_genotype %in% EXCLUDED.GENOTYPES)),
               aes(gene, fill = as.factor(apoe_genotype), color = as.factor(apoe_genotype))) +
            geom_density(alpha = 0.2) +
            ggtitle(p.title) +
            labs(fill = "APOE Genotype", color = "APOE Genotype") +
            scale_fill_discrete(labels = genotype.counts$label) +
            scale_color_discrete(labels = genotype.counts$label) +
            DENSITY.THEME +
            labs(x = ifelse(data.type == "prot", "Protein Abundance (Normalized)", "Bulk Gene Expression (Normalized)"), 
                 y = "Density")
        
    } else if (data.type == "bulk") {
        p.title <- paste0(gene, " Bulk Expression",
                          "\np-value = ", 
                          format(p.val, scientific = TRUE, digits = 3))
        
        ggplot(df %>%
                   remove_missing(na.rm = T) %>%
                   filter(!(apoe_genotype %in% EXCLUDED.GENOTYPES)),
               aes_string(gene, fill = "as.factor(apoe_genotype)", color = "as.factor(apoe_genotype)")) +
            # geom_bar(position = "identity", alpha = 0.2, bins = 30, color = "black") +
            geom_density(alpha = 0.2) +
            ggtitle(p.title) +
            labs(fill = "APOE Genotype", color = "APOE Genotype") +
            scale_fill_discrete(labels = genotype.counts$label) +
            scale_color_discrete(labels = genotype.counts$label) +
            DENSITY.THEME +
            labs(x = ifelse(data.type == "prot", "Protein Abundance (Normalized)", "Bulk Gene Expression (Normalized)"), 
                 y = "Density")
    }
    
    
    
        
}

create.density.plot.split.by <- function(df, gene, p.val = NULL, data.type = "prot", 
                                         uniprot = NULL, split.by = NULL, pval.test.opt = "two.sided") {
    if (data.type == "prot" && is.null(uniprot)) {
        stop("Error: uniprot cannot be NULL when data.type is 'prot'")
    }
    
    if (data.type == "bulk" && !is.null(uniprot)) {
        warning("Warning: uniprot is not NULL when data.type is 'bulk'")
    }
    
    # Filter excluded genotypes
    df <- df %>% filter(!(apoe_genotype %in% EXCLUDED.GENOTYPES))
    
    # If no split.by, revert to regular plot
    if (is.null(split.by)) {
        return(create.density.plot(df, gene, p.val, data.type = data.type, uniprot = uniprot))
    }
    
    # Calculate counts and p-values per level of split.by
    split.levels <- unique(df[[split.by]])
    
    p.vals <- sapply(split.levels, function(lvl) {
        df.sub <- df[df[[split.by]] == lvl, ]
        if (nrow(df.sub) == 0) return(NA)
        test <- calc.genotype.pval(df.sub, opt = pval.test.opt)
        test$p.value
    }, simplify = TRUE, USE.NAMES = TRUE)
    
    # Prepare title line with flexible format
    pval.str <- paste0(names(p.vals), " p = ", format(p.vals, digits = 3), collapse = ", ")
    
    p.title <- if (data.type == "prot") {
        paste0(gene, " (", uniprot, ")", "\n", pval.str)
    } else {
        paste0(gene, " Bulk Expression\n", pval.str)
    }
    
    # Create genotype counts for legend
    genotype.counts <- df %>%
        group_by(apoe_genotype, .data[[split.by]]) %>%
        summarise(count = n(), .groups = 'drop') %>%
        pivot_wider(names_from = split.by, values_from = count, values_fill = 0) %>%
        mutate(label = paste0(apoe_genotype, " (", 
                              paste(paste0(names(.)[3:ncol(.)], "=", .[3:ncol(.)]), collapse = ", "), 
                              ")"))
    
    label.map <- setNames(genotype.counts$label, genotype.counts$apoe_genotype)
    
    # Create plot
    ggplot(df, aes(gene, fill = as.factor(apoe_genotype), color = as.factor(apoe_genotype))) +
        geom_density(alpha = 0.2) +
        facet_wrap(as.formula(paste("~", split.by)), 
                   labeller = "label_both") +
        ggtitle(p.title) +
        labs(fill = "APOE Genotype", color = "APOE Genotype") +
        scale_fill_discrete(labels = label.map[levels(as.factor(df$apoe_genotype))]) +
        scale_color_discrete(labels = label.map[levels(as.factor(df$apoe_genotype))]) +
        DENSITY.THEME +
        labs(x = ifelse(data.type == "prot", "Protein Abundance (Normalized)", "Bulk Gene Expression (Normalized)"),
             y = "Density") + 
}

create.density.plot.split.by <- function(df, gene, p.val = NULL, data.type = "prot", 
                                         uniprot = NULL, split.by = NULL, 
                                         pval.test.opt = "two.sided") {
    # Error checks
    if (data.type == "prot" && is.null(uniprot)) {
        stop("Error: uniprot cannot be NULL when data.type is 'prot'")
    }
    if (data.type == "bulk" && !is.null(uniprot)) {
        warning("Warning: uniprot is not NULL when data.type is 'bulk'")
    }
    
    # Filter out excluded genotypes
    df <- df %>% filter(!(apoe_genotype %in% EXCLUDED.GENOTYPES))
    
    # No split.by: revert to regular density plot
    if (is.null(split.by)) {
        return(create.density.plot(df, gene, p.val, data.type = data.type, uniprot = uniprot))
    }
    
    # Calculate p-values by split level (e.g., sex)
    split.levels <- unique(df[[split.by]])
    p.vals <- sapply(split.levels, function(lvl) {
        df.sub <- df[df[[split.by]] == lvl, ]
        if (nrow(df.sub) == 0) return(NA)
        test <- calc.genotype.pval(df.sub, opt = pval.test.opt)
        test$p.value
    }, simplify = TRUE, USE.NAMES = TRUE)
    
    # Format p-value string
    pval.str <- paste0(names(p.vals), " p = ", format(p.vals, digits = 3), collapse = ", ")
    p.title <- if (data.type == "prot") {
        paste0(gene, " (", uniprot, ")\n", pval.str)
    } else {
        paste0(gene, " Bulk Expression\n", pval.str)
    }
    
    # Clean legend formatting (e.g., "23 (Females = 45, Males = 12)")
    genotype.counts <- df %>%
        group_by(apoe_genotype, !!sym(split.by)) %>%
        summarise(count = n(), .groups = "drop") %>%
        mutate(split_label = case_when(
            !!sym(split.by) == 0 ~ "Females",
            !!sym(split.by) == 1 ~ "Males",
            TRUE ~ paste0(split.by, " = ", !!sym(split.by))
        )) %>%
        group_by(apoe_genotype) %>%
        summarise(label = paste0(apoe_genotype, " (", 
                                 paste0(split_label, " = ", count, collapse = ", "), 
                                 ")"))
    
    # Map for plot labels
    label.map <- setNames(genotype.counts$label, genotype.counts$apoe_genotype)
    
    # Create density plot
    ggplot(df, aes(gene, fill = as.factor(apoe_genotype), color = as.factor(apoe_genotype))) +
        geom_density(alpha = 0.2) +
        facet_wrap(as.formula(paste("~", split.by)), labeller = "label_both") +
        ggtitle(p.title) +
        labs(fill = "APOE Genotype", color = "APOE Genotype") +
        scale_fill_discrete(labels = label.map) +
        scale_color_discrete(labels = label.map) +
        DENSITY.THEME +
        labs(x = ifelse(data.type == "prot", "Protein Abundance (Normalized)", 
                        "Bulk Gene Expression (Normalized)"), 
             y = "Density") +
        theme(legend.position = "right")
}# -------------------------------------------------------------------------- #

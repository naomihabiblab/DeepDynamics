##############################################################################
#                  Functions for Data Presentation Figures                   #
##############################################################################
# NOTE - All the functions here assume utils.R was sourced before.

##############################################################################
# Pie Chart Options #
##############################################################################

create.donut.chart <- function(dataframe, 
                               variable,
                               border.col = border.col,
                               pallete.col = apoe4.pallete.col,
                               title = NULL) {
    # Ensure the variable exists in the dataframe
    if (!variable %in% names(dataframe)) {
        stop(paste("Error: Variable", variable, "not found in the dataframe"))
    }
    
    # Create frequency table
    freq.table <- table(dataframe[[variable]])
    total <- sum(freq.table)
    
    # Create dataframe for plotting
    df <- data.frame(
        category = names(freq.table),
        value = as.numeric(freq.table),
        percentage = as.numeric(freq.table) / total
    )
    
    # Calculate the positions for the labels
    df$ypos <- cumsum(df$value) - 0.5 * df$value
    
    # Create labels for legend
    df$legend.label <- sprintf("%s (%d, %.1f%%)", 
                               df$category, 
                               df$value, 
                               df$percentage * 100)
    
    # Create the plot
    p <- ggplot(df, aes(x = 2, y = value, fill = legend.label)) +
        geom_bar(stat = "identity", width = 1, color='grey25') +
        coord_polar("y", start = 0) +
        theme_void() +
        scale_fill_brewer(palette = "Pastel1") +
        ggtitle(ifelse(is.null(title), paste("Donut Chart of", variable), title)) +
        theme(legend.position = "right") +
        xlim(0.5, 2.5) + # This creates the donut hole
        guides(fill = guide_legend(title = "ApoE4 Genotype"))
    
    return(p)
}


##############################################################################
#                             Distribution Plots                             #
##############################################################################

createDistributionPlots <- function(data, variable, groupBy = "cohort") {
    # Calculate the counts and percentages
    plot_data <- data %>%
        group_by(across(all_of(c(groupBy, variable)))) %>%
        summarise(count = n(), .groups = "drop") %>%
        group_by(across(all_of(groupBy))) %>%
        mutate(percentage = count / sum(count) * 100,
               label = paste0(count, "\n(", round(percentage, 1), "%)"))
    
    # Reverse the order of levels
    plot_data[[variable]] <- factor(plot_data[[variable]], 
                                    levels = rev(levels(factor(plot_data[[variable]]))))
    
    # Calculate total counts for each group
    total_counts <- plot_data %>%
        group_by(across(all_of(groupBy))) %>%
        summarise(total = sum(count), .groups = "keep")
    
    # Define pastel colors (reversed order)
    pastel_colors <- rev(c("#FFB3BA", "#BAE1FF"))
    
    # Create the plot
    p <- ggplot() +
        geom_bar(data = plot_data, 
                 aes(x = .data[[groupBy]], y = count, fill = .data[[variable]]),
                 stat = "identity", position = "stack", color = border.col) +
        geom_text(data = plot_data,
                  aes(x = .data[[groupBy]], y = count, label = label),
                  position = position_stack(vjust = 0.5), size = 3) +
        geom_text(data = total_counts,
                  aes(x = .data[[groupBy]], y = total, label = paste("Total:", total)),
                  vjust = -0.5, size = 4, fontface = "bold") +
        scale_fill_manual(values = pastel_colors) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "right") +
        labs(title = paste("Distribution of", variable),
             x = groupBy,
             y = "Count",
             fill = variable)
    
    return(p)
}



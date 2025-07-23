#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#                                  Figure 2                               ####
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# Load utilities:
source("src/utils/utils.R")
source("src/utils/utils.dynamics.R")
source("src/visualization/fig.dynamics.funs.R")

PATH.F2 <- file.path(PATH.FIGS, "F2")

# -------------------------------------------------------------------------- #
##                        Panel A: Heatmap of Dynamics                    ####
# -------------------------------------------------------------------------- #
Map(function(traj){
    h.ls <- create.snuc.vs.bulk(dynamics.no.split.snuc, 
                                dynamics.no.split.bulk,
                                traj = traj)
    h <- h.ls[[1]] + h.ls[[2]]
    file.name <- file.path(PATH.F2, paste0("p.A.heatmap.", traj, ".pdf"))
    pdf(file.name, 
        width = COMPARE.HEATMAP.SIZE$w, 
        height = COMPARE.HEATMAP.SIZE$h)
    print(h)
    dev.off()
    
    message('Successfully saved heatmap of dynamics as', file.name)
}, c('prAD', 'ABA'))

# Panel B can be reproduced using the script in explainability folder.

# -------------------------------------------------------------------------- #
##                 Panels C - F: SHAP Communities of Dynamics             ####
# -------------------------------------------------------------------------- #
# Define the SHAP communities:
names <- list(
    p.C = c('Ast.10', 'Oli.7', 'Mic.12'),
    p.D = c('Ast.5', 'OPC.3', 'Mic.7', 'Mic.6', 'Peri.1', 'Mic.8'),
    p.E1 = c('Ast.1', 'Ast.2', 'Mic.2'),
    p.E2 = c('OPC.2', 'Mic.3', 'Peri.2', 'Mic.4'),
    p.F = c('Oli.3', 'Oli.4')
)

# Define colors:
colors <- list(
    p.C = c("darkslategrey", "cyan4", "darkorange3"),
    p.D = c("darkslategrey", "cyan4", "darkorange3", 
            "goldenrod", "darkolivegreen", "brown4"),
    p.E1 = c("brown", "goldenrod4", "darkolivegreen"),
    p.E2 = c("darkslategrey", "cyan4", "darkorange3", 'goldenrod'),
    p.F = c("darkslategrey", "cyan4")
)

Map(function(names.ls, colors.ls, panel.name) {
    # Generate the plot
    plot <- plot.nosplit.bulk(names.ls, colors.ls)
    
    # Define the file path for saving
    file.path <- file.path(PATH.F2, paste0(panel.name, ".pdf"))
    
    # Save the plot as a PDF
    pdf(file.path, width = DYN.SIZE$w, height = DYN.SIZE$h)
    print(plot)  # Print the plot to the PDF device
    dev.off()    # Close the PDF device
}, 
names, colors, names(names)) 

# -------------------------------------------------------------------------- #
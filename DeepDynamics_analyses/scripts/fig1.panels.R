#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#                                  Figure 1                               ####
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# Load utilities:
source("src/utils/utils.R")
source("src/utils/utils.dynamics.R")
source("src/visualization/fig.dynamics.funs.R")

PATH.F1 <- file.path(PATH.FIGS, "F1")

# -------------------------------------------------------------------------- #
##                     Panel C: Distribution of APOE and sex              ####
# -------------------------------------------------------------------------- #
#TODO tbd
# -------------------------------------------------------------------------- #
##                           Panel F: sNuc vs. Bulk                       ####
# -------------------------------------------------------------------------- #
names <- c('Ast.1')
colors.snuc <- c("dodgerblue4")
colors.bulk <- c("darkslategrey")

# Generate the plots
p.snuc <- plot.nosplit.snuc(names, colors.snuc, pts.flag = T)
p.bulk <- plot.nosplit.bulk(names, colors.bulk, pts.flag = T)

# Define the file paths for saving
file.path.snuc <- file.path(PATH.F1, "p.F.snuc.pdf")
file.path.bulk <- file.path(PATH.F1, "p.F.bulk.pdf")

# Save the plots as PDFs
pdf(file.path.snuc, width = DYN.SIZE$w, height = DYN.SIZE$h)
print(p.snuc)
dev.off()

pdf(file.path.bulk, width = DYN.SIZE$w, height = DYN.SIZE$h)
print(p.bulk)
dev.off()

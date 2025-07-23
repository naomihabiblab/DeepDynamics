#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#                     Supplementary Figures R Script                      ####
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# Load utilities:
source("src/utils/utils.R")
source("src/utils/utils.dynamics.R")
source("src/visualization/fig.dynamics.funs.R")

# Define working variables:
PATH.SUPP <- file.path(PATH.SUPP)
# -------------------------------------------------------------------------- #
##                       Pathologies & Cells in ABA                       ####
# -------------------------------------------------------------------------- #
names <- c('sqrt.amyloid_mf', 'sqrt.tangles_mf', 'cogng_demog_slope')

p.path.e4 <- plot.2g.split(names,
                           trajs = c("ABA"),
                           feature.name.flag = F,
                           split = "E4")
p.path.sex <- plot.2g.split(names,
                           trajs = c("ABA"),
                           feature.name.flag = F,
                           split = "Sex")

file.path.e4 <- file.path(PATH.SUPP, "p.supp.3a.new.aba.pdf")
file.path.sex <- file.path(PATH.SUPP, "p.supp.3b.new.aba.pdf")

pdf(file.path.e4, 
    width = SINGLE.DYN.SIZE$w * 3.6, 
    height = SINGLE.DYN.SIZE$h)
print(p.path.e4)
dev.off()

pdf(file.path.sex, 
    width = SINGLE.DYN.SIZE$w * 3.6, 
    height = SINGLE.DYN.SIZE$h)
print(p.path.sex)
dev.off()

names <- c("Oli.11", "Oli.4", "Mic.1",
           "Exc.1", "Inh.6", "Ast.4",
           "Ast.2", "Exc.2", "Mic.7")

file.cells.e4 <- file.path(PATH.F4, "p.supp.5c.rm42.new.pdf")

pdf(file.cells.e4, 
    width = SINGLE.DYN.SIZE$w * length(names) * 1.2, 
    height = SINGLE.DYN.SIZE$h * 4)
print(p.cells.e4)
dev.off()

# -------------------------------------------------------------------------- #
##              Pathologies & Cells in ABA, 4-Groups Split                ####
# -------------------------------------------------------------------------- #
names <- c('sqrt.amyloid_mf', 'sqrt.tangles_mf', 'cogng_demog_slope')

# Use Map to generate plots each for each name in names, and then save pdf files named "p.C.<name>.pdf"
Map(function(name) {
    p.path <- plot.4g.split(
        name,
        trajs = c("ABA"),
        feature.name.flag = F
    )
    
    file.path <- file.path(PATH.SUPP, paste0("p.supp.new.3C.aba.", name, ".pdf"))
    
    pdf(file.path, 
        width = DYN.SIZE$w, 
        height = DYN.SIZE$h * 3)
    print(p.path)
    dev.off()
    message("Saved plot to ", file.path)
}, names)


# -------------------------------------------------------------------------- #
##                              Sex Dynamics                              ####
# -------------------------------------------------------------------------- #
#' This is the split by sex dynamics for the cellular states presented in
#'  Figure 4B

names <- c('Mic.12', 'Mic.13', 'Ast.10', 'Oli.7', 'Inh.6')

p.cells.sex <- plot.2g.split(names,
                            trajs = c("prAD"),
                            feature.name.flag = F,
                            split = "Sex")

file.cells.supp <- file.path(PATH.SUPP, "p.supp.4b.sex.new.pdf")

pdf(file.cells.supp, 
    width = SINGLE.DYN.SIZE$w * length(names) * 1.2, 
    height = SINGLE.DYN.SIZE$h)
print(p.cells.sex)
dev.off()

# -------------------------------------------------------------------------- #
##                           Same States in ABA                           ####
# -------------------------------------------------------------------------- #
names <- c('Mic.12', 'Mic.13', 'Ast.10', 'Oli.7', 'Inh.6')

p.cells.ABA <- plot.2g.split(names,
                             trajs = c("ABA"),
                             feature.name.flag = F)

file.cells.ABA <- file.path(PATH.SUPP, "p.supp.4b.new.aba.pdf")

pdf(file.cells.ABA, 
    width = SINGLE.DYN.SIZE$w * length(names) * 1.2, 
    height = SINGLE.DYN.SIZE$h)
print(p.cells.ABA)
dev.off()


p.path.e4 <- plot.2g.split(names,
                           trajs = c("ABA"),
                           feature.name.flag = F,
                           split = "E4")
p.path.sex <- plot.2g.split(names,
                            trajs = c("ABA"),
                            feature.name.flag = F,
                            split = "Sex")

# -------------------------------------------------------------------------- #
##                          Same States in ABA 4g                         ####
# -------------------------------------------------------------------------- #
names <- c('Mic.12', 'Mic.13', 'Ast.10', 'Oli.7', 'Inh.6')

Map(function(name) {
    p.path <- plot.4g.split(
        name,
        trajs = c("ABA"),
        feature.name.flag = F
    )
    
    file.path <- file.path(PATH.SUPP, paste0("p.supp.4D.new.aba.", name, ".pdf"))
    
    pdf(file.path, 
        width = DYN.SIZE$w, 
        height = DYN.SIZE$h * 3)
    print(p.path)
    dev.off()
    message("Saved plot to ", file.path)
}, names)

#TODO Move to other better place 
name <- 'Inh.16'

p.state <- plot.4g.split(
    name,
    trajs = c("prAD"),
    feature.name.flag = F
)
file.path <- file.path(PATH.SUPP, paste0("p.supp.4D.new.prAD.", name, ".pdf"))

pdf(file.path, 
    width = DYN.SIZE$w, 
    height = DYN.SIZE$h * 3)
print(p.state)
dev.off()
message("Saved plot to ", file.path)

# -------------------------------------------------------------------------- #
##        sNuc vs. Bulk without split in prAD and ABA, pathologies        ####
# -------------------------------------------------------------------------- #
names <- c('sqrt.amyloid_mf', 'sqrt.tangles_mf', 'cogng_demog_slope')
colors.snuc <- c("dodgerblue4")
colors.bulk <- c("darkslategrey")

all.p.snuc.vs.bulk <- Map(function(state, color.s, color.b) {
    p.bulk <- plot.nosplit.bulk(state, color.b, pts.flag = F) +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    p.snuc <- plot.nosplit.snuc(state, color.s, pts.flag = F) +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    
    p.list <- plot_grid(p.snuc, p.bulk, nrow = 2) +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    return(p.list)
}, names, colors.snuc, colors.bulk)

# Define the file paths for saving
file.path.snuc <- file.path(PATH.SUPP, "for.supp.path.snuc.vs.bulk.nopts.pdf")
# Save the plots as PDFs
pdf(file.path.snuc, width = DYN.SIZE$w, height = DYN.SIZE$h)
print(all.p)
dev.off()

pdf(file.path.snuc, h=8, w=5.5)
plot_grid(plotlist=all.p.snuc.vs.bulk, 
          labels = names(all.p.snuc.vs.bulk),
          label_size = 8,
          label_fontface = "plain",
          align = 'hv',  # Align both horizontally and vertically
          axis = 'tblr',  # Align all sides
          hjust = -1)
dev.off()

# -------------------------------------------------------------------------- #
##              sNuc vs. Bulk without split in prAD and ABA               ####
# -------------------------------------------------------------------------- #
names <- SIG.CLUSTERS
colors.snuc <- c("dodgerblue4")
colors.bulk <- c("darkslategrey")

all.p.snuc.vs.bulk <- Map(function(state, color.s, color.b) {
    p.bulk <- plot.nosplit.bulk(state, color.b, pts.flag = T) +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    p.snuc <- plot.nosplit.snuc(state, color.s, pts.flag = T) +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    
    p.list <- plot_grid(p.snuc, p.bulk, nrow = 2) +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    return(p.list)
}, names, colors.snuc, colors.bulk)

# Define the file paths for saving
file.path <- file.path(PATH.SUPP, "for.supp.all.snuc.vs.bulk.pdf")

# Save the plots as PDFs
pdf(file.path.snuc, width = DYN.SIZE$w, height = DYN.SIZE$h)
print(p.snuc)
dev.off()

pdf(file.path, h=34.69291, w=24.08661)
plot_grid(plotlist=all.p.snuc.vs.bulk, 
          labels = names(all.p.snuc.vs.bulk),
          label_size = 8,
          label_fontface = "plain",
          align = 'hv',  # Align both horizontally and vertically
          axis = 'tblr',  # Align all sides
          hjust = -1)
dev.off()


#TODO remove and add the new sup
# TA:
names <- c("Oli.11", "Oli.4", "Mic.1",
           "Exc.1", "Inh.6", "Ast.4", 
           "Ast.2", "Exc.2", "Mic.7")

p.cells.e4 <- plot.2g.split(names,
                            trajs = c("prAD"),
                            feature.name.flag = F,
                            split = "E4",
                            nrow = 3)

file.cells.e4 <- file.path("results", "p.e4.TA.pdf")
pdf(file.cells.e4, 
    width = SINGLE.DYN.SIZE$w * length(names) * 1.2 / 3, 
    height = SINGLE.DYN.SIZE$h * 3)
print(p.cells.e4)
dev.off()


Map(function(traj, name) {
    print(traj); print(name)
    p.snuc <- plot.nosplit.snuc(name, colors.snuc, pts.flag = T, trajs = c(traj))
    p.bulk <- plot.nosplit.bulk(name, colors.bulk, pts.flag = T, trajs = c(traj))
    
    p.both <- plot_grid(p.snuc, p.bulk, nrow = 1) +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    
    file.path <- file.path(PATH.SUPP, paste("p.supp.signals.", name, "pdf", sep = "."))
    p<-p.snuc / p.bulk
    pdf(file.path, width = 3*SINGLE.DYN.SIZE$w, height = 2*SINGLE.DYN.SIZE$h)
    print(p)
    dev.off()
}, names.list, names(names.list))




# -------------------------------------------------------------------------- #
## Supplemantry Spider Graphs ####
# -------------------------------------------------------------------------- #
source("src/utils/utils.R")
source("src/utils/utils.DEGs.R")
source("src/visualization/fig.DEGs.funs.R")

# prAD spiders:
cell.types = c('Astrocytes', 'Microglia', 'Oligodendrocytes', 'SST Neurons', 'PV', 'Exc.1')

Map(function(file, ind, top) {
    # Create spider plot:
    file.name <- file.path(PATH.SUPP, paste0("p.supp.5D.rm24.", cell.types[[ind]], 
                                           ".prad.top.", top, ".pdf"))
    p.spider <- create.spider.graph(fname.rd = file,
                                    top.degs.num = top,
                                    add.title = paste('prAD Probability > 0.5,', 
                                                      cell.types[[ind]]))
    # Maybe sizes will be changed later, now the good size is (7, 7.5)
    pdf(file.name, 
        width = SPIDER.SIZE$w, 
        height = SPIDER.SIZE$h)
    print(p.spider)
    dev.off()
    

    message('Successfully saved ', cell.types[[ind]], ' spider plot as', file.name)
}, rep(c(AST.PRAD.FILE, MIC.PRAD.FILE, OLI.PRAD.FILE, SST.PRAD.FILE, PV.PRAD.FILE, EXC.PRAD.FILE), each=4), rep(1:6, each=4), c(5,7,10,12))

# ABA spiders:
cell.types = c('Astrocytes', 'Microglia', 'Oligodendrocytes', 'SST Neurons', 'PV', 'Exc.1')

Map(function(file, ind, top) {
    # Create spider plot:
    file.name <- file.path(PATH.SUPP, paste0("p.supp.5D.rm24.", cell.types[[ind]], 
                                           ".aba.top.", top, ".pdf"))
    p.spider <- create.spider.graph(fname.rd = file,
                                    top.degs.num = top,
                                    add.title = paste('ABA Probability > 0.5,', 
                                                      cell.types[[ind]]))
    # Maybe sizes will be changed later, now the good size is (7, 7.5)
    pdf(file.name, 
        width = SPIDER.SIZE$w, 
        height = SPIDER.SIZE$h)
    print(p.spider)
    dev.off()

    message('Successfully saved ', cell.types[[ind]], ' spider plot as', file.name)
}, rep(c(AST.ABA.FILE, MIC.ABA.FILE, OLI.ABA.FILE, SST.ABA.FILE, PV.ABA.FILE, EXC.ABA.FILE), each=4), rep(1:6, each=4), c(5,7,10,12))


# All people spiders:
cell.types = c('Astrocytes', 'Microglia', 'Oligodendrocytes', 'SST Neurons', 'PV', 'Exc.1')

Map(function(file, ind, top) {
    # Create spider plot:
    file.name <- file.path(PATH.SUPP, paste0("p.supp.5D.rm12.", cell.types[[ind]], 
                                           ".all.top.", top, ".pdf"))
    p.spider <- create.spider.graph(fname.rd = file,
                                    top.degs.num = top,
                                    add.title = paste('All People,', 
                                                      cell.types[[ind]]))
    # Maybe sizes will be changed later, now the good size is (7, 7.5)
    pdf(file.name, 
        width = SPIDER.SIZE$w, 
        height = SPIDER.SIZE$h)
    print(p.spider)
    dev.off()

    message('Successfully saved ', cell.types[[ind]], ' spider plot as', file.name)
}, rep(c(AST.ALL.FILE, MIC.ALL.FILE, OLI.ALL.FILE, SST.ALL.FILE, PV.ALL.FILE, EXC.ALL.FILE), each=4), rep(1:6, each=4), c(5,7,10,12))


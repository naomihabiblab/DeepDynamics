#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#                                  Figure 5                               ####
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# Load utilities:
source("src/utils/utils.R")
source("src/utils/utils.DEGs.R")
source("src/visualization/fig.DEGs.funs.R")

PATH.F5 <- file.path(PATH.FIGS, "F5")

# -------------------------------------------------------------------------- #
##                       Panel A: Intro to Spider Graphs                  ####
# -------------------------------------------------------------------------- #
# Prepared in Adobe Illustrator
# -------------------------------------------------------------------------- #
##                          Panels B-F: Spider Graphs                     ####
# -------------------------------------------------------------------------- #
cell.types = c('Astrocytes', 'Microglia', 'Oligodendrocytes', 'SST Neurons', 'PV', 'Exc.1')

Map(function(file, ind, top) {
    # Create spider plot:
    file.name <- file.path(PATH.F5, paste0("p.spider.", cell.types[[ind]], 
                                           ".top.", top, ".pdf"))
    print(file)
    print(cell.types[[ind]])
    p.spider <- create.spider.graph(fname.rd = file,
                                    top.degs.num = top,
                                    add.title = paste('Healthy Pseudotime,', 
                                                     cell.types[[ind]]))
    # Maybe sizes will be changed later, now the good size is (7, 7.5)
    pdf(file.name,
        width = SPIDER.SIZE$w,
        height = SPIDER.SIZE$h)
    print(p.spider)
    dev.off()

    message('Successfully saved', cell.types[[ind]], 'spider plot as', file.name)
}, rep(c(AST.FILE.HEALTHY, MIC.FILE.HEALTHY, OLI.FILE.HEALTHY, SST.FILE.HEALTHY, PV.FILE.HEALTHY, EXC.FILE.HEALTHY), each=4), rep(1:6, each=4), c(5,7,10,12))

# -------------------------------------------------------------------------- #
##                   Panel G: Heatmap of Heat Shock Genes                 ####
# -------------------------------------------------------------------------- #
# Create heatmap:
file.name <- file.path(PATH.F5, "p.G.heatmap.pdf")
p.heatmap <- create.genes.heatmap(dir.path = HSP.DIR,
                                  split.by = 'celltype',
                                  filter.genes = NULL)

pdf(file.name, 
    width = 10, 
    height = 9)
draw(p.heatmap)
dev.off()


message('Successfully saved heatmap of Heat Shock Genes as ', file.name)

# -------------------------------------------------------------------------- #
##        Panel H: Comparison heatmap in microglia- SEA-AD, ROSMAP        ####
# -------------------------------------------------------------------------- #
# Mic comparison
dir.path <- "../degs.data/heatmap.mic.ros-sea"
hm <- create.simple.genes.heatmap(dir.path,
                                  gene.pattern = "^(HSP|DNAJ)",
                                  file.pattern = "",
                                  title = "Mic vs ROSMAP vs SEA-AD",
                                  color.sat =TRUE,
                                  cluster.columns = TRUE,
                                  cluster.rows = TRUE)

hma <- create.simple.genes.heatmap.annotated(dir.path,
                                             gene.pattern = "^(HSP|DNAJ)",
                                             file.pattern = "",
                                             title = "Mic vs ROSMAP vs SEA-AD",
                                             color.sat =TRUE,
                                             cluster.columns = TRUE,
                                             cluster.rows = FALSE
                                             
)

pdf(file.path(PATH.F5, "p.H.heatmap.mic.ros-sea.simple.annot.pdf"),
    width = 5, 
    height = 5
)
draw(hma)
dev.off()

# Panel I was prepared via STRING DB (ref in paper).

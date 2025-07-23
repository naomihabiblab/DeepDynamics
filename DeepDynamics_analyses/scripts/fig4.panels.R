#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#                                  Figure 4                               ####
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# Load utilities:
source("src/utils/utils.R")
source("src/utils/utils.TA.R")
source("src/utils/utils.prot.R")
source("src/utils/utils.dynamics.R")
source("src/visualization/fig.TA.funs.R")
source("src/visualization/fig.prot.funs.R")
source("src/visualization/fig.dynamics.funs.R")

PATH.F4 <- file.path(PATH.FIGS, "F4")

# -------------------------------------------------------------------------- #
##                              Panel A                                   ####
# -------------------------------------------------------------------------- #
# -------------------------------------------------------------------------- #
###                   Panel A: TA with Pathologies                        ####
# -------------------------------------------------------------------------- #
# Load bulk df
load.celmod.bulk()

bulk.df <- bulk.df.min %>% mutate(
    Female = ifelse(msex, 0, 1)
)

p.all.sex.apoe.pathology.sig <- states.pathologies(bulk.df,
                                                   trait.names = c(DEFAULT.PATHOLOGIES, 'apoe_4', 'Female'),
                                                   add2title = '\nAll Individuals, Sex and APOE4 as Pathology',
                                                   sigonly = T,
                                                   cluster.cols.dist = 'pearson')

file.path.TA <- file.path(PATH.F4, "p.A.1107.TA.pdf")

pdf(file.path.TA, 
    width = TA.SIZE$w, 
    height = TA.SIZE$h)
grid.arrange(p.all.sex.apoe.pathology.sig)
dev.off()


# -------------------------------------------------------------------------- #
###                Panel A: States with SHAP direction                    ####
# -------------------------------------------------------------------------- #
# Create DF:
SHAP.prAD.ls <- list(
    pos = c('Ast.10', 'Ast.3', 'Ast.7',
            'Oli.7', 'OPC.2', 'Mic.12', 
            'Peri.2', 'OPC.1', 'Ast.6',
            'End.2', 'Mic.8', 'Mic.3',
            'Exc.4', 'Mic.1', 'Mic.5',
            'Mic.14', 'Exc.1'),
    neg = c('Ast.5', 'Mic.7', 'OPC.3')
)

SHAP.ABA.ls <- list(
    pos = c('OPC.3', 'Ast.5', 'Mic.7',
            'Exc.3', 'Oli.3', 'Ast.2',
            'Ast.1', 'Mic.10', 'Ast.4',
            'Mic.6', 'Oli.4', 'Peri.1',
            'End.2', 'Inh.6', 'End.1',
            'Ast.9', 'Mic.2', 'Mic.4'),
    neg = c('Ast.10', 'Ast.7')
)

cells.ls <- c('End.3', 'Inh12', 'Inh.5', 'Exc.12', 
              'Inh.7', 'Mic.1','Oli.5', 'Exc.8',
              'Mic.12', 'Mic.13', 'Inh.16', 'Ast.10',
              'Oli.7', 'End.1', 'Ast.6', 'OPC.1',
              'Oli.11', 'Ast.2', 'OPC.2', 'Mic.14',
              'Ast.1', 'Inh.15', 'Oli.3', 'Ast.7',
              'Oli.9', 'Inh.6', 'Ast.4', 'Oli.4',
              'Exc.1',  'Exc.3')



hm <- create.SHAP.states.HM(cells.ls, SHAP.prAD.ls, SHAP.ABA.ls)


file.path.hm <- file.path(PATH.F4, "p.A.shap.new.hm.pdf")

pdf(file.path.hm, 
    width = 2.4, 
    height = 10)
print(hm)
dev.off()

# -------------------------------------------------------------------------- #
##                     Panel B: Proteomics Density Plots                  ####
# -------------------------------------------------------------------------- #
names <- c(AST.MARKERS, MIC.MARKERS, NEUR.MARKERS)

Map(function(name) {
    var.lst <- prepare.prot.for.density(proteomics, name)
    
    p.density <- create.density.plot(
        df = var.lst$df,
        gene = name,
        uniprot = var.lst$uni.name,
        p.val = var.lst$p.val
    )
    
    file.path <- file.path(PATH.F4, paste("p.B", name, "pdf", sep = "."))
    
    pdf(file.path, 
        width = PROT.SINGLE.SIZE$w, 
        height = PROT.SINGLE.SIZE$h)
    print(p.density)
    dev.off()
    message("Saved plot to ", file.path)
}, names)



# -------------------------------------------------------------------------- #
##               Panel C: Dynamics of Cellular Sub-populations            ####
# -------------------------------------------------------------------------- #
names <- c('Mic.13', 'Ast.10', 'Oli.7', 'Inh.16')

p.cells.e4 <- plot.2g.split(names,
                           trajs = c("prAD"),
                           feature.name.flag = F,
                           split = "E4")
p.cells.e4 <- plot.2g.split(names,
                            trajs = c("prAD"),
                            feature.name.flag = F,
                            split = "E4",
                            nrow = 3
                            )

file.cells.e4 <- file.path(PATH.F4, "p.C.e4.pdf")

pdf(file.cells.e4, 
    width = SINGLE.DYN.SIZE$w * length(names) * 1.2, 
    height = SINGLE.DYN.SIZE$h)
print(p.cells.e4)
dev.off()

# -------------------------------------------------------------------------- #
##                 Panel D: 4-Group Split of Selected States              ####
# -------------------------------------------------------------------------- #
names <- c('Ast.10','Mic.13', 'Mic.12', 'Oli.7')

# Vertical orientation
Map(function(name) {
    or <- "col"
    p.cells <- plot.4g.split(
        name,
        trajs = c("prAD"),
        feature.name.flag = F,
        orientation = or
    )
    
    file.path <- file.path(PATH.F4, paste("p.D", name, or, "pdf", sep = "."))
    
    pdf(file.path, 
        width = SINGLE.DYN.SIZE$w, 
        height = DYN.SIZE$h * 2.4)
    print(p.cells)
    dev.off()
    message("Saved plot to ", file.path)
}, names)

# Horizontal orientation
Map(function(name) {
    or <- "row"
    p.cells <- plot.4g.split(
        name,
        trajs = c("prAD"),
        feature.name.flag = F,
        orientation = or
    )
    
    file.path <- file.path(PATH.F4, paste("p.D", name, or, "pdf", sep = "."))
    
    pdf(file.path, 
        width = DYN.SIZE$w, 
        height = DYN.SIZE$h)
    print(p.cells)
    dev.off()
    message("Saved plot to ", file.path)
}, names)

# -------------------------------------------------------------------------- #
##                       Panel E: Dynamics of Pathways                    ####
# -------------------------------------------------------------------------- #
data.slot.ls <- c('psbulk.mic.lipid.catabolic', 'psbulk.ast.met.ion.trans.synaptic',
                  'psbulk.ast.react.oxy.metabolic', 'psbulk.oli.morphogenesis') # TODO check data slots

Map(function(data.slot) {
    p.pathway <- plot.2g.split(
        features = c('mean_scaled'),
        trajs = c("prAD"),
        feature.name.flag = F,
        split = "E4",
        fit.section.e4 = data.slot
    )
    
    file.path <- file.path(PATH.F4, paste("p.C", data.slot, "pdf", sep = "."))
    
    pdf(file.path, 
        width = SINGLE.DYN.SIZE$w, 
        height = SINGLE.DYN.SIZE$h)
    print(p.pathway)
    dev.off()
    message("Saved plot to ", file.path)
}, data.slot.ls)


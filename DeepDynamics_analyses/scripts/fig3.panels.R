#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#                                  Figure 3                               ####
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# Load utilities:
source("src/utils/utils.R")
source("src/utils/utils.dynamics.R")
source("src/utils/utils.TA.R")
source("src/visualization/fig.TA.funs.R")
source("src/visualization/fig.dynamics.funs.R")

# Define working variables:
PATH.F3 <- file.path(PATH.FIGS, "F3")

# -------------------------------------------------------------------------- #
##                  Panel A-B: Pathologies APOE4 linear model               ####
# -------------------------------------------------------------------------- #
p.bulk <- pathology.sex.apoe(sex.title = "Sex and Pathologies Association",
                             apoe.title = "APOE4 and Pathologies Association", 
                             control.names = c("age_death")
)

# TODO check bug
# -------------------------------------------------------------------------- #
##                        Panel C: Pathologies APOE4 dynamics             ####
# -------------------------------------------------------------------------- #
names <- c('sqrt.amyloid_mf', 'sqrt.tangles_mf', 'cogng_demog_slope')

p.path.e4 <- plot.2g.split(names,
                            trajs = c("prAD"),
                            feature.name.flag = F,
                            split = "E4")

file.path.e4 <- file.path(PATH.F3, "p.C.e4.pdf")

pdf(file.path.e4, 
    width = SINGLE.DYN.SIZE$w * 3.6, 
    height = SINGLE.DYN.SIZE$h)
print(p.path.e4)
dev.off()

# -------------------------------------------------------------------------- #
##                        Panel D: Pathologies sex dynamics               ####
# -------------------------------------------------------------------------- #
names <- c('sqrt.amyloid_mf', 'sqrt.tangles_mf', 'cogng_demog_slope')

p.path.sex <- plot.2g.split(names,
                            trajs = c("prAD"),
                            feature.name.flag = F,
                            split = "Sex")

file.path.sex <- file.path(PATH.F3, "p.D.sex.pdf")

pdf(file.path.sex, 
    width = SINGLE.DYN.SIZE$w * 3.6, 
    height = SINGLE.DYN.SIZE$h)
print(p.path.sex)
dev.off()

# -------------------------------------------------------------------------- #
##                       Panel E: Pathologies ApoE4 & Sex                 ####
# -------------------------------------------------------------------------- #
names <- c('sqrt.amyloid_mf', 'sqrt.tangles_mf', 'cogng_demog_slope')

# Use Map to generate plots each for each name in names, and then save pdf files named "p.C.<name>.pdf"
Map(function(name) {
    p.path <- plot.4g.split(
        name,
        trajs = c("prAD"),
        feature.name.flag = F
    )
    
    file.path <- file.path(PATH.F3, paste0("p.E.", name, ".pdf"))
    
    pdf(file.path, 
        width = DYN.SIZE$w, 
        height = DYN.SIZE$h * 3)
    print(p.path)
    dev.off()
    message("Saved plot to ", file.path)
}, names)

# -------------------------------------------------------------------------- #
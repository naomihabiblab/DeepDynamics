#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#             Differentially Expressed Genes (DEGs) Functions             ####
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# -------------------------------------------------------------------------- #
##                                Constants                               ####
# -------------------------------------------------------------------------- #


# File paths for DEGs results (4-groups or 2-groups)
SST.FILE <- "../data/async/traj.thr/inhibitoryInh.6.traj.thrs.All people.pscount.1e-6.poisson.logfc.0.25.pct.0.1.4g.apoe4sex.2g.RData"  # 2 groups
AST.FILE <- "../degs.data/four.g/astrocytes.no24.bulk.prob.All.\'healthy\'.-.before.0.1.pstime.poisson.logfc.0.25.pct.0.1.4g.apoe4sex.traj.meanfxn.RData"
MIC.FILE <- "../degs.data/four.g/microglia.no24.bulk.prob.All.\'healthy\'.-.before.0.1.pstime.poisson.logfc.0.25.pct.0.1.4g.apoe4sex.traj.meanfxn.RData"
OLI.FILE <- "../degs.data/four.g/oligodendrocytes.no24.bulk.prob.All.\'healthy\'.-.before.0.1.pstime.poisson.logfc.0.25.pct.0.1.4g.apoe4sex.traj.meanfxn.RData"

SST.FILE.HEALTHY = "../degs.data/after12rm.four.g/rm24/inhibitory.inh6.rm24.bulk.prob.All.\'healthy\'.-.before.0.1.pstime.poisson.logfc.0.25.pct.0.1.4g.apoe4sex.traj.meanfxn.RData"
PV.FILE.HEALTHY = "../degs.data/after12rm.four.g/rm24/inhibitory.inh16.rm24.bulk.prob.All.\'healthy\'.-.before.0.1.pstime.poisson.logfc.0.25.pct.0.1.4g.apoe4sex.traj.meanfxn.RData"
EXC.FILE.HEALTHY = "../degs.data/after12rm.four.g/rm24/cux2+.exc1.rm24.bulk.prob.All.\'healthy\'.-.before.0.1.pstime.poisson.logfc.0.25.pct.0.1.4g.apoe4sex.traj.meanfxn.RData"
AST.FILE.HEALTHY = "../degs.data/after12rm.four.g/rm24/astrocytes.rm24.bulk.prob.All.\'healthy\'.-.before.0.1.pstime.poisson.logfc.0.25.pct.0.1.4g.apoe4sex.traj.meanfxn.RData"
MIC.FILE.HEALTHY = "../degs.data/after12rm.four.g/rm24/microglia.rm24.bulk.prob.All.\'healthy\'.-.before.0.1.pstime.poisson.logfc.0.25.pct.0.1.4g.apoe4sex.traj.meanfxn.RData"
OLI.FILE.HEALTHY = "../degs.data/after12rm.four.g/rm24/oligodendrocytes.rm24.bulk.prob.All.\'healthy\'.-.before.0.1.pstime.poisson.logfc.0.25.pct.0.1.4g.apoe4sex.traj.meanfxn.RData"

# Directory path of files for DEGs heatmap:
HSP.DIR <- "../degs.data/heatmap.4g"


# For SUPP figures, prAD & ABA spider files:
AST.PRAD.FILE <- "../degs.data/after12rm.four.g/rm24/astrocytes.rm24.bulk.prob.All.after.0.1.pstime.and.prad.higher.0.5.poisson.logfc.0.25.pct.0.1.4g.apoe4sex.traj.meanfxn.RData"
MIC.PRAD.FILE <- "../degs.data/after12rm.four.g/rm24/microglia.rm24.bulk.prob.All.after.0.1.pstime.and.prad.higher.0.5.poisson.logfc.0.25.pct.0.1.4g.apoe4sex.traj.meanfxn.RData"
OLI.PRAD.FILE <- "../degs.data/after12rm.four.g/rm24/oligodendrocytes.rm24.bulk.prob.All.after.0.1.pstime.and.prad.higher.0.5.poisson.logfc.0.25.pct.0.1.4g.apoe4sex.traj.meanfxn.RData"
SST.PRAD.FILE <- "../degs.data/after12rm.four.g/rm24/inhibitory.inh6.rm24.bulk.prob.All.after.0.1.pstime.and.prad.higher.0.5.poisson.logfc.0.25.pct.0.1.4g.apoe4sex.traj.meanfxn.RData"
PV.PRAD.FILE <- "../degs.data/after12rm.four.g/rm24/inhibitory.inh16.rm24.bulk.prob.All.after.0.1.pstime.and.prad.higher.0.5.poisson.logfc.0.25.pct.0.1.4g.apoe4sex.traj.meanfxn.RData"
EXC.PRAD.FILE <- "../degs.data/after12rm.four.g/rm24/cux2+.exc1.rm24.bulk.prob.All.after.0.1.pstime.and.prad.higher.0.5.poisson.logfc.0.25.pct.0.1.4g.apoe4sex.traj.meanfxn.RData"

AST.ABA.FILE <- "../degs.data/after12rm.four.g/rm24/astrocytes.rm24.bulk.prob.All.after.0.1.pstime.and.aba.higher.0.5.poisson.logfc.0.25.pct.0.1.4g.apoe4sex.traj.meanfxn.RData"
MIC.ABA.FILE <- "../degs.data/after12rm.four.g/rm24/microglia.rm24.bulk.prob.All.after.0.1.pstime.and.aba.higher.0.5.poisson.logfc.0.25.pct.0.1.4g.apoe4sex.traj.meanfxn.RData"
OLI.ABA.FILE <- "../degs.data/after12rm.four.g/rm24/oligodendrocytes.rm24.bulk.prob.All.after.0.1.pstime.and.aba.higher.0.5.poisson.logfc.0.25.pct.0.1.4g.apoe4sex.traj.meanfxn.RData"
SST.ABA.FILE <- "../degs.data/after12rm.four.g/rm24/inhibitory.inh6.rm24.bulk.prob.All.after.0.1.pstime.and.aba.higher.0.5.poisson.logfc.0.25.pct.0.1.4g.apoe4sex.traj.meanfxn.RData"
PV.ABA.FILE <- "../degs.data/after12rm.four.g/rm24/inhibitory.inh16.rm24.bulk.prob.All.after.0.1.pstime.and.aba.higher.0.5.poisson.logfc.0.25.pct.0.1.4g.apoe4sex.traj.meanfxn.RData"
EXC.ABA.FILE <- "../degs.data/after12rm.four.g/rm24/cux2+.exc1.rm24.bulk.prob.All.after.0.1.pstime.and.aba.higher.0.5.poisson.logfc.0.25.pct.0.1.4g.apoe4sex.traj.meanfxn.RData"

# -------------------------------------------------------------------------- #
## DEGs Analysis Functions ####
# -------------------------------------------------------------------------- #
#' This function adds apoe4 to metadata if there is any 4 in the apoe genotype
#'@param obj - Seurat object with all the data
#'
#'@return obj - same object with the addition of the ApoE4.
add_apoe4_to_metadata <- function(obj) {
    tmp.data <- obj@meta.data
    tmp.data <- tmp.data %>%
        mutate(apoe_4 = case_when(
            apoe_genotype %in% c("34", "44") ~ 1,
            apoe_genotype %in% c("22", "23", "33") ~ 0,
            TRUE ~ NA_real_
        ))
    obj <- Seurat::AddMetaData(object = obj,
                               metadata = tmp.data$apoe_4,
                               col.name = "apoe_4")
    rm(tmp.data)
    # obj.apoe.list <- SplitObject(obj, split.by = 'apoe_4')
    return(obj)
}


add.traj.to.metadata <- function(obj, traj.path = PATH.snuc.df) {
    # load the traj.path
    load(traj.path)
    
    # take projid from obj:
    tmp.data <- obj@meta.data
    
    tryCatch({
        # merge into it the prAD, pstime, ABA
        tmp.data <- merge(
            tmp.data,
            snuc.traj.df %>%
                mutate(projid = ID),
            by = 'projid'
        )
    } , error = function(e) {
        print('it loaded bulk.df')
        tmp.data <<- merge(
            tmp.data,
            bulk.df %>%
                dplyr::select(ID, prAD, pseudotime, ABA) %>%
                mutate(projid = ID),
            by = 'projid'
        )
        
    })
    
    
    # Add to metadata
    obj <- Seurat::AddMetaData(object = obj,
                               metadata = tmp.data$pseudotime,
                               col.name = "pstime")
    obj <- Seurat::AddMetaData(object = obj,
                               metadata = tmp.data$prAD,
                               col.name = "prAD")
    obj <- Seurat::AddMetaData(object = obj,
                               metadata = tmp.data$ABA,
                               col.name = "ABA")
    # return object
    return(obj)
}

# This funcition is very similar to find.markers.2.groups but different with the trajectory thresholds for the object.
find.markers.traj.thrs <- function(obj, test = 'poisson', meanfxn.flag = F, 
                                   fc.thr = 0.25, pct.thr = 0.1, prad.min = -1,
                                   prad.max = 2, pstime.min = -1, pstime.max = 2,
                                   aba.min = -1, aba.max = 2, add.title = '',
                                   latent.vars = NULL, prob.choose = 'snuc') {
    if (meanfxn.flag) {
        print('with different mean.fxn')
        
        mean.fxn <- function(x) {
            return(log(x = rowMeans(x = x) + 1e-6, base = 2))
        }
        
        message('Subset object according to trajectory thresholds')
        obj <- add.traj.to.metadata(obj, traj.path = ifelse(prob.choose == 'snuc', PATH.snuc.df, PATH.bulk.df))
        obj <- subset(x = obj, subset = 
                          (pstime >= pstime.min) & 
                          (pstime <= pstime.max) & 
                          (prAD >= prad.min) & 
                          (prAD <= prad.max) &
                          (ABA >= aba.min) &
                          (ABA <= aba.max))
        
        message(paste('In', test, 'check with mean.fxn'))
        if (is.null(latent.vars) == F) message(paste('With latent variables:', latent.vars))
        
        message("APOE4")
        Idents(obj) <- 'apoe_4'
        
        # Find markers in all interseting variables
        print('1 vs 0')
        e4.10 <- FindMarkers(obj, test.use = test, mean.fxn = mean.fxn,
                             logfc.threshold = fc.thr, min.pct = pct.thr,
                             ident.1 = "1", ident.2 = "0",
                             latent.vars = latent.vars)
        print('Done! :)')
        
        print('Sex')
        Idents(obj) <- 'sex'
        
        print('Females vs. Males')
        sex.fm <- FindMarkers(obj, test.use = test, mean.fxn = mean.fxn,
                              logfc.threshold = fc.thr, min.pct = pct.thr,
                              ident.1 = "Female", ident.2 = "Male",
                              latent.vars = latent.vars)
        print('Done! :)')
        
        if (is.null(latent.vars) == F) {
            latent.add <- paste0('.', paste(latent.vars, collapse = '.'), '.')
        } else {
            latent.add <- ''
        }
        
        # Design the name of the rds file
        rd.fname <- paste0(obj@project.name,
                           add.title,
                           '.pscount.1e-6.', 
                           test,
                           '.logfc.',
                           fc.thr,
                           '.pct.',
                           pct.thr,
                           latent.add,
                           '.4g.apoe4sex.2g.RData')
        
    } else {
        # Find markers in all intersting variables
        print('Currently implementred only for with mean.fxn')
    }
    
    # Save file containing the interesting comparisons 
    save(e4.10, sex.fm, file = rd.fname)
    message('Succesfuly saved the file: ', add.title)
    
}

# This function runs differentially expressed genes on the sex.apoe ident 4 groups.
find.markers.4.groups <- function(obj, test = 'poisson', meanfxn.flag = F, 
                                  fc.thr = 0.25, pct.thr = 0.1) {
    # Change idents:
    Idents(obj) <- 'sex.apoe4'
    
    if (test == 'DEseq2') { DefaultAssay(obj) <- 'RNA'}
    
    if (meanfxn.flag) {
        print('with different mean.fxn')
        
        mean.fxn <- function(x) {
            return(log(x = rowMeans(x = x) + 1e-6, base = 2))
        }
        
        if (test == 'MAST') {
            print('In MAST no need for mean.fxn')
            # Find markers in all intersting variables
            print('Females.0 vs Females.1')
            f0f1 <- FindMarkers(obj, test.use = test,
                                logfc.threshold = fc.thr, min.pct = pct.thr,
                                ident.1 = "Female.0", ident.2 = "Female.1")
            print('Done! :)')
            
            print('Males.0 vs Males.1')
            m0m1 <- FindMarkers(obj, test.use = test,
                                logfc.threshold = fc.thr, min.pct = pct.thr,
                                ident.1 = "Male.0", ident.2 = "Male.1")
            print('Done! :)')
            
            print('Males.0 vs Females.0')
            m0f0 <- FindMarkers(obj, test.use = test,
                                logfc.threshold = fc.thr, min.pct = pct.thr,
                                ident.1 = "Male.0", ident.2 = "Female.0")
            print('Done! :)')
            
            print('Males.1 vs Females.1')
            m1f1 <- FindMarkers(obj, test.use = test,
                                logfc.threshold = fc.thr, min.pct = pct.thr,
                                ident.1 = "Male.1", ident.2 = "Female.1")
            print('Done! :)')
            
            # Design the name of the rds file
            rd.fname <- paste0(obj@project.name,
                               '.', 
                               test,
                               '.logfc.',
                               fc.thr, 
                               '.pct.',
                               pct.thr,
                               '.4g.apoe4sex.updatedthrs.RData')
        } else {
            print(paste('In', test, 'check with mean.fxn'))
            # Find markers in all interseting variables
            print('Females.0 vs Females.1')
            f0f1 <- FindMarkers(obj, test.use = test, mean.fxn = mean.fxn,
                                logfc.threshold = fc.thr, min.pct = pct.thr,
                                ident.1 = "Female.0", ident.2 = "Female.1")
            print('Done! :)')
            
            print('Males.0 vs Males.1')
            m0m1 <- FindMarkers(obj, test.use = test, mean.fxn = mean.fxn,
                                logfc.threshold = fc.thr, min.pct = pct.thr,
                                ident.1 = "Male.0", ident.2 = "Male.1")
            print('Done! :)')
            
            print('Males.0 vs Females.0')
            m0f0 <- FindMarkers(obj, test.use = test, mean.fxn = mean.fxn,
                                logfc.threshold = fc.thr, min.pct = pct.thr,
                                ident.1 = "Male.0", ident.2 = "Female.0")
            print('Done! :)')
            
            print('Males.1 vs Females.1')
            m1f1 <- FindMarkers(obj, test.use = test, mean.fxn = mean.fxn,
                                logfc.threshold = fc.thr, min.pct = pct.thr,
                                ident.1 = "Male.1", ident.2 = "Female.1")
            print('Done! :)')
            
            # Design the name of the rds file
            rd.fname <- paste0(obj@project.name,
                               '.pscount.1e-6.', 
                               test,
                               '.logfc.',
                               fc.thr,
                               '.pct.',
                               pct.thr,
                               '.4g.apoe4sex.thrs.meanfxn.RData')
        }
        
    } else {
        # Find markers in all intersting variables
        print('Females.0 vs Females.1')
        f0f1 <- FindMarkers(obj, test.use = test,
                            logfc.threshold = fc.thr, min.pct = pct.thr,
                            ident.1 = "Female.0", ident.2 = "Female.1")
        print('Done! :)')
        
        print('Males.0 vs Males.1')
        m0m1 <- FindMarkers(obj, test.use = test,
                            logfc.threshold = fc.thr, min.pct = pct.thr,
                            ident.1 = "Male.0", ident.2 = "Male.1")
        print('Done! :)')
        
        print('Males.0 vs Females.0')
        m0f0 <- FindMarkers(obj, test.use = test,
                            logfc.threshold = fc.thr, min.pct = pct.thr,
                            ident.1 = "Male.0", ident.2 = "Female.0")
        print('Done! :)')
        
        print('Males.1 vs Females.1')
        m1f1 <- FindMarkers(obj, test.use = test,
                            logfc.threshold = fc.thr, min.pct = pct.thr,
                            ident.1 = "Male.1", ident.2 = "Female.1")
        print('Done! :)')
        
        # Design the name of the rds file
        rd.fname <- paste0(obj@project.name,
                           '.', 
                           test,
                           '.logfc.',
                           fc.thr,
                           '.pct.',
                           pct.thr,
                           '.4g.apoe4sex.updatedthrs.RData')
    }
    
    # Save file containing the interesting comparisons 
    save(f0f1, m0m1, m0f0, m1f1, file = rd.fname)
} 



# This function runs differentially expressed genes on the sex.apoe ident 4 groups.
find.markers.traj.4.groups <- function(obj, test = 'poisson', fc.thr = 0.25, 
                                       pct.thr = 0.1, prad.min = -1, prad.max = 2, 
                                       pstime.min = -1, pstime.max = 2, aba.min = -1, 
                                       aba.max = 2, add.fname = '', latent.vars = NULL,
                                       prob.choose = 'snuc', meanfxn.flag = T) {
    
    obj <- add.traj.to.metadata(obj, traj.path = ifelse(prob.choose == 'snuc', PATH.snuc.df, PATH.bulk.df))
    obj <- subset(x = obj, subset = 
                      (pstime >= pstime.min) & 
                      (pstime <= pstime.max) & 
                      (prAD >= prad.min) & 
                      (prAD <= prad.max) &
                      (ABA >= aba.min) &
                      (ABA <= aba.max))
    # Change idents:
    Idents(obj) <- 'sex.apoe4'
    
    if (test == 'DEseq2') { DefaultAssay(obj) <- 'RNA'}
    
    if (meanfxn.flag) {
        print('with different mean.fxn')
        
        mean.fxn <- function(x) {
            return(log(x = rowMeans(x = x) + 1e-6, base = 2))
        }
        
        if (test == 'MAST') {
            print('In MAST no need for mean.fxn')
            # Find markers in all intersting variables
            print('Females.0 vs Females.1')
            f0f1 <- FindMarkers(obj, test.use = test,
                                logfc.threshold = fc.thr, min.pct = pct.thr,
                                ident.1 = "Female.0", ident.2 = "Female.1")
            print('Done! :)')
            
            print('Males.0 vs Males.1')
            m0m1 <- FindMarkers(obj, test.use = test,
                                logfc.threshold = fc.thr, min.pct = pct.thr,
                                ident.1 = "Male.0", ident.2 = "Male.1")
            print('Done! :)')
            
            print('Males.0 vs Females.0')
            m0f0 <- FindMarkers(obj, test.use = test,
                                logfc.threshold = fc.thr, min.pct = pct.thr,
                                ident.1 = "Male.0", ident.2 = "Female.0")
            print('Done! :)')
            
            print('Males.1 vs Females.1')
            m1f1 <- FindMarkers(obj, test.use = test,
                                logfc.threshold = fc.thr, min.pct = pct.thr,
                                ident.1 = "Male.1", ident.2 = "Female.1")
            print('Done! :)')
            
            # Design the name of the rds file
            rd.fname <- paste0(obj@project.name,
                               '.', 
                               add.fname,
                               test,
                               '.logfc.',
                               fc.thr, 
                               '.4g.apoe4sex.RData')
        } else {
            print(paste('In', test, 'check with mean.fxn'))
            # Find markers in all interseting variables
            print('Females.0 vs Females.1')
            f0f1 <- FindMarkers(obj, test.use = test, mean.fxn = mean.fxn,
                                logfc.threshold = fc.thr, min.pct = pct.thr,
                                ident.1 = "Female.0", ident.2 = "Female.1")
            print('Done! :)')
            
            print('Males.0 vs Males.1')
            m0m1 <- FindMarkers(obj, test.use = test, mean.fxn = mean.fxn,
                                logfc.threshold = fc.thr, min.pct = pct.thr,
                                ident.1 = "Male.0", ident.2 = "Male.1")
            print('Done! :)')
            
            print('Males.0 vs Females.0')
            m0f0 <- FindMarkers(obj, test.use = test, mean.fxn = mean.fxn,
                                logfc.threshold = fc.thr, min.pct = pct.thr,
                                ident.1 = "Male.0", ident.2 = "Female.0")
            print('Done! :)')
            
            print('Males.1 vs Females.1')
            m1f1 <- FindMarkers(obj, test.use = test, mean.fxn = mean.fxn,
                                logfc.threshold = fc.thr, min.pct = pct.thr,
                                ident.1 = "Male.1", ident.2 = "Female.1")
            print('Done! :)')
            
            # Design the name of the rds file
            rd.fname <- paste(obj@project.name,
                              add.fname,
                              test,
                              'logfc',
                              fc.thr,
                              'pct',
                              pct.thr,
                              '4g.apoe4sex.traj.meanfxn.RData',
                              sep = '.')
        }
        
    } else {
        message('Are you sure? this is without mean.fxn so currently not implemented.')
    }
    
    # Save file containing the interesting comparisons 
    save(f0f1, m0m1, m0f0, m1f1, file = rd.fname)
} 

# This function runs differentially expressed genes on the sex.apoe idents 2 groups.
find.markers.2.groups <- function(obj, test = 'poisson', meanfxn.flag = F, 
                                  fc.thr = 0.25, pct.thr = 0.1,
                                  add.title = '') {
    # Change idents:
    Idents(obj) <- 'sex.apoe4'
    
    if (test == 'DEseq2') { DefaultAssay(obj) <- 'RNA'}
    
    if (meanfxn.flag) { 
        print('with different mean.fxn')
        
        mean.fxn <- function(x) {
            return(log(x = rowMeans(x = x) + 1e-6, base = 2))
        }
        
        print(paste('In', test, 'check with mean.fxn'))
        print("ApoE4")
        Idents(obj) <- 'apoe_4'
        
        # Find markers in all interseting variables
        print('1 vs 0')
        e4.10 <- FindMarkers(obj, test.use = test, mean.fxn = mean.fxn,
                             logfc.threshold = fc.thr, min.pct = pct.thr,
                             ident.1 = "1", ident.2 = "0")
        print('Done! :)')
        
        print('Sex')
        Idents(obj) <- 'sex'
        
        print('Females vs. Males')
        sex.fm <- FindMarkers(obj, test.use = test, mean.fxn = mean.fxn,
                              logfc.threshold = fc.thr, min.pct = pct.thr,
                              ident.1 = "Female", ident.2 = "Male")
        print('Done! :)')
        
        # Design the name of the rds file
        rd.fname <- paste0(obj@project.name,
                           add.title,
                           '.pscount.1e-6.', 
                           test,
                           '.logfc.',
                           fc.thr,
                           '.pct.',
                           pct.thr,
                           '.4g.apoe4sex.2g.RData')
        
    } else {
        # Find markers in all intersting variables
        print('Currently implementred only for with mean.fxn')
    }
    
    # Save file containing the interesting comparisons 
    save(e4.10, sex.fm, file = rd.fname)
    
    # Create volcano plots
    
    # Install package if needed:
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("EnhancedVolcano")
    
    library(EnhancedVolcano)
    
    pdf(paste0(obj@project.name, add.title, '.sex.apoe.2g.degs.pdf'))
    print(
        EnhancedVolcano(e4.10,
                        lab = rownames(e4.10),
                        x = 'avg_log2FC',
                        y = 'p_val_adj',
                        title = 'ApoE4 Carriers vs Non-Carriers',
                        pCutoff = 0.05,
                        FCcutoff = 0.3,
                        pointSize = 3.0,
                        labSize = 4.0,
                        legendLabels = c("NS", expression(Log[2] ~ FC), 
                                         "p-value (adjusted)", 
                                         expression(p - value ~ and ~ log[2] ~ FC)),
                        drawConnectors = TRUE)
    )
    
    print(
        EnhancedVolcano(sex.fm,
                        lab = rownames(sex.fm),
                        x = 'avg_log2FC',
                        y = 'p_val_adj',
                        title = 'Females vs. Males',
                        pCutoff = 0.05,
                        FCcutoff = 0.3,
                        pointSize = 3.0,
                        labSize = 4.0,
                        legendLabels = c("NS", expression(Log[2] ~ FC), 
                                         "p-value (adjusted)", 
                                         expression(p - value ~ and ~ log[2] ~ FC)),
                        drawConnectors = TRUE)
    )
    
} 


#' This function prep an object to the comparison of sex and ApoE4.
pre.process.apoe4.sex <- function(obj) {
    # If object is microglia, remove macrophages and monocytes
    
    if (obj@project.name == 'microglia') {
        Idents(obj) <- 'cell.type'
        obj <- subset(obj, idents = 'Microglia')
        Idents(obj) <- 'sex'
    } else if (obj@project.name == 'opcs') {
        Idents(obj) <- 'cell.type'
        obj <- subset(obj, idents = 'OPCs')
    } else {
        message('No need in subsetting the object')
    }
    
    # Add apoe4 column
    obj <- add_apoe4_to_metadata(obj)
    print('added apoe4')
    # Add feature of interaction of ApoE4 and sex:
    obj <- Seurat::AddMetaData(object = obj,
                               metadata = interaction(obj@meta.data$sex, obj@meta.data$apoe_4),
                               col.name = "sex.apoe4")
    print('added to meta')
    
    # Remove people without specified sex
    Idents(obj) <- 'sex'
    obj <- subset(obj, idents = c('Male', 'Female'))
    print('take subset')
    
    # Remove individuals with bad QC:
    rm.id.QC <- c("10298957", "20201927", "20683921", "20834164", 
                  "3283241",  "44749170", "45115248", "50100932", 
                  "50406057", "51815338", "85666580", "90447310")
    # Remove QC samples
    obj <- obj[, !(obj@meta.data$projid %in% rm.id.QC)]
    print('remove 12 indsividuals from QC')
    
    return(obj)
}






# -------------------------------------------------------------------------- #
##                         Process Results Functions                      ####
# -------------------------------------------------------------------------- #
# Function to extract info from filename
# NOTE - very specific to the filename format
extract.info <- function(filename, pattern = "(?<=prob\\.).*(?=\\.poisson)") {
    parts <- str_split(basename(filename), "\\.")[[1]]
    list(
        celltype = parts[1],
        datatype = parts[2],
        description = str_extract(filename, pattern)
    )
} 


# Function to process a single sex-apoe4 DEGs file:
process.single.file <- function(file.path) {
    load(file.path)
    info <- extract.info(file.path)
    
    process.data <- function(data, sex) {
        data %>%
            filter(p_val_adj <= 0.05) %>%
            rownames_to_column('gene') %>%
            dplyr::select(gene, avg_log2FC) %>%
            mutate(
                sex = sex,
                celltype = info$celltype,
                datatype = info$datatype,
                description = info$description
            )
    }
    
    rbind(
        process.data(f0f1, "female"),
        process.data(m0m1, "male")
    )
}


# Function to create a matrix of all DEGs:
create.degs.matrix <- function(dir.path) {
    files <- list.files(dir.path, full.names = TRUE)
    degs <- map(files, safely(process.single.file)) %>%
        map("result") %>%
        compact() %>% 
        bind_rows()
    
    all.genes <- degs$gene %>% unique()
    
    degs <- degs %>%
        complete(gene = all.genes,
                 nesting(celltype, description, datatype, sex),
                 fill = list(avg_log2FC = NA)
                 ) %>%
        mutate(avg_log2FC = replace_na(avg_log2FC, 0)) %>%   # replace NAs with 0 and make positive LogFC to be unregulated in APOE4
        pivot_wider(
            id_cols = gene,
            names_from = c(celltype, description, datatype, sex),
            values_from = avg_log2FC,
            values_fill = 0  # fill NAs with 0
        )
    
    # Prepare as matrix for heatmap:
    mat <- as.matrix(degs[, -1])
    rownames(mat) <- degs$gene
    return(mat)
}



# -------------------------------------------------------------------------- #
##                           Preprocess for Plots                         ####
# -------------------------------------------------------------------------- #
# Define the function to get top genes based on given column and conditions
get.top.degs <- function(df, column, threshold, top.n, direction = "up") {
    if (direction == "up") {
        df %>%
            filter(!!sym(column) > threshold) %>%
            arrange(desc(!!sym(column))) %>%
            head(n = top.n) %>%
            pull(gene)
    } else {
        df %>%
            filter(!!sym(column) < -threshold) %>%
            arrange(!!sym(column)) %>%
            head(n = top.n) %>%
            pull(gene)
    }
}






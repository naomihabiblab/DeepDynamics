#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#                         Scrips to Run Analyses                          ####
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# ------------------------------------------------------------------------ #
##                    Null model vs. split by APOE4                     ####
# ------------------------------------------------------------------------ #
source("src/utils/utils.R")
source("src/utils/utils.TA.R") #TODO remove or change code completely, this is just for sig.clusters and default patholgoies
if (!exists("bulk.df", envir = .GlobalEnv)) {
    load(PATH.bulk.df)
} else {
    message("Note: bulk.df already exists in the global environment")
}

# make apoe_4 and sex as factors:
bulk.df <- bulk.df %>%
    mutate(
        across(
            c(sex, apoe_4_new), as.factor
        )
    )

fit.split.GAM <- function(
        feature, split = "apoe_4",
        test.use = "Chisq"
) {
    null.model <- gam(as.formula(
        paste0(feature, " ~ s(pseudotime)")
    ), 
    weights = prAD, data = bulk.df %>% filter(pseudotime <= 0.34))
    split.model <- gam(as.formula(
        paste0(feature, " ~ ", split, " + s(pseudotime, by = ", split, ")")
    ), 
    weights = prAD, data = bulk.df %>% filter(pseudotime <= 0.34))

    # Compare models
    anova.res <- anova.gam(null.model, split.model, test = test.use)
    # Return all components
    list(
        anova = anova.res,
        null.model = null.model,
        split.model = split.model
    )
    
}

fit.single.GAM <- function(
    feature, split = "apoe_4",
    test.use = "Chisq"
) {
    null.model <- gam(as.formula(
        paste0(feature, " ~ s(pseudotime)")
        ), 
        weights = prAD, data = bulk.df %>% filter(pseudotime <= 0.34))
    
    # Compare models
    anova(null.model, test = test.use)
}

fit.split.GAM.ABA <- function(
        feature, split = "apoe_4",
        test.use = "Chisq"
) {
    null.model <- gam(as.formula(
        paste0(feature, " ~ s(pseudotime)")
    ), 
    weights = ABA, data = bulk.df)
    split.model <- gam(as.formula(
        paste0(feature, " ~ ", split, " + s(pseudotime, by = ", split, ")")
    ), 
    weights = ABA, data = bulk.df)
    
    # Compare models
    anova.res <- anova(null.model, split.model, test = test.use)
    # Return all components
    list(
        anova = anova.res,
        null.model = null.model,
        split.model = split.model
    )
    
}


results <- c(SIG.CLUSTERS, DEFAULT.PATHOLOGIES) %>%
    set_names() %>%
    map_dfr( ~ {
        anova_res <- fit.split.GAM(feature = .x, test.use = "Chisq", split = 'apoe_4')
        # Extract p-value and test statistic from the second row (split.model)
        tibble(
            p.value = anova_res$`Pr(>Chi)`[2],
            test.stat = anova_res$Chi[2],
            deviance = anova_res$Deviance[2]
        )
    }, .id = "cluster")

results.s <- c(SIG.CLUSTERS, DEFAULT.PATHOLOGIES) %>%
    set_names() %>%
    map_dfr( ~ {
        res <- fit.split.GAM(feature = .x, test.use = "Chisq", split='sex')
        
        tibble(
            p.value = res$anova$`Pr(>Chi)`[2],
            test.stat = res$anova$Chi[2],
            deviance = res$anova$Deviance[2],
            AIC.diff = AIC(res$null.model) - AIC(res$split.model),
            BIC.diff = BIC(res$null.model) - BIC(res$split.model)
        )
    }, .id = "cluster")

results <- c(SIG.CLUSTERS, DEFAULT.PATHOLOGIES) %>%
    set_names() %>%
    map_dfr( ~ {
        res <- fit.split.GAM(feature = .x, test.use = "Chisq", split='apoe_4')
        
        tibble(
            p.value = res$anova$`Pr(>Chi)`[2],
            test.stat = res$anova$Chi[2],
            deviance = res$anova$Deviance[2],
            AIC.diff = AIC(res$null.model) - AIC(res$split.model),
            BIC.diff = BIC(res$null.model) - BIC(res$split.model)
        )
    }, .id = "cluster")

results.ABA.e4 <- c(SIG.CLUSTERS, DEFAULT.PATHOLOGIES) %>%
    set_names() %>%
    map_dfr( ~ {
        res <- fit.split.GAM.ABA(feature = .x, test.use = "Chisq", split='apoe_4')
        
        tibble(
            p.value = res$anova$`Pr(>Chi)`[2],
            test.stat = res$anova$Chi[2],
            deviance = res$anova$Deviance[2],
            AIC.diff = AIC(res$null.model) - AIC(res$split.model),
            BIC.diff = BIC(res$null.model) - BIC(res$split.model)
        )
    }, .id = "cluster")

results.ABA.sex <- c(SIG.CLUSTERS, DEFAULT.PATHOLOGIES) %>%
    set_names() %>%
    map_dfr( ~ {
        res <- fit.split.GAM.ABA(feature = .x, test.use = "Chisq", split='sex')
        
        tibble(
            p.value = res$anova$`Pr(>Chi)`[2],
            test.stat = res$anova$Chi[2],
            deviance = res$anova$Deviance[2],
            AIC.diff = AIC(res$null.model) - AIC(res$split.model),
            BIC.diff = BIC(res$null.model) - BIC(res$split.model)
        )
    }, .id = "cluster")

results.E3E4 <- c(SIG.CLUSTERS, DEFAULT.PATHOLOGIES) %>%
    set_names() %>%
    map_dfr( ~ {
        res <- fit.split.GAM.apoe4(feature = .x, test.use = "Chisq", split='sex')
        
        tibble(
            p.value.E3 = res$anova.E3$`Pr(>Chi)`[2],
            p.value.E4 = res$anova.E4$`Pr(>Chi)`[2]
        )
    }, .id = "cluster")

results.E3E4.ABA <- c(SIG.CLUSTERS, DEFAULT.PATHOLOGIES) %>%
    set_names() %>%
    map_dfr( ~ {
        res <- fit.split.GAM.apoe4.ABA(feature = .x, test.use = "Chisq", split='sex')
        
        tibble(
            p.value.E3 = res$anova.E3$`Pr(>Chi)`[2],
            p.value.E4 = res$anova.E4$`Pr(>Chi)`[2]
        )
    }, .id = "cluster")



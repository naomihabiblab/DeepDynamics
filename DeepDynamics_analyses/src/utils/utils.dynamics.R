#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#                           Dynamics Functions                            ####
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# -------------------------------------------------------------------------- #
##                               Constants                                ####
# -------------------------------------------------------------------------- #

# Colors for 2-group comparison
risk.col <- "#ef6c00"
norisk.col <- "#4db6ac"

# Colors for 4-group comparison]
risk.risk.col <- ""
risk.norisk.col <- ""
norisk.risk.col <- ""
norisk.norisk.col <- ""

# Colors for split ranking
high.col <- "#ef6c00"
mid.col <- "white"
low.col <- "#4db6ac"

# KS temporal score:
start.col <- "#A8D5BA"  # Light green
middle.col <- "#FFE1A8"    # Light orange
end.col <- "#FFBFA8"    # Light red

# Colors for the heatmap
border.col <- "grey26"   # Border color


# Consts (used as default values):
DEFAULT.PATHOLOGIES <- c('sqrt.amyloid_mf', 'sqrt.tangles_mf', 'cogng_demog_slope')


# To define SIG.CLUSTERS, we need to first load the 500.
# if (!exists("data", envir = .GlobalEnv)) {
#     data <- anndata::read_h5ad(PATH.500had)
#     message('Loaded 500 data as needed for dynamics')
# } else {
#     message('500 data already exists in the global environment. No changes were made.')
# }
data <- anndata::read_h5ad(PATH.500had)
tryCatch({
    SIG.CLUSTERS <- data$uns$celmod$test.corrs %>% #old.data$uns$celmod$cor
        py_to_r() %>% filter(adj.pval < .005 & corr > 0) %>%
        rownames()
}, error = function(e) {
    SIG.CLUSTERS <- c('Oli.11', 'Mic.12', 'Oli.7', 'Ast.10', 'Mic.13',
                      'Ast.7', 'Oli.3', 'Oli.4', 'Ast.4', 'Inh.6',
                      'Exc.1', 'Exc.3', 'Inh.15', 'OPC.1', 'Inh.16',
                      'Inh.5', 'Exc.8', 'End.3', 'Inh.12', 'Oli.5',
                      'Mic.1', 'Exc.12', 'Inh.7', 'Ast.6', 'End.1',
                      'Ast.2', 'OPC.2', 'Mic.14', 'Ast.1', 'Oli.9') # (only SIG for TA)
})



# -------------------------------------------------------------------------- #
##                            Fit Dynamics Model                          ####
# -------------------------------------------------------------------------- #
fit.dynamics <- function(pseudotime,
                         features,
                         trajectory.probs,
                         trajectory.terminal.pseudotime = setNames(rep(1, ncol(trajectory.probs)), colnames(trajectory.probs)),
                         min.prob.clip = 0,
                         min.prob.for.fit = 0,
                         evaluate.fit = T,
                         gam_method = "GCV.Cp",
                         ...) {
    # Equally spaced pseudotime values to predict dynamics for
    xs.pred <- seq(min(pseudotime, na.rm = T),
                   max(pseudotime, na.rm = T),
                   length.out = 50)
    if ("Time" %in% colnames(pseudotime)) {
        colnames(pseudotime)[colnames(pseudotime) == "Time"] <- "ps"
    }
    
    pb <- progress_bar$new(
        format = "Fitting dynamics along trajectories :current/:total [:bar] :percent in :elapsed. ETA :eta",
        total = ncol(features) * ncol(trajectory.probs),
        clear = F,
        width = 100,
        force = T
    )
    df <- data.frame(ps = pseudotime, trajectory.probs, features)
    
    res <- lapply(colnames(trajectory.probs), function(t) {
        # Set trajectory range between pseudotime of first sample with weight in trajectory over `min.prob.clip`
        # and pseudotime specified by trajectory's terminal pseudotime
        beg <- min(df[df[, t] >= min.prob.clip, ]$ps, na.rm = T) %>% ifelse(is.infinite(.), 0, .)
        end <- trajectory.terminal.pseudotime[t] %>% ifelse(is.infinite(.), 1, .)
        
        # Specify trajectory specific pseudotime values to predict dynamics for
        xs <- unique(c(beg, xs.pred[(beg <= xs.pred) &
                                        (xs.pred <= end)], end))
        
        # Fit- and predict dynamics for features along specific trajectory
        res <- lapply(colnames(features), function(p) {
            # subset dataframe for specified trajectory and parameter to regress
            .df <- data.frame(
                x = df$ps,
                y = df[, p],
                w = df[, t],
                row.names = rownames(df)
            ) %>%
                tidyr::drop_na() %>%
                filter(beg <= x & x <= end) %>%
                filter(w >= min.prob.for.fit)
            
            res <- list(
                fit = data.frame(
                    x = NA,
                    fit = NA,
                    se = NA,
                    fit_sd = NA,
                    se.fit_sd = NA,
                    feature = p,
                    trajectory = t
                )[0, ],
                prd = data.frame(
                    x = xs,
                    fit = NA,
                    se = NA,
                    fit_sd = NA,
                    se.fit_sd = NA,
                    feature = p,
                    trajectory = t
                )
            )
            tryCatch({
                res <- .fit.dynamics(.df, xs, evaluate.fit = evaluate.fit, ...)
                res <- lapply(res, function(.df)
                    .df %>% mutate(feature = p, trajectory = t))
            }, error = function(e)
                warning(
                    "Failed computing dynamics for feature [",
                    p,
                    "] in trajectory [",
                    t,
                    "]. Error: ",
                    e
                ))
            
            pb$tick()
            return(res)
        })
        
        # Merge results of fit- and predict from different features
        return(lapply(1:length(res[[1]]), function(i)
            do.call(rbind, lapply(res, "[[", i))))
    })
    
    # Merge results of fit- and predict from different trajectories
    res <- lapply(1:length(res[[1]]), function(i)
        do.call(rbind, lapply(res, "[[", i)))

    
    names <- c("fitted.vals", "pred.vals")
    if (evaluate.fit)
        names <- c(names, "evaluations")
    return(res %>% `names<-`(names))
}

.fit.dynamics <- function(df,
                          pred.x = seq(min(df$x, na.rm = T), max(df$x, na.rm = T), length.out =
                                           50),
                          bootstrap = T,
                          bootstrap.proportion = 1,
                          bootstrap.iterations = 100,
                          evaluate.fit = T,
                          gam_method = "GCV.Cp",
                          ...) {
    if (!bootstrap)
        bootstrap.proportion = bootstrap.iterations = 1
    
    res <- lapply(1:bootstrap.iterations, function(i) {
        # Set subset of ids to use: sample/sub-sample with/without repetition
        ids <- sample(
            rownames(df),
            size = ceiling(nrow(df) * bootstrap.proportion),
            replace = bootstrap
        )
        
        # Fit model [y~s(x), weighted by w]
        args <- modifyList(list(
            formula = y ~ s(x),
            data = df[ids, ],
            weights = df[ids, ]$w
        ) ,
        list(...))
        fit  <- do.call(mgcv::gam, args)
        
        fit.vals <- data.frame(fit$model, predict(fit, df[ids, ], se.fit =
                                                      T)) %>% mutate(i = i)
        prd.vals <- data.frame(x = pred.x, i, predict(fit, data.frame(x =
                                                                          pred.x), se.fit = T))
        
        
        if (!evaluate.fit)
            return(list(fit = fit.vals, prd = prd.vals))
        
        # Test fitted model
        models <- list(
            # against null model y~1
            null = do.call(mgcv::gam, modifyList(
                args, list(formula = y ~ 1, weights = NULL)
            )),
            # against un-weighted y~s(x) (i.e. not trajectory specific)
            unweighted = do.call(mgcv::gam, modifyList(args, list(weights =
                                                                      NULL)))
        )
        
        evaluations <- do.call(rbind, lapply(names(models), function(m)
            data.frame(
                comparison = m,
                anova(fit, models[[m]], test = "F")[2, c("Deviance", "F", "Pr(>F)")],
                i = i
            )))
        
        return(list(
            fit = fit.vals,
            prd = prd.vals,
            evaluations = evaluations
        ))
    })
    
    # Merge results of different subsampling iterations
    res <- lapply(1:length(res[[1]]), function(i)
        do.call(rbind, lapply(res, "[[", i)))
    
    # Summarise subsampling results for fitted- and predicted values
    for (j in 1:2)
        res[[j]] <- res[[j]] %>%
        dplyr::select(-i) %>%
        group_by_at(vars(-fit, -se.fit)) %>%
        summarise_all(.funs = list(mean = mean, sd = sd)) %>%
        `colnames<-`(gsub("_mean", "", colnames(.)))
    
    return(res)
}

.fit.dynamics.model <- function(df, ...) {
    args <- modifyList(list(
        formula = y ~ s(x),
        data = df[ids, ],
        weights = df[ids, ]$w
    ) , list(...))
    fit  <- do.call(mgcv::gam, args)
    
    return(fit)
}


# -------------------------------------------------------------------------- #
##                        Functions for Ranking Splits                    ####
# -------------------------------------------------------------------------- #
#' This function assumes dynamics is the OG structure that includes fitted
#'  and predicted values.
process.feature.trajectory <- function(dynamics, feature, trajectory, n_samples = 1000) {
    print(paste("Processing feature", feature, "and trajectory", trajectory))
    dyn <- dynamics$fitted.vals %>%
        py_to_r() %>% 
        filter(str_starts(feature, .env$feature) & trajectory == .env$trajectory)
    
    matrices <- concat.split.dynamics(dyn)
    bootstrap_results <- bootstrap.sample(matrices, n_samples)
    
    distance_metrics <- c("ks.temporal")
    
    
    
    if ("ks.temporal" %in% distance_metrics) {
        temp.boot <- split.temp.res(bootstrap_results)
        results.start <- map_dfr(
            distance_metrics,
            ~ tibble(
                feature = feature,
                trajectory = trajectory,
                split = "ApoE4",
                distance_metric = .x,
                result = rank.splits.at.once(temp.boot$start, .x)$statistic %>%
                    as.numeric,
                p.value = rank.splits.at.once(temp.boot$start, .x)$p.value %>%
                    as.numeric
            )
        ) %>% mutate(split = 'start')
        
        results.mid <- map_dfr(
            distance_metrics,
            ~ tibble(
                feature = feature,
                trajectory = trajectory,
                split = "ApoE4",
                distance_metric = .x,
                result = rank.splits.at.once(temp.boot$mid, .x)$statistic %>%
                    as.numeric,
                p.value = rank.splits.at.once(temp.boot$mid, .x)$p.value %>%
                    as.numeric
            )
        )  %>% mutate(split = 'mid')
        
        results.end <- map_dfr(
            distance_metrics,
            ~ tibble(
                feature = feature,
                trajectory = trajectory,
                split = "ApoE4",
                distance_metric = .x,
                result = rank.splits.at.once(temp.boot$end, .x)$statistic %>%
                    as.numeric,
                p.value = rank.splits.at.once(temp.boot$end, .x)$p.value %>%
                    as.numeric
            )
        ) %>% mutate(split = 'end')
        
        results <- bind_rows(results.start, results.mid, results.end)
        
    } else {
        stop("Invalid distance metric")
    }
    
    return(results)
}

#' This function assumes to accept only fitted values of a specific trajectory.
concat.split.dynamics <- function(dyn) {
    # Function to summarize duplicates
    summarize_duplicates <- function(x) {
        if (length(x) == 1) return(x)
        if (all(is.na(x))) return(NA)
        mean(x, na.rm = TRUE)
    }
    
    # Reshape for fit values
    fit_df <- dyn %>%
        dplyr::select(feature, x, split, fit) %>%
        group_by(split, x) %>%
        summarise(fit = summarize_duplicates(fit), .groups = "drop") %>%
        pivot_wider(
            id_cols = split,
            names_from = x,
            values_from = fit
        ) %>%
        arrange(split)
    
    # Reshape for se.fit values
    se_fit_df <- dyn %>%
        dplyr::select(feature, x, split, se.fit) %>%
        group_by(split, x) %>%
        summarise(se.fit = summarize_duplicates(se.fit), .groups = "drop") %>%
        pivot_wider(
            id_cols = split,
            names_from = x,
            values_from = se.fit
        ) %>%
        arrange(split)
    
    # Convert to matrices
    fit_matrix <- as.matrix(fit_df[, -1])  # Remove the trajectory column
    se_fit_matrix <- as.matrix(se_fit_df[, -1])  # Remove the trajectory column
    
    # Ensure column names are sorted
    col_order <- order(as.numeric(colnames(fit_matrix)))
    fit_matrix <- fit_matrix[, col_order, drop = FALSE]
    se_fit_matrix <- se_fit_matrix[, col_order, drop = FALSE]
    
    # Function to merge complementary columns
    merge_complementary_columns <- function(matrix) {
        if (nrow(matrix) < 2) return(matrix)  # Return if there's only one row
        
        col <- 1
        while (col < ncol(matrix)) {
            if (is.na(matrix[1, col]) && !is.na(matrix[2, col]) &&
                col + 1 <= ncol(matrix) && !is.na(matrix[1, col + 1]) && is.na(matrix[2, col + 1])) {
                matrix[1, col] <- matrix[1, col + 1]
                matrix[2, col] <- matrix[2, col]
                matrix <- matrix[, -c(col + 1), drop = FALSE]
            } else {
                col <- col + 1
            }
        }
        return(matrix)
    }
    
    # Merge complementary columns
    fit_matrix <- merge_complementary_columns(fit_matrix)
    se_fit_matrix <- merge_complementary_columns(se_fit_matrix)
    
    # Remove columns with NA values
    na_cols <- apply(fit_matrix, 2, function(x) any(is.na(x))) | 
        apply(se_fit_matrix, 2, function(x) any(is.na(x)))
    fit_matrix <- fit_matrix[, !na_cols, drop = FALSE]
    se_fit_matrix <- se_fit_matrix[, !na_cols, drop = FALSE]
    
    return(list(fit = fit_matrix, se.fit = se_fit_matrix))
}

# Function to bootstrap sample from each dynamics split
bootstrap.sample <- function(dynamics_matrices, n = 1000) {
    fit_matrix <- dynamics_matrices$fit
    se_fit_matrix <- dynamics_matrices$se.fit
    
    n_splits <- nrow(fit_matrix)
    n_timepoints <- ncol(fit_matrix)
    
    bootstrap_samples <- list()
    
    for (split in 1:n_splits) {
        bootstrap_matrix <- matrix(0, nrow = n, ncol = n_timepoints)
        
        for (timepoint in 1:n_timepoints) {
            fit <- fit_matrix[split, timepoint]
            se <- se_fit_matrix[split, timepoint]
            
            # Use 1.96 * se for the 95% confidence interval
            lower_ci <- fit - 1.96 * se
            upper_ci <- fit + 1.96 * se
            
            # # Sample from a truncated normal distribution
            # bootstrap_matrix[, timepoint] <- truncnorm::rtruncnorm(
            #     n, 
            #     a = lower_ci, 
            #     b = upper_ci, 
            #     mean = fit, 
            #     sd = se
            # )
            # 
            # TODO delete one of them
            bootstrap_matrix[, timepoint] <- rnorm(
                n,
                mean = fit, 
                sd = se
            )
        }
        
        bootstrap_samples[[paste0("split", split)]] <- bootstrap_matrix
    }
    
    return(bootstrap_samples)
}

# Function to easily split bootstrap results for temporal tests:
split.temp.res <- function(bootstrap_results) {
    # Function to split a matrix into thirds by columns
    split_matrix <- function(matrix) {
        cols <- ncol(matrix)
        third <- cols %/% 3  # Integer division
        remainder <- cols %% 3
        
        # Determine split points
        split1 <- third + if(remainder > 0) 1 else 0
        split2 <- 2 * third + if(remainder > 1) 2 else if(remainder > 0) 1 else 0
        
        list(
            start = matrix[, 1:split1],
            mid = matrix[, (split1+1):split2],
            end = matrix[, (split2+1):cols]
        )
    }
    
    # Apply the split function to both elements of bootstrap_results
    bootstrap_results_split <- lapply(bootstrap_results, split_matrix)
    
    # Create the three new lists
    bootstrap_results_start <- list(
        split1 = bootstrap_results_split$split1$start,
        split2 = bootstrap_results_split$split2$start
    )
    
    bootstrap_results_mid <- list(
        split1 = bootstrap_results_split$split1$mid,
        split2 = bootstrap_results_split$split2$mid
    )
    
    bootstrap_results_end <- list(
        split1 = bootstrap_results_split$split1$end,
        split2 = bootstrap_results_split$split2$end
    )
    
    return(list(
        start = bootstrap_results_start,
        mid = bootstrap_results_mid,
        end = bootstrap_results_end
    ))
}

# Function to rank the splits.
rank.splits.at.once <- function(bootstrap_samples, distance_metric = "ks.temporal"){
    # Separate the two splits
    split1_samples <- bootstrap_samples$split1
    split2_samples <- bootstrap_samples$split2
    
    # Calculate the distance metric
    if (distance_metric == "ks.temporal") {
        ks.res <- ks.test(split1_samples, split2_samples)
        return(list("statistic" = ks.res$statistic, "p.value" = ks.res$p.value))
    } else {
        stop("Invalid distance metric")
    }
    
}



# -------------------------------------------------------------------------- #
##                       Functions for Visualizations                     ####
# -------------------------------------------------------------------------- #
plot.dynamics <- function(fits,
                          preds,
                          cols = NULL,
                          facet.by = c("trajectory", "feature"),
                          overlap.pseudotime = NULL,
                          include.points = T,
                          min.point.alpha = .1,
                          label = T,
                          points.above.thr = F,
                          no.legend = T,
                          legend.position = "none",
                          ymin = NA, ymax = NA,
                          xmin = NA, xmax = NA, ...) {
    facet.by = match.arg(facet.by)
    group.by = setdiff(c("trajectory", "feature"), facet.by)
    
    groups <- unique(fits[, group.by]) %>% as.character()
    if (is.null(cols)) {
        cols <- setNames(scales::hue_pal()(length(groups)), groups)
    }
    
    p <- ggplot(fits) +
        facet_wrap(paste0("~", facet.by), ...) +
        scale_x_continuous(expand = c(0, 0)) +
        theme_classic() +
        theme(strip.background = element_blank(),
              legend.position = legend.position) +
        guides(colour = "none", shape = "none")
    
    
    if (is.na(ymax)) ymax <- max(preds$fit+2*preds$se.fit, na.rm = T)
    if (is.na(ymin)) ymin <- min(preds$fit-2*preds$se.fit, na.rm = T)
    
    if (is.na(xmin)) xmin <- min(preds$x, na.rm = T)
    
    if (include.points) {
        ymax <- max(c(fits$y, ymax), na.rm = T) * 1.05
        ymin <- min(c(fits$y, ymin), na.rm = T) * .95
    }
    
    p <- p + scale_y_continuous(expand = c(0, 0), limits = c(ymin, ymax))
    
    
    if (!is.null(overlap.pseudotime)) {
        p <- p +
            geom_ribbon(
                aes(x, ymax = ymax, ymin = ymin),
                data.frame(
                    x = c(xmin, overlap.pseudotime),
                    ymax = rep(ymax, 2),
                    ymin = rep(ymin, 2)
                ),
                fill = "black",
                alpha = .075,
                show.legend = F
            ) +
            geom_vline(
                xintercept = overlap.pseudotime,
                linetype = "dashed",
                show.legend = F,
                size = .25
            )
    }
    
    
    if (include.points) {
        for (g in groups) {
            p <- p +
                geom_point(
                    aes(x, y, color = X.weights., alpha = X.weights.),
                    fits[fits[, group.by] == g, ],
                    inherit.aes = F,
                    show.legend = F,
                    size = 0.25
                ) +
                scale_color_gradientn(colors = colorRampPalette(c("white", cols[[g]]))(10)[-1]) +
                scale_alpha(range = c(min.point.alpha, 1)) +
                labs(color = g, alpha = g) +
                new_scale("alpha") +
                new_scale("color") + #+
                guides(color_new = guide_legend(g),
                       alpha_new = guide_legend(g))
        }
    } else {
        if (points.above.thr) {
            for (g in groups) {
                p <- p +
                    geom_point(aes(x, y, colour = split, alpha = X.weights.),
                               fits[fits[, group.by] == g, ],
                               inherit.aes = F) +
                    scale_color_manual(values = setNames(c(cols[[g]], 'white'), c(T, F))) +
                    labs(color = g, alpha = g)
            }
        }
    }
    for (g in groups) {
        message(paste("g:", g))
        p <- p +
            geom_ribbon(
                aes_string(
                    "x",
                    "fit",
                    ymin = "fit-se.fit",
                    ymax = "fit+se.fit",
                    fill = group.by,
                    color = group.by
                ),
                preds[preds[, group.by] == g, ],
                alpha = .2,
                linetype = "dashed",
                size = 0.25,
                inherit.aes = F,
                show.legend = T
            ) +
            geom_line(
                aes_string("x", "fit", color = group.by),
                preds[preds[, group.by] == g, ],
                size = 0.5,
                inherit.aes = F,
                show.legend = T
            ) +
            scale_color_manual(values = cols[[g]]) +
            scale_fill_manual(values = cols[[g]], guide = F) +
            
            new_scale("color") +
            new_scale("fill") +
            guides(color_new = guide_legend(g), fill_new = guide_legend(g))
    }
    
    if (label) {
        p <- p + ggrepel::geom_label_repel(
            aes_string(
                "x",
                "fit",
                label = group.by,
                color = group.by,
                position = c(Inf, Inf, Inf, Inf)
            ),
            fits %>%
                group_by(feature, trajectory) %>%
                slice_max(x, n = 1),
            min.segment.length = unit(0, "pt")
        )
        
    }
    
    return(p)
}


plot.dynamics.wrapper <- function(dynamics, features, trajectories = c('prAD', 'ABA'), ...) {
    if ("pandas.core.frame.DataFrame" %in% class(dynamics$fitted.vals))
        dynamics$fitted.vals <- py_to_r(dynamics$fitted.vals)
    
    if ("pandas.core.frame.DataFrame" %in% class(dynamics$pred.vals))
        dynamics$pred.vals <- py_to_r(dynamics$pred.vals)
    
    fits <- dynamics$fitted.vals %>% filter(feature %in% features & trajectory %in% trajectories)
    preds <- dynamics$pred.vals %>% filter(feature %in% features & trajectory %in% trajectories)
    
    # #to remove
    # fits <- fits %>% filter(trajectory == 'prAD')
    # pred <- pred %>% filter(trajectory == 'prAD')
    return(plot.dynamics(fits, preds, ...))
}



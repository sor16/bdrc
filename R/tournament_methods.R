#' @importFrom ggplot2 ggplot geom_boxplot stat_boxplot geom_line geom_point labs
plot_tournament_fun <- function(x, type = 'boxplot'){
    if(type == 'boxplot'){
        # DmC: Deviance minus the constant C
        DmC_post_dat <- lapply(x$contestants, function(m){
            data.frame(model = class(m), DmC = c(-2 * m$posterior_log_likelihood))
        })
        DmC_post_dat <- do.call(rbind, DmC_post_dat)
        WAIC_dat <- lapply(x$contestants, function(m){
            data.frame(model = class(m), WAIC = c(m$WAIC))
        })
        WAIC_dat <- do.call(rbind, WAIC_dat)
        p <- ggplot(data = DmC_post_dat, aes(x = .data$model, y = .data$DmC)) +
            geom_boxplot(size = 0.4, width = 0.4, color = "black", outlier.size = 0.2, outlier.shape = 21, outlier.fill = "gray90", fill = "gray90") +
            stat_boxplot(geom = 'errorbar', width = 0.2) +
            geom_line(data = WAIC_dat, aes(x = .data$model, y = .data$WAIC, group = 1), color = 'gray30') +
            geom_point(data = WAIC_dat, aes(x = .data$model, y = .data$WAIC), size = 3, shape = 23, fill = 'red2', color = 'black') +
            theme_bdrc() +
            labs(y = '', x = '') +
            ggtitle('-2( Posterior log-likelihood ) & WAIC')
    }
    return(p)
}


#' @importFrom ggplot2 autoplot geom_text geom_label geom_segment scale_colour_manual theme_classic theme
#' @importFrom gridExtra arrangeGrob
plot_tournament_grob <- function(x, type = 'panel', transformed = FALSE){
    ylim_dat <- lapply(x$contestants, function(m){
        if(grepl('0' ,class(m))){
            sig_ylim <- c(min = 0, max = 1.1 * m$param_summary['sigma_eps', 'upper'])
        }else{
            sig_ylim <- c(min = 0, max = 1.1 * max(m$sigma_eps_summary$upper))
        }
        if(grepl('g', class(m))){
            f_ylim <- c(min = min(1, 0.9 * min(m$f_summary$lower)), max = max(3.5, 1.1 * max(m$f_summary$upper)))
        }else{
            f_ylim <- c(min = min(1, 0.9 * m$param_summary['b', 'lower']), max = max(3.5, 1.1 * m$param_summary['b', 'upper']))
        }
        max_res <- 1.1 * max(abs(get_residuals_dat(m)[, c('m_upper', 'm_lower', 'r_median', 'r_lower', 'r_upper')]), na.rm = TRUE)
        data.frame(sigma_eps_min = sig_ylim[1], sigma_eps_max = sig_ylim[2], f_min = f_ylim[1], f_max = f_ylim[2], residuals_min = -max_res, residuals_max = max_res)
    })
    ylim_dat <- do.call('rbind', ylim_dat)
    ylim_dat <- sapply(colnames(ylim_dat), function(col){
        if(grepl('min', col)){
            min(ylim_dat[, col])
        }else{
            max(ylim_dat[, col])
        }
    })
    ylim_dat <- c(ylim_dat, rating_curve_min = NA, rating_curve_max = NA)
    if(type == "residuals"){
        plot_list <- lapply(x$contestants, function(m){
            autoplot(m, type = type, title = class(m), ylim = ylim_dat[c('residuals_min', 'residuals_max')])
        })
        p <- do.call(arrangeGrob, c(plot_list, ncol=2))
    }else if(type == "sigma_eps"){
        plot_list <- lapply(x$contestants, function(m){
            autoplot(m, type = type, title = class(m), ylim = ylim_dat[c('sigma_eps_min', 'sigma_eps_max')])
        })
        p <- do.call(arrangeGrob, c(plot_list, ncol = 2))
    }else if(type == "f"){
        plot_list <- lapply(x$contestants, function(m){
            autoplot(m, type = type, title = class(m), ylim = ylim_dat[c('f_min', 'f_max')])
        })
        p <- do.call(arrangeGrob, c(plot_list, ncol = 2))
    }else if(type %in% c("rating_curve", "rating_curve_mean")){
        plot_list <- lapply(x$contestants, function(m){
            autoplot(m, type = type, transformed = transformed, title = class(m))
        })
        p <- do.call(arrangeGrob, c(plot_list, ncol = 2))
    }else if(type == 'convergence_diagnostics'){
        plot_list <- lapply(x$contestants, function(m){
            plot_grob(m, type = type)
        })
        p <- do.call(arrangeGrob, c(plot_list, nrow = 4))
    }else if(type == 'panel'){
        panel_types <- c('rating_curve', 'residuals', 'f', 'sigma_eps')
        p <- lapply(x$contestants, function(m){
            plot_list <- lapply(panel_types, function(ty){
                ylim_vec <- ylim_dat[c(paste0(ty, '_min'), paste0(ty, '_max'))]
                plot_fun(m, type = ty, transformed = transformed, param = NULL, ylim = ylim_vec)
            })
            do.call('arrangeGrob', c(plot_list, ncol = round(sqrt(length(panel_types)))))
        })
        names(p) <- sapply(x$contestants, class)
    }else if(type == 'tournament_results'){
        loc_pts <- data.frame(x = c(seq(0, 3), 0.5, 2.5, 1.5),
                              y = c(rep(0, 4), 1, 1, 2),
                              xend = c(0.5, 0.5, 2.5, 2.5, 1.5, 1.5, NA),
                              yend = c(rep(1, 4), 2, 2, NA),
                              model = c(sapply(x$contestants, function(m)class(m)),
                                        x$summary$model[5:6],
                                        paste0('Tournament winner  =>  ', class(x$winner),
                                               paste0(rep(' ', 20), collapse = ' '))))
        method <- if(grepl("WAIC", x$info$method)) "WAIC" else if(grepl("DIC", x$info$method)) "DIC" else "PMP"
        n_digits <- if(method == "PMP") 2 else 1
        result_dat <- data.frame(mc_stat = paste0(ifelse(method == 'PMP', paste0("P="), paste0(method, "=\n")), format(round(x$summary[[method]], digits = n_digits), nsmall = n_digits)),
                                 winner = x$summary$winner,
                                 x = c(loc_pts$x[1:4], 0.8 * (loc_pts$x[5:6] - 1.5) + 1.5),
                                 y = loc_pts$y[1:6] + 0.5)
        comparison_results <- ggplot() +
            geom_segment(data = loc_pts[1:6, ], aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend)) +
            geom_text(data = result_dat, aes(x = .data$x, y = .data$y, label = .data$mc_stat, color = .data$winner, size = 7)) +
            geom_label(data = loc_pts[5:7, ], aes(x = .data$x, y = .data$y, label = .data$model), label.padding = unit(0.5, "lines"), label.size = 0, color = "Black", fill = "white", size = 6) +
            scale_colour_manual(values = c("red", "green3")) +
            theme_classic() +
            theme(line = element_blank(),
                  text = element_blank(),
                  plot.margin = unit(c(0, 1, -0.5, 3), "cm"),
                  legend.position = "none")
        residual_plots <- plot_tournament_grob(x, type = 'residuals')
        p <- arrangeGrob(comparison_results,arrangeGrob(grobs = residual_plots$grobs, ncol = 4), nrow = 2, heights = c(1, 1))
    }else{
        stop('type is not recognized.')
    }
    return(p)
}

#' Print method for discharge rating curve tournament
#'
#' Print the results of a tournament of discharge rating curve model comparisons
#' @param x An object of class "tournament"
#' @param ... Not used in this function
#' @seealso  \code{\link{tournament}} to run a discharge rating curve tournament, \code{\link{summary.tournament}} for summaries and \code{\link{plot.tournament}} for visualizing the mode comparison.
#' @export
print.tournament <- function(x, ...){
    cat(paste0('Tournament winner: ', class(x$winner)))
}

#' Summary method for a discharge rating curve tournament
#'
#' Print the summary of a tournament of model comparisons. This function allows for an efficient and fast re-run of the tournament with different methods or winning criteria.
#'
#' @param object An object of class "tournament"
#' @param method Optional; a string specifying the method to use for the summary. If NULL, uses the method from the original tournament. Options are "WAIC", "DIC", or "PMP".
#' @param winning_criteria Optional; specifies new winning criteria for the summary. If NULL, uses the criteria from the original tournament. See Details in \code{\link{tournament}} for proper formatting.
#' @param ... Not used in this function
#'
#' @details
#' If either \code{method} or \code{winning_criteria} is provided, the function re-runs the tournament with the new parameters using the fitted models.
#'
#' @return Prints the summary to the console.
#'
#' @seealso \code{\link{tournament}} to run a discharge rating curve tournament and \code{\link{plot.tournament}} for visualizing the model comparison
#'
#' @examples
#' \donttest{
#' data(krokfors)
#' set.seed(1)
#' t_obj <- tournament(Q ~ W, krokfors, num_cores = 2)
#' summary(t_obj)
#'
#' # Re-run summary with different method
#' summary(t_obj, method = "DIC")
#'
#' # Re-run summary with different winning criteria
#' summary(t_obj, winning_criteria = "Delta_WAIC > 3")
#' }
#'
#' @export
summary.tournament <- function(object, method = NULL, winning_criteria = NULL, ...) {
    print_warning <- FALSE
    if(is.null(method) & is.null(winning_criteria)){
        out <- tournament_summary_output(object$summary, object$info$method, object$info$winning_criteria)
    } else {
        if(!is.null(method) & is.null(winning_criteria)){
            new_t_obj <- tournament(object$contestants, method = method, verbose = FALSE)
        } else if(is.null(method) & !is.null(winning_criteria)){
            new_t_obj <- tournament(object$contestants, method = object$info$method, winning_criteria = winning_criteria, verbose = FALSE)
        } else if(!is.null(method) & !is.null(winning_criteria)){
            new_t_obj <- tournament(object$contestants, method = method, winning_criteria = winning_criteria, verbose = FALSE)
        }
        out <- tournament_summary_output(new_t_obj$summary, new_t_obj$info$method, new_t_obj$info$winning_criteria)
        if(new_t_obj$info$method=="PMP") {
            print_warning <- TRUE
        }
    }
    invisible(out)
    if(print_warning) cat( '\n\u26A0 Warning: The Harmonic Mean Estimator (HME) is used to estimate the Bayes Factor for the posterior model probability (PMP), which is known to be unstable and potentially unreliable. We recommend using method "WAIC" (Widely Applicable Information Criterion) for model comparison instead.\n' )
}

#' Autoplot method for discharge rating curve tournament
#'
#' Compare the four discharge rating curves from the tournament object in different ways
#'
#' @param object An object of class "tournament"
#' @param type A character denoting what type of plot should be drawn. Possible types are
#' @param ... Not used in this function
#' \describe{
#'   \item{\code{boxplot}}{Creates a boxplot of the posterior log-likelihood values transformed to the deviance scale.}
#' }
#' @return Returns an object of class "ggplot2".
#' @seealso \code{\link{tournament}} to run a discharge rating curve tournament and \code{\link{summary.tournament}} for summaries.
#' @examples
#' \donttest{
#' library(ggplot2)
#' data(krokfors)
#' set.seed(1)
#' t_obj <- tournament(formula = Q ~ W, data = krokfors, num_cores = 2)
#' autoplot(t_obj)
#' }
#' @importFrom ggplot2 ggplot geom_boxplot stat_boxplot geom_line geom_point xlab ylab
#' @importFrom rlang .data
#' @export
autoplot.tournament <- function(object, type = 'boxplot', ...){
    legal_types <- c('boxplot')
    if(!(type %in% legal_types)){
        stop(paste('Type argument not recognized. Possible types are:\n - ', paste(legal_types, collapse = '\n - ')))
    }else if(type == "boxplot"){
        p <- plot_tournament_fun(object, type = type)
    }
    return(p)
}


#' Plot method for a discharge rating curve tournament
#'
#' Compare the four models from the tournament object in multiple ways
#'
#' @param x An object of class "tournament"
#' @param type A character denoting what type of plot should be drawn. Possible types are:
#' \describe{
#'   \item{\code{boxplot}}{Creates a boxplot of the posterior log-likelihood values, on the deviance scale.}
#'   \item{\code{rating_curve}}{Plots the rating curve.}
#'   \item{\code{rating_curve_mean}}{Plots the posterior mean of the rating curve.}
#'   \item{\code{f}}{Plots the power-law exponent.}
#'   \item{\code{sigma_eps}}{Plots the standard deviation on the data level.}
#'   \item{\code{residuals}}{Plots the log residuals.}
#'   \item{\code{tournament_results}}{Plots a diagram showing the tournament results.}
#' }
#' @param transformed A logical value indicating whether the quantity should be plotted on a transformed scale used during the Bayesian inference. Defaults to FALSE.
#' @param ... Not used in this function
#' @return No return value, called for side effects
#' @seealso \code{\link{tournament}} to run a discharge rating curve tournament and \code{\link{summary.tournament}} for summaries.
#' @examples
#' \donttest{
#' data(krokfors)
#' set.seed(1)
#' t_obj <- tournament(formula = Q ~ W, data = krokfors, num_cores = 2)
#' plot(t_obj)
#' plot(t_obj, transformed = TRUE)
#' plot(t_obj, type = 'boxplot')
#' plot(t_obj, type = 'f')
#' plot(t_obj, type = 'sigma_eps')
#' plot(t_obj, type = 'residuals')
#' plot(t_obj, type = 'tournament_results')
#' }
#' @importFrom grid grid.draw
#' @importFrom gridExtra grid.arrange
#' @export
plot.tournament <- function(x, type = 'tournament_results', transformed = FALSE, ...){
    legal_types <- c("boxplot", "tournament_results", "rating_curve", "rating_curve_mean", "sigma_eps", "f", "residuals", "convergence_diagnostics", "panel", "tournament_results")
    error_msg <- paste0('Type not recognized. Possible types are:', paste(legal_types, collapse = '\n - '))
    if( is.null(type) ){
        stop(error_msg)
    }else if(type == 'boxplot'){
        p <- autoplot(x, type = type)
    }else if(type %in% legal_types){
        p <- plot_tournament_grob(x, type = type, transformed = transformed, ...)
    }else{
        stop(error_msg)
    }
    if('ggplot' %in% class(p)){
        print(p)
    }else{
        if(type == 'panel'){
            grid.draw(p[[class(x$winner)]])
        }else{
            grid.draw(p)
        }
    }
}

#' Internal function to generate a summary output for a discharge rating curve tournament
#'
#' This function takes the summary results of a tournament object and produces a formatted
#' console output displaying the results of model comparisons. It supports different
#' model selection criteria: WAIC, DIC, and PMP.
#'
#' @param results A data.frame containing the summary results of a tournament.
#'   The structure of this data.frame determines which model-selection criterion
#'   was used (WAIC, DIC, or PMP).
#' @param method A string indicating the method used for model comparison ("WAIC", "DIC", or "PMP").
#' @param winning_criteria The criteria used to determine the winning model.
#'
#' @details
#' The function automatically detects the model-selection criterion used based on
#' the columns present in the input data frame. It then formats and prints a
#' summary of the tournament results, including the overall winner and detailed
#' results for each comparison.
#'
#' @keywords internal
#'
#' @return This function does not return a value; it prints the formatted summary
#' to the console.
tournament_summary_output <- function(results, method, winning_criteria) {
    if(!(method %in% c("WAIC", "DIC", "PMP"))) stop("Unknown method for tournament model comparison. Methods allowed are: 'WAIC' (default), 'DIC', and 'PMP'.")
    # Determine the method used and set appropriate columns
    if ("WAIC" %in% names(results)) {
        criteria_cols <- c("lppd", "eff_num_param", "WAIC", "SE_WAIC", "Delta_WAIC", "SE_Delta_WAIC")
        line_length <- 98
    } else if ("DIC" %in% names(results)) {
        criteria_cols <- c("D_hat", "eff_num_param", "DIC", "Delta_DIC")
        line_length <- 74
    } else if ("log_marg_lik" %in% names(results)) {
        criteria_cols <- c("log_marg_lik", "PMP")
        line_length <- 55
    }

    # Round numeric columns
    results[criteria_cols] <- lapply(results[criteria_cols], function(x) round(x, digits = 4))

    # Find the overall winner (winner of last comparison)
    overall_winner <- results$model[results$round == max(results$round) & results$winner]

    # Print summary header
    cat("\n=== Tournament Model Comparison Summary ===\n\n")
    cat("Method:", method, "\n")
    if (method == "PMP") {
        cat("Winning Criteria: PMP of more complex model >", winning_criteria, "\n")
    } else if (method == "DIC") {
        cat("Winning Criteria: Delta_DIC >", winning_criteria, "\n")
    } else if (method == "WAIC") {
        if(is.numeric(winning_criteria)) {
            cat("Winning Criteria: Delta_WAIC >", winning_criteria, "\n")
        } else {
            cat("Winning Criteria:", winning_criteria, "\n")
        }
    }
    cat("Overall Winner:", overall_winner, "\n\n")

    # Function to print a single comparison's results
    print_comparison <- function(comparison_data) {
        cat("Comparison", unique(comparison_data$comparison), "Results:\n")
        cat(paste(rep("-", line_length), collapse = ""), "\n")

        # Print column headers based on the method
        if (method == "WAIC") {
            cat(sprintf("%-12s %-6s %-6s %-8s %-14s %-10s %-10s %-10s %-12s\n",
                        "complexity", "model", "winner", "lppd", "eff_num_param", "WAIC", "SE_WAIC", "Delta_WAIC", "SE_Delta_WAIC"))
        } else if (method == "DIC") {
            cat(sprintf("%-12s %-6s %-6s %-10s %-14s %-10s %-10s\n",
                        "complexity", "model", "winner", "D_hat", "eff_num_param", "DIC", "Delta_DIC"))
        } else if (method == "PMP") {
            cat(sprintf("%-12s %-6s %-6s %-20s %-15s\n",
                        "complexity", "model", "winner", "log_marg_lik", "PMP"))
        }

        # Print each row
        for (i in 1:nrow(comparison_data)) {
            row <- comparison_data[i, ]
            if (method == "WAIC") {
                cat(sprintf("%-12s %-6s %-6s %-8.4f %-14.4f %-10.4f %-10.4f %-10s %-12s\n",
                            row$complexity,
                            row$model,
                            ifelse(row$winner, "<---", ""),
                            row$lppd,
                            row$eff_num_param,
                            row$WAIC,
                            row$SE_WAIC,
                            ifelse(is.na(row$Delta_WAIC), "", sprintf("%.4f", row$Delta_WAIC)),
                            ifelse(is.na(row$SE_Delta_WAIC), "", sprintf("%.4f", row$SE_Delta_WAIC))))
            } else if (method == "DIC") {
                cat(sprintf("%-12s %-6s %-6s %-10.4f %-14.4f %-10.4f %-10s\n",
                            row$complexity,
                            row$model,
                            ifelse(row$winner, "<---", ""),
                            row$D_hat,
                            row$eff_num_param,
                            row$DIC,
                            ifelse(is.na(row$Delta_DIC), "", sprintf("%.4f", row$Delta_DIC))))
            } else if (method == "PMP") {
                cat(sprintf("%-12s %-6s %-6s %-20.4f %-15.4f\n",
                            row$complexity,
                            row$model,
                            ifelse(row$winner, "<---", ""),
                            row$log_marg_lik,
                            row$PMP))
            }
        }
        cat("\n")
    }

    # Print results for each comparison
    for (comparison in unique(results$comparison)) {
        print_comparison(results[results$comparison == comparison, ])
    }

    # Print footer
    cat("=== End of Summary ===\n")
}

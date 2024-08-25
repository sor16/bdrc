#' Compare two models using a specified model-selection criteria
#'
#' evaluate_comparison uses the Widely Applicable Information Criterion (WAIC), the Deviance Information Criterion (DIC), or the posterior model probabilities (PMP), calculated with Bayes factor, to determine whether one model is more appropriate than the other given the data at hand.
#'
#' @param m A list of two model objects fit on the same dataset. The allowed model objects are "gplm", "gplm0", "plm" and "plm0"
#' @param method A string specifying the method used to estimate the predictive performance of the models. The allowed methods are "WAIC", "DIC" and "PMP".
#' @param winning_criteria A numerical value which sets the threshold which the first model in the list must exceed for it to be declared the more appropriate model. This value defaults to 2 for methods "WAIC" and "DIC", but defaults to 0.75 for method "PMP".
#' @return
#' A data.frame with the summary of the results of each comparison
#' @references Hrafnkelsson, B., Sigurdarson, H., Rögnvaldsson, S., Jansson, A. Ö., Vias, R. D., and Gardarsson, S. M. (2022). Generalization of the power-law rating curve using hydrodynamic theory and Bayesian hierarchical modeling, Environmetrics, 33(2):e2711. doi: https://doi.org/10.1002/env.2711
#'
#' @seealso \code{\link{tournament}}
#' @keywords internal
evaluate_comparison <- function(m, method, winning_criteria){

    if(method == "PMP"){

        ml <- sapply(m, function(x) exp(log_ml_harmonic_mean_est(x, x$d)))
        PR_m1 <- post_model_prob_m1(m[[2]], m[[1]], m[[1]]$data)
        df <- data.frame("marg_lik" = ml, "Post_prob" = c(PR_m1, 1 - PR_m1))

    }else if(method == "DIC"){

        D_hat <- sapply(m, function(x) x$D_hat)
        eff_num_param_DIC <- sapply(m, function(x) x$effective_num_param_DIC)
        DIC_vec <- sapply(m, function(x) x$DIC)
        DDIC <- DIC_vec[2] - DIC_vec[1]
        df <- data.frame("D_hat" = D_hat,
                         "eff_num_param" = eff_num_param_DIC,
                         "DIC" = DIC_vec,
                         "Delta_DIC" = c(DDIC, NA))

    }else if(method == "WAIC"){

        WAIC_i <- lapply(m, function(x) calc_waic(x, x$d)$waic_i)
        SE_WAIC <- sapply(WAIC_i, function(x) sqrt(length(x) * var(x)))
        SE_DWAIC <- SE_Delta_WAIC(m[[1]], m[[2]])
        lppd <- sapply(m, function(x) x$lppd)
        eff_num_param_WAIC <- sapply(m, function(x) x$effective_num_param_WAIC)
        WAIC_vec <- sapply(m, function(x) x$WAIC)
        DWAIC <- WAIC_vec[2] - WAIC_vec[1]
        df <- data.frame("lppd" = lppd,
                         "eff_num_param" = eff_num_param_WAIC,
                         "WAIC" = WAIC_vec,
                         "SE_WAIC" = SE_WAIC,
                         "Delta_WAIC" = c(DWAIC, NA),
                         "SE_Delta_WAIC" = c(SE_DWAIC, NA))
    }

    winner <- if(df[1, ncol(df)] >= winning_criteria) 1 else 2
    data.frame(model = sapply(m, class), df, winner = (1:2 == winner))

}


#' Tournament - Model comparison
#'
#' tournament compares four rating curve models of different complexities and determines the model that provides the best fit of the data at hand.
#'
#' @param formula An object of class "formula", with discharge column name as response and stage column name as a covariate.
#' @param data A data.frame containing the variables specified in formula.
#' @param model_list A list of exactly four model objects of types "plm0","plm","gplm0" and "gplm" to be used in the tournament. Note that all of the model objects are required to be run with the same data and same c_param.
#' @param method A string specifying the method used to estimate the predictive performance of the models. The allowed methods are "WAIC", "DIC" and "PMP".
#' @param winning_criteria A numerical value which sets a threshold the more complex model in each model comparison must exceed to be deemed the more appropriate model. See the Details section.
#' @param verbose A logical value indicating whether to print progress and diagnostic information. If `TRUE`, the function will print messages as it runs. If `FALSE`, the function will run silently. Default is `TRUE`.
#' @param ... Optional arguments passed to the model functions.
#' @details Tournament is a model comparison method that uses WAIC (default method) to estimate the expected prediction error of the four models and select the most appropriate model given the data. The first round of model comparisons sets up model types, "gplm" vs. "gplm0" and "plm" vs. "plm0". The two comparisons are conducted such that if the WAIC of the more complex model ("gplm" and "plm", respectively) is smaller than the WAIC of the simpler models ("gplm0" and "plm0", respectively) by an input argument called the \code{winning_criteria} (default value = 2), then it is chosen as the more appropriate model. If not, the simpler model is chosen. The more appropriate models move on to the second round and are compared in the same way. The winner of the second round is chosen as the overall tournament winner and deemed the most appropriate model given the data.
#'
#' The default method "WAIC", or the Widely Applicable Information Criterion (see Watanabe (2010)), is used to estimate the predictive performance of the models. This method is a fully Bayesian method that uses the full set of posterior draws to estimate of the expected log pointwise predictive density.
#'
#' Method "DIC", or Deviance Information Criterion (see Spiegelhalter (2002)), is similar to the "WAIC" but instead of using the full set of posterior draws to compute the estimate of the expected log pointwise predictive density, it uses a point estimate of the posterior distribution.
#'
#' Method "PMP" uses the posterior model probabilities, calculated with Bayes factor (see Jeffreys (1961) and Kass and Raftery (1995)), to compare the models, where all the models are assumed a priori to be equally likely. This method is not chosen as the default method because the Bayes factor calculations can be quite unstable.
#'
#' When methods "WAIC" or "DIC" are used, the \code{winning_criteria} should be a real number. The winning criteria is a threshold value which the more complex model in each model comparison must exceed for it to be declared the more appropriate model. Setting the winning criteria slightly above 0 (default value = 2 for both "WAIC" and "DIC") gives the less complex model in each comparison a slight advantage. When method "PMP" is used, the winning criteria should be a real value between 0 and 1 (default value = 0.75). This sets the threshold value for which the posterior probability of the more complex model, given the data, in each model comparison must exceed for it to be declared the more appropriate model. In all three cases, the default value is selected so as to give the less complex models a slight advantage, and should give more or less consistent results when applying the tournament to real world data.
#' @return
#' An object of type "tournament" with the following elements:
#' \describe{
#'   \item{\code{contestants}}{The model objects of types "plm0", "plm", "gplm0" and "gplm" being compared.}
#'   \item{\code{winner}}{The model object of the tournament winner.}
#'   \item{\code{info}}{The specifics about the tournament; the overall winner; the method used; and the winning criteria.}
#'   \item{\code{summary}}{A data frame with information on results of the different comparisons in the power-law tournament. The contents of this data frame depend on the method used:
#'     \itemize{
#'       \item For method "WAIC":
#'         \itemize{
#'           \item round: The tournament round
#'           \item comparison: The comparison number
#'           \item model: The model type
#'           \item lppd: Log pointwise predictive density
#'           \item eff_num_param: Effective number of parameters (WAIC)
#'           \item WAIC: Widely Applicable Information Criterion
#'           \item SE_WAIC: Standard error of WAIC
#'           \item Delta_WAIC: Difference in WAIC
#'           \item SE_Delta_WAIC: Standard error of the difference in WAIC
#'           \item winner: Logical value indicating if the model was selected in the corresponding comparison
#'         }
#'       \item For method "DIC":
#'         \itemize{
#'           \item round: The tournament round
#'           \item comparison: The comparison number
#'           \item model: The model type
#'           \item D_hat: Minus two times the log-likelihood evaluated at the median of the posterior samples
#'           \item p_D: Effective number of parameters (DIC)
#'           \item DIC: Deviance Information Criterion
#'           \item Delta_DIC: Difference in DIC
#'           \item winner: Logical value indicating if the model was selected in the corresponding comparison
#'         }
#'       \item For method "PMP":
#'         \itemize{
#'           \item round: The tournament round
#'           \item comparison: The comparison number
#'           \item model: The model type
#'           \item log_marg_lik: Log marginal likelihood estimated with the harmonic mean
#'           \item post_model_prob: Posterior model probability computed with Bayes factor
#'           \item winner: Logical value indicating if the model was selected in the corresponding comparison
#'         }
#'     }
#'   }
#' }
#'
#' @references Hrafnkelsson, B., Sigurdarson, H., Rögnvaldsson, S., Jansson, A. Ö., Vias, R. D., and Gardarsson, S. M. (2022). Generalization of the power-law rating curve using hydrodynamic theory and Bayesian hierarchical modeling, Environmetrics, 33(2):e2711. doi: https://doi.org/10.1002/env.2711
#' @references Jeffreys, H. (1961). Theory of Probability, Third Edition. Oxford University Press.
#' @references Kass, R., and A. Raftery, A. (1995). Bayes Factors. Journal of the American Statistical Association, 90, 773-795. doi: https://doi.org/10.1080/01621459.1995.10476572
#' @references Spiegelhalter, D., Best, N., Carlin, B., Van Der Linde, A. (2002). Bayesian measures of model complexity and fit. Journal of the Royal Statistical Society: Series B (Statistical Methodology) 64(4), 583–639. doi: https://doi.org/10.1111/1467-9868.00353
#' @references Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and widely applicable information criterion in singular learning theory. J. Mach. Learn. Res. 11, 3571–3594.
#'
#' @seealso \code{\link{plm0}} \code{\link{plm}}, \code{\link{gplm0}},\code{\link{gplm}} \code{\link{summary.tournament}} and \code{\link{plot.tournament}}
#' @examples
#' \donttest{
#' data(krokfors)
#' set.seed(1)
#' t_obj <- tournament(formula = Q ~ W, data = krokfors, num_cores = 2)
#' t_obj
#' summary(t_obj)
#' }
#' @export
tournament <- function(formula = NULL, data = NULL, model_list = NULL, method = 'WAIC', winning_criteria = NULL, verbose = TRUE, ...) {
    default_win_crit <- c('WAIC' = 2, 'DIC' = 2, 'PMP' = 0.75)
    error_msg <- "The method input must contain a string indicating the method to be used for comparing the models. The methods are 'WAIC' (default), 'DIC' and 'PMP'."
    if(is.null(method)){
        stop(error_msg)
    }else{
        if(length(method) > 1){
            stop(error_msg)
        }
        if(!(method%in%c('WAIC', 'DIC', 'PMP'))){
            stop(error_msg)
        }
    }
    if(grepl("IC", method)){
        error_msg <- "The winning_criteria when the method is set to 'WAIC' or 'DIC' must be a numerical real value. This is the threshold value which a model comparison statistic of the more complex model must surpass to be declared the more appropriate model when compared with a less complex model. This value defaults to 2 when the models are compared by their 'WAIC' or 'DIC' values."
    }else{
        error_msg <- "The winning_criteria when the method is set to 'PMP' must be a numerical value between 0 and 1. This is the threshold for which the posterior probability of a more complex model, which is calculated using the Bayes factor, needs to surpass to be declared the appropriate model when compared with a less complex model. This value defaults to 0.75 when the method is set to 'PMP'."
    }
    if(!is.null(winning_criteria)){
        if(length(method) > 1){
            stop(error_msg)
        }
        if(!inherits(winning_criteria, "numeric")){
            stop(error_msg)
        }else if( method == 'PMP' & abs(winning_criteria - 0.5) > 0.5){
            stop(error_msg)
        }
    }
    if(is.null(winning_criteria)){
        winning_criteria <- default_win_crit[method]
    }
    error_msg <- 'Please provide; formula and data (name arguments explicitly) or model_list with four model objects of types gplm, gplm0, plm and plm0.'
    models <- list()
    if(!is.null(model_list) | (is.null(model_list) & inherits(formula, 'list'))){
        if(is.null(model_list) & inherits(formula, 'list')){
            model_list = formula
        }
        if(length(model_list) != 4){
            stop(error_msg)
        }
        models_class <- unlist(lapply(model_list, class))
        if(!all(sort(models_class) == c('gplm', 'gplm0', 'plm', 'plm0'))){
            stop(error_msg)
        }
        models <- model_list
        names(models) <- models_class
        models <- models[order(names(models))]
        #TODO: make sure data argument matches models in model_list
        if(length(unique(lapply(models, function(x) x$data))) != 1){
            stop('The four models added have to be fit on the same data set')
        }
        if(length(unique(lapply(models, function(x) x$c_param))) != 1){
            stop('The four models added have to be fit either all with the same stage of zero discharge (c), or all with unknown c')
        }
    }else{
        if(verbose) cat('Running tournament [                                                ] 0%\n\n')
        models$gplm <- gplm(formula = formula, data = data, verbose = verbose, ...)
        if(verbose) cat('\n- gplm finished  [============                                    ] 25%\n\n')
        models$gplm0 <- gplm0(formula = formula, data = data, verbose = verbose, ...)
        if(verbose) cat('\n- gplm0 finished [========================                        ] 50%\n\n')
        models$plm <- plm(formula = formula, data = data, verbose = verbose, ...)
        if(verbose) cat('\n- plm finished   [====================================            ] 75%\n\n')
        models$plm0 <- plm0(formula = formula, data = data, verbose = verbose, ...)
        if(verbose) cat('\n- plm0 finished  [================================================] 100%\n')
    }
    round1 <- list(list(models$gplm, models$gplm0), list(models$plm, models$plm0))
    round1_res <- lapply(1:length(round1), function(i){
        comparison_df <- evaluate_comparison(round1[[i]], method, winning_criteria)
        round_df <- data.frame(round = 1, comparison = i)
        cbind(round_df, comparison_df)
    })
    round1_res <- do.call(rbind, round1_res)
    round1_winners <- round1_res$model[round1_res$winner]

    round2 <- lapply(1:length(round1), function(i){
        round1[[i]][[which(round1_res$winner[round1_res$comparison == i])]]
    })
    round2_res <- cbind(data.frame(round = 2, comparison = 3), evaluate_comparison(round2, method, winning_criteria))
    round2_winner <- round2_res$model[round2_res$winner]
    out_obj <- list()
    attr(out_obj, "class") <- "tournament"
    out_obj$contestants <- models
    out_obj$winner <- round2[[which(round2_res$winner)]]
    out_obj$summary <- rbind(round1_res, round2_res)
    out_obj$info <- list("winner" = class(round2[[which(round2_res$winner)]]), "method" = method, "winning_criteria" = winning_criteria)
    if(method == 'PMP' & verbose) cat( '\u26A0 Warning: The Harmonic Mean Estimator (HME) is used to estimate the Bayes Factor for the posterior model probability (PMP), which is known to be unstable and potentially unreliable. We recommend using method "WAIC" (Widely Applicable Information Criterion) for model comparison instead.' )
    if(verbose) cat(sprintf("\nTournament winner: %s\n", class(out_obj$winner)))
    return(out_obj)
}



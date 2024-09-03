#' Compare two models using a specified model-selection criteria
#'
#' evaluate_comparison uses the Widely Applicable Information Criterion (WAIC), the Deviance Information Criterion (DIC), or the posterior model probabilities (PMP), calculated with Bayes factor, to determine whether one model is more appropriate than the other given the data at hand.
#'
#' @param m A list of two model objects fit on the same dataset. The allowed model objects are "gplm", "gplm0", "plm" and "plm0"
#' @param method A string specifying the method used to estimate the predictive performance of the models. The allowed methods are "WAIC", "DIC" and "PMP".
#' @param winning_criteria For "WAIC", it can be either a numeric value or a string expression. For "DIC", it must be a numeric value. For "PMP", it must be a numeric value between 0 and 1. This sets the threshold for determining the more appropriate model. See Details for more information.
#' @return
#' A data.frame with the summary of the results of each comparison, including:
#' \itemize{
#'   \item complexity: Indicates whether a model is the "more" or "less" complex model in a comparison
#'   \item model: The type of model (gplm, gplm0, plm, or plm0)
#'   \item Method-specific columns (see Details)
#'   \item winner: Logical value indicating if the model was selected
#' }
#' @details
#' For "WAIC" method:
#' \itemize{
#'   \item If winning_criteria is numeric, the more complex model wins if Delta_WAIC > winning_criteria
#'   \item If winning_criteria is a string, it must be a valid R expression using Delta_WAIC and/or SE_Delta_WAIC
#'   \item Returns columns: lppd, eff_num_param, WAIC, SE_WAIC, Delta_WAIC, SE_Delta_WAIC
#' }
#' For "DIC" method:
#' \itemize{
#'   \item winning_criteria must be numeric
#'   \item The more complex model wins if Delta_DIC > winning_criteria
#'   \item Returns columns: D_hat, eff_num_param, DIC, Delta_DIC
#' }
#' For "PMP" method:
#' \itemize{
#'   \item winning_criteria must be a numeric value between 0 and 1
#'   \item The more complex model wins if its PMP > winning_criteria
#'   \item Returns columns: log_marg_lik, PMP
#' }
#' @references Hrafnkelsson, B., Sigurdarson, H., Rögnvaldsson, S., Jansson, A. Ö., Vias, R. D., and Gardarsson, S. M. (2022). Generalization of the power-law rating curve using hydrodynamic theory and Bayesian hierarchical modeling, Environmetrics, 33(2):e2711. doi: https://doi.org/10.1002/env.2711
#'
#' @seealso \code{\link{tournament}}
#' @keywords internal
evaluate_comparison <- function(m, method, winning_criteria){
    if(method == "PMP"){
        log_ml <- sapply(m, function(x) log_ml_harmonic_mean_est(x, x$d))
        PR_m1 <- post_model_prob_m1(m[[2]], m[[1]], m[[1]]$data)
        df <- data.frame("complexity" = c("more", "less"),
                         "log_marg_lik" = log_ml,
                         "PMP" = c(PR_m1, 1 - PR_m1))
        winner <- PR_m1 > winning_criteria
    } else if(method == "DIC"){
        D_hat <- sapply(m, function(x) x$D_hat)
        eff_num_param_DIC <- sapply(m, function(x) x$effective_num_param_DIC)
        DIC_vec <- sapply(m, function(x) x$DIC)
        DDIC <- DIC_vec[2] - DIC_vec[1]
        df <- data.frame("complexity" = c("more", "less"),
                         "D_hat" = D_hat,
                         "eff_num_param" = eff_num_param_DIC,
                         "DIC" = DIC_vec,
                         "Delta_DIC" = c(DDIC, NA))
        winner <- DDIC > winning_criteria
    } else if(method == "WAIC"){
        WAIC_i <- lapply(m, function(x) calc_waic(x, x$d)$waic_i)
        SE_WAIC <- sapply(WAIC_i, function(x) sqrt(length(x) * var(x)))
        SE_DWAIC <- SE_Delta_WAIC(m[[1]], m[[2]])
        lppd <- sapply(m, function(x) x$lppd)
        eff_num_param_WAIC <- sapply(m, function(x) x$effective_num_param_WAIC)
        WAIC_vec <- sapply(m, function(x) x$WAIC)
        DWAIC <- WAIC_vec[2] - WAIC_vec[1]
        df <- data.frame("complexity" = c("more", "less"),
                         "lppd" = lppd,
                         "eff_num_param" = eff_num_param_WAIC,
                         "WAIC" = WAIC_vec,
                         "SE_WAIC" = SE_WAIC,
                         "Delta_WAIC" = c(DWAIC, NA),
                         "SE_Delta_WAIC" = c(SE_DWAIC, NA))
        if(is.numeric(winning_criteria)){
            winner <- DWAIC > winning_criteria
        } else {
            winner <- with(df[1,], eval(parse(text = winning_criteria)))
        }
    } else {
        stop("Unknown method")
    }

    return(data.frame(complexity = df$complexity, model = sapply(m, class), df[,-1], winner = c(winner, !winner)))
}


#' Tournament - Model comparison
#'
#' tournament compares four rating curve models of different complexities and determines the model that provides the best fit of the data at hand.
#'
#' @param formula An object of class "formula", with discharge column name as response and stage column name as a covariate.
#' @param data A data.frame containing the variables specified in formula.
#' @param model_list A list of exactly four model objects of types "plm0","plm","gplm0" and "gplm" to be used in the tournament. Note that all of the model objects are required to be run with the same data and same c_param.
#' @param method A string specifying the method used to estimate the predictive performance of the models. The allowed methods are "WAIC", "DIC" and "PMP".
#' @param winning_criteria Specifies the criteria for model selection. For "WAIC", it can be a numeric value or a string expression. For "DIC", it must be a numeric value. For "PMP", it must be a numeric value between 0 and 1. See Details section.
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
#' When method "WAIC" is used, the \code{winning_criteria} can be either a numeric value or a string expression. If numeric, it sets the threshold which the more complex model must exceed to be declared the more appropriate model. If a string, it must be a valid R expression using Delta_WAIC and/or SE_Delta_WAIC (e.g., "Delta_WAIC > 2 & Delta_WAIC - SE_Delta_WAIC > 0"). For method "DIC", \code{winning_criteria} must be a numeric value. For method "PMP", the winning criteria should be a numeric value between 0 and 1 (default value = 0.75). This sets the threshold value for which the posterior probability of the more complex model, given the data, in each model comparison must exceed to be declared the more appropriate model. In all cases, the default values are selected to give the less complex models a slight advantage, which should give more or less consistent results when applying the tournament to real world data.
#' @return
#' An object of type "tournament" with the following elements:
#' \describe{
#'   \item{\code{contestants}}{The model objects of types "plm0", "plm", "gplm0" and "gplm" being compared.}
#'   \item{\code{winner}}{The model object of the tournament winner.}
#'   \item{\code{info}}{The specifics about the tournament; the overall winner; the method used; and the winning criteria.}
#'   \item{\code{summary}}{A data frame with information on results of the different comparisons in the power-law tournament. The contents of this data frame depend on the method used:
#'     \itemize{
#'       \item For all methods:
#'         \itemize{
#'           \item round: The tournament round
#'           \item comparison: The comparison number
#'           \item complexity: Indicates whether a model is the "more" or "less" complex model in a comparison
#'           \item model: The model type
#'           \item winner: Logical value indicating if the model was selected in the corresponding comparison
#'         }
#'       \item Additional columns for method "WAIC":
#'         \itemize{
#'           \item lppd: Log pointwise predictive density
#'           \item eff_num_param: Effective number of parameters (WAIC)
#'           \item WAIC: Widely Applicable Information Criterion
#'           \item SE_WAIC: Standard error of WAIC
#'           \item Delta_WAIC: Difference in WAIC
#'           \item SE_Delta_WAIC: Standard error of the difference in WAIC
#'         }
#'       \item Additional columns for method "DIC":
#'         \itemize{
#'           \item D_hat: Minus two times the log-likelihood evaluated at the median of the posterior samples
#'           \item eff_num_param: Effective number of parameters (DIC)
#'           \item DIC: Deviance Information Criterion
#'           \item Delta_DIC: Difference in DIC
#'         }
#'       \item Additional columns for method "PMP":
#'         \itemize{
#'           \item log_marg_lik: Logarithm of the marginal likelihood estimated, computed with the harmonic-mean estimator
#'           \item PMP: Posterior model probability computed with Bayes factor
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
#'
#' # Using different methods and winning criteria
#' t_obj_dic <- tournament(Q ~ W,
#'                         krokfors,
#'                         num_cores = 2,
#'                         method = "DIC",
#'                         winning_criteria = 3)
#' t_obj_pmp <- tournament(Q ~ W,
#'                         krokfors,
#'                         num_cores = 2,
#'                         method = "PMP",
#'                         winning_criteria = 0.8)
#' t_obj_waic_expr <- tournament(Q ~ W,
#'                               krokfors,
#'                               num_cores = 2,
#'                               winning_criteria = "Delta_WAIC > 2 & Delta_WAIC - SE_Delta_WAIC > 0")
#' }
#' @export
tournament <- function(formula = NULL, data = NULL, model_list = NULL, method = 'WAIC', winning_criteria = NULL, verbose = TRUE, ...) {
    # Set default winning_criteria
    default_win_crit <- list(
        'WAIC' = "Delta_WAIC > 2",
        'DIC' = 2,
        'PMP' = 0.75
    )

    # Send warning if an input is not recognized and is ignored
    all_args <- list(...)
    expected_args <- c("model_list", "method", "winning_criteria", "formula", "data", "c_param", "h_max", "parallel", "num_cores", "forcepoint")
    unexpected_args <- setdiff(names(all_args), expected_args)
    if (length(unexpected_args) > 0) {
        stop("The following argument(s) are not recognized: ",paste(unexpected_args, collapse = ", "))
    }

    # Check for erroneous "method" and "winning_criteria" inputs
    error_msg <- "The method input must contain a string indicating the method to be used for comparing the models. The methods are 'WAIC' (default), 'DIC' and 'PMP'."
    if(is.null(method)){
        stop(error_msg)
    } else if(length(method) > 1 || !(method %in% c('WAIC', 'DIC', 'PMP'))){
        stop(error_msg)
    }
    if(is.null(winning_criteria)){
        winning_criteria <- default_win_crit[[method]]
    } else {
        if(method == 'PMP') {
            if(!is.numeric(winning_criteria) || length(winning_criteria) != 1 ||
               is.nan(winning_criteria) || is.na(winning_criteria) ||
               is.infinite(winning_criteria) || winning_criteria < 0 || winning_criteria > 1) {
                stop("For method 'PMP', winning_criteria must be a single number between 0 and 1. Default is 0.75, meaning the more complex model needs >75% posterior probability to be chosen. See ?tournament for details.")
            }
        } else if(method == 'WAIC') {
            if(is.numeric(winning_criteria)) {
                if(length(winning_criteria) != 1 || is.nan(winning_criteria) ||
                   is.na(winning_criteria) || is.infinite(winning_criteria)) {
                    stop("For method 'WAIC', when numeric, winning_criteria must be a single finite number. Default is 2, equivalent to requiring 'Delta_WAIC > 2' for the more complex model to be chosen. See ?tournament for details on numeric and expression inputs.")
                }
                winning_criteria <- paste("Delta_WAIC >", winning_criteria)
            } else if(is.character(winning_criteria)) {
                allowed_vars <- c("Delta_WAIC", "SE_Delta_WAIC")
                expr <- try(parse(text = winning_criteria), silent = TRUE)
                cond1 <- grepl("NA",winning_criteria) | grepl("NULL",winning_criteria)
                cond2 <- !(grepl("Delta_WAIC", winning_criteria) | grepl("SE_Delta_WAIC", winning_criteria))
                cond3 <- inherits(expr, "try-error") || !all(all.vars(expr) %in% allowed_vars)
                if(cond1 | cond2 | cond3) {
                    stop("For method 'WAIC', when a string, winning_criteria must be a valid R expression using only Delta_WAIC and/or SE_Delta_WAIC. Example: 'Delta_WAIC > 2 & Delta_WAIC - SE_Delta_WAIC > 0'. See ?tournament for more examples and details.")
                }
                # Check for single '=' that is not part of '==', '<=', or '>='
                if(grepl("(?<![=<>])=(?!=)", winning_criteria, perl = TRUE)) {
                    stop("Use '==' for equality in winning_criteria expressions, not '='. See ?tournament for more examples and details.")
                }
            } else {
                stop("For method 'WAIC', winning_criteria must be either a numeric value or a string expression. See ?tournament for more examples and details.")
            }
        } else if(method == 'DIC') {
            if(!is.numeric(winning_criteria) || length(winning_criteria) != 1 ||
               is.nan(winning_criteria) || is.na(winning_criteria) ||
               is.infinite(winning_criteria)) {
                stop('For method "DIC", winning_criteria must be a single real number. Default is 2, equivalent to requiring "Delta_DIC > 2" for the more complex model to be chosen. See ?tournament for more details on choosing this value.')
            }
        }
    }

    # Either gather fitted models into list "models", or fit the models and input into the list.
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
        if(length(unique(lapply(models, function(x) x$data))) != 1){
            stop('The four models added have to be fit on the same data set')
        }
    }else{
        if(verbose)         cat('Running tournament  [                                                ] 0%\n\n')
        models$gplm <- gplm(formula = formula, data = data, verbose = verbose, ...)
        if(verbose) cat('\n \u2714  gplm finished   [============                                    ] 25%\n\n')
        models$gplm0 <- gplm0(formula = formula, data = data, verbose = verbose, ...)
        if(verbose) cat('\n \u2714  gplm0 finished  [========================                        ] 50%\n\n')
        models$plm <- plm(formula = formula, data = data, verbose = verbose, ...)
        if(verbose) cat('\n \u2714  plm finished    [====================================            ] 75%\n\n')
        models$plm0 <- plm0(formula = formula, data = data, verbose = verbose, ...)
        if(verbose) cat('\n \u2714  plm0 finished   [================================================] 100%\n')
    }

    # Compute the outcome of the tournament model comparison
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

    # Prepare the output object
    out_obj <- list()
    attr(out_obj, "class") <- "tournament"
    out_obj$contestants <- models
    out_obj$winner <- round2[[which(round2_res$winner)]]
    out_obj$summary <- rbind(round1_res, round2_res)
    out_obj$info <- list("winner" = class(round2[[which(round2_res$winner)]]),
                         "method" = method,
                         "winning_criteria" = winning_criteria)
    if(method == 'PMP' & verbose) cat( '\u26A0 Warning: The Harmonic Mean Estimator (HME) is used to estimate the Bayes Factor for the posterior model probability (PMP), which is known to be unstable and potentially unreliable. We recommend using method "WAIC" (Widely Applicable Information Criterion) for model comparison instead.\n' )
    return(out_obj)
}

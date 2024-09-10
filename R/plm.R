#' Power-law model with variance that varies with stage.
#'
#' plm is used to fit a discharge rating curve for paired measurements of stage and discharge using a power-law model with variance that varies with stage as described in Hrafnkelsson et al. (2022). See "Details" for a more elaborate description of the model.
#'
#' @param formula An object of class "formula", with discharge column name as response and stage column name as a covariate, i.e. of the form \code{y}~\code{x} where \code{y} is discharge in m\eqn{^3/}s and \code{x} is stage in m (it is very important that the data is in the correct units).
#' @param data A data.frame containing the variables specified in formula.
#' @param c_param The largest stage value for which there is zero discharge. If NULL, it is treated as unknown in the model and inferred from the data.
#' @param h_max The maximum stage to which the rating curve should extrapolate to. If NULL, the maximum stage value in the data is selected as an upper bound.
#' @param parallel A logical value indicating whether to run the MCMC in parallel or not. Defaults to TRUE.
#' @param num_cores An integer between 1 and 4 (number of MCMC chains) indicating how many cores to use. Only used if parallel=TRUE. If NULL, the number of cores available on the device is detected automatically.
#' @param forcepoint A logical vector of the same length as the number of rows in data. If an element at index \eqn{i} is TRUE it indicates that the rating curve should be forced through the \eqn{i}-th measurement. Use with care, as this will strongly influence the resulting rating curve.
#' @param verbose A logical value indicating whether to print progress and diagnostic information. If `TRUE`, the function will print messages as it runs. If `FALSE`, the function will run silently. Default is `TRUE`.
#'
#' @details The power-law model, which is commonly used in hydraulic practice, is of the form
#' \deqn{Q=a(h-c)^{b}}
#' where \eqn{Q} is discharge, \eqn{h} is stage and \eqn{a}, \eqn{b} and \eqn{c} are unknown constants.\cr\cr
#' The power-law model is here inferred by using a Bayesian hierarchical model. The model is on a logarithmic scale
#' \deqn{\log(Q_i) = \log(a) + b \log(h_i - c) + \varepsilon_i,     i = 1,...,n}
#' where \eqn{\varepsilon_i} follows a normal distribution with mean zero and variance \eqn{\sigma_\varepsilon(h_i)^2} that varies with stage. The error variance, \eqn{\sigma_\varepsilon^2(h)}, of the log-discharge data is modeled as an exponential of a B-spline curve, that is, a linear combination of six B-spline basis functions that are defined over the range of the stage observations. An efficient posterior simulation is achieved by sampling from the joint posterior density of the hyperparameters of the model, and then sampling from the density of the latent parameters conditional on the hyperparameters.\cr\cr
#' Bayesian inference is based on the posterior density and summary statistics such as the posterior mean and 95\% posterior intervals are based on the posterior density. Analytical formulas for these summary statistics are intractable in most cases and thus they are computed by generating samples from the posterior density using a Markov chain Monte Carlo simulation.
#' @return plm returns an object of class "plm". An object of class "plm" is a list containing the following components:
#' \describe{
#'   \item{\code{rating_curve}}{A data frame with 2.5\%, 50\% and 97.5\% percentiles of the posterior predictive distribution of the rating curve.}
#'   \item{\code{rating_curve_mean}}{A data frame with 2.5\%, 50\% and 97.5\% percentiles of the posterior distribution of the mean of the rating curve. Additionally contains columns with r_hat and the effective number of samples for each parameter as defined in Gelman et al. (2013).}
#'   \item{\code{param_summary}}{A data frame with 2.5\%, 50\% and 97.5\% percentiles of the posterior distribution of latent- and hyperparameters.}
#'   \item{\code{sigma_eps_summary}}{A data frame with 2.5\%, 50\% and 97.5\% percentiles of the posterior of \eqn{\sigma_{\varepsilon}}.}
#'   \item{\code{posterior_log_likelihood_summary}}{A data frame with 2.5\%, 50\% and 97.5\% percentiles of the posterior log-likelihood values.}
#'   \item{\code{rating_curve_posterior}}{A matrix containing the full thinned posterior samples of the posterior predictive distribution of the rating curve (excluding burn-in).}
#'   \item{\code{rating_curve_mean_posterior}}{A matrix containing the full thinned posterior samples of the posterior distribution of the mean of the rating curve (excluding burn-in).}
#'   \item{\code{a_posterior}}{A numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{a}.}
#'   \item{\code{b_posterior}}{A numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{b}.}
#'   \item{\code{c_posterior}}{A numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{c}.}
#'   \item{\code{sigma_eps_posterior}}{A numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{\sigma_{\varepsilon}}.}
#'   \item{\code{eta_1_posterior}}{A numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{\eta_1}.}
#'   \item{\code{eta_2_posterior}}{A numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{\eta_2}.}
#'   \item{\code{eta_3_posterior}}{A numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{\eta_3}.}
#'   \item{\code{eta_4_posterior}}{A numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{\eta_4}.}
#'   \item{\code{eta_5_posterior}}{A numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{\eta_5}.}
#'   \item{\code{eta_6_posterior}}{A numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{\eta_6}.}
#'   \item{\code{posterior_log_likelihood}}{A numeric vector containing the full thinned posterior log-likelihood values, excluding burn-in samples.}
#'   \item{\code{D_hat}}{A statistic defined as -2 times the log-likelihood evaluated at the median value of the parameters.}
#'   \item{\code{effective_num_param_DIC}}{The effective number of parameters, which is calculated as median(-2*posterior_log_likelihood) minus D_hat.}
#'   \item{\code{DIC}}{The Deviance Information Criterion for the model, calculated as D_hat plus 2*effective_num_parameters_DIC.}
#'   \item{\code{lppd}}{The log pointwise predictive density of the observed data under the model.}
#'   \item{\code{WAIC}}{The Widely Applicable Information Criterion for the model, defined as -2*( lppd - effective_num_param_WAIC ).}
#'   \item{\code{WAIC_i}}{The pointwise WAIC values, where WAIC := sum(WAIC_i).}
#'   \item{\code{effective_num_param_WAIC}}{The effective number of parameters, which is calculated by summing up the posterior variance of the log predictive density for each data point.}
#'   \item{\code{autocorrelation}}{A data frame with the autocorrelation of each parameter for different lags.}
#'   \item{\code{acceptance_rate}}{The proportion of accepted samples in the thinned MCMC chain (excluding burn-in).}
#'   \item{\code{formula}}{An object of type "formula" provided by the user.}
#'   \item{\code{data}}{The data provided by the user, ordered by stage.}
#'   \item{\code{run_info}}{The information about the input arguments and the specific parameters used in the MCMC chain.}
#' }
#' @references Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., and Rubin, D. B. (2013). Bayesian Data Analysis, Third Edition. Chapman & Hall/CRC Texts in Statistical Science. Taylor & Francis. doi: https://doi.org/10.1201/b16018
#' @references Hrafnkelsson, B., Sigurdarson, H., Rögnvaldsson, S., Jansson, A. Ö., Vias, R. D., and Gardarsson, S. M. (2022). Generalization of the power-law rating curve using hydrodynamic theory and Bayesian hierarchical modeling, Environmetrics, 33(2):e2711. doi: https://doi.org/10.1002/env.2711
#' @references Spiegelhalter, D., Best, N., Carlin, B., Van Der Linde, A. (2002). Bayesian measures of model complexity and fit. Journal of the Royal Statistical Society: Series B (Statistical Methodology) 64(4), 583–639. doi: https://doi.org/10.1111/1467-9868.00353
#' @references Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and widely applicable information criterion in singular learning theory. J. Mach. Learn. Res. 11, 3571–3594.
#' @seealso \code{\link{summary.plm}} for summaries, \code{\link{predict.plm}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.plm}} to help visualize the full posterior distributions.
#' @examples
#' \donttest{
#' data(spanga)
#' set.seed(1)
#' plm.fit <- plm(formula=Q~W,data=spanga,num_cores=2)
#' summary(plm.fit)
#' }
#' @export
plm <- function(formula, data, c_param = NULL, h_max = NULL, parallel = TRUE, num_cores = NULL, forcepoint = rep(FALSE, nrow(data)), verbose = TRUE){
    # argument checking
    check_arguments(formula, data, c_param, h_max, parallel, num_cores, forcepoint, verbose)

    # Parse the extended formula
    parsed_formula <- parse_extended_formula(formula)
    discharge_var <- parsed_formula$discharge
    error_var <- parsed_formula$discharge_error
    stage_var <- parsed_formula$stage

    # Check if variables exist in the data
    if (!discharge_var %in% names(data)) stop(paste("Discharge variable", discharge_var, "not found in data."))
    if (!is.null(error_var) && !(error_var %in% names(data))) stop(paste("Discharge measurment error variable", error_var, "not found in data."))
    if (!stage_var %in% names(data)) stop(paste("Water level (stage) variable", stage_var, "not found in data."))

    # Prepare data
    model_dat <- as.data.frame(data[, all.vars(formula)])
    forcepoint <- forcepoint[order(model_dat[, stage_var, drop = TRUE])]
    model_dat <- model_dat[order(model_dat[, stage_var, drop = TRUE]), ]
    Q <- model_dat[, discharge_var, drop = TRUE]
    h <- model_dat[, stage_var, drop = TRUE]

    # check data for errors
    if(dim(model_dat)[1] < 2) stop('At least two paired observations of stage and discharge are required to fit a rating curve')
    if(!is.null(c_param) && min(h) < c_param) stop('c_param must be lower than the minimum stage value in the data')
    if(any(Q <= 0)) stop('All discharge measurements must be strictly greater than zero. If you know the stage of zero discharge, use c_param.')

    # check if measurement errors are correct and create Q_me
    if(!is.null(error_var)){
        check_Q_me_for_errors(model_dat[, error_var, drop = TRUE])
        Q_me <- model_dat[, error_var, drop = TRUE]
    }else{
        Q_me <- NULL
    }

    # Bayesian inference
    MCMC_output_list <- plm.inference(y = log(Q), Q_me = Q_me, h = h, c_param = c_param, h_max = h_max, parallel = parallel, forcepoint = forcepoint, num_cores = num_cores, verbose = verbose)
    param_names <- get_param_names('plm', c_param)

    result_obj <- list(posterior = list(), log_likelihood = list(), summary = list(), diagnostics = list())
    attr(result_obj, "class") <- "plm"
    result_obj$posterior$a <- MCMC_output_list$x[1, ]
    result_obj$posterior$b <- MCMC_output_list$x[2, ]
    if(is.null(c_param)){
        result_obj$posterior$c <- MCMC_output_list$theta[1, ]
        result_obj$posterior$sigma_eta <- MCMC_output_list$theta[2, ]
        for(i in 3:nrow(MCMC_output_list$theta)){
            result_obj$posterior[[paste0('eta_', i - 2)]] <- MCMC_output_list$theta[i, ]
        }
    }else{
        result_obj$posterior$c <- NULL
        result_obj$posterior$sigma_eta <- MCMC_output_list$theta[1, ]
        for(i in 2:nrow(MCMC_output_list$theta)){
            result_obj$posterior[[paste0('eta_', i - 1)]] <- MCMC_output_list$theta[i, ]
        }
    }
    unique_h_idx <- !duplicated(MCMC_output_list$h)
    h_unique <- unique(MCMC_output_list$h)
    h_unique_order <- order(h_unique)
    h_unique_sorted <- h_unique[h_unique_order]
    h_idx_data <- match(h, h_unique_sorted)
    result_obj$posterior$theta <- MCMC_output_list$theta
    result_obj$posterior$rating_curve <- exp(MCMC_output_list$y_true_post_pred[unique_h_idx, ][h_unique_order, ])
    result_obj$posterior$rating_curve_mean <- exp(MCMC_output_list$mu_post[unique_h_idx, ][h_unique_order, ])
    result_obj$posterior$sigma_eps <- sqrt(MCMC_output_list$sigma_eps[unique_h_idx, ][h_unique_order, ])
    if(!is.null(Q_me)){
        result_obj$posterior$Q_true <- exp(MCMC_output_list$y_true)
    }
    #summary objects
    result_obj$summary$rating_curve <- get_MCMC_summary(result_obj$posterior$rating_curve, h = h_unique_sorted)
    result_obj$summary$rating_curve_mean <- get_MCMC_summary(result_obj$posterior$rating_curve_mean, h = h_unique_sorted)
    result_obj$summary$sigma_eps <- get_MCMC_summary(result_obj$posterior$sigma_eps, h = h_unique_sorted)
    result_obj$summary$parameters <- get_MCMC_summary(rbind(MCMC_output_list$x[1, ], MCMC_output_list$x[2, ], MCMC_output_list$theta))
    result_obj$summary$parameters$eff_n_samples <- MCMC_output_list$effective_num_samples
    result_obj$summary$parameters$r_hat <- MCMC_output_list$r_hat
    row.names(result_obj$summary$parameters) <- param_names
    if(!is.null(Q_me)){
        result_obj$summary$Q_true <- get_MCMC_summary(exp(MCMC_output_list$y_true), h = h)
    }
    result_obj$summary$log_likelihood <- get_MCMC_summary(MCMC_output_list$log_lik)
    # log_likelihood
    result_obj$log_likelihood$log_likelihood <- c(MCMC_output_list$log_lik)
    waic_list <- calc_waic(result_obj, model_dat, Q_me)
    result_obj$log_likelihood$WAIC <- waic_list$waic
    result_obj$log_likelihood$WAIC_i <- waic_list$waic_i
    result_obj$log_likelihood$lppd <- waic_list$lppd
    result_obj$log_likelihood$effective_num_param_WAIC <- waic_list$p_waic
    result_obj$log_likelihood$D_hat <- MCMC_output_list$D_hat
    result_obj$log_likelihood$effective_num_param_DIC <- -2 * result_obj$summary$log_likelihood[, 'median'] - result_obj$log_likelihood$D_hat
    result_obj$log_likelihood$DIC <- result_obj$log_likelihood$D_hat + 2 * result_obj$log_likelihood$effective_num_param_DIC
    # diagnostics
    autocorrelation_df <- as.data.frame(t(MCMC_output_list$autocorrelation))
    names(autocorrelation_df) <- param_names
    autocorrelation_df$lag <- 1:dim(autocorrelation_df)[1]
    result_obj$diagnostics$autocorrelation <- autocorrelation_df[, c('lag', param_names)]
    result_obj$diagnostics$acceptance_rate <- MCMC_output_list[['acceptance_rate']]
    # store other information
    result_obj$formula <- formula
    result_obj$data <- model_dat
    result_obj$run_info <- MCMC_output_list$run_info
    if(verbose){
        cat("\nMCMC sampling completed!\n")
        cat(sprintf("\nDiagnostics:\nAcceptance rate: %.2f%%.\n", result_obj$diagnostics$acceptance_rate * 100))
        convergence_diagnostics_warnings(result_obj$summary$parameters)
    }
    return(result_obj)
}

#' @importFrom stats optim
plm.inference <- function(y, Q_me, h, c_param = NULL, h_max = NULL, parallel = TRUE, forcepoint = rep(FALSE, length(h)), num_cores = NULL, num_chains = 4, nr_iter = 20000, burnin = 2000, thin = 5, verbose){
  if(verbose) cat("Progress:\nInitializing Metropolis MCMC algorithm...\n")
  RC <- get_model_components('plm', y, Q_me, h, c_param, h_max, forcepoint, h_min = NULL)
  output_list <- get_MCMC_output_list(theta_m = RC$theta_m, RC = RC, density_fun = RC$density_fun,
                                      unobserved_prediction_fun = RC$unobserved_prediction_fun,
                                      parallel = parallel, num_cores = num_cores, num_chains = num_chains, nr_iter = nr_iter,
                                      burnin = burnin, thin = thin, verbose = verbose)
  output_list$D_hat <- plm.calc_Dhat(output_list$theta, RC)
  if(is.null(RC$c)){
    output_list$theta[1, ] <- RC$h_min - exp(output_list$theta[1, ])
    output_list$theta[2, ] <- exp(output_list$theta[2, ])
    output_list$theta[3:8, ] <- RC$P %*% output_list$theta[3:8, ]
  }else{
    output_list$theta[1,] <- exp(output_list$theta[1, ])
    output_list$theta[2:7,] <- RC$P %*% output_list$theta[2:7, ]
  }
  output_list$x[1, ] <- exp(output_list$x[1, ])
  output_list[['h']] <- c(RC$h, RC$h_u)
  output_list[['run_info']] <- list('c_param' = c_param, 'h_max' = h_max, 'forcepoint' = forcepoint, 'nr_iter' = nr_iter, 'num_chains' = num_chains, 'burnin' = burnin, 'thin' = thin)
  return(output_list)
}

plm.density_evaluation_unknown_c <- function(theta, RC) {
    if(is.null(RC$Q_me)){
        return(plm_density_evaluation_unknown_c_cpp(theta = theta,
                                                    P = RC$P,
                                                    h = RC$h,
                                                    B = RC$B,
                                                    y = RC$y,
                                                    epsilon = RC$epsilon,
                                                    Sig_ab = RC$Sig_ab,
                                                    mu_x = RC$mu_x,
                                                    h_min = RC$h_min,
                                                    nugget = RC$nugget,
                                                    lambda_c = RC$lambda_c,
                                                    lambda_eta_1 = RC$lambda_eta_1,
                                                    lambda_seta = RC$lambda_seta
        ))
    } else {
        return(plm_me_density_evaluation_unknown_c_cpp(theta = theta,
                                                       P = RC$P,
                                                       h = RC$h,
                                                       B = RC$B,
                                                       y = RC$y,
                                                       tau = RC$tau,
                                                       epsilon = RC$epsilon,
                                                       Sig_ab = RC$Sig_ab,
                                                       mu_x = RC$mu_x,
                                                       h_min = RC$h_min,
                                                       nugget = RC$nugget,
                                                       lambda_c = RC$lambda_c,
                                                       lambda_eta_1 = RC$lambda_eta_1,
                                                       lambda_seta = RC$lambda_seta
        ))
    }
}

plm.density_evaluation_known_c <- function(theta, RC) {
    if(is.null(RC$Q_me)){
        return(plm_density_evaluation_known_c_cpp(theta = theta,
                                                  P = RC$P,
                                                  h = RC$h,
                                                  B = RC$B,
                                                  y = RC$y,
                                                  epsilon = RC$epsilon,
                                                  Sig_ab = RC$Sig_ab,
                                                  mu_x = RC$mu_x,
                                                  c = RC$c,
                                                  nugget = RC$nugget,
                                                  lambda_eta_1 = RC$lambda_eta_1,
                                                  lambda_seta = RC$lambda_seta
        ))
    } else {
        return(plm_me_density_evaluation_known_c_cpp(theta = theta,
                                                     P = RC$P,
                                                     h = RC$h,
                                                     B = RC$B,
                                                     y = RC$y,
                                                     tau = RC$tau,
                                                     epsilon = RC$epsilon,
                                                     Sig_ab = RC$Sig_ab,
                                                     mu_x = RC$mu_x,
                                                     c = RC$c,
                                                     nugget = RC$nugget,
                                                     lambda_eta_1 = RC$lambda_eta_1,
                                                     lambda_seta = RC$lambda_seta
        ))
    }
}

plm.predict_u_unknown_c <- function(theta, x, RC) {
    return(plm_predict_u_unknown_c_cpp(theta = theta,
                                       x = x,
                                       P = RC$P,
                                       B_u = RC$B_u,
                                       h_u = RC$h_u,
                                       h_min = RC$h_min,
                                       n_u = length(RC$h_u)
    ))
}

plm.predict_u_known_c <- function(theta, x, RC) {
    return(plm_predict_u_known_c_cpp(theta = theta,
                                     x = x,
                                     P = RC$P,
                                     B_u = RC$B_u,
                                     h_u = RC$h_u,
                                     c = RC$c
    ))
}

#' @importFrom stats dnorm
plm.calc_Dhat <- function(theta, RC){
    theta_median <- sapply(1:dim(theta)[1], function(x) median(theta[x, ]))
    if(!is.null(RC$c)){
        theta_median <- c(log(RC$h_min - RC$c), theta_median)
    }
    zeta <- theta_median[1]
    log_sig_eta <- theta_median[2]
    eta_1 <- theta_median[3]
    z <- theta_median[4:8]
    lambda <- c(RC$P %*% as.matrix(c(eta_1, exp(log_sig_eta) * z)))

    l <- c(log(RC$h - RC$h_min + exp(zeta)))
    varr <- c(RC$epsilon * exp(RC$B %*% lambda))
    if(any(varr > 10^2)) return(list(p = -1e9)) # to avoid numerical instability

    if (is.null(RC$Q_me)) {
        # Case without measurement error (original version)
        Sig_x <- RC$Sig_ab
        X <- cbind(1, l)
        Sig_eps <- diag(varr)
    } else {
        # Case with measurement error
        Sig_u2 <- diag(varr)
        Sig_x <- rbind(
            cbind(RC$Sig_ab, matrix(0, 2, RC$n)),
            cbind(matrix(0, RC$n, 2), Sig_u2)
        )
        X <- cbind(1, l, diag(RC$n))
        Sig_eps <- diag(RC$tau^2)
    }

    L <- compute_L(X, Sig_x, Sig_eps, RC$nugget)
    w <- compute_w(L, RC$y, X, RC$mu_x)
    x <- RC$mu_x + (Sig_x %*% (t(X) %*% solveArma(t(L), w)))

    if (is.null(RC$Q_me)) {
        # Original calculation for non-ME case
        mu_post <- (X %*% x)[1:RC$n, ]
        Dhat <- -2 * sum(dnorm(RC$y[1:RC$n, ], mu_post, sqrt(varr), log = TRUE))
    } else {
        # Calculation for ME case
        beta <- x[1:2]  # a_0 and b
        mu_post <- beta[1] + beta[2] * l
        Dhat <- -2 * sum(dnorm(RC$y[1:RC$n, ], mu_post, sqrt(varr + RC$tau^2), log = TRUE))
    }

    return(Dhat)
}


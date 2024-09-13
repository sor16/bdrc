#' Generalized power-law model with a constant log-error variance
#'
#' gplm0 is used to fit a discharge rating curve for paired measurements of stage and discharge using a generalized power-law model with a constant log-error variance as described in Hrafnkelsson et al. (2022). It also supports discharge measurement-error data. See "Details" for a more in depth description of the model.
#'
#' @param formula An object of class "formula", with discharge column name as response and stage column name as a covariate, i.e. of the form \code{y} ~ \code{x} where \code{y} is discharge in m\eqn{^3/}s and \code{x} is stage in m (the data must be in the correct units). To include discharge measurement-error data, use the form \code{y | y_error ~ x}, where \code{y_error} is the column name for discharge measurement errors.
#' @param data A data.frame containing the variables specified in formula.
#' @param c_param The largest stage value for which there is zero discharge. If NULL (default), it is treated as unknown in the model and inferred from the data.
#' @param h_max The maximum stage to which the rating curve should extrapolate to. If NULL (default), the maximum stage value in the data is selected as an upper bound.
#' @param parallel A logical value indicating whether to run the MCMC in parallel or not. Defaults to TRUE.
#' @param num_cores An integer between 1 and 4 (number of MCMC chains) indicating how many cores to use. Only used if parallel=TRUE. If NULL, the number of cores available on the device is detected automatically.
#' @param forcepoint A logical vector of the same length as the number of rows in data. If an element at index \eqn{i} is TRUE it indicates that the rating curve should be forced through the \eqn{i}-th measurement. Use with care, as this will strongly influence the resulting rating curve.
#' @param verbose A logical value indicating whether to print progress and diagnostic information. If `TRUE`, the function will print messages as it runs. If `FALSE`, the function will run silently. Default is `TRUE`.
#'
#' @details The generalized power-law model is of the form
#' \deqn{Q=a(h-c)^{f(h)}}
#' where \eqn{Q} is discharge, \eqn{h} is stage, \eqn{a} and \eqn{c} are unknown constants and \eqn{f} is a function of \eqn{h}, referred to as the generalized power-law exponent.\cr\cr
#' The generalized power-law model is here inferred by using a Bayesian hierarchical model. The function \eqn{f} is modeled at the latent level as a fixed constant \eqn{b} plus a continuous stochastic process, \eqn{\beta(h)}, which is assumed to be twice differentiable. The model is on a logarithmic scale
#' \deqn{\log(Q_i) = \mu_i + \varepsilon_i\hspace{20mm}}
#' \deqn{\hspace{24mm}\mu_i = \log(a) + (b + \beta(h_i)) \log(h_i - c)}
#' for \eqn{i = 1,...,n}, where \eqn{\varepsilon_i} follows a normal distribution with mean zero and constant variance \eqn{\sigma_\varepsilon^2}. The stochastic process \eqn{\beta(h)} is assumed a priori to be a Gaussian process governed by a Matern covariance function with smoothness parameter \eqn{\nu = 2.5}.\cr\cr
#' When measurement error is included, the model accounts for the uncertainty in the discharge measurements, and \eqn{\sigma_\varepsilon^2} captures the remaining structural uncertainty. The measurement-error datum \eqn{Q_{\text{SE},i}} corresponding to an observed discharge \eqn{Q_{\text{OBS},i}} is assumed to be the standard deviation of a log-normal distribution with median value at the true observations, \eqn{Q_{\text{TRUE},i}}. Then the model can be summarized as
#' \deqn{\log(Q_{\text{OBS},i})\sim \mathcal{N}(\log(Q_{\text{TRUE},i}),\tau^2_i),}
#' \deqn{\log(Q_{\text{TRUE},i})\sim \mathcal{N}(\mu_i,\sigma_\varepsilon^2),\hspace{10mm}}
#' where \eqn{\tau^2_i} is the variance of the normal distribution corresponding to the measurement-error datum \eqn{Q_{\text{SE},i}}, which can be estimated as \eqn{\hat{\tau}^2_i=\log(1+(Q_{\text{SE},i}/Q_{\text{OBS},i})^2)}.
#' Both numerical and categorical measurement-error data are supported:\cr
#' \itemize{
#'   \item Numerical data: Direct standard deviation values of the measurement error, \eqn{Q_{\text{SE}}}.
#'   \item Categorical data: USGS discharge measurement quality codes, which are automatically transformed into numerical \eqn{Q_{\text{SE}}} values as follows:
#'     \itemize{
#'       \item E (Excellent): Within 2\% of the actual flow
#'       \item G (Good): Within 5\% of the actual flow
#'       \item F (Fair): Within 8\% of the actual flow
#'       \item P (Poor): Greater than 8\% of the actual flow
#'     }
#' }
#' An efficient posterior simulation is achieved by sampling from the joint posterior density of the hyperparameters of the model, and then sampling from the density of the latent parameters conditional on the hyperparameters.\cr\cr
#' Bayesian inference is based on the posterior density and summary statistics such as the posterior mean and 95\% posterior intervals are based on the posterior density. Analytical formulas for these summary statistics are intractable in most cases and thus they are computed by generating samples from the posterior density using a Markov chain Monte Carlo simulation.
#' @return gplm0 returns an object of class "gplm0". An object of class "gplm0" is a list containing the following components:
#' \describe{
#'   \item{\code{posterior}}{A list containing the full posterior samples for various parameters and derived quantities.}
#'   \item{\code{summary}}{A list containing summary statistics (2.5\%, 50\%, 97.5\% percentiles) for various parameters and derived quantities.}
#'   \item{\code{log_likelihood}}{A list containing log-likelihood related quantities including WAIC and DIC.}
#'   \item{\code{diagnostics}}{A list containing MCMC diagnostics including autocorrelation and acceptance rate.}
#'   \item{\code{formula}}{The formula provided by the user.}
#'   \item{\code{data}}{The data provided by the user, ordered by stage.}
#'   \item{\code{run_info}}{Information about the input arguments and specific parameters used in the MCMC chain.}
#' }
#' When discharge measurement-error data is included, the posterior and summary components also include Q_true (the estimated true discharge observations).
#' @references Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., and Rubin, D. B. (2013). Bayesian Data Analysis, Third Edition. Chapman & Hall/CRC Texts in Statistical Science. Taylor & Francis. doi: https://doi.org/10.1201/b16018
#' @references Hrafnkelsson, B., Sigurdarson, H., Rögnvaldsson, S., Jansson, A. Ö., Vias, R. D., and Gardarsson, S. M. (2022). Generalization of the power-law rating curve using hydrodynamic theory and Bayesian hierarchical modeling, Environmetrics, 33(2):e2711. doi: https://doi.org/10.1002/env.2711
#' @references Spiegelhalter, D., Best, N., Carlin, B., Van Der Linde, A. (2002). Bayesian measures of model complexity and fit. Journal of the Royal Statistical Society: Series B (Statistical Methodology) 64(4), 583–639. doi: https://doi.org/10.1111/1467-9868.00353
#' @references Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and widely applicable information criterion in singular learning theory. J. Mach. Learn. Res. 11, 3571–3594.
#' @seealso \code{\link{summary.gplm0}} for summaries, \code{\link{predict.gplm0}} for prediction and \code{\link{plot.gplm0}} for plots. \code{\link{spread_draws}} and \code{\link{gather_draws}} are also useful to aid further visualization of the full posterior distributions.
#'
#' @examples
#' \donttest{
#' data(krokfors)
#' data(provo)
#' data(mokelumne)
#' set.seed(1)
#'
#' # Without measurement error data
#' gplm0.fit <- gplm0(formula = Q ~ W, data = krokfors, num_cores = 2)
#' summary(gplm0.fit)
#'
#' # With numerical measurement error data
#' gplm0_me.fit <- gplm0(formula = Q | Q_sigma ~ W, data = provo, num_cores = 2)
#' summary(gplm0_me.fit)
#'
#' # With categorical measurement error data
#' gplm0_me_cat.fit <- gplm0(formula = Q | Q_quality ~ W, data = mokelumne, num_cores = 2)
#' summary(gplm0_me_cat.fit)
#' }
#' @export
gplm0 <- function(formula, data, c_param = NULL, h_max = NULL, parallel = TRUE, num_cores = NULL, forcepoint = rep(FALSE, nrow(data)), verbose = TRUE){
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
    MCMC_output_list <- gplm0.inference(y = log(Q), Q_me = Q_me, h = h, c_param = c_param, h_max = h_max, parallel = parallel, forcepoint = forcepoint, num_cores = num_cores, verbose = verbose)
    param_names <- get_param_names('gplm0', c_param)

    #prepare S3 model object to be returned
    result_obj <- list(posterior = list(), log_likelihood = list(), summary = list(), diagnostics = list())
    attr(result_obj, "class") <- "gplm0"
    result_obj$posterior$a <- MCMC_output_list$x[1, ]
    result_obj$posterior$b <- MCMC_output_list$x[2, ]
    if(is.null(c_param)){
        result_obj$posterior$c <- MCMC_output_list$theta[1, ]
        result_obj$posterior$sigma_eps <- MCMC_output_list$theta[2, ]
        result_obj$posterior$sigma_beta <- MCMC_output_list$theta[3, ]
        result_obj$posterior$phi_beta <- MCMC_output_list$theta[4, ]
    }else{
        result_obj$posterior$c <- NULL
        result_obj$posterior$sigma_eps <- MCMC_output_list$theta[1, ]
        result_obj$posterior$sigma_beta <- MCMC_output_list$theta[2, ]
        result_obj$posterior$phi_beta <- MCMC_output_list$theta[3, ]
    }
    unique_h_idx <- !duplicated(MCMC_output_list$h)
    h_unique <- unique(MCMC_output_list$h)
    h_unique_order <- order(h_unique)
    h_unique_sorted <- h_unique[h_unique_order]
    h_idx_data <- match(h, h_unique_sorted)
    result_obj$posterior$theta <- MCMC_output_list$theta
    result_obj$posterior$rating_curve <- exp(MCMC_output_list$y_true_post_pred[unique_h_idx, ][h_unique_order, ])
    result_obj$posterior$rating_curve_mean <- exp(MCMC_output_list$mu_post[unique_h_idx, ][h_unique_order, ])
    result_obj$posterior$beta <- MCMC_output_list$x[3:nrow(MCMC_output_list$x), ][h_unique_order, ]
    result_obj$posterior$f <- matrix(rep(result_obj$posterior$b, nrow(result_obj$posterior$beta)), nrow = nrow(result_obj$posterior$beta), byrow = TRUE) + result_obj$posterior$beta
    if(!is.null(Q_me)){
        result_obj$posterior$Q_true <- exp(MCMC_output_list$y_true)
    }
    #summary objects
    result_obj$summary$rating_curve <- get_MCMC_summary(result_obj$posterior$rating_curve, h = h_unique_sorted)
    result_obj$summary$rating_curve_mean <- get_MCMC_summary(result_obj$posterior$rating_curve_mean, h = h_unique_sorted)
    result_obj$summary$beta <- get_MCMC_summary(result_obj$posterior$beta, h = h_unique_sorted)
    result_obj$summary$f <- get_MCMC_summary(result_obj$posterior$f, h = h_unique_sorted)
    result_obj$summary$parameters <- get_MCMC_summary(rbind(MCMC_output_list$x[1, ], MCMC_output_list$x[2, ], MCMC_output_list$theta))
    result_obj$summary$parameters$eff_n_samples <- MCMC_output_list$effective_num_samples
    result_obj$summary$parameters$r_hat <- MCMC_output_list$r_hat
    row.names(result_obj$summary$parameters) <- param_names
    if(!is.null(Q_me)){
        result_obj$summary$Q_true <- get_MCMC_summary(exp(MCMC_output_list$y_true), h = h)
    }
    result_obj$summary$log_likelihood <- get_MCMC_summary(MCMC_output_list$log_lik)
    # log_likelihood objects
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

#' @importFrom stats dist optim
gplm0.inference <- function(y, Q_me, h, c_param = NULL, h_max = NULL, parallel = TRUE, forcepoint = rep(FALSE, length(h)), num_cores = NULL, num_chains = 4, nr_iter = 20000, burnin = 2000, thin = 5, verbose){
    if(verbose) cat("Progress:\nInitializing Metropolis MCMC algorithm...\n")
    c_upper <- NULL
    if(is.null(c_param)){
        RC_plm0 <- get_model_components('plm0', y, Q_me, h, c_param, h_max, forcepoint, h_min = NULL)
        lhmc_sd <- sqrt(diag(matInverse(RC_plm0$H)))[1]
        lhmc_mode <- RC_plm0$theta_m[1]
        if(exp(lhmc_mode - 1.96 * lhmc_sd) > 2){
            warning('Dataset lacks measurements near point of zero flow and thus the model infers its upper bound (see c_upper in run_info).')
            c_upper <- min(h) - exp(lhmc_mode - 1.96 * lhmc_sd)
        }
    }
    RC <- get_model_components('gplm0', y, Q_me, h, c_param, h_max, forcepoint, h_min = c_upper)
    output_list <- get_MCMC_output_list(theta_m = RC$theta_m, RC = RC, density_fun = RC$density_fun,
                                        unobserved_prediction_fun = RC$unobserved_prediction_fun,
                                        parallel = parallel, num_cores = num_cores, num_chains = num_chains, nr_iter = nr_iter,
                                        burnin = burnin, thin = thin, verbose = verbose)
    output_list$D_hat <- gplm0.calc_Dhat(output_list$theta, RC)
    #refinement of list elements
    if(is.null(RC$c)){
        output_list$theta[1, ] <- RC$h_min - exp(output_list$theta[1, ])
        output_list$theta[2, ] <- sqrt(exp(output_list$theta[2, ]))
        output_list$theta[3, ] <- sqrt(exp(output_list$theta[3, ]))
        output_list$theta[4, ] <- exp(output_list$theta[4, ])
    }else{
        output_list$theta[1, ] <- sqrt(exp(output_list$theta[1, ]))
        output_list$theta[2, ] <- sqrt(exp(output_list$theta[2, ]))
        output_list$theta[3, ] <- exp(output_list$theta[3, ])
    }
    output_list$x[1,] <- exp(output_list$x[1, ])
    output_list[['h']] <- c(RC$h, RC$h_u)
    output_list[['run_info']] <- list('c_param' = c_param, 'h_max' = h_max, 'forcepoint' = forcepoint, 'nr_iter' = nr_iter, 'num_chains' = num_chains, 'burnin' = burnin, 'thin' = thin, 'c_upper' = c_upper)
    return(output_list)
}

gplm0.density_evaluation_unknown_c <- function(theta, RC) {
    if(is.null(RC$Q_me)){
        return(gplm0_density_evaluation_unknown_c_cpp(theta = theta,
                                                      h = RC$h,
                                                      y = RC$y,
                                                      A = RC$A,
                                                      dist = RC$dist,
                                                      epsilon = RC$epsilon,
                                                      h_min = RC$h_min,
                                                      nugget = RC$nugget,
                                                      n_unique = RC$n_unique,
                                                      mu_x = RC$mu_x,
                                                      Sig_ab = RC$Sig_ab,
                                                      Z = RC$Z,
                                                      lambda_c = RC$lambda_c,
                                                      lambda_se = RC$lambda_se,
                                                      lambda_sb = RC$lambda_sb,
                                                      lambda_pb = RC$lambda_pb
        ))
    } else {
        return(gplm0_me_density_evaluation_unknown_c_cpp(theta = theta,
                                                         h = RC$h,
                                                         y = RC$y,
                                                         A = RC$A,
                                                         dist = RC$dist,
                                                         epsilon = RC$epsilon,
                                                         tau = RC$tau,
                                                         h_min = RC$h_min,
                                                         nugget = RC$nugget,
                                                         n_unique = RC$n_unique,
                                                         mu_x = RC$mu_x,
                                                         Sig_ab = RC$Sig_ab,
                                                         Z = RC$Z,
                                                         lambda_c = RC$lambda_c,
                                                         lambda_se = RC$lambda_se,
                                                         lambda_sb = RC$lambda_sb,
                                                         lambda_pb = RC$lambda_pb
        ))
    }
}

gplm0.density_evaluation_known_c <- function(theta, RC) {
    if(is.null(RC$Q_me)){
        return(gplm0_density_evaluation_known_c_cpp(theta = theta,
                                                    h = RC$h,
                                                    y = RC$y,
                                                    A = RC$A,
                                                    dist = RC$dist,
                                                    epsilon = RC$epsilon,
                                                    c = RC$c,
                                                    nugget = RC$nugget,
                                                    n_unique = RC$n_unique,
                                                    mu_x = RC$mu_x,
                                                    Sig_ab = RC$Sig_ab,
                                                    Z = RC$Z,
                                                    lambda_se = RC$lambda_se,
                                                    lambda_sb = RC$lambda_sb,
                                                    lambda_pb = RC$lambda_pb
        ))
    } else {
        return(gplm0_me_density_evaluation_known_c_cpp(theta = theta,
                                                       h = RC$h,
                                                       y = RC$y,
                                                       A = RC$A,
                                                       dist = RC$dist,
                                                       epsilon = RC$epsilon,
                                                       tau = RC$tau,
                                                       c = RC$c,
                                                       nugget = RC$nugget,
                                                       n_unique = RC$n_unique,
                                                       mu_x = RC$mu_x,
                                                       Sig_ab = RC$Sig_ab,
                                                       Z = RC$Z,
                                                       lambda_se = RC$lambda_se,
                                                       lambda_sb = RC$lambda_sb,
                                                       lambda_pb = RC$lambda_pb
        ))
    }
}

gplm0.predict_u_unknown_c <- function(theta, x, RC) {
    return(gplm0_predict_u_unknown_c_cpp(theta = theta,
                                         x = x,
                                         h_unique = RC$h_unique,
                                         h_u = RC$h_u,
                                         dist_all = RC$dist_all,
                                         h_min = RC$h_min,
                                         nugget = RC$nugget,
                                         n_unique = RC$n_unique,
                                         n_u = RC$n_u
    ))
}

gplm0.predict_u_known_c <- function(theta, x, RC) {
    return(gplm0_predict_u_known_c_cpp(theta = theta,
                                       x = x,
                                       h_unique = RC$h_unique,
                                       h_u = RC$h_u,
                                       dist_all = RC$dist_all,
                                       c = RC$c,
                                       nugget = RC$nugget,
                                       n_unique = RC$n_unique,
                                       n_u = RC$n_u
    ))
}

#' @importFrom stats dnorm
gplm0.calc_Dhat <- function(theta, RC){
    theta_median <- sapply(1:dim(theta)[1], function(x) median(theta[x, ]))
    if(!is.null(RC$c)){
        theta_median <- c(log(RC$h_min - RC$c), theta_median)
    }
    zeta <- theta_median[1]
    log_sig_eps2 <- theta_median[2]
    log_sig_b <- theta_median[3]
    log_phi_b <- theta_median[4]

    l <- c(log(RC$h - RC$h_min + exp(zeta)))
    varr <- RC$epsilon * exp(log_sig_eps2)
    if(any(varr > 10^2)) return(list(p = -1e9)) # to avoid numerical instability

    # repeated calculations
    phi_b <- exp(log_phi_b)
    var_b <- exp(2 * log_sig_b)
    sqrt_5 <- sqrt(5)

    # Matern covariance
    R_Beta <- var_b * ((1 + sqrt_5 * RC$dist / phi_b + 5 * RC$dist^2 / (3 * phi_b^2)) * exp(-sqrt_5 * RC$dist / phi_b) + diag(RC$n_unique) * RC$nugget)

    if (is.null(RC$Q_me)) {
        # Case without measurement error (original version)
        Sig_x <- rbind(cbind(RC$Sig_ab, RC$m1), cbind(RC$m2, R_Beta))
        X <- rbind(cbind(1, l, matMult(diag(l), RC$A)), RC$Z)
        Sig_eps <- diag(c(varr, 0))
    } else {
        # Case with measurement error
        Sig_u2 <- diag(varr)
        Sig_x <- rbind(
            cbind(RC$Sig_ab, RC$m1, matrix(0, 2, RC$n)),
            cbind(RC$m2, R_Beta, matrix(0, RC$n_unique, RC$n)),
            cbind(matrix(0, RC$n, 2 + RC$n_unique), Sig_u2)
        )
        X <- rbind(
            cbind(matrix(1, RC$n, 1), l, matMult(diag(l), RC$A), diag(RC$n)),
            RC$Z
        )
        Sig_eps <- diag(c(RC$tau^2, 0))
    }

    L <- compute_L(X, Sig_x, diag(c(varr, 0)), RC$nugget)
    w <- compute_w(L, RC$y, X, RC$mu_x)
    x <- RC$mu_x + (Sig_x %*% (t(X) %*% solveArma(t(L), w)))

    if (is.null(RC$Q_me)) {
        # Original calculation for non-ME case
        mu_post <- (matMult(X, x))[1:RC$n, ]
        Dhat <- -2 * sum(dnorm(RC$y[1:RC$n, ], mu_post, sqrt(varr), log = TRUE))
    } else {
        # Calculation for ME case
        beta <- x[1:2]  # a_0 and b
        u1 <- x[3:(2 + RC$n_unique)]  # beta(h)
        mu_post <- beta[1] + (beta[2] + RC$A %*% u1) * l
        Dhat <- -2 * sum(dnorm(RC$y[1:RC$n, ], mu_post, sqrt(varr + RC$tau^2), log = TRUE))
    }

    return(Dhat)
}

#' Generalized power-law model with variance that varies with stage.
#'
#' gplm is used to fit a discharge rating curve for paired measurements of stage and discharge using a generalized power-law model with variance that varies with stage as described in Hrafnkelsson et al. (2022).  See "Details" for a more elaborate description of the model.
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
#' @details The generalized power-law model is of the form
#' \deqn{Q=a(h-c)^{f(h)}}
#' where \eqn{Q} is discharge, \eqn{h} is stage, \eqn{a} and \eqn{c} are unknown constants and \eqn{f} is a function of \eqn{h}, referred to as the generalized power-law exponent.\cr\cr
#' The generalized power-law model is here inferred by using a Bayesian hierarchical model. The function \eqn{f} is modeled at the latent level as a fixed constant \eqn{b} plus a continuous stochastic process, \eqn{\beta(h)}, which is assumed to be twice differentiable. The model is on a logarithmic scale
#' \deqn{\log(Q_i) = \log(a) + (b + \beta(h_i)) \log(h_i - c) + \varepsilon_i,     i = 1,...,n}
#' where \eqn{\varepsilon_i} follows a normal distribution with mean zero and variance \eqn{\sigma_\varepsilon(h_i)^2} that varies with stage. The stochastic process \eqn{\beta(h)} is assumed a priori to be a Gaussian process governed by a Matern covariance function with smoothness parameter \eqn{\nu = 2.5}. The error variance, \eqn{\sigma_\varepsilon^2(h)}, of the log-discharge data  is modeled as an exponential of a B-spline curve, that is, a linear combination of six B-spline basis functions that are defined over the range of the stage observations. An efficient posterior simulation is achieved by sampling from the joint posterior density of the hyperparameters of the model, and then sampling from the density of the latent parameters conditional on the hyperparameters.\cr\cr
#' Bayesian inference is based on the posterior density and summary statistics such as the posterior mean and 95\% posterior intervals are based on the posterior density. Analytical formulas for these summary statistics are intractable in most cases and thus they are computed by generating samples from the posterior density using a Markov chain Monte Carlo simulation.
#' @return gplm returns an object of class "gplm". An object of class "gplm" is a list containing the following components:
#' \describe{
#'   \item{\code{rating_curve}}{A data frame with 2.5\%, 50\% and 97.5\% percentiles of the posterior predictive distribution of the rating curve.}
#'   \item{\code{rating_curve_mean}}{A data frame with 2.5\%, 50\% and 97.5\% percentiles of the posterior distribution of the mean of the rating curve.}
#'   \item{\code{param_summary}}{A data frame with 2.5\%, 50\% and 97.5\% percentiles of the posterior distribution of latent- and hyperparameters. Additionally contains columns with r_hat and the effective number of samples for each parameter as defined in Gelman et al. (2013).}
#'   \item{\code{f_summary}}{A data frame with 2.5\%, 50\% and 97.5\% percentiles of the posterior distribution of \eqn{f(h)}.}
#'   \item{\code{beta_summary}}{A data frame with 2.5\%, 50\% and 97.5\% percentiles of the posterior distribution of \eqn{\beta(h)}.}
#'   \item{\code{sigma_eps_summary}}{A data frame with 2.5\%, 50\% and 97.5\% percentiles of the posterior distribution of \eqn{\sigma_\varepsilon(h)}.}
#'   \item{\code{posterior_log_likelihood_summary}}{A data frame with 2.5\%, 50\% and 97.5\% percentiles of the posterior log-likelihood values.}
#'   \item{\code{rating_curve_posterior}}{A matrix containing the full thinned posterior samples of the posterior predictive distribution of the rating curve excluding burn-in samples.}
#'   \item{\code{rating_curve_mean_posterior}}{A matrix containing the full thinned posterior samples of the posterior distribution of the mean of the rating curve excluding burn-in samples.}
#'   \item{\code{a_posterior}}{A numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{a} excluding burn-in samples.}
#'   \item{\code{b_posterior}}{A numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{b} excluding burn-in samples.}
#'   \item{\code{c_posterior}}{A numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{c} excluding burn-in samples.}
#'   \item{\code{sigma_beta_posterior}}{A numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{\sigma_\beta} excluding burn-in samples.}
#'   \item{\code{phi_beta_posterior}}{A numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{\phi_\beta} excluding burn-in samples.}
#'   \item{\code{sigma_eta_posterior}}{A numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{\sigma_\eta} excluding burn-in samples.}
#'   \item{\code{eta_1_posterior}}{A numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{\eta_1} excluding burn-in samples.}
#'   \item{\code{eta_2_posterior}}{A numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{\eta_2} excluding burn-in samples.}
#'   \item{\code{eta_3_posterior}}{A numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{\eta_3} excluding burn-in samples.}
#'   \item{\code{eta_4_posterior}}{A numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{\eta_4} excluding burn-in samples.}
#'   \item{\code{eta_5_posterior}}{A numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{\eta_5} excluding burn-in samples.}
#'   \item{\code{eta_6_posterior}}{A numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{\eta_6} excluding burn-in samples.}
#'   \item{\code{f_posterior}}{A numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{f(h)} excluding burn-in samples.}
#'   \item{\code{beta_posterior}}{A numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{\beta(h)} excluding burn-in samples.}
#'   \item{\code{sigma_eps_posterior}}{A numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{\sigma_\varepsilon(h)} excluding burn-in samples.}
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
#' @seealso \code{\link{summary.gplm}} for summaries, \code{\link{predict.gplm}} for prediction and \code{\link{plot.gplm}} for plots. \code{\link{spread_draws}} and \code{\link{gather_draws}} are also useful to aid further visualization of the full posterior distributions.
#'
#' @examples
#' \donttest{
#' data(norn)
#' set.seed(1)
#' gplm.fit <- gplm(formula=Q~W,data=norn,num_cores=2)
#' summary(gplm.fit)
#' }
#' @export
gplm <- function(formula, data, c_param = NULL, h_max = NULL, parallel = TRUE, num_cores = NULL, forcepoint = rep(FALSE, nrow(data)), verbose = TRUE){
    #argument checking
    stopifnot(inherits(formula, 'formula'))
    stopifnot(inherits(data, 'data.frame'))
    stopifnot(is.null(c_param) | is.double(c_param))
    stopifnot(is.null(h_max) | is.double(h_max))
    stopifnot(is.null(num_cores) | is.numeric(num_cores))
    stopifnot(length(forcepoint) == dim(data)[1] & is.logical(forcepoint))
    formula_args <- all.vars(formula)
    stopifnot(length(formula_args) == 2 & all(formula_args %in% names(data)))
    model_dat <- as.data.frame(data[,all.vars(formula)])
    forcepoint <- forcepoint[order(model_dat[, 2,drop = TRUE])]
    model_dat <- model_dat[order(model_dat[, 2, drop = TRUE]),]
    if(dim(model_dat)[1] < 2) stop('At least two paired observations of stage and discharge are required to fit a rating curve')
    Q <- model_dat[, 1, drop = TRUE]
    h <- model_dat[, 2, drop = TRUE]
    if(!is.null(c_param) && min(h) < c_param) stop('c_param must be lower than the minimum stage value in the data')
    if(any(Q <= 0)) stop('All discharge measurements must but strictly greater than zero. If you know the stage of zero discharge, use c_param.')
    MCMC_output_list <- gplm.inference(y = log(Q), h = h, c_param = c_param, h_max = h_max, parallel = parallel, forcepoint = forcepoint, num_cores = num_cores, verbose = verbose)
    param_names <- get_param_names('gplm', c_param)
    #prepare S3 model object to be returned
    result_obj=list()
    attr(result_obj, "class") <- "gplm"
    result_obj$a_posterior = MCMC_output_list$x[1, ]
    result_obj$b_posterior = MCMC_output_list$x[2, ]
    if(is.null(c_param)){
        result_obj$c_posterior <- MCMC_output_list$theta[1, ]
        result_obj$sigma_beta_posterior <- MCMC_output_list$theta[2, ]
        result_obj$phi_beta_posterior <- MCMC_output_list$theta[3, ]
        result_obj$sigma_eta_posterior <- MCMC_output_list$theta[4, ]
        for(i in 5:dim(MCMC_output_list$theta)[1]){
            result_obj[[paste0('eta_', i - 4, '_posterior')]] <- MCMC_output_list$theta[i, ]
        }
    }else{
        result_obj$c_posterior <- NULL
        result_obj$sigma_beta_posterior <- MCMC_output_list$theta[1, ]
        result_obj$phi_beta_posterior <- MCMC_output_list$theta[2, ]
        result_obj$sigma_eta_posterior <- MCMC_output_list$theta[3, ]
        for(i in 4:dim(MCMC_output_list$theta)[1]){
            result_obj[[paste0('eta_', i - 3, '_posterior')]] <- MCMC_output_list$theta[i, ]
        }
    }
    unique_h_idx <- !duplicated(MCMC_output_list$h)
    h_unique <- unique(MCMC_output_list$h)
    h_unique_order <- order(h_unique)
    h_unique_sorted <- h_unique[h_unique_order]
    h_idx_data <- match(h,h_unique_sorted)
    result_obj$theta <- MCMC_output_list$theta
    result_obj$rating_curve_posterior <- exp(MCMC_output_list$y_post_pred[unique_h_idx, ][h_unique_order, ])
    result_obj$rating_curve_mean_posterior <- exp(MCMC_output_list$y_post[unique_h_idx, ][h_unique_order, ])
    result_obj$beta_posterior <- MCMC_output_list$x[3:nrow(MCMC_output_list$x), ][h_unique_order, ]
    result_obj$f_posterior <- matrix(rep(result_obj$b_posterior, nrow(result_obj$beta_posterior)), nrow = nrow(result_obj$beta_posterior), byrow = TRUE) + result_obj$beta_posterior
    result_obj$sigma_eps_posterior <- sqrt(MCMC_output_list$sigma_eps[unique_h_idx, ][h_unique_order, ])
    result_obj$posterior_log_likelihood <- c(MCMC_output_list$log_lik)
    #summary objects
    result_obj$rating_curve <- get_MCMC_summary(result_obj$rating_curve_posterior, h = h_unique_sorted)
    result_obj$rating_curve_mean <- get_MCMC_summary(result_obj$rating_curve_mean_posterior, h = h_unique_sorted)
    result_obj$beta_summary <- get_MCMC_summary(result_obj$beta_posterior, h = h_unique_sorted)
    result_obj$f_summary <- get_MCMC_summary(result_obj$f_posterior, h = h_unique_sorted)
    result_obj$sigma_eps_summary <- get_MCMC_summary(result_obj$sigma_eps_posterior, h = h_unique_sorted)
    result_obj$param_summary <- get_MCMC_summary(rbind(MCMC_output_list$x[1, ], MCMC_output_list$x[2, ], MCMC_output_list$theta))
    result_obj$param_summary$eff_n_samples <- MCMC_output_list$effective_num_samples
    result_obj$param_summary$r_hat <- MCMC_output_list$r_hat
    row.names(result_obj$param_summary) <- param_names
    result_obj$posterior_log_likelihood_summary <- get_MCMC_summary(MCMC_output_list$log_lik)
    # DIC calculations
    result_obj$D_hat <- MCMC_output_list$D_hat
    result_obj$effective_num_param_DIC <- -2 * result_obj$posterior_log_likelihood_summary[, 'median'] - result_obj$D_hat
    result_obj$DIC <- result_obj$D_hat + 2 * result_obj$effective_num_param_DIC
    #WAIC calculations
    waic_list <- calc_waic(result_obj, model_dat)
    result_obj$lppd <- waic_list$lppd
    result_obj$effective_num_param_WAIC <- waic_list$p_waic
    result_obj$WAIC <- waic_list$waic
    result_obj$WAIC_i <- waic_list$waic_i
    #Rhat and autocorrelation
    autocorrelation_df <- as.data.frame(t(MCMC_output_list$autocorrelation))
    names(autocorrelation_df) <- param_names
    autocorrelation_df$lag <- 1:dim(autocorrelation_df)[1]
    result_obj$autocorrelation <- autocorrelation_df[, c('lag', param_names)]
    # store other information
    result_obj$acceptance_rate <- MCMC_output_list[['acceptance_rate']]
    result_obj$formula <- formula
    result_obj$data <- model_dat
    result_obj$run_info <- MCMC_output_list$run_info
    if(verbose){
        cat("\nMCMC sampling completed!\n")
        cat(sprintf("\nDiagnostics:\nAcceptance rate: %.2f%%.\n", result_obj$acceptance_rate * 100))
        convergence_diagnostics_warnings(result_obj$param_summary)
    }
    return(result_obj)
}

#' @importFrom stats dist optim
gplm.inference <- function(y, h, c_param = NULL, h_max = NULL, parallel = TRUE, forcepoint = rep(FALSE, length(h)), num_cores = NULL, num_chains = 4, nr_iter = 20000, burnin = 2000, thin = 5, verbose){
    if(verbose) cat("Progress:\nInitializing Metropolis MCMC algorithm...\n")
    c_upper <- NULL
    if(is.null(c_param)){
        RC_plm0 <- get_model_components('plm0', y, h, c_param, h_max, forcepoint, h_min = NULL)
        lhmc_sd <- sqrt(diag(matInverse(RC_plm0$H)))[1]
        lhmc_mode <- RC_plm0$theta_m[1]
        if(exp(lhmc_mode - 1.96 * lhmc_sd) > 2){
            warning('Dataset lacks measurements near point of zero flow and thus the model infers its upper bound (see c_upper in run_info).')
            c_upper <- min(h) - exp(lhmc_mode - 1.96 * lhmc_sd)
        }
    }
    RC <- get_model_components('gplm', y, h, c_param, h_max, forcepoint, h_min = c_upper)
    output_list <- get_MCMC_output_list(theta_m = RC$theta_m, RC = RC, density_fun = RC$density_fun,
                                        unobserved_prediction_fun = RC$unobserved_prediction_fun,
                                        parallel = parallel, num_cores = num_cores, num_chains = num_chains, nr_iter = nr_iter,
                                        burnin = burnin, thin = thin, verbose = verbose)
    #Calculate Dhat
    output_list$D_hat <- gplm.calc_Dhat(output_list$theta,RC)
    #refinement of list elements
    if(is.null(RC$c)){
        output_list$theta[1, ] <- RC$h_min-exp(output_list$theta[1, ])
        output_list$theta[2, ] <- exp(output_list$theta[2, ])
        output_list$theta[3, ] <- exp(output_list$theta[3, ])
        output_list$theta[4, ] <- exp(output_list$theta[4, ])
        output_list$theta[6:10, ] <- output_list$theta[6:10, ] * rep(output_list$theta[4, ], each = 5)
        output_list$theta[5:10, ] <- RC$P%*%output_list$theta[5:10, ]
    }else{
        output_list$theta[1, ] <- exp(output_list$theta[1, ])
        output_list$theta[2, ] <- exp(output_list$theta[2, ])
        output_list$theta[3, ] <- exp(output_list$theta[3, ])
        output_list$theta[5:9, ] <- output_list$theta[5:9, ] * rep(output_list$theta[3, ], each = 5)
        output_list$theta[4:9, ] <- RC$P%*%output_list$theta[4:9, ]
    }
    output_list$x[1, ] <- exp(output_list$x[1, ])
    output_list[['h']] <- c(RC$h, RC$h_u)
    output_list[['run_info']] <- list('c_param' = c_param, 'h_max' = h_max, 'forcepoint' = forcepoint, 'nr_iter' = nr_iter, 'num_chains' = num_chains, 'burnin' = burnin, 'thin' = thin, 'c_upper' = c_upper)
    return(output_list)
}

gplm.density_evaluation_unknown_c <- function(theta, RC) {
    return( gplm_density_evaluation_unknown_c_cpp( theta = theta,
                                                   P = RC$P,
                                                   h = RC$h,
                                                   B = RC$B,
                                                   dist = RC$dist,
                                                   A = RC$A,
                                                   y = RC$y,
                                                   epsilon = RC$epsilon,
                                                   h_min = RC$h_min,
                                                   nugget = RC$nugget,
                                                   n_unique = RC$n_unique,
                                                   mu_x = RC$mu_x,
                                                   Sig_ab = RC$Sig_ab,
                                                   Z = RC$Z,
                                                   lambda_c = RC$lambda_c,
                                                   lambda_sb = RC$lambda_sb,
                                                   lambda_pb = RC$lambda_pb,
                                                   lambda_eta_1 = RC$lambda_eta_1,
                                                   lambda_seta = RC$lambda_seta) )
}

gplm.density_evaluation_known_c <- function(theta, RC) {
    return(gplm_density_evaluation_known_c_cpp( theta = theta,
                                                P = RC$P,
                                                h = RC$h,
                                                B = RC$B,
                                                dist = RC$dist,
                                                A = RC$A,
                                                y = RC$y,
                                                epsilon = RC$epsilon,
                                                nugget = RC$nugget,
                                                n_unique = RC$n_unique,
                                                mu_x = RC$mu_x,
                                                Sig_ab = RC$Sig_ab,
                                                Z = RC$Z,
                                                c = RC$c,
                                                lambda_sb = RC$lambda_sb,
                                                lambda_pb = RC$lambda_pb,
                                                lambda_eta_1 = RC$lambda_eta_1,
                                                lambda_seta = RC$lambda_seta))
}

gplm.predict_u_unknown_c <- function(theta, x, RC) {
    return( gplm_predict_u_unknown_c_cpp( theta = theta,
                                          x = x,
                                          P = RC$P,
                                          B_u = RC$B_u,
                                          h_unique = RC$h_unique,
                                          h_u = RC$h_u,
                                          dist_all = RC$dist_all,
                                          h_min = RC$h_min,
                                          nugget = RC$nugget,
                                          n_unique = RC$n_unique,
                                          n_u = RC$n_u ) )
}

gplm.predict_u_known_c <- function(theta, x, RC) {
    return( gplm_predict_u_known_c_cpp( theta = theta,
                                        x = x,
                                        P = RC$P,
                                        B_u = RC$B_u,
                                        h_unique = RC$h_unique,
                                        h_u = RC$h_u,
                                        dist_all = RC$dist_all,
                                        c = RC$c,
                                        nugget = RC$nugget,
                                        n_unique = RC$n_unique,
                                        n_u = RC$n_u ) )
}

#' @importFrom stats dnorm
gplm.calc_Dhat <- function(theta, RC){
    theta_median <- sapply(1:dim(theta)[1], function(x) median(theta[x, ]))
    if(!is.null(RC$c)){
        theta_median <- c(log(RC$h_min - RC$c), theta_median)
    }
    zeta <- theta_median[1]
    log_sig_b <- theta_median[2]
    log_phi_b <- theta_median[3]
    log_sig_eta <- theta_median[4]
    eta_1 <- theta_median[5]
    z <- theta_median[6:10]

    eta=c(RC$P %*% as.matrix(c(eta_1, exp(log_sig_eta) * z)))

    l=c(log(RC$h - RC$h_min + exp(zeta)))
    varr=c(RC$epsilon * exp(RC$B %*% eta))
    if(any(varr > 10^2)) return(list(p = -1e9)) # to avoid numerical instability

    # repeated calculations
    phi_b <- exp(log_phi_b)
    var_b <- exp(2 * log_sig_b)
    sqrt_5 <- sqrt(5)

    #Matern covariance
    Sig_x <- rbind(cbind(RC$Sig_ab, RC$m1), cbind(RC$m2, var_b * ((1 + sqrt_5 * RC$dist / phi_b + 5 * RC$dist^2 / (3 * phi_b^2)) * exp(-sqrt_5 * RC$dist / phi_b) + diag(RC$n_unique) * RC$nugget)))

    X <- rbind(cbind(1, l, matMult(diag(l), RC$A)), RC$Z)
    L <- compute_L(X, Sig_x, diag(c(varr, 0)), RC$nugget)
    w <- compute_w(L, RC$y, X, RC$mu_x)
    x <- RC$mu_x + (Sig_x %*% (t(X) %*% solveArma(t(L), w)))

    yp=(matMult(X, x))[1:RC$n, ]

    Dhat <- -2 * sum(dnorm(RC$y[1:RC$n, ], yp, sqrt(varr), log = T))

    return(Dhat)
}

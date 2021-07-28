#' Power-law model with a constant variance
#'
#' plm0 is used to fit a discharge rating curve for paired measurements of stage and discharge using a power-law model with a constant variance as described in Hrafnkelsson et al. (2020). See "Details" for a more elaborate description of the model.
#' @param formula an object of class "formula", with discharge column name as response and stage column name as a covariate, i.e. of the form \code{y}~\code{x} where \code{y} is discharge in m\eqn{^3/}s and \code{x} is stage in m (it is very important that the data is in the correct units).
#' @param data data.frame containing the variables specified in formula.
#' @param c_param stage for which there is zero discharge. If NULL, it is treated as unknown in the model and inferred from the data.
#' @param h_max maximum stage to which the rating curve should extrapolate to. If NULL, the maximum stage value in the data is selected as an upper bound.
#' @param parallel logical value indicating whether to run the MCMC in parallel or not. Defaults to TRUE.
#' @param num_cores integer between 1 and 4 (number of MCMC chains) indicating how many cores to use. Only used if parallel=TRUE. If NULL, the number of cores available on the device is detected automatically.
#' @param forcepoint logical vector of the same length as the number of rows in data. If an element at index \eqn{i} is TRUE it indicates that the rating curve should be forced through the \eqn{i}-th measurement. Use with care, as this will strongly influence the resulting rating curve.
#'
#' @details The power-law model, which is commonly used in hydraulic practice, is of the form
#' \deqn{Q=a(h-c)^{b}}
#' where \eqn{Q} is discharge, \eqn{h} is stage and \eqn{a}, \eqn{b} and \eqn{c} are unknown constants.\cr\cr
#' The power-law model is here inferred by using a Bayesian hierarchical model. The model is on a logarithmic scale
#' \deqn{\log(Q_i) = \log(a) + b \log(h_i - c) + \varepsilon,     i = 1,...,n}
#' where \eqn{\varepsilon} follows a normal distribution with mean zero and variance \eqn{\sigma_\varepsilon^2}, independent of stage. An efficient posterior simulation is achieved by sampling from the joint posterior density of the hyperparameters of the model, and then sampling from the density of the latent parameters conditional on the hyperparameters.\cr\cr
#' Bayesian inference is based on the posterior density and summary statistics such as the posterior mean and 95\% posterior intervals are based on the posterior density. Analytical formulas for these summary statistics are intractable in most cases and thus they are computed by generating samples from the posterior density using a Markov chain Monte Carlo simulation.
#' @return plm0 returns an object of class "plm0". An object of class "plm0" is a list containing the following components: \cr
#' \item{\code{rating_curve}}{a data frame with 2.5\%, 50\% and 97.5\% percentiles of the posterior predictive distribution of the rating curve.}
#' \item{\code{rating_curve_mean}}{a data frame with 2.5\%, 50\% and 97.5\% percentiles of the posterior distribution of the mean of the rating curve.}
#' \item{\code{param_summary}}{a data frame with 2.5\%, 50\% and 97.5\% percentiles of the posterior distribution of latent- and hyperparameters. Additionally contains columns with r_hat and the effective number of samples for each parameter as defined in Gelman et al. (2013).}
#' \item{\code{Deviance_summary}}{a data frame with 2.5\%, 50\% and 97.5\% percentiles of the posterior distribution of the deviance.}
#' \item{\code{rating_curve_posterior}}{a matrix containing the full thinned posterior samples of the posterior predictive distribution of the rating curve (excluding burn-in).}
#' \item{\code{rating_curve_mean_posterior}}{a matrix containing the full thinned posterior samples of the posterior distribution of the mean of the rating curve (excluding burn-in).}
#' \item{\code{a_posterior}}{a numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{a}.}
#' \item{\code{b_posterior}}{a numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{b}.}
#' \item{\code{c_posterior}}{a numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{c}.}
#' \item{\code{sigma_eps_posterior}}{a numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{\sigma_{\varepsilon}}.}
#' \item{\code{Deviance_posterior}}{a numeric vector containing the full thinned posterior samples of the posterior distribution of the deviance excluding burn-in samples.}
#' \item{\code{D_hat}}{deviance at the median value of the parameters}
#' \item{\code{num_effective_param}}{number of effective parameters, which is calculated as median(Deviance_posterior) minus D_hat.}
#' \item{\code{DIC}}{Deviance Information Criterion for the model, calculated as D_hat plus 2*num_effective_parameters.}
#' \item{\code{autocorrelation}}{a data frame with the autocorrelation of each parameter for different lags.}
#' \item{\code{acceptance_rate}}{proportion of accepted samples in the thinned MCMC chain (excluding burn-in).}
#' \item{\code{formula}}{object of type "formula" provided by the user.}
#' \item{\code{data}}{data provided by the user, ordered by stage.}
#' \item{\code{run_info}}{information about the input arguments and the specific parameters used in the MCMC chain.}
#' @references Hrafnkelsson, B., Sigurdarson, H., and Gardarsson, S. M. (2020). Generalization of the power-law rating curve using hydrodynamic theory and Bayesian hierarchical modeling. arXiv preprint 2010.04769.
#' @references Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., and Rubin, D. B. (2013). Bayesian Data Analysis, Third Edition. Chapman & Hall/CRC Texts in Statistical Science. Taylor & Francis.
#' @seealso \code{\link{summary.plm0}} for summaries, \code{\link{predict.plm0}} for prediction. It is also useful to look at \code{\link{spread_draws}} and \code{\link{plot.plm0}} to help visualize the full posterior distributions.
#' @examples
#' \donttest{
#' data(skogsliden)
#' set.seed(1)
#' plm0.fit <- plm0(formula=Q~W,data=skogsliden,num_cores=2)
#' summary(plm0.fit)
#' }
#' @export
plm0 <- function(formula,data,c_param=NULL,h_max=NULL,parallel=TRUE,num_cores=NULL,forcepoint=rep(FALSE,nrow(data))){
    #argument checking
    stopifnot('formula' %in% class(formula))
    stopifnot('data.frame' %in% class(data))
    stopifnot(is.null(c_param) | is.double(c_param))
    stopifnot(is.null(h_max) | is.double(h_max))
    stopifnot(is.null(num_cores) | is.numeric(num_cores))
    stopifnot(length(forcepoint)==nrow(data) & is.logical(forcepoint))
    formula_args <- all.vars(formula)
    stopifnot(length(formula_args)==2 & all(formula_args %in% names(data)))
    model_dat <- as.data.frame(data[,formula_args])
    forcepoint <- forcepoint[order(model_dat[,2,drop=TRUE])]
    model_dat <- model_dat[order(model_dat[,2,drop=TRUE]),]
    Q <- model_dat[,1,drop=TRUE]
    h <- model_dat[,2,drop=TRUE]
    if(!is.null(c_param) && min(h)<c_param) stop('c_param must be lower than the minimum stage value in the data')
    if(any(Q<=0)) stop('All discharge measurements must but strictly greater than zero. If you know the stage of zero discharge, use c_param.')
    MCMC_output_list <- plm0.inference(y=log(Q),h=h,c_param=c_param,h_max=h_max,parallel=parallel,forcepoint=forcepoint,num_cores=num_cores)
    param_names <- get_param_names('plm0',c_param)
    result_obj=list()
    attr(result_obj, "class") <- "plm0"
    result_obj$a_posterior = MCMC_output_list$x[1,]
    result_obj$b_posterior = MCMC_output_list$x[2,]
    #refinement of list elements
    if(is.null(c_param)){
      result_obj$c_posterior <- MCMC_output_list$theta[1,]
      result_obj$sigma_eps_posterior <- MCMC_output_list$theta[2,]
    }else{
      result_obj$c_posterior <- NULL
      result_obj$sigma_eps_posterior <- MCMC_output_list$theta[1,]
    }
    unique_h_idx <- !duplicated(MCMC_output_list$h)
    h_unique <- unique(MCMC_output_list$h)
    h_unique_order <- order(h_unique)
    h_unique_sorted <- h_unique[h_unique_order]
    h_idx_data <- match(h,h_unique_sorted)
    result_obj$rating_curve_posterior <- exp(MCMC_output_list$y_post_pred[unique_h_idx,][h_unique_order,])
    result_obj$rating_curve_mean_posterior <- exp(MCMC_output_list$y_post[unique_h_idx,][h_unique_order,])
    result_obj$Deviance_posterior <- c(MCMC_output_list$D)
    #summary objects
    result_obj$rating_curve <- get_MCMC_summary(result_obj$rating_curve_posterior,h=h_unique_sorted)
    result_obj$rating_curve_mean <- get_MCMC_summary(result_obj$rating_curve_mean_posterior,h=h_unique_sorted)
    result_obj$param_summary <- get_MCMC_summary(rbind(MCMC_output_list$x[1,],MCMC_output_list$x[2,],MCMC_output_list$theta))
    result_obj$param_summary$n_eff_samples <- MCMC_output_list$num_effective_samples
    result_obj$param_summary$r_hat <- MCMC_output_list$r_hat
    row.names(result_obj$param_summary) <- param_names
    result_obj$Deviance_summary <- get_MCMC_summary(MCMC_output_list$D)
    #Deviance calculations
    result_obj$D_hat <- MCMC_output_list$D_hat
    result_obj$num_effective_param <- result_obj$Deviance_summary[,'median']-result_obj$D_hat
    result_obj$DIC <- result_obj$D_hat + 2*result_obj$num_effective_param
    #Rhat and autocorrelation
    autocorrelation_df <- as.data.frame(t(MCMC_output_list$autocorrelation))
    names(autocorrelation_df) <- param_names
    autocorrelation_df$lag <- 1:nrow(autocorrelation_df)
    result_obj$autocorrelation <- autocorrelation_df[,c('lag',param_names)]
    # store other information
    result_obj$acceptance_rate <- MCMC_output_list[['acceptance_rate']]
    result_obj$formula <- formula
    result_obj$data <- model_dat
    result_obj$run_info <- MCMC_output_list$run_info
    return(result_obj)
}

#' @importFrom stats optim
plm0.inference <- function(y,h,c_param=NULL,h_max=NULL,parallel=TRUE,forcepoint=rep(FALSE,length(h)),num_cores=NULL,num_chains=4,nr_iter=20000,burnin=2000,thin=5){
    RC <- get_model_components('plm0',y,h,c_param,h_max,forcepoint,h_min=NULL)
    output_list <- get_MCMC_output_list(theta_m=RC$theta_m,RC=RC,density_fun=RC$density_fun,
                                        unobserved_prediction_fun=RC$unobserved_prediction_fun,
                                        parallel=parallel,num_cores=num_cores,num_chains=num_chains,nr_iter=nr_iter,
                                        burnin=burnin,thin=thin)
    output_list$D_hat <- plm0.calc_Dhat(output_list$theta,RC)
    #refinement of list elements
    if(is.null(RC$c)){
      output_list$theta[1,] <- RC$h_min-exp(output_list$theta[1,])
      output_list$theta[2,] <- sqrt(exp(output_list$theta[2,]))
    }else{
      output_list$theta[1,] <- sqrt(exp(output_list$theta[1,]))
    }
    output_list$x[1,] <- exp(output_list$x[1,])
    output_list$h <- c(RC$h,RC$h_u)
    output_list$run_info <- list('c_param'=c_param,'h_max'=h_max,'forcepoint'=forcepoint,'nr_iter'=nr_iter,'num_chains'=num_chains,'burnin'=burnin,'thin'=thin)
    return(output_list)
}

#' @importFrom stats rnorm dlnorm
plm0.density_evaluation_known_c <- function(theta,RC){
    log_sig_eps2 <- theta[1]
    l=c(log(RC$h-RC$c))
    varr=RC$epsilon*exp(log_sig_eps2)
    if(any(varr>10^2)) return(list(p=-1e9)) # to avoid numerical instability
    Sig_eps=diag(varr)
    Sig_x=RC$Sig_ab

    X=cbind(rep(1,length(l)),l)
    L=t(chol(X%*%Sig_x%*%t(X)+Sig_eps+diag(nrow(Sig_eps))*RC$nugget))
    w=solve(L,RC$y-X%*%RC$mu_x)
    p=-0.5%*%t(w)%*%w-sum(log(diag(L)))+pri('sigma_eps2',log_sig_eps2 = log_sig_eps2,lambda_se=RC$lambda_se)

    W=solve(L,X%*%Sig_x)
    x_u=RC$mu_x+t(chol(Sig_x))%*%rnorm(length(RC$mu_x))
    sss=(X%*%x_u)-RC$y+sqrt(varr)*as.matrix(rnorm(RC$n))
    x=as.matrix(x_u-t(W)%*%solve(L,sss))
    yp=(X %*% x)[1:RC$n,]
    #posterior predictive draw
    ypo=yp+as.matrix(rnorm(RC$n))*sqrt(varr)
    D=-2*sum(log(dlnorm(exp(RC$y[1:RC$n,]),yp,sqrt(varr))))
    return(list("p"=p,"x"=x,"y_post"=yp,"y_post_pred"=ypo,"D"=D))
}

#' @importFrom stats rnorm dlnorm
plm0.density_evaluation_unknown_c <- function(theta,RC){
    zeta <- theta[1]
    log_sig_eps2 <- theta[2]
    l=c(log(RC$h-RC$h_min+exp(zeta)))
    varr=RC$epsilon*exp(log_sig_eps2)
    if(any(varr>10^2)) return(list(p=-1e9)) # to avoid numerical instability
    Sig_eps=diag(varr)
    Sig_x=RC$Sig_ab
    X=cbind(rep(1,length(l)),l)
    L=t(chol(X%*%Sig_x%*%t(X)+Sig_eps+diag(nrow(Sig_eps))*RC$nugget))
    w=solve(L,RC$y-X%*%RC$mu_x)
    p=-0.5%*%t(w)%*%w-sum(log(diag(L)))+
    pri('c',zeta = zeta,lambda_c = RC$lambda_c) +
    pri('sigma_eps2',log_sig_eps2 = log_sig_eps2,lambda_se=RC$lambda_se)

    W=solve(L,X%*%Sig_x)
    x_u=RC$mu_x+t(chol(Sig_x))%*%rnorm(length(RC$mu_x))
    sss=(X%*%x_u)-RC$y+sqrt(varr)*as.matrix(rnorm(RC$n))
    x=as.matrix(x_u-t(W)%*%solve(L,sss))
    yp=(X %*% x)[1:RC$n,]
    #posterior predictive draw
    ypo=yp+as.matrix(rnorm(RC$n))*sqrt(varr)
    D=-2*sum(log(dlnorm(exp(RC$y[1:RC$n,]),yp,sqrt(varr))))

    return(list("p"=p,"x"=x,"y_post"=yp,"y_post_pred"=ypo,"D"=D))
}

#' @importFrom stats dlnorm
plm0.calc_Dhat <- function(theta,RC){
  theta_median <- apply(theta,1,median)
  if(!is.null(RC$c)){
    theta_median <- c(log(RC$h_min-RC$c),theta_median)
  }
  zeta <- theta_median[1]
  log_sig_eps2 <- theta_median[2]
  l=c(log(RC$h-RC$h_min+exp(zeta)))
  varr=RC$epsilon*exp(log_sig_eps2)
  if(any(varr>10^2)) return(list(p=-1e9)) # to avoid numerical instability
  Sig_eps=diag(varr)
  Sig_x=RC$Sig_ab
  X=cbind(rep(1,length(l)),l)
  L=t(chol(X%*%Sig_x%*%t(X)+Sig_eps+diag(nrow(Sig_eps))*RC$nugget))
  w=solve(L,RC$y-X%*%RC$mu_x)
  x=RC$mu_x+Sig_x%*%(t(X)%*%solve(t(L),w))
  yp=(X %*% x)[1:RC$n,]
  D=-2*sum(log(dlnorm(exp(RC$y[1:RC$n,]),yp,sqrt(varr))))
  return(D)
}

#' @importFrom stats rnorm
plm0.predict_u_known_c <- function(theta,x,RC){
  #collecting parameters from the MCMC sample
  log_sig_eps2 = theta[1]
  m = length(RC$h_u)
  #building blocks of the explanatory matrix X calculated
  l=log(RC$h_u-RC$c)
  X=cbind(rep(1,m),l)
  #sample from the posterior predictive distr of y
  yp_u <- c(X%*%x)
  ypo_u <- yp_u + rnorm(m) * sqrt(exp(log_sig_eps2))
  return(list('y_post'=yp_u,'y_post_pred'=ypo_u))
}

#' @importFrom stats rnorm
plm0.predict_u_unknown_c <- function(theta,x,RC){
    #collecting parameters from the MCMC sample
    zeta=theta[1]
    log_sig_eps2 = theta[2]
    m=length(RC$h_u)
    above_c <- RC$h_min-exp(zeta) < RC$h_u
    m_above_c <- sum(above_c)
    #building blocks of the explanatory matrix X calculated
    l=log(RC$h_u[above_c]-RC$h_min+exp(zeta))
    X=cbind(rep(1,m_above_c),l)
    #sample from the posterior predictive distr of y
    yp_u <- c(X%*%x)
    ypo_u = yp_u + rnorm(m_above_c) * sqrt(exp(log_sig_eps2))
    return(list('y_post'=c(rep(-Inf,m-m_above_c),yp_u),'y_post_pred'=c(rep(-Inf,m-m_above_c),ypo_u)))
}




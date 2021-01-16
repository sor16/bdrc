#' Bayesian Generalized Power Law Model
#'
#' bgplm is used to fit a rating curve for paired measurements of stage and discharge using a Bayesian Power Law Model as described in Hrafnkelsson et al.
#'
#' @param formula an object of class "formula", with discharge column name as response and stage column name as a covariate.The details of model specification are given under "Details".
#' @param data data.frame containing the variables specified in formula
#' @param c_param stage for which there is zero discharge. If NULL, it is treated as unknown in the model and inferred from the data
#' @param h_max maximum stage to which the rating curve should extrapolate for. If NULL, the maximum stage value in data is selected as an upper bound.
#' @param forcepoint A boolean vector of the same length as the number of rows in data. If an element at index i is TRUE it indicates that the rating curve should be forced through the i-th measurement. Use with care, as this will strongly influence the resulting rating curve.
#'
#' @details Fits a discharge rating curve using a Bayesian generalized power-law model described in x. The generalized power-law model is of the form
#' \deqn{Q=a(h-c)^{f(h)}}
#' where \eqn{Q} is discharge, \eqn{h} is stage, \eqn{a} and \eqn{c} are unknown constants and \eqn{f} is a function of  \eqn{h} referred to as the generalized power-law exponent.\cr\cr
#' The Bayesian generalized power-law model is presented as a Bayesian hierarchical model. The function \eqn{f} is modelled at the latent level as a fixed constant b plus a continuous stochastic process which is assumed to be twice differentiable. The model is presented on a logarithmic scale
#' \deqn{log(Q_i) = log(a) + (b + \beta(h_i)) log(h_i - c) + \epsilon_i,     i = 1,...,n}
#' where \eqn{\epsilon_i} follows a normal distribution with mean zero and variance \eqn{\sigma_\epsilon(h_i)^{2}} that can vary with stage.
#' The stochastic process \eqn{\beta(h)} is assumed a priori to be a Gaussian process governed by a Matern covariance function with smoothness parameter \eqn{\nu = 2.5}, see Matern (1960).
#' An efficient posterior simulation is achieved by sampling from the joint posterior density of the hyperparameters of the model, and then sampling from the conditional density the latent parameters conditioned on the hyperparameters.\cr\cr
#' Bayesian inference is based on the posterior density and summary statistics such as the posterior mean and 95\% posterior intervals are based on the posterior density.
#' Analytical formulas for these summary statistics are intractable in most cases and thus they are computed by generating samples from the posterior density using Markov chain Monte Carlo simulation.
#'
#' @return
#' bgplm returns an object of class "bgplm". An object of class "bgplm" is a list containing the following components:
#' \describe{
#'  \item{\code{rating_curve}}{a data frame with 2.5\%, 50\% and 97.5\% quantiles of the posterior distribution of the rating curve.}
#'  \item{\code{rating_curve_mean}}{a data frame with 2.5\%, 50\% and 97.5\% quantiles of the posterior distribution of the mean of the rating curve.}
#'  \item{\code{param_summary}}{a data frame with 2.5\%, 50\% and 97.5\% quantiles of the posterior distribution of latent- and hyperparameters.}
#'  \item{\code{DIC_summary}}{a data frame with 2.5\%, 50\% and 97.5\% quantiles of the posterior distribution of the Deviance Information Criterion.}
#'  \item{\code{rating_curve_posterior}}{a matrix containing the full thinned posterior samples of the posterior distribution of the rating curve (excluding burn-in).}
#'  \item{\code{rating_curve_mean_posterior}}{a matrix containing the full thinned posterior samples of the posterior distribution of the mean of the rating curve (excluding burn-in).}
#'  \item{\code{a_posterior}}{a numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{a}.}
#'  \item{\code{b_posterior}}{a numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{b}.}
#'  \item{\code{c_posterior}}{a numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{c}.}
#'  \item{\code{sigma_beta_posterior}}{a numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{\sigma_\beta}.}
#'  \item{\code{phi_beta_posterior}}{a numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{\phi_\beta}.}
#'  \item{\code{sigma_eta_posterior}}{a numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{\sigma_\eta}.}
#'  \item{\code{eta_1_posterior}}{a numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{\eta_1}.}
#'  \item{\code{eta_2_posterior}}{a numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{\eta_2}.}
#'  \item{\code{eta_3_posterior}}{a numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{\eta_3}.}
#'  \item{\code{eta_4_posterior}}{a numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{\eta_4}.}
#'  \item{\code{eta_5_posterior}}{a numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{\eta_5}.}
#'  \item{\code{eta_6_posterior}}{a numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{\eta_6}.}
#'  \item{\code{f_posterior}}{a numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{f(h)}.}
#'  \item{\code{beta_posterior}}{a numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{\beta(h)}.}
#'  \item{\code{sigma_eps_posterior}}{a numeric vector containing the full thinned posterior samples of the posterior distribution of \eqn{\sigma_\epsilon(h)}.}
#'  \item{\code{DIC_posterior}}{a numeric vector containing the full thinned posterior samples of the posterior distribution of the Deviance Information Criterion.}
#'  \item{\code{f_summary}}{a data frame with 2.5\%, 50\% and 97.5\% quantiles of the posterior distribution of \eqn{f(h)}.}
#'  \item{\code{beta_summary}}{a data frame with 2.5\%, 50\% and 97.5\% quantiles of the posterior distribution of \eqn{\beta(h)}.}
#'  \item{\code{sigma_eps_summary}}{a data frame with 2.5\%, 50\% and 97.5\% quantiles of the posterior distribution of \eqn{\sigma_\epsilon(h)}.}
#'  \item{\code{Bayes_factor}}{a numeric value containing the rating curves Bayes factor.}
#'  \item{\code{formula}}{object of type "formula" provided by the user.}
#'  \item{\code{data}}{data provided by the user.}
#'  \item{\code{run_info}}{Information about the specific parameters used in the MCMC chain.}
#' }
#'
#' @references B. Hrafnkelsson, H. Sigurdarson, S.M. Gardarsson, 2020, Generalization of the power-law rating curve using hydrodynamic theory and Bayesian hierarchical modeling. arXiv preprint 2010.04769.
#'
#' @seealso \code{\link{summary.bgplm}} for summaries, \code{\link{predict.bgplm}} for prediction and \code{\link{plot.bgplm}} for plots. \code{\link{spread_draws}} and \code{\link{gather_draws}} are also useful to aid further visualization of the full posterior distributions.
#'
#' @examples
#' data(V316_river)
#' formula <- Q~W
#' bgplm.fit <- bgplm(formula,V316_river)
#' summary(bgplm.fit)
#' plot(bgplm.fit)
#' bgplm.fit_known_c <- bgplm(f,V316_river,c_param=0.75,h_max=2)
#' summary(bgplm.fit_known_c)
#' plot(bgplm.fit_known_c)
#' @export
bgplm <- function(formula,data,c_param=NULL,h_max=NULL,forcepoint=rep(FALSE,nrow(data)),...){
  #TODO:argument checking
  model_dat <- data[,all.vars(formula)]
  model_dat <- model_dat[order(model_dat[,2,drop=T]),]
  Q <- model_dat[,1,drop=T]
  h <- model_dat[,2,drop=T]
  MCMC_output_list <- bgplm.inference(y=log(Q),h,c_param,h_max,forcepoint,...)
  #prepare S3 model object to be returned
  result_obj=list()
  attr(result_obj, "class") <- "bgplm"
  result_obj$formula <- formula
  result_obj$data <- model_dat
  result_obj$a_posterior = MCMC_output_list$x[1,]
  result_obj$b_posterior = MCMC_output_list$x[2,]
  if(is.null(c_param)){
    result_obj$c_posterior <- MCMC_output_list$theta[1,]
    result_obj$sigma_beta_posterior <- MCMC_output_list$theta[2,]
    result_obj$phi_beta_posterior <- MCMC_output_list$theta[3,]
    result_obj$sigma_eta_posterior <- MCMC_output_list$theta[4,]
    for(i in 5:nrow(MCMC_output_list$theta)){
      result_obj[[paste0('eta_',i-4,'_posterior')]] <- MCMC_output_list$theta[i,]
    }
  }else{
    result_obj$c_posterior <- NULL
    result_obj$sigma_beta_posterior <- MCMC_output_list$theta[1,]
    result_obj$phi_beta_posterior <- MCMC_output_list$theta[2,]
    result_obj$sigma_eta_posterior <- MCMC_output_list$theta[3,]
    for(i in 4:nrow(MCMC_output_list$theta)){
      result_obj[[paste0('eta_',i-3,'_posterior')]] <- MCMC_output_list$theta[i,]
    }
  }
  unique_h_idx <- !duplicated(MCMC_output_list$h)
  h_unique <- unique(MCMC_output_list$h)
  h_unique_order <- order(h_unique)
  h_unique_sorted <- h_unique[h_unique_order]
  result_obj$rating_curve_posterior <- exp(MCMC_output_list$y_post_pred[unique_h_idx,][h_unique_order,])
  result_obj$rating_curve_mean_posterior <- exp(MCMC_output_list$y_post[unique_h_idx,][h_unique_order,])
  result_obj$beta_posterior <- MCMC_output_list$x[3:nrow(MCMC_output_list$x),][h_unique_order,]
  result_obj$f_posterior <- matrix(rep(result_obj$b_posterior,nrow(MCMC_output_list$x)-2),nrow=nrow(MCMC_output_list$x)-2,byrow=T)+MCMC_output_list$x[3:nrow(MCMC_output_list$x),]
  result_obj$sigma_eps_posterior <- sqrt(MCMC_output_list$sigma_eps[unique_h_idx,][h_unique_order,])
  result_obj$DIC_posterior <- MCMC_output_list$DIC
  #summary objects
  result_obj$rating_curve <- get_MCMC_summary(result_obj$rating_curve_posterior,h=h_unique_sorted)
  result_obj$rating_curve_mean <- get_MCMC_summary(result_obj$rating_curve_mean_posterior,h=h_unique_sorted)
  result_obj$beta_summary <- get_MCMC_summary(result_obj$beta_posterior,h=h_unique_sorted)
  result_obj$f_summary <- get_MCMC_summary(result_obj$f_posterior,h=h_unique_sorted)
  result_obj$sigma_eps_summary <- get_MCMC_summary(result_obj$sigma_eps_posterior,h=h_unique_sorted)
  result_obj$param_summary <- get_MCMC_summary(rbind(MCMC_output_list$x[1,],MCMC_output_list$x[2,],MCMC_output_list$theta))
  row.names(result_obj$param_summary) <- get_param_names('bgplm',c_param)
  result_obj$DIC_summary <- get_MCMC_summary(result_obj$DIC_posterior)
  result_obj$Bayes_factor <- 1/(exp(0.5*mean(result_obj$DIC_posterior)))
  result_obj$run_info <- MCMC_output_list$run_info
  return(result_obj)
}

bgplm.inference <- function(y,h,c_param=NULL,h_max=NULL,forcepoint=rep(FALSE,length(h)),num_chains=4,nr_iter=20000,burnin=2000,thin=5){
  #TODO: add error message if length(formula)! = 3 or if it contains more than one covariate. Also make sure that names in formula exist in data
  RC <- priors('bgplm',c_param)
  RC$y <- rbind(as.matrix(y),RC$mu_b)
  RC$h <- as.matrix(h)
  RC$h_min <- min(RC$h)
  RC$h_max <- max(RC$h)
  RC$h_tild <- RC$h-RC$h_min
  RC$h_unique <- unique(RC$h)
  RC$n <- length(RC$h)
  RC$n_unique <- length(RC$h_unique)
  RC$A <- create_A(RC$h)
  RC$dist <- as.matrix(stats::dist(c(RC$h_unique)))

  RC$mu_x <- as.matrix(c(RC$mu_a,RC$mu_b, rep(0,RC$n_unique)))
  RC$P <- lower.tri(matrix(rep(1,36),6,6),diag=T)*1
  RC$B <- B_splines(t(RC$h_tild)/RC$h_tild[length(RC$h_tild)])
  RC$epsilon <- rep(1,RC$n)

  RC$epsilon[forcepoint] <- 1/RC$n

  RC$Z <- cbind(t(c(0,1)),t(rep(0,RC$n_unique)))
  RC$m1 <- matrix(0,nrow=2,ncol=RC$n_unique)
  RC$m2 <- matrix(0,nrow=RC$n_unique,ncol=2)
  if(!is.null(RC$c)){
    if(RC$c<=RC$h_min){
      density_fun <- bgplm.density_evaluation_known_c
      unobserved_prediction_fun <- bgplm.predict_u_known_c
    }else{
      stop(paste0('the given c must be less than the lowest stage measurement, which is ',RC$h_min,' m'))
    }
  }else{
    density_fun <- bgplm.density_evaluation_unknown_c
    unobserved_prediction_fun <- bgplm.predict_u_unknown_c
  }
  #determine proposal density
  RC$theta_length <- if(is.null(RC$c)) 10 else 9
  theta_init <- rep(0,RC$theta_length)
  loss_fun  <-  function(theta) {-density_fun(theta,RC)$p}
  optim_obj <- stats::optim(par=theta_init,loss_fun,method="L-BFGS-B",hessian=TRUE)
  theta_m <- optim_obj$par
  H <- optim_obj$hessian
  RC$LH <- t(chol(H))/0.8
  h_min <- ifelse(is.null(RC$c),min(RC$h)-exp(theta_m[1]),RC$c)
  if(is.null(h_max)){
    h_max <- RC$h_max
  }
  if(h_max<RC$h_max){
    stop(paste0('maximum stage value must be larger than the maximum stage value in the data, which is ', RC$h_max,' m'))
  }
  RC$h_u <- h_unobserved(RC,h_min,h_max)
  RC$n_u <- length(RC$h_u)
  h_u_std <- ifelse(RC$h_u < RC$h_min,0.0,ifelse(RC$h_u>RC$h_max,1.0,(RC$h_u-RC$h_min)/(RC$h_max-RC$h_min)))
  RC$B_u <- B_splines(h_u_std)
  #determine length of each part of the output, in addition to theta
  RC$desired_output <- get_desired_output('bgplm',RC)
  #MCMC parameters added, number of iterations,burnin and thin
  if(num_chains>4){
    stop('Max number of chains is 4. Please pick a lower number of chains')
  }
  MCMC_output_list <- parallel::mclapply(1:num_chains,mc.cores=num_chains,FUN=function(i){
    run_MCMC(theta_m,RC,density_fun,unobserved_prediction_fun,nr_iter,num_chains,burnin,thin)
  })
  output_list <- list()
  for(elem in names(MCMC_output_list[[1]])){
    output_list[[elem]] <- do.call(cbind,lapply(1:num_chains,function(i) MCMC_output_list[[i]][[elem]]))
  }
  # #refinement of list elements
  if(is.null(RC$c)){
    output_list$theta[1,] <- RC$h_min-exp(output_list$theta[1,])
    output_list$theta[2,] <- exp(output_list$theta[2,])
    output_list$theta[3,] <- exp(output_list$theta[3,])
    output_list$theta[4,] <- exp(output_list$theta[4,])
    output_list$theta[5:10,] <- RC$P%*%output_list$theta[5:10,]
  }else{
    output_list$theta[1,] <- exp(output_list$theta[1,])
    output_list$theta[2,] <- exp(output_list$theta[2,])
    output_list$theta[3,] <- exp(output_list$theta[3,])
    output_list$theta[4:9,] <- RC$P%*%output_list$theta[4:9,]
  }
  output_list$x[1,] <- exp(output_list$x[1,])
  output_list[['h']] <- c(RC$h,RC$h_u)
  output_list[['run_info']] <- list('c_param'=c_param,'h_max'=h_max,'forcepoint'=forcepoint,'nr_iter'=nr_iter,'num_chains'=num_chains,'burnin'=burnin,'thin'=thin)
  return(output_list)
}

bgplm.density_evaluation_known_c <- function(theta,RC){
  log_sig_b <- theta[1]
  log_phi_b <- theta[2]
  log_sig_eta <- theta[3]
  eta_1 <- theta[4]
  z <- theta[5:9]

  eta=c(RC$P%*%as.matrix(c(eta_1,exp(log_sig_eta)*z)))

  l=c(log(RC$h-RC$c))

  varr=c(RC$epsilon*exp(RC$B%*%eta))
  Sig_eps=diag(c(varr,0))
  #Matern covariance
  R_Beta=(1+sqrt(5)*RC$dist/exp(log_phi_b)+5*RC$dist^2/(3*exp(log_phi_b)^2))*exp(-sqrt(5)*RC$dist/exp(log_phi_b))+diag(RC$n_unique)*RC$nugget
  Sig_x=rbind(cbind(RC$Sig_ab,RC$m1),cbind(RC$m2,exp(2*log_sig_b)*R_Beta))

  X=rbind(cbind(1,l,diag(l)%*%RC$A),RC$Z)
  L=t(chol(X%*%Sig_x%*%t(X)+Sig_eps+diag(nrow(Sig_eps))*RC$nugget))
  w=solve(L,RC$y-X%*%RC$mu_x)
  p=-0.5%*%t(w)%*%w-sum(log(diag(L))) +
    pri('sigma_b',log_sig_b = log_sig_b, lambda_sb = RC$lambda_sb) +
    pri('phi_b', log_phi_b = log_phi_b, lambda_pb = RC$lambda_pb) +
    pri('eta_1',eta_1=eta_1,lambda_eta_1=RC$lambda_eta_1) +
    pri('eta_minus1',z=z) +
    pri('sigma_eta',log_sig_eta=log_sig_eta,lambda_seta=RC$lambda_seta)


  W=solve(L,X%*%Sig_x)
  x_u=RC$mu_x+t(chol(Sig_x))%*%stats::rnorm(RC$n_unique+2)
  sss=(X%*%x_u)-RC$y+rbind(sqrt(varr)*as.matrix(stats::rnorm(RC$n)),0)
  x=as.matrix(x_u-t(W)%*%solve(L,sss))
  yp=(X %*% x)[1:RC$n,]
  #posterior predictive draw
  ypo=yp+as.matrix(stats::rnorm(RC$n))*sqrt(varr)
  D=-2*sum(log(stats::dlnorm(exp(RC$y[1:RC$n,]),yp,sqrt(varr))))

  return(list("p"=p,"x"=x,"y_post"=yp,"y_post_pred"=ypo,"sigma_eps"=varr,"DIC"=D))
}

bgplm.density_evaluation_unknown_c <- function(theta,RC){
  zeta <- theta[1]
  log_sig_b <- theta[2]
  log_phi_b <- theta[3]
  log_sig_eta <- theta[4]
  eta_1 <- theta[5]
  z <- theta[6:10]

  eta=c(RC$P%*%as.matrix(c(eta_1,exp(log_sig_eta)*z)))

  l=c(log(RC$h_tild+exp(zeta)))

  varr=c(RC$epsilon*exp(RC$B%*%eta))
  Sig_eps=diag(c(varr,0))
  #Matern covariance
  R_Beta=(1+sqrt(5)*RC$dist/exp(log_phi_b)+5*RC$dist^2/(3*exp(log_phi_b)^2))*exp(-sqrt(5)*RC$dist/exp(log_phi_b))+diag(RC$n_unique)*RC$nugget
  Sig_x=rbind(cbind(RC$Sig_ab,RC$m1),cbind(RC$m2,exp(2*log_sig_b)*R_Beta))

  X=rbind(cbind(1,l,diag(l)%*%RC$A),RC$Z)
  L=t(chol(X%*%Sig_x%*%t(X)+Sig_eps+diag(nrow(Sig_eps))*RC$nugget))
  w=solve(L,RC$y-X%*%RC$mu_x)

  p=-0.5%*%t(w)%*%w-sum(log(diag(L)))+
    pri('c',zeta = zeta,lambda_c = RC$lambda_c) +
    pri('sigma_b',log_sig_b = log_sig_b, lambda_sb = RC$lambda_sb) +
    pri('phi_b', log_phi_b = log_phi_b, lambda_pb = RC$lambda_pb) +
    pri('eta_1',eta_1=eta_1,lambda_eta_1=RC$lambda_eta_1) +
    pri('eta_minus1',z=z) +
    pri('sigma_eta',log_sig_eta=log_sig_eta,lambda_seta=RC$lambda_seta)

  W=solve(L,X%*%Sig_x)
  x_u=RC$mu_x+t(chol(Sig_x))%*%stats::rnorm(RC$n_unique+2)
  sss=(X%*%x_u)-RC$y+rbind(sqrt(varr)*as.matrix(stats::rnorm(RC$n)),0)
  x=as.matrix(x_u-t(W)%*%solve(L,sss))
  yp=(X %*% x)[1:RC$n,]
  #posterior predictive draw
  ypo=yp+as.matrix(stats::rnorm(RC$n))*sqrt(varr)

  D=-2*sum(log(stats::dlnorm(exp(RC$y[1:RC$n,]),yp,sqrt(varr))))

  return(list("p"=p,"x"=x,"y_post"=yp,"y_post_pred"=ypo,"sigma_eps"=varr,"DIC"=D))
}

bgplm.predict_u_known_c <- function(theta,x,RC){
  #store particular hyperparameter values
  sig_b <- exp(theta[1])
  phi_b <- exp(theta[2])
  log_sig_eta <- theta[3]
  eta_1 <- theta[4]
  z <- theta[5:9]
  eta <- c(RC$P%*%as.matrix(c(eta_1,exp(log_sig_eta)*z)))

  n=RC$n_unique
  m=RC$n_u
  #get sample of data variance using splines
  varr_u = c(exp(RC$B_u %*% eta))
  #combine stages from data with unobserved stages
  h_all=c(RC$h_unique,RC$h_u)
  #calculating distance matrix for h_all
  dist_mat=as.matrix(stats::dist(h_all))
  #Covariance of the joint prior for betas from data and beta unobserved.
  #Matern covariance formula used for v=5/2
  sigma_all=sig_b^2*(1 + sqrt(5)*dist_mat/phi_b+(5*dist_mat^2)/(3*phi_b^2))*exp(-sqrt(5)*dist_mat/phi_b) + diag(length(h_all))*RC$nugget
  sigma_11=sigma_all[1:n,1:n]
  sigma_22=sigma_all[(n+1):(m+n),(n+1):(m+n)]
  sigma_12=sigma_all[1:n,(n+1):(n+m)]
  sigma_21=sigma_all[(n+1):(n+m),1:n]
  #parameters for the posterior of beta_u
  mu_x_u=sigma_21%*%solve(sigma_11,x[3:length(x)])
  Sigma_x_u=(sigma_22-sigma_21%*%solve(sigma_11,sigma_12))
  #a sample from posterior of beta_u drawn
  beta_u=as.numeric(mu_x_u) + stats::rnorm(ncol(Sigma_x_u)) %*% chol(Sigma_x_u)
  #buidling blocks of the explanatory matrix X calculated
  l=log(RC$h_u-RC$c)
  X=cbind(rep(1,m),l,diag(l))
  x_u=c(x[1:2],beta_u)
  #sample from the posterior of discharge y
  yp_u <- c(X%*%x_u)
  #make sure the log discharge at point of zero discharge is -Inf
  #yp_u[1] <- -Inf
  ypo_u = yp_u + stats::rnorm(m) * sqrt(varr_u)
  return(list('x'=beta_u,'sigma_eps'=varr_u,'y_post'=yp_u,'y_post_pred'=ypo_u))
}

bgplm.predict_u_unknown_c <- function(theta,x,RC){
  #store particular hyperparameter values
  zeta <- theta[1]
  sig_b <- exp(theta[2])
  phi_b <- exp(theta[3])
  log_sig_eta <- theta[4]
  eta_1 <- theta[5]
  z <- theta[6:10]
  eta <- c(RC$P%*%as.matrix(c(eta_1,exp(log_sig_eta)*z)))

  n=RC$n_unique
  m=RC$n_u
  #get sample of data variance using splines
  varr_u = c(exp(RC$B_u %*% eta))
  #combine stages from data with unobserved stages
  h_all=c(RC$h_unique,RC$h_u)
  #calculating distance matrix for h_all
  dist_mat=as.matrix(stats::dist(h_all))
  #Covariance of the joint prior for betas from data and beta unobserved.
  #Matern covariance formula used for v=5/2
  sigma_all=sig_b^2*(1 + sqrt(5)*dist_mat/phi_b+(5*dist_mat^2)/(3*phi_b^2))*exp(-sqrt(5)*dist_mat/phi_b) + diag(length(h_all))*RC$nugget
  sigma_11=sigma_all[1:n,1:n]
  sigma_22=sigma_all[(n+1):(m+n),(n+1):(m+n)]
  sigma_12=sigma_all[1:n,(n+1):(n+m)]
  sigma_21=sigma_all[(n+1):(n+m),1:n]
  #parameters for the posterior of beta_u
  mu_x_u=sigma_21%*%solve(sigma_11,x[3:length(x)])
  Sigma_x_u=(sigma_22-sigma_21%*%solve(sigma_11,sigma_12))
  #a sample from posterior of beta_u drawn
  beta_u=as.numeric(mu_x_u) + stats::rnorm(ncol(Sigma_x_u)) %*% chol(Sigma_x_u)
  above_c <- -(exp(zeta)-RC$h_min) < RC$h_u
  m_above_c <- sum(above_c)
  #buidling blocks of the explanatory matrix X calculated
  l=log(RC$h_u[above_c]-RC$h_min+exp(zeta))
  X=cbind(rep(1,m_above_c),l,diag(l))
  #vector of parameters
  x_u=c(x[1:2],beta_u[above_c])
  #sample from the posterior of discharge y
  yp_u <- c(X%*%x_u)
  ypo_u = yp_u + stats::rnorm(m_above_c) * sqrt(varr_u[above_c])
  return(list('x'=beta_u,'sigma_eps'=varr_u,'y_post'=c(rep(-Inf,m-m_above_c),yp_u),'y_post_pred'=c(rep(-Inf,m-m_above_c),ypo_u)))
}

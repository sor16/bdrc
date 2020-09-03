#' Bayesian Power Law Model with constant variance
#'
#' Infers a rating curve for paired measurements of stage and discharge using a  power law model described in Hrafnkelsson et al.
#'@param formula formula with name of discharge column in data as response and name of stage column in data as the single covariate.
#'@param data data.frame containing the columns in formula
#'@param w_limits vector of length 2 setting the lower and upper bound of stage values at which a rating curve should be predicted. If NULL, the known value of c or the mle of c will be used as lower bound (depending on the value of the input parameter c) and maximum stage value in data as upper bound.
#'@param country Name of the country the prior parameters should be defined for, default value is "Iceland".
#'@param Wmin Positive numeric value for the lowest stage the user wants to calculate a rating curve. If input is an empty string (default) Wmin will
#'automatically be set to c_hat.
#'@param Wmax Positive numeric value for the highest stage the user wants to calculate a rating curve. If input is an empty string (default) Wmax will
#'automatically be set to the maximum stage of the data.
#'@return List containing information on the calculated rating curve,
#'the data frames observedData, betaData, completePrediction, observedPrediction, TableOfData, FitTable, LowerTable, UpperTable, plotTable.
#'@references Birgir Hrafnkelsson, Helgi Sigurdarson and Sigurdur M. Gardarson (2015) \emph{Bayesian Generalized Rating Curves}
#'@seealso \code{\link{clean}}

bplm0 <- function(formula,data,c_param=NULL,w_limits=NULL,forcepoint=rep(FALSE,nrow(data)),...){
    #argument checking
    model_dat <- data[,all.vars(formula)]
    model_dat <- model_dat[order(model_dat[,2,drop=T]),]
    Q <- model_dat[,1,drop=T]
    w <- model_dat[,2,drop=T]
    MCMC_output_list <- bplm0.inference(y=log(Q),w=w,c_param,w_limits,forcepoint,...)
    result_obj=list()
    attr(result_obj, "class") <- "bplm0"
    result_obj$formula <- formula
    result_obj$data <- model_dat
    result_obj$a_posterior = MCMC_output_list$x[1,]
    result_obj$b_posterior = MCMC_output_list$x[2,]
    if(is.null(c_param)){
      result_obj$c_posterior <- MCMC_output_list$theta[1,]
      result_obj$sig_eps_posterior <- MCMC_output_list$theta[2,]

    }else{
      result_obj$c_posterior <- NULL
      result_obj$sig_eps_posterior <- MCMC_output_list$theta[1,]
    }
    result_obj$Q_posterior_predictive <- exp(MCMC_output_list$y_post_pred)
    result_obj$Q_posterior <- exp(MCMC_output_list$y_post)
    result_obj$DIC_posterior <- MCMC_output_list$DIC
    #summary objects
    result_obj$rating_curve <- get_MCMC_summary(result_obj$Q_posterior_predictive,w=MCMC_output_list$w)
    result_obj$rating_curve_mean <- get_MCMC_summary(result_obj$Q_posterior,w=MCMC_output_list$w)
    result_obj$param_summary <- get_MCMC_summary(rbind(MCMC_output_list$x[1,],MCMC_output_list$x[2,],MCMC_output_list$theta))
    row.names(result_obj$param_summary) <- get_param_names('bplm0',c_param)
    result_obj$DIC_summary <- get_MCMC_summary(result_obj$DIC_posterior)
    return(result_obj)
}

bplm0.inference <- function(y,w,c_param=NULL,w_limits=NULL,forcepoint=rep(FALSE,nrow(data)),num_chains=4,nr_iter=20000,burnin=2000,thin=5){
    require(parallel)
    RC <- priors('bplm0',c_param)
    RC$y <- as.matrix(y)
    RC$w <- w
    RC$w_min <- min(RC$w)
    RC$w_max <- max(RC$w)
    RC$w_tild <- RC$w-RC$w_min
    RC$n <- length(w)
    forcepoint_dat=model_dat[forcepoint,]
    RC$epsilon[forcepoint]=1/RC$n
    if(!is.null(RC$c)){
        density_fun <- bplm0.density_evaluation_known_c
        unobserved_prediction_fun <- bplm0.predict_u_known_c
    }else{
        density_fun <- bplm0.density_evaluation_unknown_c
        unobserved_prediction_fun <- bplm0.predict_u_unknown_c
    }
    #determine proposal density
    RC$theta_length <- if(is.null(RC$c)) 2 else 1
    theta_init <- rep(0,RC$theta_length)
    loss_fun = function(theta) {-density_fun(theta,RC)$p}
    optim_obj=optim(par=theta_init,loss_fun,method="L-BFGS-B",hessian=TRUE)
    theta_m =optim_obj$par
    H=optim_obj$hessian
    RC$LH=t(chol(H))/(2.38/sqrt(2))
    #create equally spaced grid of stages
    if(is.null(w_limits)){ # if not user defined, round w_min and w_max to nearest decimeter
      w_max <- ceiling(RC$w_max*10)/10
      w_min <- ceiling(10*ifelse(is.null(RC$c),RC$w_min-exp(theta_m[1]),RC$c))/10
    }else{
      w_min <- w_limits[1]
      w_max <- w_limits[2]
    }
    RC$w_u <- W_unobserved(RC,w_min,w_max)
    RC$n_u <- length(RC$w_u)
    #determine length of each part of the output, in addition to theta
    RC$desired_output <- get_desired_output('bplm0',RC)
    #MCMC parameters added, number of iterations,burnin and thin
    if(num_chains>4){
      stop('Max number of chains is 4. Please pick a lower number of chains')
    }
    MCMC_output_list <- mclapply(1:num_chains,mc.cores=num_chains,FUN=function(i){
      run_MCMC(theta_m,RC,density_fun,unobserved_prediction_fun,nr_iter,num_chains,burnin,thin)
    })
    output_list <- list()
    for(elem in names(MCMC_output_list[[1]])){
      output_list[[elem]] <- do.call(cbind,lapply(1:num_chains,function(i) MCMC_output_list[[i]][[elem]]))
    }
    #refinement of list elements
    if(is.null(RC$c)){
      output_list$theta[1,] <- RC$w_min-exp(output_list$theta[1,])
      output_list$theta[2,] <- sqrt(exp(output_list$theta[2,]))
    }else{
      output_list$theta[1,] <- sqrt(exp(output_list$theta[1,]))
    }
    output_list$x[1,] <- exp(output_list$x[1,])
    output_list[['w']] <- c(RC$w,RC$w_u)
    return(output_list)
}

#'Density evaluation for model2
#'
#'Evaluates the log density of the posterior distribution of the parameters of .
#'@param theta A vector with length 9 containing parameters
#'@param RC A list containing prior parameters, matrices and the data.
#'@return Returns a list containing predictive values of the parameters drawn out of the evaluated density.
#'@references Birgir Hrafnkelsson, Helgi Sigurdarson and Sigurdur M. Gardarson (2015) \emph{Bayesian Generalized Rating Curves}
bplm0.density_evaluation_known_c <- function(theta,RC){
    log_sig_eps2 <- theta[1]
    l=c(log(RC$w-RC$c))
    X=cbind(rep(1,length(l)),l)
    L=t(chol(RC$Sig_xinv + (t(X) %*% X)/exp(log_sig_eps2)))
    q=solve(L,(RC$Sinvmu+t(X)%*%RC$y/exp(log_sig_eps2)))
    p=0.5*sum(q^2)+log(L[1,1])+log(L[2,2])-
     0.5*sum(RC$y^2)/exp(log_sig_eps2)-RC$n*log_sig_eps2/2 - log_sig_eps2
    x=solve(t(L),(q+as.matrix(rnorm(2))))
    yp=X %*% x
    ypo=yp+as.matrix(rnorm(RC$n,sd=sqrt(exp(log_sig_eps2))))
    D=-2*sum(log(dlnorm(exp(RC$y),X%*%x,sqrt(exp(log_sig_eps2)))))

    return(list("p"=p,"x"=x,"y_post"=yp,"y_post_pred"=ypo,"DIC"=D))
}

#'Density evaluation for model2
#'
#'Evaluates the log density of the posterior distribution of the parameters of model2BH.
#'@param theta A vector with length 9 containing parameters
#'@param RC A list containing prior parameters, matrices and the data.
#'@return Returns a list containing predictive values of the parameters drawn out of the evaluated density.
#'@references Birgir Hrafnkelsson, Helgi Sigurdarson and Sigurdur M. Gardarson (2015) \emph{Bayesian Generalized Rating Curves}
bplm0.density_evaluation_unknown_c <- function(theta,RC){
    zeta <- theta[1]
    log_sig_eps2 <- theta[2]
    l=c(log(RC$w-RC$w_min+exp(zeta)))
    X=cbind(rep(1,length(l)),l)
    L=t(chol(RC$Sig_xinv + (t(X) %*% X)/exp(log_sig_eps2)))
    q=solve(L,(RC$Sinvmu+t(X)%*%RC$y/exp(log_sig_eps2)))

    p=0.5*sum(q^2)+log(L[1,1])+log(L[2,2])-
    0.5*sum(RC$y^2)/exp(log_sig_eps2)-RC$n*log_sig_eps2/2 +
    zeta-exp(zeta)*RC$mu_c

    x=solve(t(L),(q+as.matrix(rnorm(2))))
    yp=X %*% x
    ypo=yp+as.matrix(rnorm(RC$n,sd=sqrt(exp(log_sig_eps2))))
    D=-2*sum(log(dlnorm(exp(RC$y),X%*%x,sqrt(exp(log_sig_eps2)))))

    return(list("p"=p,"x"=x,"y_post"=yp,"y_post_pred"=ypo,"DIC"=D))
}

#' Predictive values for unoberved stages
#'
#'Calculates predictive values for unobserved stages
#'@param param A vector of samples of theta and samples of betas from MCMC. Theta is a vector containing c (stage at which discharge is zero), two hyperparameters sig_b^2 and phi_b
#'and six lambda parameters that affect the variance through the Bspline functions.
#'@param RC A list containing prior parameters, matrices and the data which are calculated in \code{\link{model2BH}}
#'
#'@return
#'\itemize{
#'\item Vector containing predictive values ypo and values of beta for every stage measurement.
#'}
#'@references Birgir Hrafnkelsson, Helgi Sigurdarson and Sigurdur M. Gardarson (2015) \emph{Bayesian Generalized Rating Curves}
bplm0.predict_u_known_c <- function(theta,x,RC){
  #collecting parameters from the MCMC sample
  log_sig_eps2 = theta[1]
  m = length(RC$w_u)
  #building blocks of the explanatory matrix X calculated
  l=log(RC$w_u-RC$c)
  X=cbind(rep(1,m),l)
  #sample from the posterior predictive distr of y
  yp_u <- X%*%x
  ypo_u <- yp_u + as.matrix(rnorm(m)) * sqrt(exp(log_sig_eps2))
  return(list('y_post'=yp_u,'y_post_pred'=ypo_u))
}

#' Predictive values for unoberved stages
#'
#'Calculates predictive values for unobserved stages
#'@param param A vector of samples of theta and samples of betas from MCMC. Theta is a vector containing c (stage at which discharge is zero), two hyperparameters sig_b^2 and phi_b
#'and six lambda parameters that affect the variance through the Bspline functions.
#'@param RC A list containing prior parameters, matrices and the data which are calculated in \code{\link{model2BH}}
#'
#'@return
#'\itemize{
#'\item Vector containing predictive values ypo and values of beta for every stage measurement.
#'}
#'@references Birgir Hrafnkelsson, Helgi Sigurdarson and Sigurdur M. Gardarson (2015) \emph{Bayesian Generalized Rating Curves}
bplm0.predict_u_unknown_c <- function(theta,x,RC){
    #collecting parameters from the MCMC sample
    zeta=theta[1]
    log_sig_eps2 = theta[2]
    m=length(RC$w_u)
    above_c <- RC$w_min-exp(zeta) < RC$w_u
    m_above_c <- sum(above_c)
    #building blocks of the explanatory matrix X calculated
    l=log(RC$w_u[above_c]-RC$w_min+exp(zeta))
    X=cbind(rep(1,m_above_c),l)
    #sample from the posterior predictive distr of y
    yp_u <- X%*%x
    ypo_u = yp_u + as.matrix(rnorm(m_above_c)) * sqrt(exp(log_sig_eps2))
    return(list('y_post'=c(rep(-Inf,m-m_above_c),yp_u),'y_post_pred'=c(rep(-Inf,m-m_above_c),ypo_u)))
}




#' Generalized Bayesian Power Law Model
#'
#' Infers a rating curve for paired measurements of stage and discharge using a generalized power law model described in Hrafnkelsson et al.
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

bgplm <- function(formula,data,c_param=NULL,w_limits=NULL,forcepoint=rep(FALSE,nrow(data)),...){
  #TODO:argument checking
  model_dat <- data[,all.vars(formula)]
  model_dat <- model_dat[order(model_dat[,2,drop=T]),]
  Q <- model_dat[,1,drop=T]
  w <- model_dat[,2,drop=T]
  MCMC_output_list <- bgplm.inference(y=log(Q),w,c_param,w_limits,forcepoint,...)
  #prepare S3 model object to be returned
  result_obj=list()
  attr(result_obj, "class") <- "bgplm"
  result_obj$formula <- formula
  result_obj$data <- model_dat
  result_obj$a_posterior = MCMC_output_list$x[1,]
  result_obj$b_posterior = MCMC_output_list$x[2,]
  if(is.null(c_param)){
    result_obj$c_posterior <- MCMC_output_list$theta[1,]
    result_obj$sig_beta_posterior <- MCMC_output_list$theta[2,]
    result_obj$phi_beta_posterior <- MCMC_output_list$theta[3,]
    result_obj$lambda_posterior <- MCMC_output_list$theta[4:nrow(MCMC_output_list$theta),]
  }else{
    result_obj$c_posterior <- NULL
    result_obj$sig_beta_posterior <- MCMC_output_list$theta[1,]
    result_obj$phi_beta_posterior <- MCMC_output_list$theta[2,]
    result_obj$lambda_posterior <- MCMC_output_list$theta[3:nrow(MCMC_output_list$theta),]
  }
  result_obj$Q_posterior_predictive <- exp(MCMC_output_list$y_post_pred)
  result_obj$Q_posterior <- exp(MCMC_output_list$y_post)
  result_obj$beta_posterior <- matrix(rep(result_obj$b_posterior,nrow(MCMC_output_list$x)-2),nrow=nrow(MCMC_output_list$x)-2,byrow=T)+MCMC_output_list$x[3:nrow(MCMC_output_list$x),]
  result_obj$sigma_eps_posterior <- sqrt(exp(MCMC_output_list$sigma_eps))
  result_obj$DIC_posterior <- MCMC_output_list$DIC
  #summary objects
  result_obj$rating_curve <- get_MCMC_summary(result_obj$Q_posterior_predictive,w=MCMC_output_list$w)
  result_obj$rating_curve_mean <- get_MCMC_summary(result_obj$Q_posterior,w=MCMC_output_list$w)
  result_obj$beta_summary <- get_MCMC_summary(result_obj$beta_posterior,w=unique(MCMC_output_list$w))
  result_obj$sigma_eps_summary <- get_MCMC_summary(result_obj$sigma_eps_posterior,w=MCMC_output_list$w)
  result_obj$param_summary <- get_MCMC_summary(rbind(MCMC_output_list$x[1,],MCMC_output_list$x[2,],MCMC_output_list$theta))
  row.names(result_obj$param_summary) <- get_param_names('bgplm',c_param)
  result_obj$DIC_summary <- get_MCMC_summary(result_obj$DIC_posterior)
  return(result_obj)
}

bgplm.inference <- function(y,w,c_param=NULL,w_limits=NULL,forcepoint=rep(FALSE,nrow(data)),num_chains=4,nr_iter=20000,burnin=2000,thin=5){
  suppressPackageStartupMessages(require(doParallel))
  #TODO: add error message if length(formula)! = 3 or if it contains more than one covariate. Also make sure that names in formula exist in data
  RC <- priors('bgplm',c_param)
  RC$y <- rbind(as.matrix(y),0)
  RC$w <- as.matrix(w)
  RC$w_min <- min(RC$w)
  RC$w_max <- max(RC$w)
  RC$w_tild <- RC$w-RC$w_min
  RC$w_unique <- unique(RC$w)
  RC$n <- length(RC$w)
  RC$n_unique <- length(RC$w_unique)
  RC$A <- create_A(RC$w)
  RC$dist <- as.matrix(dist(c(RC$w_unique)))

  RC$mu_x <- as.matrix(c(RC$mu_a,RC$mu_b, rep(0,RC$n_unique)))
  RC$P <- diag(nrow=5,ncol=5,6)-matrix(nrow=5,ncol=5,1)
  RC$B <- B_splines(t(RC$w_tild)/RC$w_tild[length(RC$w_tild)])
  RC$epsilon <- rep(1,RC$n)
  #Spyrja Bigga út í varíans hér
  RC$epsilon[forcepoint] <- 1/RC$n

  RC$Z <- cbind(t(rep(0,2)),t(rep(1,RC$n_unique)))
  RC$m1 <- matrix(0,nrow=2,ncol=RC$n_unique)
  RC$m2 <- matrix(0,nrow=RC$n_unique,ncol=2)
  if(!is.null(RC$c)){
    density_fun <- bgplm.density_evaluation_known_c
    unobserved_prediction_fun <- bgplm.predict_u_known_c
  }else{
    density_fun <- bgplm.density_evaluation_unknown_c
    unobserved_prediction_fun <- bgplm.predict_u_unknown_c
  }
  #determine proposal density
  RC$theta_length <- if(is.null(RC$c)) 9 else 8
  theta_init <- rep(0,RC$theta_length)
  loss_fun  <-  function(theta) {-density_fun(theta,RC)$p}
  optim_obj <- optim(par=theta_init,loss_fun,method="L-BFGS-B",hessian=TRUE)
  theta_m <- optim_obj$par
  H <- optim_obj$hessian
  RC$LH <- t(chol(H))/0.8

  #make Wmin and Wmax divisable by 10 up, both in order to make rctafla and so l_m is defined
  if(is.null(w_limits)){
    w_max <- ceiling(max(RC$w)*10)/10
    w_min <- ceiling(10*ifelse(is.null(RC$c),min(RC$w)-exp(theta_m[1]),RC$c))/10
  }else{
    w_min <- w_limits[1]
    w_max <- w_limits[2]
  }
  RC$w_u <- W_unobserved(RC,w_min,w_max)
  RC$n_u <- length(RC$w_u)
  w_u_std <- ifelse(RC$w_u < RC$w_min,0.0,ifelse(RC$w_u>RC$w_max,1.0,(RC$w_u-RC$w_min)/(RC$w_max-RC$w_min)))
  RC$B_u <- B_splines(w_u_std)
  #determine length of each part of the output, in addition to theta
  RC$desired_output <- get_desired_output('bgplm',RC)
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
    output_list$theta[2,] <- exp(output_list$theta[2,])
    output_list$theta[3,] <- exp(output_list$theta[3,])
  }else{
    output_list$theta[1,] <- exp(output_list$theta[1,])
    output_list$theta[2,] <- exp(output_list$theta[2,])
  }
  output_list$x[1,] <- exp(output_list$x[1,])
  output_list[['w']] <- c(RC$w,RC$w_u)
  return(output_list)
}



#create a predict method for interpolation of posterior predictive

#'Density evaluation for model2
#'
#'Evaluates the log density of the posterior distribution of the parameters of .
#'@param th A vector with length 9 containing parameters
#'@param RC A list containing prior parameters, matrices and the data.
#'@return Returns a list containing predictive values of the parameters drawn out of the evaluated density.
#'@references Birgir Hrafnkelsson, Helgi Sigurdarson and Sigurdur M. Gardarson (2015) \emph{Bayesian Generalized Rating Curves}
bgplm.density_evaluation_known_c <- function(theta,RC){
  sig_b2=theta[1]
  phi_b=theta[2]
  lambda=theta[3:8]

  f=lambda[1:5]-lambda[6]
  l=c(log(RC$w-RC$c))

  varr=c(RC$epsilon*exp(RC$B%*%lambda))
  Sig_eps=diag(c(varr,0))
  #Matern covariance
  R_Beta=(1+sqrt(5)*RC$dist/exp(phi_b)+5*RC$dist^2/(3*exp(phi_b)^2))*exp(-sqrt(5)*RC$dist/exp(phi_b))+diag(RC$n_unique)*RC$nugget
  Sig_x=rbind(cbind(RC$Sig_ab,RC$m1),cbind(RC$m2,exp(sig_b2)*R_Beta))

  X=rbind(cbind(1,l,diag(l)%*%RC$A),RC$Z)
  L=t(chol(X%*%Sig_x%*%t(X)+Sig_eps))
  w=solve(L,RC$y-X%*%RC$mu_x)
  p=-0.5%*%t(w)%*%w-sum(log(diag(L)))-
    (RC$v+5-1)/2*log(RC$v*RC$s+f%*%RC$P%*%f)+
    sig_b2-exp(sig_b2)/RC$mu_sb-0.5/RC$tau_pb2*(phi_b-RC$mu_pb)^2

  W=solve(L,X%*%Sig_x)
  x_u=RC$mu_x+t(chol(Sig_x))%*%rnorm(RC$n_unique+2)
  sss=(X%*%x_u)-RC$y+rbind(sqrt(varr)*as.matrix(rnorm(RC$n)),0)
  x=as.matrix(x_u-t(W)%*%solve(L,sss))
  yp=(X %*% x)[1:RC$n,]
  #posterior predictive draw
  ypo=yp+as.matrix(rnorm(RC$n))*sqrt(varr)
  D=-2*sum(log(dlnorm(exp(RC$y[1:RC$n,]),yp,sqrt(varr))))

  return(list("p"=p,"x"=x,"y_post"=yp,"y_post_pred"=ypo,"sigma_eps"=varr,"DIC"=D))
}

#'Density evaluation for model2
#'
#'Evaluates the log density of the posterior distribution of the parameters of model2BH.
#'@param theta A vector with length 9 containing parameters
#'@param RC A list containing prior parameters, matrices and the data.
#'@return Returns a list containing predictive values of the parameters drawn out of the evaluated density.
#'@references Birgir Hrafnkelsson, Helgi Sigurdarson and Sigurdur M. Gardarson (2015) \emph{Bayesian Generalized Rating Curves}
bgplm.density_evaluation_unknown_c <- function(theta,RC){
  phi_b=theta[3]
  sig_b2=theta[2]
  zeta=theta[1]
  lambda=theta[4:9]

  f=lambda[1:5]-lambda[6]
  l=c(log(RC$w_tild+exp(theta[1])))

  varr=c(RC$epsilon*exp(RC$B%*%lambda))
  Sig_eps=diag(c(varr,0))
  #Matern covariance
  R_Beta=(1+sqrt(5)*RC$dist/exp(phi_b)+5*RC$dist^2/(3*exp(phi_b)^2))*
          exp(-sqrt(5)*RC$dist/exp(phi_b))+diag(RC$n_unique)*RC$nugget
  Sig_x=rbind(cbind(RC$Sig_ab,RC$m1),cbind(RC$m2,exp(sig_b2)*R_Beta))

  X=rbind(cbind(1,l,diag(l)%*%RC$A),RC$Z)
  L=t(chol(X%*%Sig_x%*%t(X)+Sig_eps))
  w=solve(L,RC$y-X%*%RC$mu_x)
  p=-0.5%*%t(w)%*%w-sum(log(diag(L)))-
    (RC$v+5-1)/2*log(RC$v*RC$s+f%*%RC$P%*%f)+
    sig_b2-exp(sig_b2)/RC$mu_sb+zeta-exp(zeta)/RC$mu_c-0.5/RC$tau_pb2*(phi_b-RC$mu_pb)^2

  W=solve(L,X%*%Sig_x)
  x_u=RC$mu_x+t(chol(Sig_x))%*%rnorm(RC$n_unique+2)
  sss=(X%*%x_u)-RC$y+rbind(sqrt(varr)*as.matrix(rnorm(RC$n)),0)
  x=as.matrix(x_u-t(W)%*%solve(L,sss))
  yp=(X %*% x)[1:RC$n,]
  #posterior predictive draw
  ypo=yp+as.matrix(rnorm(RC$n))*sqrt(varr)

  D=-2*sum(log(dlnorm(exp(RC$y[1:RC$n,]),yp,sqrt(varr))))

  return(list("p"=p,"x"=x,"y_post"=yp,"y_post_pred"=ypo,"sigma_eps"=varr,"DIC"=D))
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
bgplm.predict_u_known_c <- function(theta,x,RC){
  #store particular hyperparameter values
  sig_b2=exp(theta[1])
  phi_b=exp(theta[2])
  lambda <- theta[3:8]
  n=RC$n_unique
  m=RC$n_u
  #get sample of data variance using splines
  varr_u = c(exp(RC$B_u %*% lambda))
  #combine stages from data with unobserved stages
  w_all=c(RC$w_unique,RC$w_u)
  #calculating distance matrix for W_all
  dist_mat=as.matrix(dist(w_all))
  #Covariance of the joint prior for betas from data and beta unobserved.
  #Matern covariance formula used for v=5/2
  sigma_all=sig_b2*(1 + sqrt(5)*dist_mat/phi_b+(5*dist_mat^2)/(3*phi_b^2))*exp(-sqrt(5)*dist_mat/phi_b) + diag(length(w_all))*RC$nugget
  sigma_11=sigma_all[1:n,1:n]
  sigma_22=sigma_all[(n+1):(m+n),(n+1):(m+n)]
  sigma_12=sigma_all[1:n,(n+1):(n+m)]
  sigma_21=sigma_all[(n+1):(n+m),1:n]
  #parameters for the posterior of beta_u
  mu_x_u=sigma_21%*%solve(sigma_11,x[3:length(x)])
  Sigma_x_u=(sigma_22-sigma_21%*%solve(sigma_11,sigma_12))
  #a sample from posterior of beta_u drawn
  beta_u=as.numeric(mu_x_u) + rnorm(ncol(Sigma_x_u)) %*% chol(Sigma_x_u)
  #buidling blocks of the explanatory matrix X calculated
  l=log(RC$w_u-RC$c)
  X=cbind(rep(1,m),l,diag(l))
  x_u=c(x[1:2],beta_u)
  #sample from the posterior of discharge y
  yp_u <- c(X%*%x_u)
  ypo_u = yp_u + as.matrix(rnorm(m)) * sqrt(varr_u)
  return(list('x'=beta_u,'sigma_eps'=varr_u,'y_post'=yp_u,'y_post_pred'=ypo_u))
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
bgplm.predict_u_unknown_c <- function(theta,x,RC){
  #store particular hyperparameter values
  zeta <- theta[1]
  sig_b2=exp(theta[2])
  phi_b=exp(theta[3])
  lambda <- theta[4:9]
  n=RC$n_unique
  m=RC$n_u
  #get sample of data variance using splines
  varr_u = c(exp(RC$B_u %*% lambda))
  #combine stages from data with unobserved stages
  w_all=c(RC$w_unique,RC$w_u)
  #calculating distance matrix for W_all
  dist_mat=as.matrix(dist(w_all))
  #Covariance of the joint prior for betas from data and beta unobserved.
  #Matern covariance formula used for v=5/2
  sigma_all=sig_b2*(1 + sqrt(5)*dist_mat/phi_b+(5*dist_mat^2)/(3*phi_b^2))*exp(-sqrt(5)*dist_mat/phi_b) + diag(length(w_all))*RC$nugget
  sigma_11=sigma_all[1:n,1:n]
  sigma_22=sigma_all[(n+1):(m+n),(n+1):(m+n)]
  sigma_12=sigma_all[1:n,(n+1):(n+m)]
  sigma_21=sigma_all[(n+1):(n+m),1:n]
  #parameters for the posterior of beta_u
  mu_x_u=sigma_21%*%solve(sigma_11,x[3:length(x)])
  Sigma_x_u=(sigma_22-sigma_21%*%solve(sigma_11,sigma_12))
  #a sample from posterior of beta_u drawn
  beta_u=as.numeric(mu_x_u) + rnorm(ncol(Sigma_x_u)) %*% chol(Sigma_x_u)
  above_c <- -(exp(zeta)-RC$w_min) < RC$w_u
  m_above_c <- sum(above_c)
  #buidling blocks of the explanatory matrix X calculated
  l=log(RC$w_u[above_c]-RC$w_min+exp(zeta))
  X=cbind(rep(1,m_above_c),l,diag(l))
  #vector of parameters
  x_u=c(x[1:2],beta_u[above_c])
  #sample from the posterior of discharge y
  yp_u <- c(X%*%x_u)
  ypo_u = yp_u + as.matrix(rnorm(m_above_c)) * sqrt(varr_u[above_c])
  return(list('x'=beta_u,'sigma_eps'=varr_u,'y_post'=c(rep(-Inf,m-m_above_c),yp_u),'y_post_pred'=c(rep(-Inf,m-m_above_c),ypo_u)))
}

#' Generalized Bayesian Power Law Model
#'
#' Infers a rating curve for paired measurements of stage and discharge using a generalized power law model described in Hrafnkelsson et al.
#'@param formula formula with name of discharge column in data as response and name of stage column in data as the single covariate.
#'@param data data.frame containing the columns in formula
#'@param W_limits vector of length 2 setting the lower and upper bound of stage values at which a rating curve should be predicted. If NULL, the known value of c or the mle of c will be used as lower bound (depending on the value of the input parameter c) and maximum stage value in data as upper bound.
#'@param country Name of the country the prior parameters should be defined for, default value is "Iceland".
#'@param Wmin Positive numeric value for the lowest stage the user wants to calculate a rating curve. If input is an empty string (default) Wmin will
#'automatically be set to c_hat.
#'@param Wmax Positive numeric value for the highest stage the user wants to calculate a rating curve. If input is an empty string (default) Wmax will
#'automatically be set to the maximum stage of the data.
#'@return List containing information on the calculated rating curve,
#'the data frames observedData, betaData, completePrediction, observedPrediction, TableOfData, FitTable, LowerTable, UpperTable, plotTable.
#'@references Birgir Hrafnkelsson, Helgi Sigurdarson and Sigurdur M. Gardarson (2015) \emph{Bayesian Generalized Rating Curves}
#'@seealso \code{\link{clean}}
gbplm <- function(formula,data,c_param=NULL,W_limits=NULL,country="Iceland",forcepoint=rep(FALSE,nrow(data))){
  suppressPackageStartupMessages(require(doParallel))
  #TODO: add error message if length(formula)!=3 or if it contains more than one covariate. Also make sure that names in formula exist in data
  model_dat <- data[,all.vars(formula)]
  RC=priors(country)
  RC$formula = formula
  RC$nugget=10^-8
  RC$mu_sb=0.5
  RC$mu_pb=0.5
  RC$tau_pb2=0.25^2
  RC$s=3
  RC$v=5
  RC$y=rbind(as.matrix(log(model_dat[,1])),0)
  RC$w=as.matrix(model_dat[,2])
  RC$w_tild=RC$w-min(RC$w)
  Adist1 <- Adist(RC$w)
  RC$A=Adist1$A
  RC$dist=Adist1$dist
  RC$n=Adist1$n
  RC$N=Adist1$N
  RC$O=Adist1$O

  RC$Sig_ab= rbind(c(RC$sig_a^2, RC$p_ab*RC$sig_a*RC$sig_b), c(RC$p_ab*RC$sig_a*RC$sig_b, RC$sig_b^2))
  RC$mu_x=as.matrix(c(RC$mu_a,RC$mu_b, rep(0,RC$n)))

  RC$P=diag(nrow=5,ncol=5,6)-matrix(nrow=5,ncol=5,1)
  RC$B=B_splines(t(RC$w_tild)/RC$w_tild[length(RC$w_tild)])
  RC$epsilon=rep(1,RC$N)
  #Spyrja Bigga út í varíans hér
  forcepoint_dat=model_dat[forcepoint,]
  RC$epsilon[forcepoint]=1/RC$N

  RC$Z=cbind(t(rep(0,2)),t(rep(1,RC$n)))
  RC$m1=matrix(0,nrow=2,ncol=RC$n)
  RC$m2=matrix(0,nrow=RC$n,ncol=2)
  RC$c=c_param
  if(!is.null(RC$c)){
    density_fun <- density_evaluation_known_c
    unobserved_prediction_fun <- predict_u_known_c
  }else{
    density_fun <- density_evaluation_unknown_c
    unobserved_prediction_fun <- predict_u_unknown_c
  }
  #determine proposal density
  theta_init=rep(0,9)
  loss_fun = function(th) {-density_fun(th,RC)$p}
  optim_obj=optim(par=theta_init,loss_fun,method="L-BFGS-B",hessian=TRUE)
  t_m =optim_obj$par
  H=optim_obj$hessian
  LH=t(chol(H))/0.8

  #make Wmin and Wmax divisable by 10 up, both in order to make rctafla and so l_m is defined
  if(is.null(W_limits)){
    Wmax=ceiling(max(RC$w)*10)/10
    Wmin=ceiling(10*ifelse(is.null(RC$c),min(RC$w)-exp(t_m[1]),RC$c))/10
  }else{
    Wmin=W_limits[1]
    Wmax=W_limits[2]
  }
  WFill=W_unobserved(c(RC$O),min=Wmin,max=Wmax)
  RC$W_u=WFill$W_u
  RC$W_u_tild=WFill$W_u_tild
  Bsiminput=t(RC$W_u_tild)/RC$W_u_tild[length(RC$W_u_tild)]
  Bsiminput[is.na(Bsiminput)]=0
  RC$Bsim=B_splines(Bsiminput)

  #MCMC parameters added, number of iterations,burnin and thin
  Nit=20000
  burnin=2000
  thin=5
  cl <- makeCluster(4)
  registerDoParallel(cl)
  MCMC <- foreach(i=1:4,.combine=cbind,.export=c("density_fun","unobserved_prediction_fun")) %dopar% {
    output=matrix(0,nrow=length(t_m)+RC$n+2+length(RC$W_u)+RC$n+length(RC$W_u),ncol=Nit)
    t_old=as.matrix(t_m)
    Dens<-density_fun(t_old,RC)
    p_old=Dens$p
    ypo_old=Dens$ypo
    x_old=Dens$x
    unobserved_old=unobserved_prediction_fun(c(t_old,x_old),RC)
    print(i)
    for(j in 1:Nit){
      print(j)
      t_new=t_old+solve(t(LH),rnorm(length(t_m),0,1))
      Densnew <- density_fun(t_new,RC)
      x_new=Densnew$x
      ypo_new=Densnew$ypo
      p_new=Densnew$p
      unobserved_new=unobserved_prediction_fun(c(t_old,x_old),RC)
      logR=p_new-p_old

      if (logR>log(runif(1))){
        t_old=t_new
        p_old=p_new
        ypo_old=ypo_new
        x_old=x_new
        unobserved_old=unobserved_new
      }
      output[,j]=rbind(t_old,x_old,unobserved_old[1:length(RC$W_u)],ypo_old,unobserved_old[(length(RC$W_u)+1):length(unobserved_old)])
    }
    seq=seq(burnin,Nit,thin)
    ypo_obs=ypo_obs[,seq]
    param=param[,seq]
    output=rbind(param,unobserved[1:length(RC$W_u)],ypo_obs,unobserved[(length(RC$W_u)+1):nrow(unobserved)])
    output=list('theta'=param[1:length(t_m),],'a'=exp(param[length(t_m)+1,]),
                'b'=param[length(t_m)+2,],'beta'=rbind(param[(length(t_m)+3):nrow(param)],unobserved[1:length(RC$W_u),]),
                'ypo'=exp(rbind(ypo_obs,unobserved[(length(RC$W_u)+1):nrow(unobserved)])))
    if(!is.null(RC$c)){
      output$theta[1,] <- min(RC$O)-exp(output$theta[1,])
      output$theta[2,] <- exp(output$theta[2,])
      output$theta[3,] <- exp(output$theta[3,])
    }else{
      output$theta[1,] <- exp(output$theta[1,])
      output$theta[2,] <- exp(output$theta[2,])
    }

    return(output)
  }
  #TODO: create S3 object to store results from MCMC chain
  stopCluster(cl)
  rating_curve <- as.data.frame(t(apply(MCMC$ypo,1,quantile, probs = c(0.025,0.5, 0.975),na.rm=T)))
  names(rating_curve) <- c('lower','median','upper')
  param_summary <- as.data.frame(t(apply(rbind(MCMC$a,MCMC$b,MCMC$theta),1,quantile, probs = c(0.025,0.5, 0.975),na.rm=T)))
  names(param_summary) <- c('lower','median','upper')
  beta_summary <- as.data.frame(t(apply(MCMC$ypo,1,quantile, probs = c(0.025,0.5, 0.975),na.rm=T)))
  names(beta_summary) <- c('lower','median','upper')
  W=c(RC$O,RC$W_u)
  param_names <- c('var_beta','phi_beta',paste0('lambda_',1:6))
  if(is.null(RC$c)){
    param_names <- c('c',param_names)
  }
  param_summary <- cbind(data.frame(parameter=param_names),param_summary)
  rating_curve <- rating_curve[order(W),]
  beta_summary <- beta_summary[order(W),]

  #S3 object gbplm Test
  result_obj=list()
  result_obj$formula <- formula
  result_obj$data <- data
  result_obj$W_full <- W
  result_obj$post_a = MCMC$a
  result_obj$post_b = MCMC$b
  if(!is.null(RC$c)){
    result_obj$post_c <- MCMC$theta[1,]
    result_obj$post_var_beta <- MCMC$theta[2,]
    result_obj$post_phi_beta <- MCMC$theta[3,]
  }else{
    result_obj$post_c <- NULL
    result_obj$post_var_beta <- MCMC$theta[2,]
    result_obj$post_phi_beta <- MCMC$theta[3,]
  }
  result_obj$param_summary <- param_summary
  result_obj$beta <- beta_summary
  result_obj$rating_curve <- rating_curve

  return(result_obj)
}

print.gbplm <- function(x,...){
  cat("\nCall:\n",
      paste(deparse(x$formula), sep = "\n", collapse = "\n"), "\n\n", sep = "")
}


summary.gbplm <- function(x,...){
  cat("\nFormula: \n",
      paste(deparse(x$formula), sep = "\n", collapse = "\n"),
      "\nParameters:\n",
      paste(deparse(x$post_a), sep = "\n", collapse = "\n"),
      "\na:",
      paste(deparse(x$post_a), sep = "\n", collapse = "\n"),
      "\nb:",
      paste(deparse(x$post_b), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")


  if(!is.null(RC$c)){
    cat("\n c:",
    paste(deparse(x$post_c), sep = "\n", collapse = "\n"),
    "\nPosterior beta:",
    paste(deparse(x$post_a), sep = "\n", collapse = "\n"),
    "\nPosterior phi:",
    paste(deparse(x$post_b), sep = "\n", collapse = "\n"),
    "\n\n", sep = "")

  } else{
    cat("\nPosterior beta:",
    paste(deparse(x$post_a), sep = "\n", collapse = "\n"),
    "\nPosterior phi:",
    paste(deparse(x$post_b), sep = "\n", collapse = "\n"),
    "\n\n", sep = "")
    }


}

#'Density evaluation for model2
#'
#'Evaluates the log density of the posterior distribution of the parameters of .
#'@param th A vector with length 9 containing parameters
#'@param RC A list containing prior parameters, matrices and the data.
#'@return Returns a list containing predictive values of the parameters drawn out of the evaluated density.
#'@references Birgir Hrafnkelsson, Helgi Sigurdarson and Sigurdur M. Gardarson (2015) \emph{Bayesian Generalized Rating Curves}
density_evaluation_known_c <- function(th,RC){
  phi_b=th[2]
  sig_b2=th[1]
  lambda=th[4:8]

  f=lambda[1:5]-lambda[6]
  l=log(RC$w-RC$c)

  varr=c(RC$epsilon*exp(RC$B%*%lambda))
  Sig_eps=diag(c(varr,0))
  #Matern covariance
  R_Beta=(1+sqrt(5)*RC$dist/exp(phi_b)+5*RC$dist^2/(3*exp(phi_b)^2))*exp(-sqrt(5)*RC$dist/exp(phi_b))+diag(RC$n)*RC$nugget
  Sig_x=rbind(cbind(RC$Sig_ab,RC$m1),cbind(RC$m2,exp(sig_b2)*R_Beta))

  X=rbind(cbind(1,l,diag(l)%*%RC$A),RC$Z)#1337 micro
  L=t(chol(X%*%Sig_x%*%t(X)+Sig_eps))#2521 micro
  w=solve(L,RC$y-X%*%RC$mu_x)#754 micro
  p=-0.5%*%t(w)%*%w-sum(log(diag(L)))-
    (RC$v+5-1)/2*log(RC$v*RC$s+f%*%RC$P%*%f)+
    sig_b2-exp(sig_b2)/RC$mu_sb-0.5/RC$tau_pb2*(phi_b-RC$mu_pb)^2 #63 micro

  W=solve(L,X%*%Sig_x)
  x_u=RC$mu_x+t(chol(Sig_x))%*%rnorm(RC$n+2)
  sss=(X%*%x_u)-RC$y+rbind(sqrt(varr)*as.matrix(rnorm(RC$N)),0)
  x=as.matrix(x_u-t(W)%*%solve(L,sss))
  yp=(X %*% x)[1:RC$N,]
  #posterior predictive draw
  ypo=yp+as.matrix(rnorm(RC$N))*sqrt(varr)

  return(list("p"=p,"x"=x,"yp"=yp,"ypo"=ypo,"varr"=varr))
}

#'Density evaluation for model2
#'
#'Evaluates the log density of the posterior distribution of the parameters of model2BH.
#'@param th A vector with length 9 containing parameters
#'@param RC A list containing prior parameters, matrices and the data.
#'@return Returns a list containing predictive values of the parameters drawn out of the evaluated density.
#'@references Birgir Hrafnkelsson, Helgi Sigurdarson and Sigurdur M. Gardarson (2015) \emph{Bayesian Generalized Rating Curves}
density_evaluation_unknown_c <- function(th,RC){
  phi_b=th[3]
  sig_b2=th[2]
  zeta=th[1]
  lambda=th[4:9]

  f=lambda[1:5]-lambda[6]
  l=c(log(RC$w_tild+exp(th[1])))

  varr=c(RC$epsilon*exp(RC$B%*%lambda))
  Sig_eps=diag(c(varr,0))
  #Matern covariance
  R_Beta=(1+sqrt(5)*RC$dist/exp(phi_b)+5*RC$dist^2/(3*exp(phi_b)^2))*exp(-sqrt(5)*RC$dist/exp(phi_b))+diag(RC$n)*RC$nugget
  Sig_x=rbind(cbind(RC$Sig_ab,RC$m1),cbind(RC$m2,exp(sig_b2)*R_Beta))

  X=rbind(cbind(1,l,diag(l)%*%RC$A),RC$Z)#1337 micro
  L=t(chol(X%*%Sig_x%*%t(X)+Sig_eps))#2521 micro
  w=solve(L,RC$y-X%*%RC$mu_x)#754 micro
  p=-0.5%*%t(w)%*%w-sum(log(diag(L)))-
    (RC$v+5-1)/2*log(RC$v*RC$s+f%*%RC$P%*%f)+
    sig_b2-exp(sig_b2)/RC$mu_sb+zeta-exp(zeta)/RC$mu_c-0.5/RC$tau_pb2*(phi_b-RC$mu_pb)^2 #63 micro

  W=solve(L,X%*%Sig_x)
  x_u=RC$mu_x+t(chol(Sig_x))%*%rnorm(RC$n+2)
  sss=(X%*%x_u)-RC$y+rbind(sqrt(varr)*as.matrix(rnorm(RC$N)),0)
  x=as.matrix(x_u-t(W)%*%solve(L,sss))
  yp=(X %*% x)[1:RC$N,]
  #posterior predictive draw
  ypo=yp+as.matrix(rnorm(RC$N))*sqrt(varr)

  #D=-2*sum(log(dlnorm(exp(RC$y[1:RC$N,]),yp,sqrt(varr))))#45.04

  return(list("p"=p,"x"=x,"yp"=yp,"ypo"=ypo,"varr"=varr))
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
predict_u_known_c <- function(param,RC){
  #collecting parameters from the MCMC sample
  th=param[1:8]
  x=param[9:length(param)]
  sig_b2=exp(th[1])
  phi_b=exp(th[2])
  lambda = th[3:length(th)]
  #calculate spline variance from B_splines
  varr = c(exp(RC$Bsim %*% lambda))
  m=length(RC$W_u)
  n=RC$n
  #combine stages from data with unobserved stages
  W_all=c(RC$O,RC$W_u)
  #calculating distance matrix for W_all
  dist=abs(outer(W_all,W_all,FUN="-"))
  #defining the variance of the joint prior for betas from data and beta unobserved, that is p(beta,beta_u).
  #Matern covariance formula used for v=5/2
  sigma_all=sig_b2*(1 + sqrt(5)*dist/phi_b+(5*dist^2)/(3*phi_b^2))*exp(-sqrt(5)*dist/phi_b) + diag(length(W_all))*RC$nugget
  sigma_11=sigma_all[1:n,1:n]
  sigma_22=sigma_all[(n+1):(m+n),(n+1):(m+n)]
  sigma_12=sigma_all[1:n,(n+1):(n+m)]
  sigma_21=sigma_all[(n+1):(n+m),1:n]
  #parameters for the posterior of beta_u
  mu_u=sigma_21%*%solve(sigma_11,x[3:length(x)])
  Sigma_u=(sigma_22-sigma_21%*%solve(sigma_11,sigma_12))
  #a sample from posterior of beta_u drawn
  beta_u=as.numeric(mu_u) + rnorm(ncol(Sigma_u)) %*% chol(Sigma_u)
  #buidling blocks of the explanatory matrix X calculated
  l=log(RC$W_u-RC$c)
  X=cbind(rep(1,m),l,matrix(0,m,n),diag(l))
  #vector of parameters
  x_extended=c(x,beta_u)
  #sample from the posterior of discharge y
  ypo_extended = X%*%x_extended + as.matrix(rnorm(m)) * sqrt(varr)
  return(c(beta_u,ypo_extended[(RC$n+1):length(ypo_extended)]))
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
predict_u_unknown_c <- function(param,RC){
  #collecting parameters from the MCMC sample
  th=param[1:9]
  x=param[10:length(param)]
  zeta=th[1]
  phi_b=exp(th[3])
  sig_b2=exp(th[2])
  lambda = th[4:9]
  #calculate spline variance from B_splines
  varr = c(exp(RC$Bsim %*% lambda))
  m=length(RC$W_u)
  n=RC$n
  #combine stages from data with unobserved stages
  W_all=c(RC$O,RC$W_u)
  #calculating distance matrix for W_all
  dist=abs(outer(W_all,W_all,FUN="-"))
  #defining the variance of the joint prior for betas from data and beta unobserved, that is p(beta,beta_u).
  #Matern covariance formula used for v=5/2
  sigma_all=sig_b2*(1 + sqrt(5)*dist/phi_b+(5*dist^2)/(3*phi_b^2))*exp(-sqrt(5)*dist/phi_b) + diag(length(W_all))*RC$nugget
  sigma_11=sigma_all[1:n,1:n]
  sigma_22=sigma_all[(n+1):(m+n),(n+1):(m+n)]
  sigma_12=sigma_all[1:n,(n+1):(n+m)]
  sigma_21=sigma_all[(n+1):(n+m),1:n]
  #parameters for the posterior of beta_u
  mu_u=sigma_21%*%solve(sigma_11,x[3:length(x)])
  Sigma_u=(sigma_22-sigma_21%*%solve(sigma_11,sigma_12))
  #a sample from posterior of beta_u drawn
  beta_u=as.numeric(mu_u) + rnorm(ncol(Sigma_u)) %*% chol(Sigma_u)
  #buidling blocks of the explanatory matrix X calculated
  c=min(RC$O)-exp(zeta)
  l=log(RC$W_u-c)
  X=cbind(rep(1,m),l,diag(l))
  #vector of parameters
  x_extended=c(x[1],x[2],beta_u)
  #sample from the posterior of discharge y
  ypo_extended = X%*%x_extended + as.matrix(rnorm(m)) * sqrt(varr)
  return(as.matrix(c(beta_u,ypo_extended)))
}




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

bgplm <- function(formula,data,c_param=NULL,w_limits=NULL,country="Iceland",forcepoint=rep(FALSE,nrow(data)),...){
  #argument checking
  model_dat <- data[,all.vars(formula)]
  model_dat <- model_dat[order(model_dat[,2,drop=T]),]
  Q <- log(model_dat[,1,drop=T])
  W <- model_dat[,2,drop=T]
  MCMC_output_list <- bgplm.inference(Q,W,c_param,w_limits,country,forcepoint,...)
  rating_curve <- data.frame(W=MCMC_output_list$W,as.data.frame(t(apply(MCMC_output_list$ypo,1,quantile, probs = c(0.025,0.5, 0.975),na.rm=T))),row.names=NULL)
  names(rating_curve) <- c('W','lower','median','upper')
  rating_curve <- rating_curve[order(rating_curve$W),]
  beta_summary <- data.frame(W=unique(MCMC_output_list$W),as.data.frame(t(apply(MCMC_output_list$beta,1,quantile, probs = c(0.025,0.5, 0.975),na.rm=T))),row.names=NULL)
  names(beta_summary) <- c('W','lower','median','upper')
  beta_summary <- beta_summary[order(beta_summary$W),]
  param_summary <- as.data.frame(t(apply(rbind(MCMC_output_list$a,MCMC_output_list$b,MCMC_output_list$theta),1,quantile, probs = c(0.025,0.5, 0.975),na.rm=T)))
  names(param_summary) <- c('lower','median','upper')
  param_names <- c('var_beta','phi_beta',paste0('lambda_',1:6))
  if(is.null(c_param)){
    param_names <- c('c',param_names)
  }
  param_names <- c('a','b',param_names)
  param_summary <- data.frame(param_summary,row.names=param_names)

  DIC <- quantile(MCMC_output_list$DIC,probs=c(0.025,0.05,0.975))
  names(DIC) <- c('lower','median','upper')

  #S3 object gbplm Test
  result_obj=list()
  result_obj$formula <- formula
  result_obj$data <- model_dat
  result_obj$W_full <- W
  result_obj$post_a = MCMC_output_list$a
  result_obj$post_b = MCMC_output_list$b
  if(!is.null(c_param)){
    result_obj$post_c <- MCMC_output_list$theta[1,]
    result_obj$post_var_beta <- MCMC_output_list$theta[2,]
    result_obj$post_phi_beta <- MCMC_output_list$theta[3,]
  }else{
    result_obj$post_c <- NULL
    result_obj$post_var_beta <- MCMC_output_list$theta[2,]
    result_obj$post_phi_beta <- MCMC_output_list$theta[3,]
  }
  result_obj$DIC <- DIC
  result_obj$param_summary <- param_summary
  result_obj$beta <- beta_summary
  result_obj$rating_curve <- rating_curve
  attr(result_obj, "class") <- "bgplm"
  return(result_obj)
}

bgplm.inference <- function(y,w,c_param=NULL,w_limits=NULL,country="Iceland",forcepoint=rep(FALSE,nrow(data)),num_chains=4,nr_iter=20000,burnin=2000,thin=5){
  suppressPackageStartupMessages(require(doParallel))
  #TODO: add error message if length(formula)!=3 or if it contains more than one covariate. Also make sure that names in formula exist in data
  RC=priors(country)
  RC$nugget=10^-8
  RC$mu_sb=0.5
  RC$mu_pb=0.5
  RC$tau_pb2=0.25^2
  RC$s=3
  RC$v=5
  RC$y=rbind(as.matrix(y),0)
  RC$w=as.matrix(w)
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
  theta_init=rep(0,8)
  loss_fun = function(th) {-density_fun(th,RC)$p}
  optim_obj=optim(par=theta_init,loss_fun,method="L-BFGS-B",hessian=TRUE)
  t_m =optim_obj$par
  H=optim_obj$hessian
  LH=t(chol(H))/0.8

  #make Wmin and Wmax divisable by 10 up, both in order to make rctafla and so l_m is defined
  if(is.null(w_limits)){
    Wmax=ceiling(max(RC$w)*10)/10
    Wmin=ceiling(10*ifelse(is.null(RC$c),min(RC$w)-exp(t_m[1]),RC$c))/10
  }else{
    Wmin=w_limits[1]
    Wmax=w_limits[2]
  }
  WFill=W_unobserved(c(RC$O),min=Wmin,max=Wmax)
  RC$W_u=WFill$W_u
  RC$W_u_tild=WFill$W_u_tild
  Bsiminput=t(RC$W_u_tild)/RC$W_u_tild[length(RC$W_u_tild)]
  Bsiminput[is.na(Bsiminput)]=0
  RC$Bsim=B_splines(Bsiminput)

  #MCMC parameters added, number of iterations,burnin and thin
  if(num_chains>4){
    stop('Max number of chains is 4. Please pick a lower number of chains')
  }
  cl <- makeCluster(num_chains)
  registerDoParallel(cl)
  MCMC_output_mat <- foreach(i=1:num_chains,.combine=cbind) %dopar% {
    latent_and_hyper_param <- matrix(0,nrow=length(t_m)+RC$n+2,ncol=nr_iter)
    y_pred_mat <- matrix(0,nrow=RC$N,ncol=nr_iter)
    DIC <- rep(0,nr_iter)
    t_old=as.matrix(t_m)
    Dens<-density_fun(t_old,RC)
    p_old=Dens$p
    ypo_old=Dens$ypo
    x_old=Dens$x
    DIC_old <- Dens$DIC
    for(j in 1:nr_iter){
      t_new=t_old+solve(t(LH),rnorm(length(t_m),0,1))
      Densnew <- density_fun(t_new,RC)
      x_new=Densnew$x
      ypo_new=Densnew$ypo
      p_new=Densnew$p
      DIC_new <- Densnew$DIC
      logR=p_new-p_old
      if (logR>log(runif(1))){
        t_old=t_new
        p_old=p_new
        ypo_old=ypo_new
        x_old=x_new
        DIC_old <- DIC_new
      }
      latent_and_hyper_param[,j] <- c(t_old,x_old)
      y_pred_mat[,j] <- ypo_old
      DIC[j] <- DIC_old
    }

    seq=seq(burnin,nr_iter,thin)
    latent_and_hyper_param=latent_and_hyper_param[,seq]
    y_pred_mat=y_pred_mat[,seq]
    DIC=DIC[seq]
    unobserved_mat=apply(latent_and_hyper_param,2,function(x) unobserved_prediction_fun(x,RC))
    output_mat=rbind(latent_and_hyper_param,unobserved_mat[1:length(RC$W_u),],y_pred_mat,unobserved_mat[(length(RC$W_u)+1):nrow(unobserved_mat),],t(matrix(DIC)))
    if(is.null(RC$c)){
      output_mat[1,] <- min(RC$O)-exp(output_mat[1,])
      output_mat[2,] <- exp(output_mat[2,])
      output_mat[3,] <- exp(output_mat[3,])
    }else{
      output_mat[1,] <- exp(output_mat[1,])
      output_mat[2,] <- exp(output_mat[2,])
    }
    return(output_mat)
  }
  stopCluster(cl)
  MCMC_output_list=list('W'=c(RC$w,RC$W_u),'theta'=MCMC_output_mat[1:length(t_m),],'a'=exp(MCMC_output_mat[length(t_m)+1,]),
           'b'=MCMC_output_mat[length(t_m)+2,],'beta'=MCMC_output_mat[(length(t_m)+3):(length(t_m)+2+RC$n+length(RC$W_u)),],
           'ypo'=exp(MCMC_output_mat[(length(t_m)+2+RC$n+length(RC$W_u)+1):(nrow(MCMC_output_mat)-1),]),'DIC'=MCMC_output_mat[nrow(MCMC_output_mat),])

  return(MCMC_output_list)
}



#create a predict method for interpolation of posterior predictive

#'Density evaluation for model2
#'
#'Evaluates the log density of the posterior distribution of the parameters of .
#'@param th A vector with length 9 containing parameters
#'@param RC A list containing prior parameters, matrices and the data.
#'@return Returns a list containing predictive values of the parameters drawn out of the evaluated density.
#'@references Birgir Hrafnkelsson, Helgi Sigurdarson and Sigurdur M. Gardarson (2015) \emph{Bayesian Generalized Rating Curves}
density_evaluation_known_c <- function(th,RC){
  sig_b2=th[1]
  phi_b=th[2]
  lambda=th[3:8]

  f=lambda[1:5]-lambda[6]
  l=c(log(RC$w-RC$c))

  varr=c(RC$epsilon*exp(RC$B%*%lambda))
  Sig_eps=diag(c(varr,0))
  #Matern covariance
  R_Beta=(1+sqrt(5)*RC$dist/exp(phi_b)+5*RC$dist^2/(3*exp(phi_b)^2))*exp(-sqrt(5)*RC$dist/exp(phi_b))+diag(RC$n)*RC$nugget
  Sig_x=rbind(cbind(RC$Sig_ab,RC$m1),cbind(RC$m2,exp(sig_b2)*R_Beta))

  X=rbind(cbind(1,l,diag(l)%*%RC$A),RC$Z)
  L=t(chol(X%*%Sig_x%*%t(X)+Sig_eps))
  w=solve(L,RC$y-X%*%RC$mu_x)
  p=-0.5%*%t(w)%*%w-sum(log(diag(L)))-
    (RC$v+5-1)/2*log(RC$v*RC$s+f%*%RC$P%*%f)+
    sig_b2-exp(sig_b2)/RC$mu_sb-0.5/RC$tau_pb2*(phi_b-RC$mu_pb)^2

  W=solve(L,X%*%Sig_x)
  x_u=RC$mu_x+t(chol(Sig_x))%*%rnorm(RC$n+2)
  sss=(X%*%x_u)-RC$y+rbind(sqrt(varr)*as.matrix(rnorm(RC$N)),0)
  x=as.matrix(x_u-t(W)%*%solve(L,sss))
  yp=(X %*% x)[1:RC$N,]
  #posterior predictive draw
  ypo=yp+as.matrix(rnorm(RC$N))*sqrt(varr)
  D=-2*sum(log(dlnorm(exp(RC$y[1:RC$N,]),yp,sqrt(varr))))

  return(list("p"=p,"x"=x,"yp"=yp,"ypo"=ypo,"varr"=varr,"DIC"=D))
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
  R_Beta=(1+sqrt(5)*RC$dist/exp(phi_b)+5*RC$dist^2/(3*exp(phi_b)^2))*
          exp(-sqrt(5)*RC$dist/exp(phi_b))+diag(RC$n)*RC$nugget
  Sig_x=rbind(cbind(RC$Sig_ab,RC$m1),cbind(RC$m2,exp(sig_b2)*R_Beta))

  X=rbind(cbind(1,l,diag(l)%*%RC$A),RC$Z)
  L=t(chol(X%*%Sig_x%*%t(X)+Sig_eps))
  w=solve(L,RC$y-X%*%RC$mu_x)
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

  D=-2*sum(log(dlnorm(exp(RC$y[1:RC$N,]),yp,sqrt(varr))))

  return(list("p"=p,"x"=x,"yp"=yp,"ypo"=ypo,"varr"=varr,"DIC"=D))
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
  l=log(RC$W_u_tild+exp(zeta))
  X=cbind(rep(1,m),l,matrix(0,m,n),diag(l))
  #vector of parameters
  x=c(x,beta_u)
  #sample from the posterior of discharge y
  ypo_u = X%*%x + as.matrix(rnorm(m)) * sqrt(varr)
  return(as.matrix(c(beta_u,ypo_u)))
}




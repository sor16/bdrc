#' Bayesian Power Law Model with varying variance
#'
#' Infers a rating curve for paired measurements of stage and discharge using a  power law model described in Hrafnkelsson et al.
#'@param formula formula with name of discharge column in data as response and name of stage column in data as the single covariate.
#'@param data data.frame containing the columns in formula
#'@param w_max vector of length 2 setting the lower and upper bound of stage values at which a rating curve should be predicted. If NULL, the known value of c or the mle of c will be used as lower bound (depending on the value of the input parameter c) and maximum stage value in data as upper bound.
#'@param country Name of the country the prior parameters should be defined for, default value is "Iceland".
#'@param Wmin Positive numeric value for the lowest stage the user wants to calculate a rating curve. If input is an empty string (default) Wmin will
#'automatically be set to c_hat.
#'@param Wmax Positive numeric value for the highest stage the user wants to calculate a rating curve. If input is an empty string (default) Wmax will
#'automatically be set to the maximum stage of the data.
#'@return List containing information on the calculated rating curve,
#'the data frames observedData, betaData, completePrediction, observedPrediction, TableOfData, FitTable, LowerTable, UpperTable, plotTable.
#'@references Birgir Hrafnkelsson, Helgi Sigurdarson and Sigurdur M. Gardarson (2015) \emph{Bayesian Generalized Rating Curves}
#'@seealso \code{\link{clean}}

bplm <- function(formula,data,c_param=NULL,w_max=NULL,forcepoint=rep(FALSE,nrow(data)),...){
    #TODO: argument checking
    model_dat <- data[,all.vars(formula)]
    model_dat <- model_dat[order(model_dat[,2,drop=T]),]
    Q <- model_dat[,1,drop=T]
    w <- model_dat[,2,drop=T]
    MCMC_output_list <- bplm.inference(y=log(Q),w=w,c_param,w_max,forcepoint)
    result_obj=list()
    attr(result_obj, "class") <- "bplm"
    result_obj$formula <- formula
    result_obj$data <- model_dat
    result_obj$a_posterior = MCMC_output_list$x[1,]
    result_obj$b_posterior = MCMC_output_list$x[2,]
    if(is.null(c_param)){
        result_obj$c_posterior <- MCMC_output_list$theta[1,]
        result_obj$sigma_eta_posterior <- MCMC_output_list$theta[2,]
        for(i in 3:nrow(MCMC_output_list$theta)){
          result_obj[[paste0('eta',i,'_posterior')]] <- MCMC_output_list$theta[i,]
        }
    }else{
        result_obj$c_posterior <- NULL
        result_obj$sigma_eta_posterior <- MCMC_output_list$theta[1,]
        for(i in 2:nrow(MCMC_output_list$theta)){
          result_obj[[paste0('eta',i,'_posterior')]] <- MCMC_output_list$theta[i,]
        }
    }
    result_obj$Q_posterior_predictive <- exp(MCMC_output_list$y_post_pred)
    result_obj$Q_posterior <- exp(MCMC_output_list$y_post)
    result_obj$sigma_eps_posterior <- sqrt(MCMC_output_list$sigma_eps)
    result_obj$DIC_posterior <- MCMC_output_list$DIC
    #summary objects
    result_obj$rating_curve <- get_MCMC_summary(result_obj$Q_posterior_predictive,w=MCMC_output_list$w)
    result_obj$rating_curve_mean <- get_MCMC_summary(result_obj$Q_posterior,w=MCMC_output_list$w)
    result_obj$sigma_eps_summary <- get_MCMC_summary(result_obj$sigma_eps_posterior,w=MCMC_output_list$w)
    result_obj$param_summary <- get_MCMC_summary(rbind(MCMC_output_list$x[1,],MCMC_output_list$x[2,],MCMC_output_list$theta))
    row.names(result_obj$param_summary) <- get_param_names('bplm',c_param)
    result_obj$DIC_summary <- get_MCMC_summary(result_obj$DIC_posterior)
    result_obj$run_info <- MCMC_output_list$run_info
    return(result_obj)
}

bplm.inference <- function(y,w,c_param=NULL,w_max=NULL,forcepoint=rep(FALSE,nrow(data)),num_chains=4,nr_iter=20000,burnin=2000,thin=5){
    require(parallel)
    #TODO: add error message if length(formula)!=3 or if it contains more than one covariate. Also make sure that names in formula exist in data
    RC=priors('bplm',c_param)
    RC$y <- as.matrix(y)
    RC$w <- w
    RC$w_min <- min(RC$w)
    RC$w_max <- max(RC$w)
    RC$w_tild <- RC$w-RC$w_min
    RC$n <- length(w)

    RC$P <- lower.tri(matrix(rep(1,36),6,6),diag=T)*1
    RC$B <- B_splines(t(RC$w_tild)/RC$w_tild[length(RC$w_tild)])

    RC$epsilon <- rep(1,RC$n)
    RC$epsilon[forcepoint]=1/RC$n
    if(!is.null(RC$c)){
      if(RC$c<=RC$w_min){
        density_fun <- bplm.density_evaluation_known_c
        unobserved_prediction_fun <- bplm.predict_u_known_c
      }else{
        stop(paste0('the given c must be less than the lowest stage measurement, which is ',RC$w_min,' m'))
      }
    }else{
        density_fun <- bplm.density_evaluation_unknown_c
        unobserved_prediction_fun <- bplm.predict_u_unknown_c
    }
    #determine proposal density
    RC$theta_length <- if(is.null(RC$c)) 8 else 7
    theta_init <- rep(0,RC$theta_length)
    loss_fun  <-  function(th) {-density_fun(th,RC)$p}
    optim_obj <- optim(par=theta_init,loss_fun,method="L-BFGS-B",hessian=TRUE)
    theta_m <- optim_obj$par
    H <- optim_obj$hessian
    RC$LH <- t(chol(H))/0.8

    w_min <- ifelse(is.null(RC$c),min(RC$w)-exp(theta_m[1]),RC$c)
    if(is.null(w_max)){
      w_max <- RC$w_max
    }
    if(w_max<RC$w_max){
      stop(paste0('maximum stage value must be larger than the maximum stage value in the data, which is ', RC$w_max,' m'))
    }
    RC$w_u <- W_unobserved(RC,w_min,w_max)
    RC$n_u <- length(RC$w_u)
    w_u_std <- ifelse(RC$w_u < RC$w_min,0.0,ifelse(RC$w_u>RC$w_max,1.0,(RC$w_u-RC$w_min)/(RC$w_max-RC$w_min)))
    RC$B_u <- B_splines(w_u_std)
    #determine length of each part of the output, in addition to theta
    RC$desired_output <- get_desired_output('bplm',RC)
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
    }
    output_list$x[1,] <- exp(output_list$x[1,])
    output_list[['w']] <- c(RC$w,RC$w_u)
    output_list[['run_info']] <- list('nr_iter'=nr_iter,'num_chains'=num_chains,'burnin'=burnin,'thin'=thin)
    return(output_list)
}

#'Density evaluation for model2
#'
#'Evaluates the log density of the posterior distribution of the parameters of .
#'@param th A vector with length 9 containing parameters
#'@param RC A list containing prior parameters, matrices and the data.
#'@return Returns a list containing predictive values of the parameters drawn out of the evaluated density.
#'@references Birgir Hrafnkelsson, Helgi Sigurdarson and Sigurdur M. Gardarson (2015) \emph{Bayesian Generalized Rating Curves}
bplm.density_evaluation_known_c <- function(theta,RC){
    log_sig_eta2 <- theta[1]
    eta_1 <- theta[2]
    eta_minus1 <- theta[3:7]
    lambda=c(RC$P%*%as.matrix(c(eta_1,exp(log_sig_eta2)*eta_minus1)))

    l=c(log(RC$w-RC$c))

    varr=c(RC$epsilon*exp(RC$B%*%lambda))
    Sig_eps=diag(varr)
    X=cbind(1,l)
    L=t(chol(X%*%RC$Sig_x%*%t(X)+Sig_eps))
    w=solve(L,RC$y-X%*%RC$mu_x)
    p=-0.5%*%t(w)%*%w-sum(log(diag(L)))+
      pri('eta_1',eta_1=eta_1,lambda_eta_1=RC$lambda_eta_1) +
      pri('eta_minus1',eta_minus1=eta_minus1) +
      pri('sigma_eta',log_sig_eta2=log_sig_eta2,lambda_seta=RC$lambda_seta)
    W=solve(L,X%*%RC$Sig_x)
    x_u=RC$mu_x+t(chol(RC$Sig_x))%*%rnorm(nrow(RC$mu_x))
    sss=(X%*%x_u)-RC$y+sqrt(varr)*as.matrix(rnorm(RC$n))
    x=as.matrix(x_u-t(W)%*%solve(L,sss))
    yp=X%*%x
    #posterior predictive draw
    ypo=yp+as.matrix(rnorm(RC$n))*sqrt(varr)

    D=-2*sum(log(dlnorm(exp(RC$y),yp,sqrt(varr))))

    return(list("p"=p,"x"=x,"y_post"=yp,"y_post_pred"=ypo,"sigma_eps"=varr,"DIC"=D))
}

#'Density evaluation for model2
#'
#'Evaluates the log density of the posterior distribution of the parameters of model2BH.
#'@param th A vector with length 9 containing parameters
#'@param RC A list containing prior parameters, matrices and the data.
#'@return Returns a list containing predictive values of the parameters drawn out of the evaluated density.
#'@references Birgir Hrafnkelsson, Helgi Sigurdarson and Sigurdur M. Gardarson (2015) \emph{Bayesian Generalized Rating Curves}
bplm.density_evaluation_unknown_c <- function(theta,RC){
    zeta <- theta[1]
    log_sig_eta2 <- theta[2]
    eta_1 <- theta[3]
    eta_minus1 <- theta[4:8]
    lambda=c(RC$P%*%as.matrix(c(eta_1,exp(log_sig_eta2)*eta_minus1)))

    l=c(log(RC$w_tild+exp(zeta)))

    varr=c(RC$epsilon*exp(RC$B%*%lambda))
    Sig_eps=diag(varr)
    X=cbind(1,l)
    L=t(chol(X%*%RC$Sig_x%*%t(X)+Sig_eps))
    w=solve(L,RC$y-X%*%RC$mu_x)
    p=-0.5%*%t(w)%*%w-sum(log(diag(L))) +
      pri('c', zeta = zeta, lambda_c = RC$lambda_c) +
      pri('eta_1',eta_1=eta_1,lambda_eta_1=RC$lambda_eta_1) +
      pri('eta_minus1',eta_minus1=eta_minus1) +
      pri('sigma_eta',log_sig_eta2=log_sig_eta2,lambda_seta=RC$lambda_seta)

    W=solve(L,X%*%RC$Sig_x)
    x_u=RC$mu_x+t(chol(RC$Sig_x))%*%rnorm(nrow(RC$mu_x))
    sss=(X%*%x_u)-RC$y+sqrt(varr)*as.matrix(rnorm(RC$n))
    x=as.matrix(x_u-t(W)%*%solve(L,sss))
    yp=X%*%x
    #posterior predictive draw
    ypo=yp+as.matrix(rnorm(RC$n))*sqrt(varr)

    D=-2*sum(log(dlnorm(exp(RC$y),yp,sqrt(varr))))

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
bplm.predict_u_known_c <- function(theta,x,RC){
    log_sig_eta2 <- theta[1]
    eta_1 <- theta[2]
    eta_minus1 <- theta[3:7]
    lambda=c(RC$P%*%as.matrix(c(eta_1,exp(log_sig_eta2)*eta_minus1)))
    m <- length(RC$w_u)
    #calculate spline variance from B_splines
    varr_u <- c(exp(RC$B_u %*% lambda))
    l <- log(RC$w_u-RC$c)
    X <- cbind(rep(1,m),l)
    #sample from the posterior of discharge y
    yp_u <- X%*%x
    ypo_u <- yp_u  + as.matrix(rnorm(m)) * sqrt(varr_u)
    return(list('y_post'=yp_u,'y_post_pred'=ypo_u,'sigma_eps'=varr_u))
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
bplm.predict_u_unknown_c <- function(theta,x,RC){
    #collecting parameters from the MCMC sample
    zeta=theta[1]
    log_sig_eta2 <- theta[2]
    eta_1 <- theta[3]
    eta_minus1 <- theta[4:8]
    lambda=c(RC$P%*%as.matrix(c(eta_1,exp(log_sig_eta2)*eta_minus1)))
    m=length(RC$w_u)
    #calculate spline variance from B_splines
    varr_u <- c(exp(RC$B_u %*% lambda))
    above_c <- RC$w_min-exp(zeta) < RC$w_u
    m_above_c <- sum(above_c)
    #building blocks of the explanatory matrix X calculated
    l=log(RC$w_u[above_c]-RC$w_min+exp(zeta))
    X=cbind(rep(1,m_above_c),l)
    #sample from the posterior of discharge y
    yp_u <- X%*%x
    ypo_u = X%*%x + as.matrix(rnorm(m_above_c)) * sqrt(varr_u[above_c])
    return(list('y_post'=c(rep(-Inf,m-m_above_c),yp_u),'y_post_pred'=c(rep(-Inf,m-m_above_c),ypo_u),'sigma_eps'=varr_u))
}




#' Bayesian Power Law Model with varying variance
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

bplm <- function(formula,data,c_param=NULL,w_limits=NULL,country="Iceland",forcepoint=rep(FALSE,nrow(data))){
    #TODO: argument checking
    model_dat <- data[,all.vars(formula)]
    model_dat <- model_dat[order(model_dat[,2,drop=T]),]
    Q <- model_dat[,1,drop=T]
    W <- model_dat[,2,drop=T]
    MCMC_output_list <- bplm.inference(y=log(Q),w=W,c_param,w_limits,country,forcepoint)
    rating_curve <- data.frame(W=MCMC_output_list$W,as.data.frame(t(apply(MCMC_output_list$ypo,1,quantile, probs = c(0.025,0.5, 0.975),na.rm=T))),row.names=NULL)
    names(rating_curve) <- c('W','lower','median','upper')
    rating_curve <- rating_curve[order(rating_curve$W),]
    param_summary <- as.data.frame(t(apply(rbind(MCMC_output_list$a,MCMC_output_list$b,MCMC_output_list$theta),1,quantile, probs = c(0.025,0.5, 0.975),na.rm=T)))
    names(param_summary) <- c('lower','median','upper')
    param_names <- paste0('lambda_',1:6)
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
    }else{
        result_obj$post_c <- NULL
    }
    result_obj$DIC <- DIC
    result_obj$param_summary <- param_summary
    result_obj$rating_curve <- rating_curve
    attr(result_obj, "class") <- "bplm"
    return(result_obj)
}

bplm.inference <- function(y,w,c_param=NULL,w_limits=NULL,country="Iceland",forcepoint=rep(FALSE,nrow(data)),num_chains=4,nr_iter=20000,burnin=2000,thin=5){
    suppressPackageStartupMessages(require(doParallel))
    #TODO: add error message if length(formula)!=3 or if it contains more than one covariate. Also make sure that names in formula exist in data
    RC=priors(country)
    RC$s <- 3
    RC$v <- 5
    RC$y=y
    RC$w=w
    RC$w_min <- min(RC$w)
    RC$w_max <- max(RC$w)
    RC$w_tild=RC$w-RC$w_min
    RC$n=length(RC$y)
    Adist1 <- Adist(RC$w)
    RC$A <- Adist1$A
    RC$dist <- Adist1$dist
    RC$N <- Adist1$N
    RC$O <- Adist1$O

    RC$P <- diag(nrow=5,ncol=5,6)-matrix(nrow=5,ncol=5,1)
    RC$B <- B_splines(t(RC$w_tild)/RC$w_tild[length(RC$w_tild)])

    RC$epsilon <- rep(1,RC$N)
    RC$epsilon[forcepoint]=1/RC$N
    RC$c=c_param
    if(!is.null(RC$c)){
        density_fun <- bplm.density_evaluation_known_c
        unobserved_prediction_fun <- bplm.predict_u_known_c
    }else{
        density_fun <- bplm.density_evaluation_unknown_c
        unobserved_prediction_fun <- bplm.predict_u_unknown_c
    }
    #determine proposal density
    theta_length <- if(is.null(RC$c)) 7 else 6
    theta_init <- rep(0,theta_length)
    loss_fun  <-  function(th) {-density_fun(th,RC)$p}
    optim_obj <- optim(par=theta_init,loss_fun,method="L-BFGS-B",hessian=TRUE)
    t_m <- optim_obj$par
    H <- optim_obj$hessian
    LH <- t(chol(H))/0.8

    #make Wmin and Wmax divisable by 10 up, both in order to make rctafla and so l_m is defined
    if(is.null(w_limits)){
        w_max <- ceiling(max(RC$w)*10)/10
        w_min <- ceiling(10*ifelse(is.null(RC$c),min(RC$w)-exp(t_m[1]),RC$c))/10
    }else{
        w_min <- w_limits[1]
        w_max <- w_limits[2]
    }
    RC$w_u <- W_unobserved(RC,w_min,w_max)
    w_u_std <- ifelse(RC$w_u < RC$w_min,0.0,ifelse(RC$w_u>RC$w_max,1.0,(RC$w_u-RC$w_min)/(RC$w_max-RC$w_min)))
    RC$B_u <- B_splines(w_u_std)

    #MCMC parameters added, number of iterations,burnin and thin
    if(num_chains>4){
        stop('Max number of chains is 4. Please pick a lower number of chains')
    }
    cl <- makeCluster(num_chains)
    registerDoParallel(cl)
    MCMC_output_mat <- foreach(i=1:num_chains,.combine=cbind) %dopar% {
        latent_and_hyper_param <- matrix(0,nrow=length(t_m)+2,ncol=nr_iter)
        y_pred_mat <- matrix(0,nrow=RC$N,ncol=nr_iter)
        DIC <- rep(0,nr_iter)
        t_old <- as.matrix(t_m)
        Dens <- density_fun(t_old,RC)
        p_old <- Dens$p
        ypo_old <- Dens$ypo
        x_old <- Dens$x
        DIC_old <- Dens$DIC
        for(j in 1:nr_iter){
            t_new <- t_old+solve(t(LH),rnorm(length(t_m),0,1))
            Densnew <- density_fun(t_new,RC)
            x_new <- Densnew$x
            ypo_new <- Densnew$ypo
            p_new <- Densnew$p
            DIC_new <- Densnew$DIC
            logR <- p_new-p_old
            if (logR>log(runif(1))){
                t_old <- t_new
                p_old <- p_new
                ypo_old <- ypo_new
                x_old <- x_new
                DIC_old <- DIC_new
            }
            latent_and_hyper_param[,j] <- c(t_old,x_old)
            y_pred_mat[,j] <- ypo_old
            DIC[j] <- DIC_old
        }

        seq <- seq(burnin,nr_iter,thin)
        latent_and_hyper_param <- latent_and_hyper_param[,seq]
        y_pred_mat <- y_pred_mat[,seq]
        DIC <- DIC[seq]
        unobserved_mat <- apply(latent_and_hyper_param,2,function(x) unobserved_prediction_fun(x,RC))
        output_mat <- rbind(latent_and_hyper_param,y_pred_mat,unobserved_mat,t(matrix(DIC)))
        if(is.null(RC$c)){
            output_mat[1,] <- min(RC$O)-exp(output_mat[1,])
        }
        return(output_mat)
    }
    stopCluster(cl)
    MCMC_output_list <- list('W'=c(RC$w,RC$w_u),'theta'=MCMC_output_mat[1:theta_length,],
                             'a'=exp(MCMC_output_mat[theta_length+1,]),
                             'b'=MCMC_output_mat[theta_length+2,],
                             'ypo'=exp(MCMC_output_mat[(theta_length+2+1):(nrow(MCMC_output_mat)-1),]),
                             'DIC'=MCMC_output_mat[nrow(MCMC_output_mat),])

    return(MCMC_output_list)
}

#'Density evaluation for model2
#'
#'Evaluates the log density of the posterior distribution of the parameters of .
#'@param th A vector with length 9 containing parameters
#'@param RC A list containing prior parameters, matrices and the data.
#'@return Returns a list containing predictive values of the parameters drawn out of the evaluated density.
#'@references Birgir Hrafnkelsson, Helgi Sigurdarson and Sigurdur M. Gardarson (2015) \emph{Bayesian Generalized Rating Curves}
bplm.density_evaluation_known_c <- function(th,RC){
    lambda=th[1:6]

    f=lambda[1:5]-lambda[6]
    l=c(log(RC$w-RC$c))

    varr=c(RC$epsilon*exp(RC$B%*%lambda))
    Sig_eps=diag(varr)
    X=cbind(1,l)
    L=t(chol(X%*%RC$Sig_x%*%t(X)+Sig_eps))
    w=solve(L,RC$y-X%*%RC$mu_x)
    p=-0.5%*%t(w)%*%w-sum(log(diag(L)))-
        (RC$v+5-1)/2*log(RC$v*RC$s+f%*%RC$P%*%f) #63 micro
    W=solve(L,X%*%RC$Sig_x)
    x_u=RC$mu_x+t(chol(RC$Sig_x))%*%rnorm(nrow(RC$mu_x))
    sss=(X%*%x_u)-RC$y+sqrt(varr)*as.matrix(rnorm(RC$N))
    x=as.matrix(x_u-t(W)%*%solve(L,sss))
    yp=X%*%x
    #posterior predictive draw
    ypo=yp+as.matrix(rnorm(RC$N))*sqrt(varr)

    D=-2*sum(log(dlnorm(exp(RC$y),yp,sqrt(varr))))

    return(list("p"=p,"x"=x,"yp"=yp,"ypo"=ypo,"varr"=varr,"DIC"=D))
}

#'Density evaluation for model2
#'
#'Evaluates the log density of the posterior distribution of the parameters of model2BH.
#'@param th A vector with length 9 containing parameters
#'@param RC A list containing prior parameters, matrices and the data.
#'@return Returns a list containing predictive values of the parameters drawn out of the evaluated density.
#'@references Birgir Hrafnkelsson, Helgi Sigurdarson and Sigurdur M. Gardarson (2015) \emph{Bayesian Generalized Rating Curves}
bplm.density_evaluation_unknown_c <- function(th,RC){
    zeta=th[1]
    lambda=th[2:7]

    f=lambda[1:5]-lambda[6]
    l=c(log(RC$w_tild+exp(th[1])))

    varr=c(RC$epsilon*exp(RC$B%*%lambda))
    Sig_eps=diag(varr)
    X=cbind(1,l)
    L=t(chol(X%*%RC$Sig_x%*%t(X)+Sig_eps))
    w=solve(L,RC$y-X%*%RC$mu_x)
    p=-0.5%*%t(w)%*%w-sum(log(diag(L)))-
        (RC$v+5-1)/2*log(RC$v*RC$s+f%*%RC$P%*%f)+zeta-exp(zeta)/RC$mu_c #63 micro
    W=solve(L,X%*%RC$Sig_x)
    x_u=RC$mu_x+t(chol(RC$Sig_x))%*%rnorm(nrow(RC$mu_x))
    sss=(X%*%x_u)-RC$y+sqrt(varr)*as.matrix(rnorm(RC$N))
    x=as.matrix(x_u-t(W)%*%solve(L,sss))
    yp=X%*%x
    #posterior predictive draw
    ypo=yp+as.matrix(rnorm(RC$N))*sqrt(varr)

    D=-2*sum(log(dlnorm(exp(RC$y),yp,sqrt(varr))))

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
bplm.predict_u_known_c <- function(param,RC){
    #collecting parameters from the MCMC sample
    th <- param[1:6]
    x <- param[7:8]
    lambda <- th[1:6]
    m <- length(RC$w_u)
    #calculate spline variance from B_splines
    varr <- c(exp(RC$B_u %*% lambda))
    l=log(RC$w_u-RC$c)
    X=cbind(rep(1,m),l)
    #sample from the posterior of discharge y
    ypo_u = X%*%x + as.matrix(rnorm(m)) * sqrt(varr)
    return(ypo_u)
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
bplm.predict_u_unknown_c <- function(param,RC){
    #collecting parameters from the MCMC sample
    th <- param[1:7]
    x <- param[8:9]
    zeta <- th[1]
    lambda <- th[2:7]
    #calculate spline variance from B_splines
    varr <- c(exp(RC$B_u %*% lambda))
    m <- length(RC$w_u)
    above_c <- -(exp(zeta)-RC$w_min) < RC$w_u
    m_above_c <- sum(above_c)
    #buidling blocks of the explanatory matrix X calculated
    l <- log(RC$w_u[above_c]-RC$w_min+exp(zeta))
    X <- cbind(rep(1,m_above_c),l)
    #sample from the posterior of discharge y
    ypo_u = X%*%x + as.matrix(rnorm(m_above_c)) * sqrt(varr[above_c])
    return(as.matrix(c(rep(-Inf,m-m_above_c),ypo_u)))
}




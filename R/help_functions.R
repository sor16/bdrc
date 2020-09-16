library(parallel)

#' Prior parameter specification
#'
#' Specifies the prior parameters into \code{\link{model1BH}} and \code{\link{model2BH}}
#'@param country A string with the name of the country of which the prior parameters into the models should be specified for
#'@return Returns a list of prior parameters into \code{\link{model1BH}} and \code{\link{model2BH}}.
#'The priors are based from data from rivers of a given country.If you want to add your country to this function,
#' please contact the developers at sor16@@hi.is or aoj8@@hi.is.
#'@references Birgir Hrafnkelsson, Helgi Sigurdarson and Sigurdur M. Gardarson (2015) \emph{Bayesian Generalized Rating Curves}
priors <- function(model,c_param) {
    RC=list()
    #Prior parameters for all models
    RC$mu_a <- 3.20;
    RC$mu_b <- 2.29;
    RC$sig_a <- sqrt(1.21);
    RC$sig_b <- sqrt(0.48);
    RC$p_ab <- -0.61;
    if(is.null(c_param)){
        RC$mu_c <- 1.9;
    }else{
        RC$c <- c_param
    }
    #Prior parameters depending on model
    if(model %in% c('bplm0','bplm')){
        RC$Sig_x <- rbind(c(RC$sig_a^2, RC$p_ab*RC$sig_a*RC$sig_b), c(RC$p_ab*RC$sig_a*RC$sig_b, RC$sig_b^2))
        RC$mu_x <- as.matrix(c(RC$mu_a, RC$mu_b))
        RC$Sig_xinv <- solve(RC$Sig_x)
        RC$Sinvmu <- RC$Sig_xinv%*%RC$mu_x
    }else{
        RC$nugget <- 10^-8
        RC$mu_sb <- 0.5
        RC$mu_pb <- 0.5
        RC$tau_pb2 <- 0.25^2
        RC$Sig_ab <- rbind(c(RC$sig_a^2, RC$p_ab*RC$sig_a*RC$sig_b), c(RC$p_ab*RC$sig_a*RC$sig_b, RC$sig_b^2))
    }
    if(model %in% c('bplm','bgplm')){
        RC$s <- 3
        RC$v <- 5
    }
    return(RC)
}

#'Linking unique water level measurements to actual
#'water level measurements
#'
#'Adist links unique water level measurements (\strong{w'}) to actual
#'water level measurements (w) such that \strong{w}=\strong{Aw'}.
#'from the measurements.
#'@param w full stage measurements
#'@return
#'\itemize{
#'\item A: Matrix \strong{A} linking unique water level measurements (\strong{w'}) to actual
#'water level measurements (w) such that \strong{w}=\strong{Aw'}
#'}
#'@references Birgir Hrafnkelsson, Helgi Sigurdarson and Sigurdur M. Gardarson (2015) \emph{Bayesian Generalized Rating Curves}
create_A <- function(w){
    n <- length(w)
    A=matrix(0,nrow=n,ncol=length(unique(w)))
    A[1,1]=1
    i=1
    for(ii in 2:n){
        if(w[ii]==w[ii-1]){
            A[ii,i]=1
        }else{
            i=i+1
            A[ii,i]=1
        }
    }
    return(A)
}


initiate_output_list <- function(desired_output,nr_iter){
    output_list <- list()
    for(elem in names(desired_output)){
        output_list[[elem]] <- matrix(0,nrow=desired_output[[elem]][['observed']]+desired_output[[elem]][['unobserved']],ncol=nr_iter)
    }
    return(output_list)
}

run_MCMC <- function(theta_m,RC,density_fun,unobserved_prediction_fun,nr_iter=20000,num_chains=4,burnin=2000,thin=5){
    theta_mat <- matrix(0,nrow=RC$theta_length,ncol=nr_iter)
    output_list <- initiate_output_list(RC$desired_output,nr_iter)
    density_eval_m <- density_fun(theta_m,RC)
    theta_old <- theta_m
    density_eval_old <- density_eval_m
    for(i in 1:nr_iter){
        theta_new <- theta_old+solve(t(RC$LH),rnorm(RC$theta_length,0,1))
        density_eval_new <- density_fun(theta_new,RC)
        logR <- density_eval_new[['p']]-density_eval_old[['p']]
        if (logR>log(runif(1))){
            theta_old <- theta_new
            density_eval_old <- density_eval_new
        }
        theta_mat[,i] <- theta_old
        for(elem in names(RC$desired_output)){
            output_list[[elem]][1:RC$desired_output[[elem]][['observed']],i] <- density_eval_old[[elem]]
        }
    }
    idx <- seq(burnin,nr_iter,thin)
    theta_mat <- theta_mat[,idx,drop=F]
    output_list <- sapply(output_list,FUN=function(x) x[,idx,drop=F],simplify=F,USE.NAMES=T)
    for(i in 1:ncol(theta_mat)){
        unobserved_list <- unobserved_prediction_fun(theta_mat[,i],output_list[['x']][1:(RC$desired_output[['x']][['observed']]),i],RC)
        for(elem in names(unobserved_list)){
            output_list[[elem]][(RC$desired_output[[elem]][['observed']]+1):nrow(output_list[[elem]]),i] <- unobserved_list[[elem]]
        }
    }
    output_list[['theta']] <- theta_mat
    return(output_list)
}

get_MCMC_summary <- function(X,w=NULL){
    summary_dat <- as.data.frame(t(apply(X,1,quantile, probs = c(0.025,0.5, 0.975),na.rm=T)))
    names(summary_dat) <- c('lower','median','upper')
    if(!is.null(w)){
        summary_dat <- data.frame(w=w,summary_dat,row.names=NULL)
        summary_dat <- summary_dat[order(summary_dat$w),]
    }
    return(summary_dat)
}

get_param_names <- function(model,c_param){
    if(model=='bplm0'){
        hyper_param <- 'sigma_eps'
    }else if(model=='bplm'){
        hyper_param <- paste('lambda',1:6,sep='_')
    }else if(model=='bgplm0'){
        hyper_param <- c('sigma_eps','sigma_beta','phi_beta')
    }else if(model=='bgplm'){
        hyper_param <- c('sigma_beta','phi_beta',paste('lambda',1:6,sep='_'))
    }
    if(is.null(c_param)){
        hyper_param <- c('c',hyper_param)
    }
    return(c('a','b',hyper_param))
}

get_desired_output <- function(model,RC){
    const_var <- model %in% c('bplm0','bgplm0')
    const_b <- model %in% c('bplm0','bplm')
    desired_output <- list('y_post'=list('observed'=RC$n,'unobserved'=RC$n_u),
                           'y_post_pred'=list('observed'=RC$n,'unobserved'=RC$n_u),
                           'DIC'=list('observed'=1,'unobserved'=0))
    if(!const_var){
        desired_output$sigma_eps <- list('observed'=RC$n,'unobserved'=RC$n_u)
    }
    if(!const_b){
        desired_output$x <- list('observed'=2+RC$n_unique,'unobserved'=RC$n_u)
    }else{
        desired_output$x <- list('observed'=2,'unobserved'=0)
    }
    return(desired_output)
}

pri <- function(type,...){
  args = list(...)
  if(type == 'c'){
    p = args$zeta-exp(args$zeta)/args$mu_c
    
  }else if(type == 'sig_eps2'){
    p = 0
    
  }else if(type == 'sig_b2'){
    p = args$sig_b2-exp(args$sig_b2)/args$mu_sb
    
  }else if(type == 'phi_b'){
    p = -(0.5/args$tau_pb2*(args$phi_b-args$mu_pb)^2)
    
  }else if(type == 'eta'){
    p = -(args$v+5-1)/2*log(args$v*args$s+args$f%*%args$P%*%args$f)
    
  }
  return(p)
}

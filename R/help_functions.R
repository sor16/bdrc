priors <- function(model,c_param=NULL) {
    RC=list()
    #Prior parameters for all models
    RC$mu_a <- 3;
    RC$mu_b <- 1.835;
    RC$sig_a <- 3;
    RC$p_ab <- 0;
    RC$nugget <- 10^-8
    if(is.null(c_param)){
        RC$lambda_c <- 2;
    }else{
        RC$c <- c_param
    }
    #if f(h)=b vs f(h)=b+beta(h)
    if(model %in% c('plm0','plm')){
        RC$sig_b <- 0.426;
        RC$Sig_x <- rbind(c(RC$sig_a^2, RC$p_ab*RC$sig_a*RC$sig_b), c(RC$p_ab*RC$sig_a*RC$sig_b, RC$sig_b^2))
        RC$mu_x <- as.matrix(c(RC$mu_a, RC$mu_b))
        RC$Sig_xinv <- solve(RC$Sig_x)
        RC$Sinvmu <- RC$Sig_xinv%*%RC$mu_x
    }else{
        RC$sig_b <- 0.01;
        RC$lambda_sb <- 5.405
        RC$lambda_pb <- 3.988
    }
    RC$Sig_ab <- rbind(c(RC$sig_a^2, RC$p_ab*RC$sig_a*RC$sig_b), c(RC$p_ab*RC$sig_a*RC$sig_b, RC$sig_b^2))
    #if fixed variance vs not fixed variance
    if(model %in% c('plm0','gplm0')){
      RC$lambda_se <- 28.78
    }else{
      RC$lambda_eta_1 <- 28.78
      RC$lambda_seta <- 8.62
    }
    return(RC)
}

get_model_components <- function(model,y,h,c_param,h_max,forcepoint,h_min){
  RC <- priors(model,c_param)
  RC$y <- as.matrix(y)
  RC$h <- h
  RC$h_min <- if(is.null(h_min)) min(RC$h) else h_min
  RC$h_max <- max(RC$h)
  RC$n <- length(h)
  RC$epsilon <- rep(1,RC$n)
  RC$epsilon[forcepoint] <- 1/RC$n
  if(model %in% c('plm','gplm')){
    RC$P <- lower.tri(matrix(rep(1,36),6,6),diag=TRUE)*1
    h_tilde <- RC$h-min(RC$h)
    RC$B <- B_splines(t(h_tilde)/h_tilde[RC$n])
  }
  if(model %in% c('gplm0','gplm')){
    RC$y <- rbind(RC$y,RC$mu_b)
    RC$h_unique <- unique(RC$h)
    RC$n_unique <- length(RC$h_unique)
    RC$A <- create_A(RC$h)
    RC$dist <- as.matrix(dist(c(RC$h_unique)))
    RC$mu_x <- as.matrix(c(RC$mu_a,RC$mu_b, rep(0,RC$n_unique)))
    RC$Z <- cbind(t(c(0,1)),t(rep(0,RC$n_unique)))
    RC$m1 <- matrix(0,nrow=2,ncol=RC$n_unique)
    RC$m2 <- matrix(0,nrow=RC$n_unique,ncol=2)
  }
  if(!is.null(RC$c)){
    density_fun_name <- paste0(model,'.density_evaluation_known_c')
    unobserved_prediction_fun_name <- paste0(model,'.predict_u_known_c')
  }else{
    density_fun_name <- paste0(model,'.density_evaluation_unknown_c')
    unobserved_prediction_fun_name <- paste0(model,'.predict_u_unknown_c')
  }
  RC$density_fun <- get(density_fun_name)
  RC$unobserved_prediction_fun <- get(unobserved_prediction_fun_name)
  theta_length_vec=c('plm0'=2,'plm'=8,'gplm0'=4,'gplm'=10)
  #determine proposal density
  RC$theta_length <- if(is.null(RC$c)) theta_length_vec[model] else theta_length_vec[model]-1
  theta_init <- rep(0,RC$theta_length)
  loss_fun <- function(theta) {-RC$density_fun(theta,RC)$p}
  optim_obj <- optim(par=theta_init,loss_fun,method="L-BFGS-B",hessian=TRUE)
  RC$theta_m <- optim_obj$par
  RC$H <- optim_obj$hessian
  proposal_scaling <- 2.38^2/RC$theta_length
  RC$LH <- t(chol(RC$H))/sqrt(proposal_scaling)
  h_min_pred <- ifelse(is.null(RC$c),RC$h_min-exp(RC$theta_m[1]),RC$c)
  h_max_pred <- h_max
  if(is.null(h_max_pred)){
    h_max_pred <- RC$h_max
  }else if(h_max_pred<RC$h_max){
    stop(paste0('maximum stage value must be larger than the maximum stage value in the data, which is ', RC$h_max,' m'))
  }
  RC$h_u <- h_unobserved(RC,h_min_pred,h_max_pred)
  RC$n_u <- length(RC$h_u)
  if(model %in% c('plm','gplm')){
    h_u_std <- ifelse(RC$h_u < min(RC$h),0.0,ifelse(RC$h_u>RC$h_max,1.0,(RC$h_u-min(RC$h))/(RC$h_max-min(RC$h))))
    RC$B_u <- B_splines(h_u_std)
  }
  #determine length of each part of the output, in addition to theta
  RC$desired_output <- get_desired_output(model,RC)
  return(RC)
}

create_A <- function(h){
    n <- length(h)
    A=matrix(0,nrow=n,ncol=length(unique(h)))
    A[1,1]=1
    i=1
    for(ii in 2:n){
        if(h[ii]==h[ii-1]){
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

calc_variogram <- function(i,param_mat,burnin=2000,nr_iter=20000){
  rowSums((param_mat[,(i+1):ncol(param_mat)] - param_mat[,1:(ncol(param_mat)-i)])^2)
}

#' @importFrom stats rnorm runif
run_MCMC <- function(theta_m,RC,density_fun,unobserved_prediction_fun,nr_iter=20000,burnin=2000,thin=5,T_max=50){
    theta_mat <- matrix(0,nrow=RC$theta_length,ncol=nr_iter)
    output_list <- initiate_output_list(RC$desired_output,nr_iter)
    density_eval_m <- density_fun(theta_m,RC)
    theta_old <- theta_m
    density_eval_old <- density_eval_m
    acceptance_vec <- rep(FALSE,nr_iter)
    for(i in 1:nr_iter){
        theta_new <- theta_old+solve(t(RC$LH),rnorm(RC$theta_length,0,1))
        density_eval_new <- density_fun(theta_new,RC)
        logR <- density_eval_new[['p']]-density_eval_old[['p']]
        if (logR>log(runif(1))){
            acceptance_vec[i] <- TRUE
            theta_old <- theta_new
            density_eval_old <- density_eval_new
        }
        theta_mat[,i] <- theta_old
        for(elem in names(RC$desired_output)){
            output_list[[elem]][1:RC$desired_output[[elem]][['observed']],i] <- density_eval_old[[elem]]
        }
    }
    param_mat <- rbind(output_list$x[1:2,],theta_mat)
    split_idx <- round(0.5*(nr_iter-burnin))
    param_mat1 <- param_mat[,seq(burnin,burnin+split_idx)]
    param_mat2 <- param_mat[,seq(burnin+split_idx+1,nr_iter)]
    param_mean <- cbind(rowMeans(param_mat1),rowMeans(param_mat2))
    param_var <- cbind(apply(param_mat1,1,var),apply(param_mat2,1,var))
    variogram_chain <- cbind(sapply(1:T_max,calc_variogram,param_mat1,burnin,nr_iter),
                             sapply(1:T_max,calc_variogram,param_mat2,burnin,nr_iter))
    rm(param_mat,param_mat1,param_mat2)
    idx <- seq(burnin,nr_iter,thin)
    acceptance_vec <- acceptance_vec[idx]
    theta_mat <- theta_mat[,idx,drop=FALSE]
    output_list <- sapply(output_list,FUN=function(x) x[,idx,drop=FALSE],simplify=FALSE,USE.NAMES=TRUE)
    for(i in 1:ncol(theta_mat)){
        unobserved_list <- unobserved_prediction_fun(theta_mat[,i],output_list[['x']][1:(RC$desired_output[['x']][['observed']]),i],RC)
        for(elem in names(unobserved_list)){
            output_list[[elem]][(RC$desired_output[[elem]][['observed']]+1):nrow(output_list[[elem]]),i] <- unobserved_list[[elem]]
        }
    }
    output_list$theta <- theta_mat
    output_list$param_mean <- param_mean
    output_list$param_var <- param_var
    output_list$variogram_chain <- variogram_chain
    output_list$acceptance_vec <- t(as.matrix(acceptance_vec))
    return(output_list)
}
#' @importFrom parallel detectCores makeCluster clusterSetRNGStream clusterExport parLapply stopCluster
get_MCMC_output_list <- function(theta_m,RC,density_fun,unobserved_prediction_fun,parallel,num_cores=NULL,num_chains=4,nr_iter=20000,burnin=2000,thin=5){
  #pb <- utils::txtProgressBar(min=0, max=nr_iter*(1 + (1-parallel)*(num_chains-1)),style=3)
  if(num_chains>4){
    stop('Max number of chains is 4. Please pick a lower number of chains')
  }
  T_max <- 50
  run_MCMC_wrapper <- function(i){
    run_MCMC(theta_m=theta_m,RC=RC,density_fun=density_fun,
             unobserved_prediction_fun=unobserved_prediction_fun,
             nr_iter=nr_iter,burnin=burnin,thin=thin,T_max=T_max)
  }
  if(parallel){
    num_cores_on_device <- detectCores()
    num_cores_default <- min(num_cores_on_device,num_chains)
    if (!is.null(num_cores)) {
      if(!(num_cores %in% 1:num_chains)){
        num_cores <- num_cores_default
        warning(paste0('num_cores argument used must be an integer between 1 and ',num_chains,' (the number of chains). Using ',num_cores_default,' cores instead.'))
      }
    }else{
      num_cores <- num_cores_default
    }
    cl <- makeCluster(num_cores,setup_strategy='sequential')
    clusterSetRNGStream(cl=cl) #set RNG to type L'Ecuyer
    clusterExport(cl,c('run_MCMC','initiate_output_list','pri','calc_variogram',
                       'theta_m','RC','density_fun','unobserved_prediction_fun',
                       'parallel','nr_iter','burnin','thin','T_max'),envir = environment())
    MCMC_output_list <- parLapply(cl,1:num_chains,run_MCMC_wrapper)
    stopCluster(cl)
  }else{
    MCMC_output_list <- lapply(1:num_chains,run_MCMC_wrapper)
  }
  output_list <- list()
  for(elem in names(MCMC_output_list[[1]])){
    output_list[[elem]] <- do.call(cbind,lapply(1:num_chains,function(i) MCMC_output_list[[i]][[elem]]))
  }
  m <- 2*num_chains
  n <- round(0.5*(nr_iter-burnin))
  variogram <- sapply(1:T_max,function(i) rowMeans(output_list$variogram_chain[,T_max*seq(0,m-1)+i])/(n-i))
  between_chain_var <- n*apply(output_list$param_mean,1,var)
  within_chain_var <- rowMeans(output_list$param_var)
  chain_var_hat <- ((n-1)*within_chain_var + between_chain_var)/n
  output_list$r_hat <- sqrt(chain_var_hat/within_chain_var)
  output_list$autocorrelation <- 1-variogram/(2*matrix(rep(chain_var_hat,T_max),nrow=nrow(variogram)))
  output_list$num_effective_samples <-round(m*n/(1+2*rowSums(output_list$autocorrelation)))
  output_list$acceptance_rate <- sum(output_list$acceptance_vec)/ncol(output_list$acceptance_vec)
  output_list$param_mean <- output_list$param_var <- output_list$variogram_chain <- NULL
  return(output_list)
}

#' @importFrom stats quantile
get_MCMC_summary <- function(X,h=NULL){
    summary_dat <- apply(X,1,function(x) {
                      c(quantile(x,probs=0.025,na.rm=TRUE),
                        quantile(x,probs=0.5,na.rm=TRUE),
                        quantile(x,probs=0.975,na.rm=TRUE))
                    })
    summary_dat <- as.data.frame(t(summary_dat))
    names(summary_dat) <- c('lower','median','upper')
    if(!is.null(h)){
        summary_dat <- data.frame(h=h,summary_dat,row.names=NULL)
    }
    return(summary_dat)
}

get_param_names <- function(model,c_param){
    if(model=='plm0'){
        hyper_param <- 'sigma_eps'
    }else if(model=='plm'){
        hyper_param <- c('sigma_eta',paste('eta',1:6,sep='_'))
    }else if(model=='gplm0'){
        hyper_param <- c('sigma_eps','sigma_beta','phi_beta')
    }else if(model=='gplm'){
        hyper_param <- c('sigma_beta','phi_beta','sigma_eta',paste('eta',1:6,sep='_'))
    }
    if(is.null(c_param)){
        hyper_param <- c('c',hyper_param)
    }
    return(c('a','b',hyper_param))
}

get_param_expression <- function(param){
  expr_vec <- c('a'='a','b'='b','c'='c','sigma_eps'='sigma[epsilon]',
                'sigma_beta'='sigma[beta]','phi_beta'='phi[beta]',
                'sigma_eta'='sigma[eta]','eta_1'='eta[1]','eta_2'='eta[2]',
                'eta_3'='eta[3]','eta_4'='eta[4]','eta_5'='eta[5]',
                'eta_6'='eta[6]','log(a)'='log(a)','log(h_min-c)'='log(h[min]-c)',
                '2log(sigma_eps)'='log(sigma[epsilon]^2)',
                'log(sigma_beta)'='log(sigma[beta])',
                'log(phi_beta)'='log(phi[beta])',
                'log(sigma_eta)'='log(sigma[eta])',
                'z_1'='z[1]','z_2'='z[2]','z_3'='z[3]',
                'z_4'='z[4]','z_5'='z[5]','z_6'='z[6]')
  param_expr <- expr_vec[param]
  if(any(is.na(param_expr))){
    stop('param not found')
  }
  return(param_expr)
}

get_args_rollout <- function(args,param_vec){
  rollout <- unlist(sapply(args,function(x) {
    if(x=='latent_parameters'){
      return(param_vec[1:2])
    }else if(x=='hyperparameters'){
      return(param_vec[3:length(param_vec)])
    }else{
      return(x)
    }
  }))
  return(unique(rollout))
}

get_parameter_levels <- function(param_vec){
    order_vec <- c('a'=1,'log(a)'=2,'b'=3,'c'=4,'log(h_min-c)'=5,'sigma_eps'=6,
                   '2log(sigma_eps)'=7,'sigma_beta'=8,'log(sigma_beta)'=9,
                   'phi_beta'=10,'log(phi_beta)'=11,'sigma_eta'=12,'log(sigma_eta)'=13,
                   'eta_1'=14,'eta_2'=15,'z_1'=16,'eta_3'=17,'z_2'=18,
                   'eta_4'=19,'z_3'=20,'eta_5'=21,'z_4'=22,'eta_6'=23,'z_5'=24)
  return(names(sort(rank(order_vec[param_vec]))))
}

get_transformed_param <- function(v,param_name,mod,...){
  args <- list(...)
  # fun_vec <- c('a'=log,
  #              'b'=identity,
  #              'c'=function(x,h_min) log())
  if(param_name=='a'){
    out_v <- log(v)
    names(out_v) <- rep('log(a)',length(v))
  }else if(param_name=='b'){
    out_v <- v
    names(out_v) <- rep('b',length(v))
  }else if(param_name=='c'){
    out_v <- log(args$h_min-v)
    names(out_v) <- rep('log(h_min-c)',length(v))
  }else if(param_name=='sigma_eps'){
    out_v <- 2*log(v)
    names(out_v) <- rep('2log(sigma_eps)',length(v))
  }else if(param_name=='sigma_beta'){
    out_v <- log(v)
    names(out_v) <- rep('log(sigma_beta)',length(v))
  }else if(param_name=='phi_beta'){
    out_v <- log(v)
    names(out_v) <- rep('log(phi_beta)',length(v))
  }else if(param_name=='sigma_eta'){
    out_v <- log(v)
    names(out_v) <- rep('log(sigma_eta)',length(v))
  }else if(param_name=='eta_1'){
    out_v <- v
    names(out_v) <- rep('eta_1',length(v))
  }else if(param_name %in% paste0('eta_',2:6)){
    eta_nr <- as.numeric(unlist(strsplit(param_name,split='_'))[2])
    out_v <- v-mod[[paste0('eta_',eta_nr-1,'_posterior')]]
    names(out_v) <- rep(paste0('z_',eta_nr-1),length(v))
  }else{
    stop('param not found')
  }
  return(out_v)
}

get_desired_output <- function(model,RC){
    const_var <- model %in% c('plm0','gplm0')
    const_b <- model %in% c('plm0','plm')
    desired_output <- list('y_post'=list('observed'=RC$n,'unobserved'=RC$n_u),
                           'y_post_pred'=list('observed'=RC$n,'unobserved'=RC$n_u),
                           'D'=list('observed'=1,'unobserved'=0))
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
    p <- args$zeta - exp(args$zeta)*args$lambda_c
  }else if(type == 'sigma_eps2'){
    p <- 0.5*args$log_sig_eps2 - exp(0.5*args$log_sig_eps2)*args$lambda_se
  }else if(type == 'sigma_b'){
    p <- args$log_sig_b - exp(args$log_sig_b)*args$lambda_sb
  }else if(type == 'phi_b'){
    p <- - 0.5*args$log_phi_b - args$lambda_pb*sqrt(0.5)*exp(-0.5*args$log_phi_b)
  }else if(type == 'eta_1'){
    p <- 0.5*args$eta_1 - exp(0.5*args$eta_1)*args$lambda_eta_1
  }else if(type == 'eta_minus1'){
    p <- -0.5*t(as.matrix(args$z))%*%as.matrix(args$z)
  }else if(type == 'sigma_eta'){
    p <- args$log_sig_eta - exp(args$log_sig_eta)*args$lambda_seta
  }
  return(p)
}

h_unobserved <- function(RC,h_min=NA,h_max=NA){
  h_u=NULL
  h=100*c(RC$h) #work in cm
  max_h_diff=5
  #distance between subsequent elements in vector with additional dummy point 1000
  distvect=abs(h-c(h[2:length(h)],1000))
  #add datapoints to corresponding distances to see range of distance
  distwithdata=rbind(h,distvect,c(h[2:length(h)],1000))
  distfilter=distwithdata[,distvect>max_h_diff,drop=FALSE]
  #remove dummy distance
  distfilter=as.matrix(distfilter[,-ncol(distfilter)])
  if(ncol(distfilter)!=0){
    #make sequence from the ranges with length.out equal to corresponding elelement in distvect
    h_u=0.01*unlist(apply(distfilter,2,FUN=function(x){setdiff(seq(x[1],x[3],length.out=2+ceiling(x[2]/max_h_diff)),c(x[1],x[3]))
    }))
  }
  h_before_data=setdiff(seq(h_min,RC$h_min,by=0.05),c(RC$h_min))
  h_after_data=setdiff(seq(RC$h_max,h_max,length.out=2+ceiling(20*(h_max-RC$h_max))),RC$h_max)
  h_u=c(h_before_data,h_u,h_after_data)
  return(h_u)
}

B_splines <- function(ZZ){
  #The number of equally spaced interior knots.
  kx=2
  #Delta x and Delta y.
  dx=1/(kx+1)
  #The order of the splines.
  M = 4
  #Determine the number of functions.
  Nx = kx + M
  #The epsilon-knots
  epsilon_x = dx*seq(0,kx+1,by=1)
  #the tau-knots.
  tau_x = matrix(0,nrow=1,ncol=(kx+2*M))
  tau_x[1:M] = epsilon_x[1]*matrix(1,nrow=1,ncol=M)
  tau_x[(M+1):(kx+M)]=epsilon_x[2:(kx+1)]
  tau_x[(kx+M+1):(kx+2*M)]=epsilon_x[kx+2]*matrix(1,nrow=1,ncol=M)
  #Vector with values of x and y.
  lx = length(ZZ)
  #Compute the x-splines and the y-splines.
  XX = matrix(0,nrow=(kx+M),ncol=length(ZZ))
  # i = 1
  XX[1,] = (1/dx^3)*(tau_x[M+1]-ZZ)*(tau_x[M+1]-ZZ)*(tau_x[M+1]-ZZ)*(tau_x[M]<=ZZ)*(ZZ<tau_x[M+1]);
  # i = 2
  XX[2,] = (1/dx^3)*(ZZ-tau_x[2])*(tau_x[M+1]-ZZ)*(tau_x[M+1]-ZZ)*(tau_x[M]<=ZZ)*(ZZ<tau_x[M+1])+
    (1/2/dx^3)*(tau_x[M+2]-ZZ)*(ZZ-tau_x[3])*(tau_x[M+1]-ZZ)*(tau_x[M]<=ZZ)*(ZZ<tau_x[M+1])+
    (1/4/dx^3)*(tau_x[M+2]-ZZ)*(tau_x[M+2]-ZZ)*(ZZ-tau_x[M])*(tau_x[M]<=ZZ)*(ZZ<tau_x[M+1])+
    (1/4/dx^3)*(tau_x[M+2]-ZZ)*(tau_x[M+2]-ZZ)*(tau_x[M+2]-ZZ)*(tau_x[M+1]<=ZZ)*(ZZ<tau_x[M+2])
  # i = 3
  XX[3,] = (1/2/dx^3)*(ZZ-tau_x[3])*(ZZ-tau_x[3])*(tau_x[M+1]-ZZ)*(tau_x[M]<=ZZ)*(ZZ<tau_x[M+1])+
    (1/4/dx^3)*(ZZ-tau_x[3])*(tau_x[M+2]-ZZ)*(ZZ-tau_x[M])*(tau_x[M]<=ZZ)*(ZZ<tau_x[M+1])+
    (1/4/dx^3)*(ZZ-tau_x[3])*(tau_x[M+2]-ZZ)*(tau_x[M+2]-ZZ)*(tau_x[M+1]<=ZZ)*(ZZ<tau_x[M+2])+
    (1/6/dx^3)*(tau_x[M+3]-ZZ)*(ZZ-tau_x[M])*(ZZ-tau_x[M])*(tau_x[M]<=ZZ)*(ZZ<tau_x[M+1])+
    (1/6/dx^3)*(tau_x[M+3]-ZZ)*(ZZ-tau_x[M])*(tau_x[M+2]-ZZ)*(tau_x[M+1]<=ZZ)*(ZZ<tau_x[M+2])+
    (1/6/dx^3)*(tau_x[M+3]-ZZ)*(tau_x[M+3]-ZZ)*(ZZ-tau_x[M+1])*(tau_x[M+1]<=ZZ)*(ZZ<tau_x[M+2])+
    (1/6/dx^3)*(tau_x[M+3]-ZZ)*(tau_x[M+3]-ZZ)*(tau_x[M+3]-ZZ)*(tau_x[M+2]<=ZZ)*(ZZ<tau_x[M+3])
  # i = kx + 2
  XX[kx+2,] =  -(1/6/dx^3)*(tau_x[kx+2]-ZZ)*(tau_x[kx+2]-ZZ)*(tau_x[kx+2]-ZZ)*(tau_x[kx+2]<=ZZ)*(ZZ<tau_x[kx+3])-
    (1/6/dx^3)*(tau_x[kx+2]-ZZ)*(tau_x[kx+2]-ZZ)*(ZZ-tau_x[kx+4])*(tau_x[kx+3]<=ZZ)*(ZZ<tau_x[kx+4])-
    (1/6/dx^3)*(tau_x[kx+2]-ZZ)*(ZZ-tau_x[kx+5])*(tau_x[kx+3]-ZZ)*(tau_x[kx+3]<=ZZ)*(ZZ<tau_x[kx+4])-
    (1/6/dx^3)*(tau_x[kx+2]-ZZ)*(ZZ-tau_x[kx+5])*(ZZ-tau_x[kx+5])*(tau_x[kx+4]<=ZZ)*(ZZ<tau_x[kx+5])-
    (1/4/dx^3)*(ZZ-tau_x[kx+6])*(tau_x[kx+3]-ZZ)*(tau_x[kx+3]-ZZ)*(tau_x[kx+3]<=ZZ)*(ZZ<tau_x[kx+4])-
    (1/4/dx^3)*(ZZ-tau_x[kx+6])*(tau_x[kx+3]-ZZ)*(ZZ-tau_x[kx+5])*(tau_x[kx+4]<=ZZ)*(ZZ<tau_x[kx+5])-
    (1/2/dx^3)*(ZZ-tau_x[kx+6])*(ZZ-tau_x[kx+6])*(tau_x[kx+4]-ZZ)*(tau_x[kx+4]<=ZZ)*(ZZ<tau_x[kx+5])
  # i = kx + 3
  XX[kx+3,] = - (1/4/dx^3)*(tau_x[kx+3]-ZZ)*(tau_x[kx+3]-ZZ)*(tau_x[kx+3]-ZZ)*(tau_x[kx+3]<=ZZ)*(ZZ<tau_x[kx+4])-
    (1/4/dx^3)*(tau_x[kx+3]-ZZ)*(tau_x[kx+3]-ZZ)*(ZZ-tau_x[kx+5])*(tau_x[kx+4]<=ZZ)*(ZZ<tau_x[kx+5])-
    (1/2/dx^3)*(tau_x[kx+3]-ZZ)*(ZZ-tau_x[kx+6])*(tau_x[kx+4]-ZZ)*(tau_x[kx+4]<=ZZ)*(ZZ<tau_x[kx+5])-
    (1/dx^3)*(ZZ-tau_x[kx+7])*(tau_x[kx+4]-ZZ)*(tau_x[kx+4]-ZZ)*(tau_x[kx+4]<=ZZ)*(ZZ<tau_x[kx+5])
  # i = kx + 4
  XX[kx+4,] = -(1/dx^3)*(tau_x[kx+4]-ZZ)*(tau_x[kx+4]-ZZ)*(tau_x[kx+4]-ZZ)*(tau_x[kx+4]<=ZZ)*(ZZ<=tau_x[kx+5])
  XX = t(XX)
  return(XX)
}


predict_wider <- function(p_dat){
  p_dat <- p_dat[,c('h','median')]
  str_h <- format(p_dat$h)
  p_dat$decimal <- as.numeric(substr(str_h,1,nchar(str_h)-1))
  first_decimal <- length(p_dat$decimal[p_dat$decimal==p_dat$decimal[1]])
  if(first_decimal!=10) {
    n <- 10-first_decimal
    top_rows <- data.frame(h=sapply(n:1,function(x) p_dat$h[1]-0.01*x),median=rep(0,n),decimal=rep(p_dat$decimal[1],n))
    p_dat <- rbind(top_rows,p_dat)
  }
  last_decimal <- length(p_dat$decimal[p_dat$decimal==p_dat$decimal[length(p_dat$decimal)]])
  if(last_decimal!=10){
    m <- 10-last_decimal
    bot_rows <- data.frame(h=sapply(1:m,function(x) p_dat$h[length(p_dat$h)]+0.01*x),median=rep(NA,m),decimal=rep(p_dat$decimal[length(p_dat$decimal)],m))
    p_dat <- rbind(p_dat,bot_rows)
  }
  p_mat <- lapply(unique(p_dat$decimal),function(d) p_dat$median[p_dat$decimal==d])
  p_mat <- do.call('rbind',p_mat)
  rownames(p_mat) <- unique(p_dat$decimal)
  colnames(p_mat) <- seq(0,0.09,by=0.01)
  return(p_mat)
}

#' @importFrom stats median
get_residuals_dat <- function(m){
  h_min <- min(m$data[[all.vars(m$formula)[2]]])
  rc_dat <- merge(m$rating_curve_mean[,c('h','median','lower','upper')],m$rating_curve[,c('h','median','lower','upper')],by.x='h',by.y='h')
  resid_dat <- merge(rc_dat[rc_dat$h>=h_min,],m$data,by.x='h',by.y=all.vars(m$formula)[2],all.x=TRUE)
  colnames(resid_dat) <- c('h','rcm_median','rcm_lower','rcm_upper','rc_median','rc_lower','rc_upper','Q')
  c_hat <- if(is.null(m$run_info$c_param)) median(m$c_posterior) else m$run_info$c_param
  resid_dat[,'log(h-c_hat)'] <- log(resid_dat$h-c_hat)
  resid_dat$r_median <- log(resid_dat$Q)-log(resid_dat$rc_median)
  resid_dat$m_lower <- log(resid_dat$rcm_lower)-log(resid_dat$rcm_median)
  resid_dat$m_upper <- log(resid_dat$rcm_upper)-log(resid_dat$rcm_median)
  resid_dat$r_lower <- log(resid_dat$rc_lower)-log(resid_dat$rc_median)
  resid_dat$r_upper <- log(resid_dat$rc_upper)-log(resid_dat$rc_median)
  return(resid_dat)
}


#' @importFrom stats var
chain_statistics <- function(chains){
  chains_length <- nrow(chains)
  split_idx <- round(0.5*chains_length)
  chains1 <- chains[seq(1,split_idx),]
  chains2 <- chains[seq(split_idx+1,chains_length),]
  chains_mean <- c(colMeans(chains1),colMeans(chains2))
  chains_var <- c(apply(chains1,2,var),apply(chains2,2,var))
  n <- round(0.5*chains_length)
  m <- ncol(chains)*2
  within_chain_var <- mean(chains_var)
  between_chain_var <- n*var(chains_mean)
  var_hat <- ((n-1)*within_chain_var + between_chain_var)/n
  return(list('W'=within_chain_var,'var_hat'=var_hat))
}

R_hat <- function(chains){
  staistics <- chain_statistics(chains)
  W <- staistics$W
  var_hat <- staistics$var_hat
  rhat <- sqrt(var_hat/W)
  return(rhat)
}


get_rhat_dat <- function(m,param,smoothness=20){
  thin <- m$run_info$thin
  burnin <- m$run_info$burnin
  draws_list <- lapply(param,function(x){
    draws <- gather_draws(m,x,transformed=TRUE)
    disjoint <- split(draws$value,draws$chain,drop=TRUE)
    disjoint <- do.call('cbind',disjoint)
  })
  names(draws_list) <- param
  rhat_dat <- lapply(param,function(x){
    draws  <- draws_list[[x]]
    real_iter <- seq(4*thin+burnin,((nrow(draws))*thin)+burnin,by=smoothness*thin)
    by_n <- seq(4,nrow(draws),by=smoothness)
    data.frame('iterations'=real_iter,'parameters'=rep(x,length(by_n)),'Rhat'=sapply(by_n, function(r) R_hat(draws[1:r,])))
  })
  rhat_dat <- do.call('rbind',rhat_dat)
  return(rhat_dat)
}


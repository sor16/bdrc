#' Prior parameter specification
#'
#'
#'@param country A string with the name of the country of which the prior parameters into the models should be specified for
#'@return
#'The priors are based from data from rivers of a given country.If you want to add your country to this function,
#'@references Birgir Hrafnkelsson, Helgi Sigurdarson and Sigurdur M. Gardarson (2015) \emph{Bayesian Generalized Rating Curves}
priors <- function(model,c_param) {
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
    if(model %in% c('bplm0','bplm')){
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
    if(model %in% c('bplm0','bgplm0')){
      RC$lambda_se <- 28.78
    }else{
      RC$lambda_eta_1 <- 28.78
      RC$lambda_seta <- 8.62
    }
    return(RC)
}

#'Linking unique water level measurements to actual
#'water level measurements
#'
#'Adist links unique water level measurements (\strong{h'}) to actual
#'water level measurements (h) such that \strong{h}=\strong{Ah'}.
#'from the measurements.
#'@param h numeric vector of stage measurements in meters
#'@return
#'\itemize{
#'\item A: Matrix \strong{A} linking unique water level measurements (\strong{h'}) to actual
#'water level measurements (h) such that \strong{h}=\strong{Ah'}
#'}
#'@references Birgir Hrafnkelsson, Helgi Sigurdarson and Sigurdur M. Gardarson (2015) \emph{Bayesian Generalized Rating Curves}
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

run_MCMC <- function(theta_m,RC,density_fun,unobserved_prediction_fun,nr_iter=20000,num_chains=4,burnin=2000,thin=5){
    theta_mat <- matrix(0,nrow=RC$theta_length,ncol=nr_iter)
    output_list <- initiate_output_list(RC$desired_output,nr_iter)
    density_eval_m <- density_fun(theta_m,RC)
    theta_old <- theta_m
    density_eval_old <- density_eval_m
    for(i in 1:nr_iter){
        theta_new <- theta_old+solve(t(RC$LH),stats::rnorm(RC$theta_length,0,1))
        density_eval_new <- density_fun(theta_new,RC)
        logR <- density_eval_new[['p']]-density_eval_old[['p']]
        if (logR>log(stats::runif(1))){
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

get_MCMC_summary <- function(X,h=NULL){
    summary_dat <- as.data.frame(t(apply(X,1,stats::quantile, probs = c(0.025,0.5, 0.975),na.rm=T)))
    names(summary_dat) <- c('lower','median','upper')
    if(!is.null(h)){
        summary_dat <- data.frame(h=h,summary_dat,row.names=NULL)
        summary_dat <- summary_dat[order(summary_dat$h),]
    }
    return(summary_dat)
}

get_param_names <- function(model,c_param){
    if(model=='bplm0'){
        hyper_param <- 'sigma_eps'
    }else if(model=='bplm'){
        hyper_param <- c('sigma_eta',paste('eta',1:6,sep='_'))
    }else if(model=='bgplm0'){
        hyper_param <- c('sigma_eps','sigma_beta','phi_beta')
    }else if(model=='bgplm'){
        hyper_param <- c('sigma_beta','phi_beta','sigma_eta',paste('eta',1:6,sep='_'))
    }
    if(is.null(c_param)){
        hyper_param <- c('c',hyper_param)
    }
    return(c('a','b',hyper_param))
}

# get_param_expression <- function(param){
#   expr_vec <- c('a'='a','b'='b','c'='c','sigma_eps'='sigma[epsilon]',
#                 'sigma_beta'='sigma[beta]','phi_beta'='phi[beta]',
#                 'sigma_eta'='sigma[eta]','eta_1'='eta[1]','eta_2'='eta[2]',
#                 'eta_3'='eta[3]','eta_4'='eta[4]','eta_5'='eta[5]',
#                 'eta_6'='eta[6]','log(a)'='log(a)','log(h_min-c)'='log(h[min]-c)',
#                 '2log(sigma_eps)'='log(sigma[epsilon]^2)',
#                 'log(sigma_beta)'='log(sigma[beta])',
#                 'log(phi_beta)'='log(phi[beta])',
#                 'log(sigma_eta)'='log(sigma[eta])',
#                 'z_1'='z[1]','z_2'='z[2]','z_3'='z[3]',
#                 'z_4'='z[4]','z_5'='z[5]','z_6'='z[6]')
#   param_expr <- expr_vec[param]
#   if(is.na(param_expr)){
#     stop('param not found')
#   }
#   return(param_expr)
# }


get_param_expression <- function(param){
  expr_vec <- c('a'='\\textit{a}','b'='\\textit{b}','c'='\\textit{c}','sigma_eps'='\\sigma_\\varepsilon',
                'sigma_beta'='\\sigma_\\beta','phi_beta'='\\phi_\\beta',
                'sigma_eta'='\\sigma_\\eta','eta_1'='\\eta_1','eta_2'='\\eta_2',
                'eta_3'='\\eta_3','eta_4'='\\eta_4','eta_5'='\\eta_5',
                'eta_6'='\\eta_6','log(a)'='\\log(\\textit{a})','log(h_min-c)'='\\log(\\textit{h_{min}-c})',
                '2log(sigma_eps)'='\\log(\\sigma_{\\epsilon}^2)',
                'log(sigma_beta)'='\\log(\\sigma_{\\beta})',
                'log(phi_beta)'='\\log(\\phi_{\\beta})',
                'log(sigma_eta)'='\\log(\\sigma_{\\eta})',
                'z_1'='\\textit{z_1}','z_2'='\\textit{z_2}','z_3'='\\textit{z_3}',
                'z_4'='\\textit{z_4}','z_5'='\\textit{z_5}','z_6'='\\textit{z_6}')
  param_expr <- paste0('$',expr_vec[param],'$')
  param_expr <- latex2exp::TeX(param_expr,output = 'character')
  #param_expr <- latex2exp::latex2exp(expr_vec[param],output = 'character')
  if(is.na(param_expr)){
    stop('param not found')
  }
  return(param_expr)
}

get_parameter_levels <- function(param_vec){
    order_vec <- c('a'=1,'log(a)'=2,'b'=3,'c'=4,'log(h_min-c)'=5,'sigma_eps'=6,
                   '2log(sigma_eps)'=7,'sigma_beta'=8,'log(sigma_beta)'=9,
                   'phi_beta'=10,'log(phi_beta)'=11,'sigma_eta'=12,'log(sigma_eta)'=13,
                   'eta_1'=14,'eta_2'=15,'z_1'=16,'eta_3'=17,'z_2'=18,
                   'eta_4'=19,'z_3'=20,'eta_5'=21,'z_4'=22,'eta_6'=23,'z_5'=24)
  return(param_vec[rank(order_vec[param_vec])])
}

get_transformed_param <- function(v,param_name,mod,...){
  args <- list(...)
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


#'Unobserved stages
#'
#'h_unobserved returns the stages that are needed to make an equally spaced grid of stages from data of stages.
#'
#'@param h_unique vector containing unique stages from river data.
#'@param min minimum stage of rating curve.
#'@param max maximum stage of rating curve.
#'@return h_unobserved returns a list of vectors, h_u and h_u_tild. h_u is a vector of unobserved stage values
#' needed to make an equally spaced grid of stages. h_u_tild is a vector which is calculated by h_u-min(h_unique) needed to input into B_splines.
#' The unobserved stages are lower or higher than that of the data, take the same value in h_u_tild as the minimum value and maximum value of the
#' data respectively. This is done to ensure constant variance below and above observed data.
#'@references Birgir Hrafnkelsson, Helgi Sigurdarson and Sigurdur M. Gardarson (2015) \emph{Bayesian Generalized Rating Curves}
h_unobserved <- function(RC,h_min=NA,h_max=NA){
  h_u=NULL
  h=100*c(RC$h) #work in cm
  max_h_diff=5
  #distance between subsequent elements in vector with additional dummy point 1000
  distvect=abs(h-c(h[2:length(h)],1000))
  #add datapoints to corresponding distances to see range of distance
  distwithdata=rbind(h,distvect,c(h[2:length(h)],1000))
  distfilter=distwithdata[,distvect>max_h_diff]
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

#'Bsplines in a generalized rating curve
#'
#'A function to test the B-splines in a rating curve. When calculating error variance of log discharge in a rating curve the data depends on stage. It is modeled as an exponential of a B-splines curve of order 4,
#'with 2 interior knots and 6 basis functions.
#'
#'@param ZZ A numeric matrix of dimension 1xn where n is number of osbervations. The input is calculated as follows:
#'(h-min(h)) divided by last element of the resulting vector, where h is stage observations.
#'@return The function returns a linear combination of scaled B-spline basis functions for every stage observation.
#'@references Birgir Hrafnkelsson, Helgi Sigurdarson and Sigurdur M. Gardarson (2015) \emph{Bayesian Generalized Rating Curves}
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





#' Spread MCMC chain draws to data.frame on a wide format
#'
#' Useful to convert MCMC chain draws of particular parameters or output from the model object to a wide format for further data wrangling
#'@param mod an object of class "bplm0","bplm","bgplm0" or "bgplm".
#'@param param character with the name of the parameter of interest. Also accepts "latent_parameters" and "hyperparameters".
#'@param transformed boolean value determining whether the parameter is to be represented on the transformed scale used for sampling in the MCMC chain or the original scale. Defaults to FALSE.
#'@return Data frame with columns
#'\code{chain}
#'\code{iter}
#'\code{param}
#'\code{value}
#'@references Birgir Hrafnkelsson, Helgi Sigurdarson, & Sigurdur M. Gardarsson. (2020). Generalization of the power-law rating curve using hydrodynamic theory and Bayesian hierarchical modeling.
#'@seealso \code{\link{bplm0}},\code{\link{bplm}},\code{\link{bgplm0}},\code{\link{bgplm}} for further information on parameters
#'@examples
#'data(sim_dat)
#'f <- Q~W
#'bplm0.fit <- bplm0(f,sim_dat)
#'summary(bplm9.fit)
#'plot(bplm0.fit,type='rating_curve')
#'@export
spread_draws <- function(mod,...,transformed=F){
    gathered_dat <- gather_draws(mod,...,transformed = F)
    if('h' %in% names(gathered_dat)){
        spread_dat <- expand.grid(iter=sort(unique(gathered_dat$iter)),
                                  chain=sort(unique(gathered_dat$chain)),
                                  h=sort(unique(gathered_dat$h)))
        spread_dat <- spread_dat[,c('chain','iter','h')]
    }else{
        spread_dat <- expand.grid(iter=sort(unique(gathered_dat$iter)),
                                  chain=sort(unique(gathered_dat$chain)))
        spread_dat <- spread_dat[,c('chain','iter')]
    }
    mod_res_list <- lapply(unique(gathered_dat$name),function(n){
        gathered_dat$value[gathered_dat$name==n]
    })
    mod_res <- as.data.frame(do.call('cbind',mod_res_list))
    names(mod_res) <- unique(gathered_dat$name)
    spread_dat <- cbind(spread_dat,mod_res)
}

#' Gather MCMC chain draws to data.frame on a long format
#'
#' Useful to convert MCMC chain draws of particular parameters or output from the model object to a long format for further data wrangling
#'@param mod an object of class "bplm0","bplm","bgplm0" or "bgplm".
#'@param ... character vectors of parameters or other output from the model output of interest. Also accepts "latent_parameters" and "hyperparameters".
#'@param transformed boolean value determining whether the parameter is to be represented on the transformed scale used for sampling in the MCMC chain or the original scale. Defaults to FALSE.
#'@return Data frame with columns
#'\code{chain}
#'\code{iter}
#'\code{param}
#'\code{value}
#'@references Birgir Hrafnkelsson, Helgi Sigurdarson, & Sigurdur M. Gardarsson. (2020). Generalization of the power-law rating curve using hydrodynamic theory and Bayesian hierarchical modeling.
#'@seealso \code{\link{bplm0}},\code{\link{bplm}},\code{\link{bgplm0}},\code{\link{bgplm}} for further information on parameters
#'@examples
#'data(sim_dat)
#'f <- Q~W
#'bplm0.fit <- bplm0(f,sim_dat)
#'summary(bplm9.fit)
#'plot(bplm0.fit,type='rating_curve')
#'@export
gather_draws <- function(mod,...,transformed=F){
    args <- c(...)
    ###TODO: refine class check
    if(!(class(mod) %in% c('bplm0','bplm','bgplm0','bgplm'))){
        stop('mod must be of class "bplm0","bplm","bgplm0" or "bgplm"')
    }
    mod_params <- get_param_names(class(mod),c_param=mod$run_info$c_param)
    args_rollout <- unlist(sapply(args,function(x) {
        if(x=='latent_parameters'){
            return(mod_params[1:2])
        }else if(x=='hyperparameters'){
            return(mod_params[3:length(mod_params)])
        }else{
            return(x)
        }
    }))
    args_rollout <- unique(args_rollout)
    if(all(args_rollout %in% gsub('_posterior','',names(mod)))){
        stage_dependent <- any(sapply(args_rollout,function(x) !is.null(dim(mod[[paste0(x,'_posterior')]]))))
        if(stage_dependent){
            baseline_dat <- expand.grid(iter=seq_len((mod$run_info$nr_iter-mod$run_info$burnin)/mod$run_info$thin+1),
                                        chain=1:mod$run_info$num_chains,
                                        h=mod$rating_curve$h)
            baseline_dat <- baseline_dat[,c('chain','iter','h')]
        }else{
            baseline_dat <- expand.grid(iter=seq_len((mod$run_info$nr_iter-mod$run_info$burnin)/mod$run_info$thin+1),
                                        chain=1:mod$run_info$num_chains)
            baseline_dat <- baseline_dat[,c('chain','iter')]
        }
        out_dat_list <- lapply(args_rollout,function(x){
            if(x %in% mod_params){
                dat <- gather_draws_param(mod,x,transformed=transformed,baseline_dat)
            }else{
                dat <- gather_draws_stage_dependent(mod,x,baseline_dat)
            }
            return(dat)
        })
        out_dat <- do.call('rbind',out_dat_list)
    }else{
        not_recognized <- which(!(args %in% c(paste0(names(mod),'_posterior'),'latent_parameters','hyperparameters')))
        stop(paste0('Does not recognize the following input arguments in the model object:\n',paste('\t-',args[not_recognized],collapse='\n')))
    }
    return(out_dat)
}
######## help functions
gather_draws_param <- function(mod,param,transformed,baseline_dat){
    post_param_name <- paste0(param,'_posterior')
    MCMC_output <- mod[[post_param_name]]
    param_name <- param
    if(transformed){
        if(param=='c'){
            h_min <- min(mod$data[[all.vars(mod$formula)[2]]])
            MCMC_output <- get_transformed_param(MCMC_output,param,mod,h_min=h_min)
        }else{
            MCMC_output <- get_transformed_param(MCMC_output,param,mod)
        }
        param_name <- unique(names(MCMC_output))
    }
    out_dat <- baseline_dat
    if('h' %in% names(baseline_dat)){
        param_dat <- expand.grid(name=param_name,h=mod$rating_curve$h,value=MCMC_output)
        out_dat <- cbind(out_dat,param_dat[,c('name','value')])
    }else{
        out_dat$name <- param_name
        out_dat$value <- MCMC_output
    }
    return(out_dat)
}

gather_draws_stage_dependent <- function(mod,name,baseline_dat){
    post_name <- paste0(name,'_posterior')
    MCMC_output <- mod[[post_name]]
    out_dat_list <- lapply(1:nrow(MCMC_output),
                           function(i){
                               out_dat <- data.frame(name=name,value=MCMC_output[i,,drop=T])
                               return(out_dat)
                           })
    out_dat <- cbind(baseline_dat,do.call('rbind',out_dat_list))
    return(out_dat)
}



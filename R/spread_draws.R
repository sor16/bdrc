#' Spread MCMC chain draws to data.frame
#'
#' Useful to convert MCMC chains draws of particular parameters in a model object to a tidyverse friendly data frame for further data wrangling and plotting. The MCMC samples returned lack the burnin samples and have been thinned according to run_info
#' element in the model object
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
# spread_draws <- function(mod,param,transformed=F){
#     num_chains <- mod$run_info$num_chains
#     if(paste0(param,'_posterior') %in% names(mod)){
#         list_name <- paste0(param,'_posterior')
#         MCMC_output <- mod[[list_name]]
#         if(is.null(dim(MCMC_output))){
#             if(transformed){
#                 if(param=='c'){
#                     w_min <- min(mod$data[[all.vars(mod$formula)[2]]])
#                     MCMC_output <- get_transformed_param(MCMC_output,param,mod,w_min=w_min)
#                 }else{
#                     MCMC_output <- get_transformed_param(MCMC_output,param,mod)
#                 }
#             }else{
#                 names(MCMC_output) <- rep(param,length(MCMC_output))
#             }
#             out_dat <- data.frame(chain=rep(1:num_chains,rep(length(MCMC_output)/num_chains,num_chains)),
#                                   iter=rep(1:(length(MCMC_output)/num_chains)))
#             out_dat[,unique(names(MCMC_output))] <- MCMC_output
#         }else{
#             c_hat <- if(is.null(mod$run_info_c_param)) median(mod$c_posterior) else mod$run_info$c_param
#             out_dat_list <- lapply(1:nrow(MCMC_output),
#                                    function(i){
#                                         out_dat <- data.frame(chain=rep(1:num_chains,rep(ncol(MCMC_output)/num_chains,num_chains)),
#                                                               iter=rep(1:(ncol(MCMC_output)/num_chains),num_chains),w=mod$rating_curve$w[i])
#                                         if(transformed){
#                                             out_dat[,'log(w-c_hat)'] <- suppressWarnings(log(out_dat$w-c_hat))
#                                         }
#                                         out_dat[,param] <- MCMC_output[i,,drop=T]
#                                         return(out_dat)
#                                    })
#             out_dat <- do.call('rbind',out_dat_list)
#             if(transformed){
#                 w_name <- all.vars(mod$formula)[2]
#                 out_dat <- out_dat[out_dat$w>=min(mod$data[,w_name]) & out_dat$w<=max(mod$data[,w_name]),]
#             }
#         }
#     }else if(param %in% c('latent_parameters','hyperparameters')){
#         params <- get_param_names(class(mod),c_param=mod$run_info$c_param)
#         param_idx <- if(param=='latent_parameters') 1:2 else 3:length(params)
#         out_dat_list <- lapply(params[param_idx],function(param_name){
#             list_name <- paste0(param_name,'_posterior')
#             MCMC_output <- mod[[list_name]]
#             if(transformed){
#                 if(param_name=='c'){
#                     w_min <- min(mod$data[[all.vars(mod$formula)[2]]])
#                     MCMC_output <- get_transformed_param(MCMC_output,param_name,mod,w_min=w_min)
#                 }else{
#                     MCMC_output <- get_transformed_param(MCMC_output,param_name,mod)
#                 }
#             }else{
#                 names(MCMC_output) <- rep(param_name,length(MCMC_output))
#             }
#             out_dat <- data.frame(chain=rep(1:num_chains,rep(length(MCMC_output)/num_chains,num_chains)),
#                                   iter=rep(1:(length(MCMC_output)/num_chains)),
#                                   param=names(MCMC_output),
#                                   value=MCMC_output)
#         })
#         out_dat <- do.call('rbind',out_dat_list)
#     }else{
#         stop('param not found in model object')
#     }
#     rownames(out_dat) <- NULL
#     return(out_dat)
# }

spread_draws <- function(mod,...,transformed=F){
    gathered_dat <- gather_draws(mod,...,transformed = F)
    if('w' %in% names(gathered_dat)){
        spread_dat <- expand.grid(w=sort(unique(gathered_dat$w)),
                                  iter=sort(unique(gathered_dat$iter)),
                                  chain=sort(unique(gathered_dat$chain)))
        spread_dat <- spread_dat[,c('chain','iter','w')]
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

# gather_draws <- function(mod,...,transformed=F){
#     args <- c(...)
#     ###TODO: refine class check
#     if(!(class(mod) %in% c('bplm0','bplm','bgplm0','bgplm'))){
#         stop('mod must be of class "bplm0","bplm","bgplm0" or "bgplm"')
#     }
#     mod_params <- get_param_names(class(mod),c_param=mod$run_info$c_param)
#     if(all(args %in% c(gsub('_posterior','',names(mod)),'latent_parameters','hyperparameters'))){
#         stage_dependent <- any(
#                             sapply(args,function(x){
#                                 if(x %in% c('latent_parameters','hyperparameters')){
#                                     return(FALSE)
#                                 }else{
#                                     return(!is.null(dim(mod[[paste0(x,'_posterior')]])))
#                                 }
#                             }))
#         if(stage_dependent){
#             stage_dependent_dat <- expand.grid(w=mod$rating_curve$w,
#                                                iter=seq_len((mod$run_info$nr_iter-mod$run_info$burnin)/mod$run_info$thin+1),
#                                                chain=1:mod$run_info$num_chains)
#             stage_dependent_dat <- stage_dependent_dat[,c('chain','iter','w')]
#         }
#         out_dat_list <- lapply(args,function(x){
#             if(x %in% c('hyperparameters','latent_parameters')){
#                 param_idx <- if(x=='latent_parameters') 1:2 else 3:length(mod_params)
#                 dat_list <- lapply(mod_params[param_idx],function(param_name){
#                     gather_draws_param(mod,param_name,transformed=transformed)
#                 })
#                 dat <- do.call('rbind',dat_list)
#             }else if(x %in% mod_params){
#                 dat <- gather_draws_param(mod,x,transformed=transformed)
#             }else{
#                 dat <- gather_draws_stage_dependent(mod,x)
#             }
#             if(stage_dependent & !('w' %in% names(dat))){
#                 dat <- merge(stage_dependent_dat,by=c('chain','iter'),dat)
#             }
#             return(dat)
#         })
#         out_dat <- do.call('rbind',out_dat_list)
#     }else{
#         not_recognized <- which(!(args %in% c(paste0(names(mod),'_posterior'),'latent_parameters','hyperparameters')))
#         stop(paste0('Does not recognize the following input arguments in the model object:\n',paste('\t-',args[not_recognized],collapse='\n')))
#     }
#     out_dat <- out_dat[!duplicated(out_dat), ]
#     return(out_dat)
# }

gather_draws_param <- function(mod,param,transformed,baseline_dat){
    post_param_name <- paste0(param,'_posterior')
    MCMC_output <- mod[[post_param_name]]
    param_name <- param
    if(transformed){
        if(param=='c'){
            w_min <- min(mod$data[[all.vars(mod$formula)[2]]])
            MCMC_output <- get_transformed_param(MCMC_output,param,mod,w_min=w_min)
        }else{
            MCMC_output <- get_transformed_param(MCMC_output,param,mod)
        }
        param_name <- unique(names(MCMC_output))
    }
    out_dat <- baseline_dat
    if('w' %in% names(baseline_dat)){
        param_dat <- expand.grid(name=param_name,w=mod$rating_curve$w,value=MCMC_output)
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
                                        w=mod$rating_curve$w)
            baseline_dat <- baseline_dat[,c('chain','iter','w')]
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


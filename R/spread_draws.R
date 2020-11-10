#' Spread MCMC chain draws to data.frame
#'
#' Useful to convert MCMC chains draws of particular parameters in a model object to a tidyverse friendly data frame for further data wrangling and plotting
#'@param mod an object of class "bplm0","bplm","bgplm0" or "bgplm".
#'@param param character with the name of the parameter of interest. Also accepts "latent_parameters" and "hyperparameters".
#'@param transformed boolean value determining whether the parameter is to be represented on the transformed scale used in sampling or the original scale. Defaults to FALSE.
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
spread_draws <- function(mod,param,transformed=F){
    num_chains <- mod$run_info$num_chains
    if(paste0(param,'_posterior') %in% names(mod)){
        list_name <- paste0(param,'_posterior')
        MCMC_output <- mod[[list_name]]
        if(is.null(dim(MCMC_output))){
            if(transformed){
                if(param=='c'){
                    w_min <- min(mod$data[[all.vars(mod$formula)[2]]])
                    MCMC_output <- get_transformed_param(MCMC_output,param,mod,w_min=w_min)
                }else{
                    MCMC_output <- get_transformed_param(MCMC_output,param,mod)
                }
            }else{
                names(MCMC_output) <- rep(param,length(MCMC_output))
            }
            out_dat <- data.frame(chain=rep(1:num_chains,rep(length(MCMC_output)/num_chains,num_chains)),
                                  iter=rep(1:(length(MCMC_output)/num_chains)))
            out_dat[,unique(names(MCMC_output))] <- MCMC_output
        }else{
            c_hat <- if(is.null(mod$run_info_c_param)) median(mod$c_posterior) else mod$run_info$c_param
            out_dat_list <- lapply(1:nrow(MCMC_output),
                                   function(i){
                                        out_dat <- data.frame(chain=rep(1:num_chains,rep(ncol(MCMC_output)/num_chains,num_chains)),
                                                              iter=rep(1:(ncol(MCMC_output)/num_chains),num_chains),w=mod$rating_curve$w[i])
                                        if(transformed){
                                            out_dat[,'log(w-c_hat)'] <- suppressWarnings(log(out_dat$w-c_hat))
                                        }
                                        out_dat[,param] <- MCMC_output[i,,drop=T]
                                        return(out_dat)
                                   })
            out_dat <- do.call('rbind',out_dat_list)
            if(transformed){
                w_name <- all.vars(mod$formula)[2]
                out_dat <- out_dat[out_dat$w>=min(mod$data[,w_name]) & out_dat$w<=max(mod$data[,w_name]),]
            }
        }
    }else if(param %in% c('latent_parameters','hyperparameters')){
        params <- get_param_names(class(mod),c_param=mod$run_info$c_param)
        param_idx <- if(param=='latent_parameters') 1:2 else 3:length(params)
        out_dat_list <- lapply(params[param_idx],function(param_name){
            list_name <- paste0(param_name,'_posterior')
            MCMC_output <- mod[[list_name]]
            if(transformed){
                if(param_name=='c'){
                    w_min <- min(mod$data[[all.vars(mod$formula)[2]]])
                    MCMC_output <- get_transformed_param(MCMC_output,param_name,mod,w_min=w_min)
                }else{
                    MCMC_output <- get_transformed_param(MCMC_output,param_name,mod)
                }
            }else{
                names(MCMC_output) <- rep(param_name,length(MCMC_output))
            }
            out_dat <- data.frame(chain=rep(1:num_chains,rep(length(MCMC_output)/num_chains,num_chains)),
                                  iter=rep(1:(length(MCMC_output)/num_chains)),
                                  param=names(MCMC_output),
                                  value=MCMC_output)
        })
        out_dat <- do.call('rbind',out_dat_list)
    }else{
        stop('param not found in model object')
    }
    rownames(out_dat) <- NULL
    return(out_dat)
}


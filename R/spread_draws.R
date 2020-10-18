spread_draws <- function(mod,param){
    num_chains <- mod$run_info$num_chains
    if(paste0(param,'_posterior') %in% names(mod)){
        list_name <- paste0(param,'_posterior')
        MCMC_output <- mod[[list_name]]
        if(is.null(dim(MCMC_output))){
            out_dat <- data.frame(chain=rep(1:num_chains,rep(length(MCMC_output)/num_chains,num_chains)),
                                  iter=rep(1:(length(MCMC_output)/num_chains)))
            out_dat[,param] <- MCMC_output
        }else{
            out_dat_list <- lapply(1:nrow(MCMC_output),
                                   function(i){
                                        out_dat <- data.frame(chain=rep(1:num_chains,rep(ncol(MCMC_output)/num_chains,num_chains)),
                                                              iter=rep(1:(ncol(MCMC_output)/num_chains),num_chains),W=mod$rating_curve$w[i])
                                        out_dat[,param] <- MCMC_output[i,,drop=T]
                                        return(out_dat)
                                   })
            out_dat <- do.call('rbind',out_dat_list)
        }
    }else if(param %in% c('latent_parameters','hyperparameters')){
        params <- get_param_names(class(mod),!is.null(mod$c_posterior))
        param_idx <- if(param=='latent_parameters') 1:2 else 3:length(params)
        out_dat_list <- lapply(params[param_idx],function(param_name){
            list_name <- paste0(param_name,'_posterior')
            MCMC_output <- mod[[list_name]]
            out_dat <- data.frame(chain=rep(1:num_chains,rep(length(MCMC_output)/num_chains,num_chains)),
                                  iter=rep(1:(length(MCMC_output)/num_chains)),
                                  param=param_name,
                                  value=MCMC_output)
        })
        out_dat <- do.call('rbind',out_dat_list)
    }else{
        stop('param not found in model object')
    }
    return(out_dat)
}


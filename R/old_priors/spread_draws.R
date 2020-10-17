spread_draws <- function(mod,param){
    list_name <- paste0(param,'_posterior')
    if(list_name %in% names(mod)){
        MCMC_output <- mod[[list_name]]
        num_chains <- mod$run_info$num_chains
        if(is.null(dim(MCMC_output))){
            out_dat <- data.frame(chain=rep(1:num_chains,rep(length(MCMC_output)/num_chains,num_chains)),
                                  iter=rep(1:(length(MCMC_output)/num_chains)))
            out_dat[,param] <- MCMC_output
        }else{
            out_dat_list <- lapply(1:nrow(MCMC_output),
                                   function(i){
                                        out_dat <- data.frame(chain=rep(1:num_chains,rep(ncol(MCMC_output)/num_chains,num_chains)),
                                                              iter=rep(1:(ncol(MCMC_output)/num_chains),num_chains),W=mod$rating_curve$w[i])
                                        out_dat[,param] <- MCMC_output[i,drop=T]
                                        return(out_dat)
                                   })
            out_dat <- do.call('rbind',out_dat_list)
        }
    }else{
        stop('param not found in model object')
    }
    return(out_dat)
}

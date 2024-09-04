#' Spread MCMC chain draws to data.frame on a wide format
#'
#' Useful to convert MCMC chain draws of particular parameters or output from the model object to a wide format for further data wrangling
#'
#' @param mod An object of class "plm0", "plm", "gplm0" or "gplm".
#' @param ... Any number of character vectors containing valid names of parameters in the model or "rating_curve" and "rating_curve_mean". Also accepts "latent_parameters" and "hyperparameters".
#' @param transformed A boolean value determining whether the output is to be represented on the transformed scale used for sampling in the MCMC chain or the original scale. Defaults to FALSE.
#' @return A data frame with columns:
#' \describe{
#'   \item{\code{chain}}{The chain number.}
#'   \item{\code{iter}}{The iteration number.}
#'   \item{\code{param}}{The parameter name.}
#'   \item{\code{value}}{The parameter value.}
#' }
#' @references Hrafnkelsson, B., Sigurdarson, H., Rögnvaldsson, S., Jansson, A. Ö., Vias, R. D., and Gardarsson, S. M. (2022). Generalization of the power-law rating curve using hydrodynamic theory and Bayesian hierarchical modeling, Environmetrics, 33(2):e2711. doi: https://doi.org/10.1002/env.2711
#' @seealso \code{\link{plm0}}, \code{\link{plm}}, \code{\link{gplm0}}, \code{\link{gplm}} for further information on parameters
#' @examples
#' \donttest{
#'  data(krokfors)
#'  set.seed(1)
#'  plm0.fit <- plm0(formula=Q~W,data=krokfors,num_cores=2)
#'  hyp_samples <- spread_draws(plm0.fit,'hyperparameters')
#'  head(hyp_samples)
#'  rating_curve_samples <- spread_draws(plm0.fit,'rating_curve','rating_curve_mean')
#'  head(rating_curve_samples)
#' }
#' @export
spread_draws <- function(mod, ..., transformed = FALSE){
    gathered_dat <- gather_draws(mod, ..., transformed = transformed)
    if('h' %in% names(gathered_dat)){
        spread_dat <- expand.grid(iter = sort(unique(gathered_dat$iter)),
                                  chain = sort(unique(gathered_dat$chain)),
                                  h = sort(unique(gathered_dat$h)), stringsAsFactors = FALSE)
        spread_dat <- spread_dat[, c('chain', 'iter', 'h')]
    }else{
        spread_dat <- expand.grid(iter = sort(unique(gathered_dat$iter)),
                                  chain = sort(unique(gathered_dat$chain)), stringsAsFactors = FALSE)
        spread_dat <- spread_dat[, c('chain', 'iter')]
    }
    mod_res_list <- lapply(unique(gathered_dat$name), function(n){
        gathered_dat$value[gathered_dat$name == n]
    })
    mod_res <- as.data.frame(do.call('cbind', mod_res_list))
    names(mod_res) <- unique(gathered_dat$name)
    spread_dat <- cbind(spread_dat, mod_res)
    return(spread_dat)
}

#' Gather MCMC chain draws to data.frame on a long format
#'
#' Useful to convert MCMC chain draws of particular parameters or output from the model object to a long format for further data wrangling
#'
#' @param mod An object of class "plm0", "plm", "gplm0" or "gplm".
#' @param ... Any number of character vectors containing valid names of parameters in the model or "rating_curve" and "rating_curve_mean". Also accepts "latent_parameters" and "hyperparameters".
#' @param transformed A boolean value determining whether the parameter is to be represented on the transformed scale used for sampling in the MCMC chain or the original scale. Defaults to FALSE.
#' @return A data frame with columns:
#' \describe{
#'   \item{\code{chain}}{The chain number.}
#'   \item{\code{iter}}{The iteration number.}
#'   \item{\code{param}}{The parameter name.}
#'   \item{\code{value}}{The parameter value.}
#' }
#' @references Hrafnkelsson, B., Sigurdarson, H., Rögnvaldsson, S., Jansson, A. Ö., Vias, R. D., and Gardarsson, S. M. (2022). Generalization of the power-law rating curve using hydrodynamic theory and Bayesian hierarchical modeling, Environmetrics, 33(2):e2711. doi: https://doi.org/10.1002/env.2711
#' @seealso \code{\link{plm0}}, \code{\link{plm}}, \code{\link{gplm0}}, \code{\link{gplm}} for further information on parameters
#' @examples
#' \donttest{
#'  data(krokfors)
#'  set.seed(1)
#'  plm0.fit <- plm0(formula=Q~W,data=krokfors,num_cores=2)
#'  hyp_samples <- gather_draws(plm0.fit,'hyperparameters')
#'  head(hyp_samples)
#'  rating_curve_samples <- gather_draws(plm0.fit,'rating_curve','rating_curve_mean')
#'  head(rating_curve_samples)
#' }
#' @export
gather_draws <- function(mod, ..., transformed = F){
    args <- c(...)
    if(!(class(mod) %in% c('plm0', 'plm', 'gplm0', 'gplm'))){
        stop('mod must be of class "plm0", "plm", "gplm0", or "gplm"')
    }
    mod_params <- get_param_names(class(mod), c_param = mod$run_info$c_param)
    args_rollout <- get_args_rollout(args, mod_params)
    f_not_generalized <- any(grepl('^f$', args_rollout)) & is.null(mod$f_posterior)
    if(f_not_generalized){
        args_rollout[grepl('^f$', args_rollout)] <- 'b'
    }
    if(all(args_rollout %in% gsub('_posterior', '', names(mod)))){
        stage_dependent <- any(sapply(args_rollout, function(x) !is.null(dim(mod[[paste0(x, '_posterior')]]))))
        if(stage_dependent){
            baseline_dat <- expand.grid(iter = seq_len((mod$run_info$nr_iter - mod$run_info$burnin) / mod$run_info$thin + 1),
                                        chain = 1:mod$run_info$num_chains,
                                        h = mod$rating_curve$h, stringsAsFactors = F)
            baseline_dat <- baseline_dat[, c('chain', 'iter', 'h')]
        }else{
            baseline_dat <- expand.grid(iter = seq_len((mod$run_info$nr_iter - mod$run_info$burnin) / mod$run_info$thin + 1),
                                        chain = 1:mod$run_info$num_chains, stringsAsFactors = F)
            baseline_dat <- baseline_dat[, c('chain', 'iter')]
        }
        out_dat_list <- lapply(args_rollout, function(x){
            if(x %in% mod_params){
                dat <- gather_draws_param(mod, x, transformed = transformed, baseline_dat)
            }else{
                dat <- gather_draws_stage_dependent(mod, x, baseline_dat)
            }
            return(dat)
        })
        out_dat <- do.call('rbind', out_dat_list)
    }else{
        not_recognized <- which(!(args %in% c(paste0(names(mod), '_posterior'), 'latent_parameters', 'hyperparameters')))
        stop(paste0('Does not recognize the following input arguments in the model object:\n', paste('\t-', args[not_recognized], collapse = '\n')))
    }
    if(f_not_generalized){
        out_dat$name[out_dat$name == 'b'] <- 'f'
    }
    return(out_dat)
}
######## help functions
gather_draws_param <- function(mod, param, transformed, baseline_dat){
    post_param_name <- paste0(param, '_posterior')
    MCMC_output <- mod[[post_param_name]]
    param_name <- param
    if(transformed){
        if(param == 'c'){
            h_min <- min(mod$data[[all.vars(mod$formula)[2]]])
            MCMC_output <- get_transformed_param(MCMC_output, param, mod, h_min = h_min)
        }else{
            MCMC_output <- get_transformed_param(MCMC_output, param, mod)
        }
        param_name <- unique(names(MCMC_output))
    }
    out_dat <- baseline_dat
    if('h' %in% names(baseline_dat)){
        param_dat <- expand.grid(name = param_name, value = MCMC_output, h = mod$rating_curve$h, stringsAsFactors = F)
        out_dat <- cbind(out_dat,param_dat[, c('name', 'value')])
    }else{
        out_dat$name <- param_name
        out_dat$value <- MCMC_output
    }
    return(out_dat)
}

gather_draws_stage_dependent <- function(mod, name, baseline_dat){
    post_name <- paste0(name, '_posterior')
    MCMC_output <- mod[[post_name]]
    out_dat_list <- lapply(1:nrow(MCMC_output),
                           function(i){
                               out_dat <- data.frame(name = name, value = MCMC_output[i, , drop = TRUE])
                               return(out_dat)
                           })
    out_dat <- cbind(baseline_dat, do.call('rbind', out_dat_list))
    return(out_dat)
}




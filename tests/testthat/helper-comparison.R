tol <- 1e-8

test_stage_indep_param <- function(mod, param){
    expect_true(is.double(mod[[paste0(param, '_posterior')]]))
    expect_equal(length(mod[[paste0(param, '_posterior')]]),
                 mod$run_info$num_chains * ((mod$run_info$nr_iter - mod$run_info$burnin) / mod$run_info$thin + 1))
    expect_equal(unname(unlist(mod$param_summary[param, c('lower', 'median', 'upper')])),
                 as.double(get_MCMC_summary(matrix(mod[[paste0(param, '_posterior')]], nrow = 1))),
                 tolerance = tol)
}

test_stage_dep_component <- function(mod, component){
    expect_true(is.double(mod[[paste0(component, '_posterior')]]))
    expect_equal(dim(mod[[paste0(component, '_posterior')]]),
                 c(length(unique(mod$rating_curve$h)), mod$run_info$num_chains * ((mod$run_info$nr_iter - mod$run_info$burnin) / mod$run_info$thin + 1)))
    summary_component <- if(component %in% c('rating_curve', 'rating_curve_mean')) component else paste0(component, '_summary')
    expect_equal(c(as.matrix(unname(mod[[summary_component]][, c('lower', 'median', 'upper')]))),
                 c(as.matrix(unname(get_MCMC_summary(mod[[paste0(component, '_posterior')]])))),
                 tolerance = tol)
}


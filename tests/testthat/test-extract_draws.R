context("Extract draws functions")

check_spread_draws_MCMC_mat_concordance <- function(mod, dat){
    components <- unique(names(dat)[!(names(dat) %in% c('chain', 'iter', 'h'))])
    for(x in components){
        MCMC_output <- mod[[paste0(x, '_posterior')]]
        if(is.null(dim(MCMC_output))){
            if('h' %in% names(dat)){
                expect_equal(dat[[x]][dat$h == dat$h[1]], MCMC_output)
            }else{
                expect_equal(dat[[x]], MCMC_output)
            }
        }else{
            expect_equal(dat[[x]], c(t(MCMC_output)))
        }
    }
}

check_gather_draws_MCMC_mat_concordance <- function(mod, dat){
    components <- unique(dat$name)
    for(x in components){
        MCMC_output <- mod[[paste0(x, '_posterior')]]
        curr_dat <- dat[dat$name==x, ]
        if(is.null(dim(MCMC_output))){
            if('h' %in% names(curr_dat)){
                expect_equal(curr_dat$value[curr_dat$h == curr_dat$h[1]], MCMC_output)
            }else{
                expect_equal(curr_dat$value, MCMC_output)
            }
        }else{
            expect_equal(curr_dat$value, c(t(MCMC_output)))
        }
    }
}

test_that("check different inputs to spread draws and gather draws", {
    # skip_on_cran()
    expect_equal(spread_draws(plm0.fit, 'c', 'sigma_eps'),
                 spread_draws(plm0.fit, c('c', 'sigma_eps')))
    expect_equal(spread_draws(plm.fit, 'c', paste0('eta_', 1:3)),
                 spread_draws(plm.fit, c('c', 'c', paste0('eta_', 1:3))))
    expect_equal(spread_draws(gplm0.fit, 'c', 'c'),
                 spread_draws(gplm0.fit, 'c'))
    expect_equal(spread_draws(gplm.fit, c('c', 'sigma_beta'), c('eta_1', 'sigma_eta')),
                 spread_draws(gplm.fit,'c','sigma_beta', 'eta_1', 'sigma_eta'))
    expect_error(spread_draws(gplm.fit_known_c, 'c'))

    expect_equal(gather_draws(plm0.fit, 'c', 'sigma_eps'),
                 gather_draws(plm0.fit, c('c', 'sigma_eps')))
    expect_equal(gather_draws(plm.fit, 'c', paste0('eta_', 1:3)),
                 gather_draws(plm.fit, c('c', 'c', paste0('eta_', 1:3))))
    expect_equal(gather_draws(gplm0.fit, 'c', 'c'),
                 gather_draws(gplm0.fit, 'c'))
    expect_equal(gather_draws(gplm.fit, c('c', 'sigma_beta'), c('eta_1', 'sigma_eta')),
                 gather_draws(gplm.fit, 'c', 'sigma_beta', 'eta_1', 'sigma_eta'))
    expect_error(gather_draws(gplm.fit_known_c, 'c'))
})

test_that("spread draws and gather draws are concordant with MCMC matrices", {
    # skip_on_cran()
    #spread_draws
    check_spread_draws_MCMC_mat_concordance(plm0.fit,spread_draws(plm0.fit, 'c', 'sigma_eps'))
    check_spread_draws_MCMC_mat_concordance(plm0.fit,spread_draws(plm0.fit, 'rating_curve', 'rating_curve_mean', 'c'))

    check_spread_draws_MCMC_mat_concordance(plm.fit,spread_draws(plm.fit, 'c', paste0('eta_', 1:6)))
    check_spread_draws_MCMC_mat_concordance(plm.fit,spread_draws(plm.fit, 'rating_curve', 'rating_curve_mean', 'sigma_eps', 'sigma_eta'))

    check_spread_draws_MCMC_mat_concordance(gplm0.fit,spread_draws(gplm0.fit, 'c', 'phi_beta', 'sigma_beta'))
    check_spread_draws_MCMC_mat_concordance(gplm0.fit,spread_draws(gplm0.fit, 'rating_curve', 'rating_curve_mean', 'beta', 'f', 'sigma_eps'))

    check_spread_draws_MCMC_mat_concordance(gplm.fit,spread_draws(gplm.fit, 'c', paste0('eta_', 1:6)))
    check_spread_draws_MCMC_mat_concordance(gplm.fit,spread_draws(gplm.fit, 'rating_curve', 'rating_curve_mean', 'sigma_eps', 'beta', 'f', 'sigma_eta'))
    #gather_draws
    check_gather_draws_MCMC_mat_concordance(plm0.fit,spread_draws(plm0.fit, 'c', 'sigma_eps'))
    check_gather_draws_MCMC_mat_concordance(plm0.fit,spread_draws(plm0.fit, 'rating_curve', 'rating_curve_mean', 'c'))

    check_gather_draws_MCMC_mat_concordance(plm.fit,spread_draws(plm.fit, 'c', paste0('eta_', 1:6)))
    check_gather_draws_MCMC_mat_concordance(plm.fit,spread_draws(plm.fit, 'rating_curve', 'rating_curve_mean', 'sigma_eps', 'sigma_eta'))

    check_gather_draws_MCMC_mat_concordance(gplm0.fit,spread_draws(gplm0.fit, 'c', 'phi_beta', 'sigma_beta'))
    check_gather_draws_MCMC_mat_concordance(gplm0.fit,spread_draws(gplm0.fit, 'rating_curve', 'rating_curve_mean', 'beta', 'f', 'sigma_eps'))

    check_gather_draws_MCMC_mat_concordance(gplm.fit,spread_draws(gplm.fit, 'c', paste0('eta_', 1:6)))
    check_gather_draws_MCMC_mat_concordance(gplm.fit,spread_draws(gplm.fit, 'rating_curve', 'rating_curve_mean', 'sigma_eps', 'beta', 'f', 'sigma_eta'))
})

test_that("spread_draws handles various input combinations", {
    # Test with single parameter
    single_param <- spread_draws(plm0.fit, "c")
    expect_true("c" %in% names(single_param))
    expect_equal(ncol(single_param), 3)

    # Test with multiple parameters
    multi_param <- spread_draws(plm0.fit, "c", "sigma_eps")
    expect_true(all(c("c", "sigma_eps") %in% names(multi_param)))
    expect_equal(ncol(multi_param), 4)

    # Test with transformed parameter
    transformed_param <- spread_draws(plm0.fit, "c", transformed = TRUE)
    expect_true("log(h_min-c)" %in% names(transformed_param))

    # Test with rating curve
    rc_param <- spread_draws(plm0.fit, "rating_curve")
    expect_true("rating_curve" %in% names(rc_param))
    expect_true("h" %in% names(rc_param))
})

test_that("gather_draws handles various input combinations", {
    # Test with single parameter
    single_param <- gather_draws(plm0.fit, "c")
    expect_equal(unique(single_param$name), "c")
    expect_equal(ncol(single_param), 4)

    # Test with multiple parameters
    multi_param <- gather_draws(plm0.fit, "c", "sigma_eps")
    expect_equal(sort(unique(multi_param$name)), c("c", "sigma_eps"))

    # Test with transformed parameter
    transformed_param <- gather_draws(plm0.fit, "c", transformed = TRUE)
    expect_true(any(grepl("log\\(h_min-c\\)", transformed_param$name)))

    # Test with rating curve
    rc_param <- gather_draws(plm0.fit, "rating_curve")
    expect_true("rating_curve" %in% unique(rc_param$name))
    expect_true("h" %in% names(rc_param))
})

test_that("spread_draws and gather_draws handle errors correctly", {
    # Test with invalid parameter name
    expect_error(spread_draws(plm0.fit, "invalid_param"), "Does not recognize the following input arguments")
    expect_error(gather_draws(plm0.fit, "invalid_param"), "Does not recognize the following input arguments")

    # Test with invalid model object
    invalid_model <- list(class = "invalid")
    expect_error(spread_draws(invalid_model, "c"), "mod must be of class")
    expect_error(gather_draws(invalid_model, "c"), "mod must be of class")
})

test_that("spread_draws and gather_draws produce consistent results", {
    params <- c("c", "sigma_eps", "rating_curve")
    spread_result <- spread_draws(plm0.fit, params)
    gather_result <- gather_draws(plm0.fit, params)

    # Check that the number of rows is consistent
    expect_equal(nrow(spread_result), nrow(gather_result) / length(params))

    # Check that the values are consistent
    for (param in params) {
        spread_values <- spread_result[[param]]
        gather_values <- gather_result$value[gather_result$name == param]
        expect_equal(spread_values, gather_values)
    }
})




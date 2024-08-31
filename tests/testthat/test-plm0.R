context('plm0')
tol <- 1e-8

test_that("plm0 can handle different inputs", {
    expect_error(plm0(Q ~ W, c(1, 2, 3)))
    expect_error(plm0('Q ~ W', krokfors))
    expect_error(plm0(V ~ W, krokfors))
    expect_error(plm0(Q ~ W + X, krokfors))
    expect_error(plm0(Q ~ W, krokfors, c_param = min(krokfors$W) + 0.5)) # c_param higher than lowest stage measurements
    expect_error(plm0(Q ~ W, krokfors, c_param = 1L)) # c_param not double
    expect_error(plm0(Q ~ W, krokfors, h_max = max(krokfors$W) - 0.5)) #h_max lower than highest stage measurement
    expect_error(plm0(Q ~ W, krokfors[1,]), "At least two paired observations of stage and discharge")
    expect_error(plm0(Q ~ W, -1 * krokfors), "All discharge measurements must but strictly greater than zero")
    skip_on_cran()
    krokfors_new_names <- krokfors
    names(krokfors_new_names) <- c('t1', 't2')
    set.seed(1)
    plm0.fit_new_names <- plm0(t2 ~ t1, krokfors_new_names, num_cores = 2)
    expect_equal(plm0.fit_new_names$rating_curve, plm0.fit$rating_curve, tolerance = tol)
})

test_that("the plm0 object with unknown c is in tact", {
    expect_is(plm0.fit, "plm0")
    # latent parameters
    test_stage_indep_param(plm0.fit, 'a')
    test_stage_indep_param(plm0.fit, 'b')
    # hyperparameters
    test_stage_indep_param(plm0.fit, 'c')
    test_stage_indep_param(plm0.fit, 'sigma_eps')
    # log-likelihood
    expect_true(is.double(plm0.fit$posterior_log_likelihood))
    expect_equal(length(plm0.fit$posterior_log_likelihood),
                 plm0.fit$run_info$num_chains * ((plm0.fit$run_info$nr_iter - plm0.fit$run_info$burnin) / plm0.fit$run_info$thin + 1))
    expect_equal(unname(unlist(plm0.fit$posterior_log_likelihood_summary[1, ])),
                 as.double(get_MCMC_summary(matrix(plm0.fit$posterior_log_likelihood, nrow = 1))),
                 tolerance = tol)
    # rating curve
    test_stage_dep_component(plm0.fit, 'rating_curve')
    test_stage_dep_component(plm0.fit, 'rating_curve_mean')
    # Other information
    expect_equal(plm0.fit$formula, Q ~ W)
    expect_equal(plm0.fit$data, krokfors[order(krokfors$W), c('Q', 'W')])
})

test_that("the plm0 object with known c with a maximum stage value is in tact", {
    skip_on_cran()
    set.seed(1)
    plm0.fit_known_c <- plm0(Q ~ W, krokfors, c_param = known_c, h_max = h_extrap, num_cores = 2)
    expect_is(plm0.fit_known_c, "plm0")
    # latent parameters
    test_stage_indep_param(plm0.fit_known_c, 'a')
    test_stage_indep_param(plm0.fit_known_c, 'b')
    # hyperparameters
    expect_true(is.null(plm0.fit_known_c[['c_posterior']]))
    expect_false('c' %in% row.names(plm0.fit_known_c))
    test_stage_indep_param(plm0.fit_known_c, 'sigma_eps')
    # log-likelihood
    expect_true(is.double(plm0.fit_known_c$posterior_log_likelihood))
    expect_equal(length(plm0.fit_known_c$posterior_log_likelihood),
                 plm0.fit_known_c$run_info$num_chains * ((plm0.fit_known_c$run_info$nr_iter - plm0.fit_known_c$run_info$burnin) / plm0.fit_known_c$run_info$thin + 1))
    expect_equal(unname(unlist(plm0.fit_known_c$posterior_log_likelihood_summary[1, ])),
                 as.double(get_MCMC_summary(matrix(plm0.fit_known_c$posterior_log_likelihood, nrow = 1))),
                 tolerance = tol)
    # rating curve
    test_stage_dep_component(plm0.fit_known_c, 'rating_curve')
    test_stage_dep_component(plm0.fit_known_c, 'rating_curve_mean')
    # check if maxmimum stage was in line with output
    expect_equal(max(plm0.fit_known_c$rating_curve$h), h_extrap, tolerance = tol)
    expect_true(max(diff(plm0.fit_known_c$rating_curve$h)) <= (0.05 + tol)) # added tolerance
})


# C++ functions tests


test_that("plm0.density_evaluation_unknown_c works correctly", {
    RC <- get_model_components('plm0',
                               y = y,
                               h = h,
                               c_param = NULL,
                               h_max = max(h),
                               forcepoint = rep(FALSE, length(h)),
                               h_min = min(h))

    theta <- c(log(1), log(0.1))
    result <- plm0.density_evaluation_unknown_c(theta, RC)

    expect_type(result, "list")
    expect_true(all(c("p", "x", "y_post", "y_post_pred", "log_lik") %in% names(result)))
    expect_true(all(sapply(result, is.numeric)))
})

test_that("plm0.density_evaluation_known_c works correctly", {
    RC <- get_model_components('plm0',
                               y = y,
                               h = h,
                               c_param = min(h) - 0.1,
                               h_max = max(h),
                               forcepoint = rep(FALSE, length(h)),
                               h_min = min(h))

    theta <- c(log(0.1))
    result <- plm0.density_evaluation_known_c(theta, RC)

    expect_type(result, "list")
    expect_true(all(c("p", "x", "y_post", "y_post_pred", "log_lik") %in% names(result)))
    expect_true(all(sapply(result, is.numeric)))
})

test_that("plm0.predict_u_unknown_c works correctly", {
    RC <- get_model_components('plm0',
                               y = y,
                               h = h,
                               c_param = NULL,
                               h_max = max(h),
                               forcepoint = rep(FALSE, length(h)),
                               h_min = min(h))

    theta <- c(log(1), log(0.1))
    x <- c(1, 2)

    result <- plm0.predict_u_unknown_c(theta, x, RC)

    expect_type(result, "list")
    expect_true(all(c("y_post", "y_post_pred") %in% names(result)))
    expect_true(all(sapply(result, is.numeric)))
    expect_equal(length(result$y_post), length(RC$h_u))
    expect_equal(length(result$y_post_pred), length(RC$h_u))
})

test_that("plm0.predict_u_known_c works correctly", {
    RC <- get_model_components('plm0',
                               y = y,
                               h = h,
                               c_param = min(h) - 0.1,
                               h_max = max(h),
                               forcepoint = rep(FALSE, length(h)),
                               h_min = min(h))

    theta <- c(log(0.1))
    x <- c(1, 2)

    result <- plm0.predict_u_known_c(theta, x, RC)

    expect_type(result, "list")
    expect_true(all(c("y_post", "y_post_pred") %in% names(result)))
    expect_true(all(sapply(result, is.numeric)))
    expect_equal(length(result$y_post), length(RC$h_u))
    expect_equal(length(result$y_post_pred), length(RC$h_u))
})

test_that("plm0.calc_Dhat works correctly", {
    RC <- get_model_components('plm0',
                               y = y,
                               h = h,
                               c_param = NULL,
                               h_max = max(h),
                               forcepoint = rep(FALSE, length(h)),
                               h_min = min(h))

    theta <- matrix(c(log(1), log(0.1)), nrow = 2)
    result <- plm0.calc_Dhat(theta, RC)

    expect_type(result, "double")
    expect_true(is.finite(result))
})







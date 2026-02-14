context('gplm0')
tol <- 1e-8

test_that("gplm0 can handle different inputs", {
    expect_error(gplm0(Q ~ W, c(1, 2, 3)))
    expect_error(gplm0('Q ~ W', krokfors))
    expect_error(gplm0(V ~ W, krokfors))
    expect_error(gplm0(Q ~ W + X, krokfors))
    expect_error(gplm0(Q ~ W, krokfors, c_param = min(krokfors$W) + 0.5)) # c_param higher than lowest stage measurements
    expect_error(gplm0(Q ~ W, krokfors, c_param = 1L)) # c_param not double
    expect_error(gplm0(Q ~ W, krokfors, h_max = max(krokfors$W) - 0.5)) #h_max lower than highest stage measurement
    expect_error(gplm0(Q ~ W, krokfors[1,]), "At least two paired observations of stage and discharge")
    expect_error(gplm0(Q ~ W, -1 * krokfors), "All discharge measurements must but strictly greater than zero")
    skip_on_cran()
    krokfors_new_names <- krokfors
    names(krokfors_new_names) <- c('t1', 't2')
    set.seed(1)
    gplm0.fit_new_names <- gplm0(t2 ~ t1, krokfors_new_names, num_cores = 2)
    expect_equal(gplm0.fit_new_names$rating_curve, gplm0.fit$rating_curve, tolerance = tol)
})

test_that("the gplm0 object with unknown c is in tact", {
    expect_is(gplm0.fit, "gplm0")
    # latent parameters
    test_stage_indep_param(gplm0.fit, 'a')
    test_stage_indep_param(gplm0.fit, 'b')
    # hyperparameters
    test_stage_indep_param(gplm0.fit, 'c')
    test_stage_indep_param(gplm0.fit, 'sigma_eps')
    test_stage_indep_param(gplm0.fit, 'sigma_beta')
    test_stage_indep_param(gplm0.fit, 'phi_beta')
    # log-likelihood
    expect_true(is.double(gplm0.fit$posterior_log_likelihood))
    expect_equal(length(gplm0.fit$posterior_log_likelihood),
                 gplm0.fit$run_info$num_chains * ((gplm0.fit$run_info$nr_iter - gplm0.fit$run_info$burnin) / gplm0.fit$run_info$thin + 1))
    expect_equal(unname(unlist(gplm0.fit$posterior_log_likelihood_summary[1, ])),
                 as.double(get_MCMC_summary(matrix(gplm0.fit$posterior_log_likelihood, nrow = 1))),
                 tolerance = tol)
    # rating curve and stage dependent parameters
    test_stage_dep_component(gplm0.fit, 'rating_curve')
    test_stage_dep_component(gplm0.fit, 'rating_curve_mean')
    test_stage_dep_component(gplm0.fit, 'beta')
    test_stage_dep_component(gplm0.fit, 'f')
    # Other information
    expect_equal(gplm0.fit$formula, Q ~ W)
    expect_equal(gplm0.fit$data, krokfors[order(krokfors$W), c('Q', 'W')])
})

test_that("the gplm0 object with known c with a maximum stage value is in tact", {
    skip_on_cran()
    set.seed(1)
    gplm0.fit_known_c <- gplm0(Q ~ W, krokfors, c_param = known_c, h_max = h_extrap, num_cores = 2)
    expect_is(gplm0.fit_known_c, "gplm0")
    # latent parameters
    test_stage_indep_param(gplm0.fit_known_c, 'a')
    test_stage_indep_param(gplm0.fit_known_c, 'b')
    # hyperparameters
    expect_true(is.null(gplm0.fit_known_c[['c_posterior']]))
    expect_false('c' %in% row.names(gplm0.fit_known_c))
    test_stage_indep_param(gplm0.fit_known_c, 'sigma_eps')
    test_stage_indep_param(gplm0.fit, 'sigma_beta')
    test_stage_indep_param(gplm0.fit, 'phi_beta')
    # log-likelihood
    expect_true(is.double(gplm0.fit_known_c$posterior_log_likelihood))
    expect_equal(length(gplm0.fit_known_c$posterior_log_likelihood),
                 gplm0.fit_known_c$run_info$num_chains * ((gplm0.fit_known_c$run_info$nr_iter - gplm0.fit_known_c$run_info$burnin) / gplm0.fit_known_c$run_info$thin + 1))
    expect_equal(unname(unlist(gplm0.fit_known_c$posterior_log_likelihood_summary[1, ])),
                 as.double(get_MCMC_summary(matrix(gplm0.fit_known_c$posterior_log_likelihood, nrow = 1))),
                 tolerance = tol)
    # rating curve
    test_stage_dep_component(gplm0.fit_known_c, 'rating_curve')
    test_stage_dep_component(gplm0.fit_known_c, 'rating_curve_mean')
    test_stage_dep_component(gplm0.fit, 'beta')
    test_stage_dep_component(gplm0.fit, 'f')
    # check if maxmimum stage was in line with output
    expect_equal(max(gplm0.fit_known_c$rating_curve$h), h_extrap, tolerance = tol)
    expect_true(max(diff(gplm0.fit_known_c$rating_curve$h)) <= (0.05 + tol)) # added tolerance
})

test_that("gplm0 sends a warning about the estimated c_upper parameter", {
    skip_on_cran()
    W_grid <- seq(2.1,10,0.5)
    data_far_from_c <- data.frame("W" = W_grid,
                                  "Q" = 5 + exp(rnorm(length(W_grid), 0, 0.05)) * (W_grid ^ 2.5))
    set.seed(1)
    expect_warning(gplm0(Q ~ W,data_far_from_c, verbose = FALSE, num_cores = 2), "Dataset lacks measurements near point of zero flow and thus")
})

# C++ functions tests


test_that("gplm0.density_evaluation_unknown_c works correctly", {
    RC <- get_model_components('gplm0',
                               y = y,
                               h = h,
                               c_param = NULL,
                               h_max = max(h),
                               forcepoint = rep(FALSE, length(h)),
                               h_min = min(h))

    theta <- c(log(1), log(0.1), log(0.1), log(1))
    result <- gplm0.density_evaluation_unknown_c(theta, RC)

    expect_type(result, "list")
    expect_true("p" %in% names(result))
    expect_true(is.numeric(result$p))
    expect_true("x" %in% names(result))
    expect_true("y_post" %in% names(result))
    expect_true("y_post_pred" %in% names(result))
    expect_true("log_lik" %in% names(result))
})

test_that("gplm0.density_evaluation_known_c works correctly", {
    RC <- get_model_components('gplm0',
                               y = y,
                               h = h,
                               c_param = min(h) - 0.1,
                               h_max = max(h),
                               forcepoint = rep(FALSE, length(h)),
                               h_min = min(h))

    theta <- c(log(0.1), log(0.1), log(1))
    result <- gplm0.density_evaluation_known_c(theta, RC)

    expect_type(result, "list")
    expect_true("p" %in% names(result))
    expect_true(is.numeric(result$p))
    expect_true("x" %in% names(result))
    expect_true("y_post" %in% names(result))
    expect_true("y_post_pred" %in% names(result))
    expect_true("log_lik" %in% names(result))
})

test_that("gplm0.predict_u_unknown_c works correctly", {
    RC <- get_model_components('gplm0',
                               y = y,
                               h = h,
                               c_param = NULL,
                               h_max = max(h),
                               forcepoint = rep(FALSE, length(h)),
                               h_min = min(h))

    theta <- c(log(1), log(0.1), log(0.1), log(1))
    x <- c(1, 2, rep(0, length(unique(h))))

    result <- gplm0.predict_u_unknown_c(theta, x, RC)

    expect_type(result, "list")
    expect_true(all(c("x", "y_post", "y_post_pred") %in% names(result)))
    expect_true(all(sapply(result, is.numeric)))
    expect_equal(length(result$x), RC$n_u)
    expect_equal(length(result$y_post), RC$n_u)
    expect_equal(length(result$y_post_pred), RC$n_u)
})

test_that("gplm0.predict_u_known_c works correctly", {
    RC <- get_model_components('gplm0',
                               y = y,
                               h = h,
                               c_param = min(h) - 0.1,
                               h_max = max(h),
                               forcepoint = rep(FALSE, length(h)),
                               h_min = min(h))

    theta <- c(log(0.1), log(0.1), log(1))
    x <- c(1, 2, rep(0, length(unique(h))))

    result <- gplm0.predict_u_known_c(theta, x, RC)

    expect_type(result, "list")
    expect_true(all(c("x", "y_post", "y_post_pred") %in% names(result)))
    expect_true(all(sapply(result, is.numeric)))
    expect_equal(length(result$x), RC$n_u)
    expect_equal(length(result$y_post), RC$n_u)
    expect_equal(length(result$y_post_pred), RC$n_u)
})

test_that("gplm0.calc_Dhat works correctly", {
    RC <- get_model_components('gplm0',
                               y = y,
                               h = h,
                               c_param = NULL,
                               h_max = max(h),
                               forcepoint = rep(FALSE, length(h)),
                               h_min = min(h))

    theta <- matrix(c(log(1), log(0.1), log(0.1), log(1)), nrow = 4)
    result <- gplm0.calc_Dhat(theta, RC)

    expect_type(result, "double")
    expect_true(is.finite(result))
})


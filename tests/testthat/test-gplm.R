context('gplm')
tol <- 1e-8

test_that("gplm can handle different inputs", {
    expect_error(gplm(Q ~ W, c(1, 2, 3)))
    expect_error(gplm('Q ~ W', krokfors))
    expect_error(gplm(V ~ W, krokfors))
    expect_error(gplm(Q ~ W + X, krokfors))
    expect_error(gplm(Q ~ W, krokfors, c_param = min(krokfors$W) + 0.5)) # c_param higher than lowest stage measurements
    expect_error(gplm(Q ~ W, krokfors, c_param = 1L)) # c_param not double
    expect_error(gplm(Q ~ W, krokfors, h_max = max(krokfors$W) - 0.5)) #h_max lower than highest stage measurement
    expect_error(gplm(Q ~ W, krokfors[1,]), "At least two paired observations of stage and discharge")
    expect_error(gplm(Q ~ W, -1 * krokfors), "All discharge measurements must but strictly greater than zero")
    skip_on_cran()
    krokfors_new_names <- krokfors
    names(krokfors_new_names) <- c('t1', 't2')
    set.seed(1)
    gplm.fit_new_names <- gplm(t2 ~ t1, krokfors_new_names, num_cores = 2)
    expect_equal(gplm.fit_new_names$rating_curve, gplm.fit$rating_curve, tolerance = tol)
})

test_that("the gplm object with unknown c is in tact", {
    expect_is(gplm.fit, "gplm")
    # latent parameters
    test_stage_indep_param(gplm.fit, 'a')
    test_stage_indep_param(gplm.fit, 'b')
    # hyperparameters
    test_stage_indep_param(gplm.fit, 'c')
    test_stage_indep_param(gplm.fit, 'sigma_beta')
    test_stage_indep_param(gplm.fit, 'phi_beta')
    test_stage_indep_param(gplm.fit, 'sigma_eta')
    test_stage_indep_param(gplm.fit, 'eta_1')
    test_stage_indep_param(gplm.fit, 'eta_2')
    test_stage_indep_param(gplm.fit, 'eta_3')
    test_stage_indep_param(gplm.fit, 'eta_4')
    test_stage_indep_param(gplm.fit, 'eta_5')
    test_stage_indep_param(gplm.fit, 'eta_6')
    # log-likelihood
    expect_true(is.double(gplm.fit$posterior_log_likelihood))
    expect_equal(length(gplm.fit$posterior_log_likelihood),
                 gplm.fit$run_info$num_chains * ((gplm.fit$run_info$nr_iter - gplm.fit$run_info$burnin) / gplm.fit$run_info$thin + 1))
    expect_equal(unname(unlist(gplm.fit$posterior_log_likelihood_summary[1, ])),
                 as.double(get_MCMC_summary(matrix(gplm.fit$posterior_log_likelihood, nrow = 1))),
                 tolerance = tol)
    # rating curve and stage dependent parameters
    test_stage_dep_component(gplm.fit, 'rating_curve')
    test_stage_dep_component(gplm.fit, 'rating_curve_mean')
    test_stage_dep_component(gplm.fit, 'beta')
    test_stage_dep_component(gplm.fit, 'f')
    test_stage_dep_component(gplm.fit, 'sigma_eps')
    # Other information
    expect_equal(gplm.fit$formula, Q ~ W)
    expect_equal(gplm.fit$data, krokfors[order(krokfors$W), c('Q', 'W')])
})

test_that("the gplm object with known c with a maximum stage value is in tact", {
    skip_on_cran()
    set.seed(1)
    gplm.fit_known_c <- gplm(Q ~ W, krokfors, c_param = known_c, h_max = h_extrap, num_cores = 2)
    expect_is(gplm.fit_known_c, "gplm")
    # latent parameters
    test_stage_indep_param(gplm.fit_known_c, 'a')
    test_stage_indep_param(gplm.fit_known_c, 'b')
    # hyperparameters
    expect_true(is.null(gplm.fit_known_c[['c_posterior']]))
    expect_false('c' %in% row.names(gplm.fit_known_c))
    test_stage_indep_param(gplm.fit, 'sigma_beta')
    test_stage_indep_param(gplm.fit, 'phi_beta')
    test_stage_indep_param(gplm.fit, 'sigma_eta')
    test_stage_indep_param(gplm.fit, 'eta_1')
    test_stage_indep_param(gplm.fit, 'eta_2')
    test_stage_indep_param(gplm.fit, 'eta_3')
    test_stage_indep_param(gplm.fit, 'eta_4')
    test_stage_indep_param(gplm.fit, 'eta_5')
    test_stage_indep_param(gplm.fit, 'eta_6')
    # log-likelihood
    expect_true(is.double(gplm.fit_known_c$posterior_log_likelihood))
    expect_equal(length(gplm.fit_known_c$posterior_log_likelihood),
                 gplm.fit_known_c$run_info$num_chains * ((gplm.fit_known_c$run_info$nr_iter - gplm.fit_known_c$run_info$burnin) / gplm.fit_known_c$run_info$thin + 1))
    expect_equal(unname(unlist(gplm.fit_known_c$posterior_log_likelihood_summary[1, ])),
                 as.double(get_MCMC_summary(matrix(gplm.fit_known_c$posterior_log_likelihood, nrow = 1))),
                 tolerance = tol)
    # rating curve and stage dependent parameters
    test_stage_dep_component(gplm.fit, 'rating_curve')
    test_stage_dep_component(gplm.fit, 'rating_curve_mean')
    test_stage_dep_component(gplm.fit, 'beta')
    test_stage_dep_component(gplm.fit, 'f')
    test_stage_dep_component(gplm.fit, 'sigma_eps')
    # check if maxmimum stage was in line with output
    expect_equal(max(gplm.fit_known_c$rating_curve$h) ,h_extrap, tolerance = tol)
    expect_true(max(diff(gplm.fit_known_c$rating_curve$h)) <= (0.05 + tol)) # added tolerance
})

test_that("gplm sends a warning about the estimated c_upper parameter", {
    skip_on_cran()
    W_grid <- seq(2.1,10,0.5)
    data_far_from_c <- data.frame("W" = W_grid,
                                  "Q" = 5 + exp(rnorm(length(W_grid), 0, 0.05)) * (W_grid ^ 2.5))
    set.seed(1)
    expect_warning(gplm(Q ~ W,data_far_from_c, verbose = FALSE, num_cores = 2), "Dataset lacks measurements near point of zero flow and thus")
})

# C++ functions tests

test_that("gplm.density_evaluation_unknown_c works correctly", {
    RC <- get_model_components('gplm',
                               y = y,
                               h = h,
                               c_param = NULL,
                               h_max = max(h),
                               forcepoint = rep(FALSE, length(h)),
                               h_min = min(h))

    theta <- c(log(1), log(0.1), log(1), log(0.1), 0, runif(5))
    result <- gplm.density_evaluation_unknown_c(theta, RC)

    expect_type(result, "list")
    expect_true(all(c("p", "x", "y_post", "y_post_pred", "sigma_eps", "log_lik") %in% names(result)))
    expect_true(all(sapply(result, is.numeric)))
})

test_that("gplm.density_evaluation_known_c works correctly", {
    RC <- get_model_components('gplm',
                               y = y,
                               h = h,
                               c_param = min(h) - 0.1,
                               h_max = max(h),
                               forcepoint = rep(FALSE, length(h)),
                               h_min = min(h))

    theta <- c(log(0.1), log(1), log(0.1), 0, runif(5))
    result <- gplm.density_evaluation_known_c(theta, RC)

    expect_type(result, "list")
    expect_true(all(c("p", "x", "y_post", "y_post_pred", "sigma_eps", "log_lik") %in% names(result)))
    expect_true(all(sapply(result, is.numeric)))
})

test_that("gplm.predict_u_unknown_c works correctly", {
    RC <- get_model_components('gplm',
                               y = y,
                               h = h,
                               c_param = NULL,
                               h_max = max(h),
                               forcepoint = rep(FALSE, length(h)),
                               h_min = min(h))

    theta <- c(log(1), log(0.1), log(1), log(0.1), 0, runif(5))
    x <- c(1, 2, rep(0, length(unique(h))))

    result <- gplm.predict_u_unknown_c(theta, x, RC)

    expect_type(result, "list")
    expect_true(all(c("x", "sigma_eps", "y_post", "y_post_pred") %in% names(result)))
    expect_true(all(sapply(result, is.numeric)))
    expect_equal(length(result$x), RC$n_u)
    expect_equal(length(result$y_post), RC$n_u)
    expect_equal(length(result$y_post_pred), RC$n_u)
})

test_that("gplm.predict_u_known_c works correctly", {
    RC <- get_model_components('gplm',
                               y = y,
                               h = h,
                               c_param = min(h) - 0.1,
                               h_max = max(h),
                               forcepoint = rep(FALSE, length(h)),
                               h_min = min(h))

    theta <- c(log(0.1), log(1), log(0.1), 0, runif(5))
    x <- c(1, 2, rep(0, length(unique(h))))

    result <- gplm.predict_u_known_c(theta, x, RC)

    expect_type(result, "list")
    expect_true(all(c("x", "sigma_eps", "y_post", "y_post_pred") %in% names(result)))
    expect_true(all(sapply(result, is.numeric)))
    expect_equal(length(result$x), RC$n_u)
    expect_equal(length(result$y_post), RC$n_u)
    expect_equal(length(result$y_post_pred), RC$n_u)
})

test_that("gplm.calc_Dhat works correctly", {
    RC <- get_model_components('gplm',
                               y = y,
                               h = h,
                               c_param = NULL,
                               h_max = max(h),
                               forcepoint = rep(FALSE, length(h)),
                               h_min = min(h))

    theta <- matrix(c(log(1), log(0.1), log(1), log(0.1), 0, runif(5)), nrow = 10)
    result <- gplm.calc_Dhat(theta, RC)

    expect_type(result, "double")
    expect_true(is.finite(result))
})

context('plm')
tol <- 1e-8

test_that("plm can handle different inputs", {
    expect_error(plm(Q ~ W, c(1, 2, 3)))
    expect_error(plm('Q ~ W', krokfors))
    expect_error(plm(V ~ W, krokfors))
    expect_error(plm(Q ~ W + X, krokfors))
    expect_error(plm(Q ~ W, krokfors, c_param = min(krokfors$W) + 0.5)) # c_param higher than lowest stage measurements
    expect_error(plm(Q ~ W, krokfors, c_param = 1L)) # c_param not double
    expect_error(plm(Q ~ W, krokfors, h_max = max(krokfors$W) - 0.5)) #h_max lower than highest stage measurement
    expect_error(plm(Q ~ W, krokfors[1,]), "At least two paired observations of stage and discharge")
    expect_error(plm(Q ~ W, -1 * krokfors), "All discharge measurements must but strictly greater than zero")
    skip_on_cran()
    krokfors_new_names <- krokfors
    names(krokfors_new_names) <- c('t1', 't2')
    set.seed(1)
    plm.fit_new_names <- plm(t2 ~ t1, krokfors_new_names, num_cores = 2)
    expect_equal(plm.fit_new_names$rating_curve, plm.fit$rating_curve, tolerance = tol)
})

test_that("the plm object with unknown c is in tact", {
    expect_is(plm.fit, "plm")
    # latent parameters
    test_stage_indep_param(plm.fit, 'a')
    test_stage_indep_param(plm.fit, 'b')
    # hyperparameters
    test_stage_indep_param(plm.fit, 'c')
    test_stage_indep_param(plm.fit, 'sigma_eta')
    test_stage_indep_param(plm.fit, 'eta_1')
    test_stage_indep_param(plm.fit, 'eta_2')
    test_stage_indep_param(plm.fit, 'eta_3')
    test_stage_indep_param(plm.fit, 'eta_4')
    test_stage_indep_param(plm.fit, 'eta_5')
    test_stage_indep_param(plm.fit, 'eta_6')
    # log-likelihood
    expect_true(is.double(plm.fit$posterior_log_likelihood))
    expect_equal(length(plm.fit$posterior_log_likelihood),
                 plm.fit$run_info$num_chains * ((plm.fit$run_info$nr_iter - plm.fit$run_info$burnin) / plm.fit$run_info$thin + 1))
    expect_equal(unname(unlist(plm.fit$posterior_log_likelihood_summary[1, ])),
                 as.double(get_MCMC_summary(matrix(plm.fit$posterior_log_likelihood, nrow = 1))),
                 tolerance = tol)
    # rating curve and stage dependent parameters
    test_stage_dep_component(plm.fit,'rating_curve')
    test_stage_dep_component(plm.fit,'rating_curve_mean')
    test_stage_dep_component(plm.fit,'sigma_eps')
    # Other information
    expect_equal(plm.fit$formula, Q ~ W)
    expect_equal(plm.fit$data, krokfors[order(krokfors$W), c('Q', 'W')])
})

test_that("the plm object with known c with a maximum stage value is in tact", {
    skip_on_cran()
    set.seed(1)
    plm.fit_known_c <- plm(Q ~ W, krokfors, c_param = known_c, h_max = h_extrap, num_cores = 2)
    expect_is(plm.fit_known_c, "plm")
    # latent parameters
    test_stage_indep_param(plm.fit_known_c, 'a')
    test_stage_indep_param(plm.fit_known_c, 'b')
    # hyperparameters
    expect_true(is.null(plm.fit_known_c[['c_posterior']]))
    expect_false('c' %in% row.names(plm.fit_known_c))
    test_stage_indep_param(plm.fit, 'sigma_eta')
    test_stage_indep_param(plm.fit, 'eta_1')
    test_stage_indep_param(plm.fit, 'eta_2')
    test_stage_indep_param(plm.fit, 'eta_3')
    test_stage_indep_param(plm.fit, 'eta_4')
    test_stage_indep_param(plm.fit, 'eta_5')
    test_stage_indep_param(plm.fit, 'eta_6')
    # log-likelihood
    expect_true(is.double(plm.fit_known_c$posterior_log_likelihood))
    expect_equal(length(plm.fit_known_c$posterior_log_likelihood),
                 plm.fit_known_c$run_info$num_chains * ((plm.fit_known_c$run_info$nr_iter - plm.fit_known_c$run_info$burnin) / plm.fit_known_c$run_info$thin + 1))
    expect_equal(unname(unlist(plm.fit_known_c$posterior_log_likelihood_summary[1,])),
                 as.double(get_MCMC_summary(matrix(plm.fit_known_c$posterior_log_likelihood, nrow = 1))),
                 tolerance = tol)
    # rating curve and stage dependent parameters
    test_stage_dep_component(plm.fit, 'rating_curve')
    test_stage_dep_component(plm.fit, 'rating_curve_mean')
    test_stage_dep_component(plm.fit, 'sigma_eps')
    # check if maxmimum stage was in line with output
    expect_equal(max(plm.fit_known_c$rating_curve$h), h_extrap)
    expect_true(max(diff(plm.fit_known_c$rating_curve$h)) <= (0.05 + tol)) # added tolerance
})

# C++ functions tests

test_that("plm.density_evaluation_unknown_c works correctly", {
    RC <- get_model_components('plm',
                               y = y,
                               h = h,
                               c_param = NULL,
                               h_max = max(h),
                               forcepoint = rep(FALSE, length(h)),
                               h_min = min(h))

    theta <- c(log(1), log(0.1), 0, runif(5))
    result <- plm.density_evaluation_unknown_c(theta, RC)

    expect_type(result, "list")
    expect_true(all(c("p", "x", "y_post", "y_post_pred", "sigma_eps", "log_lik") %in% names(result)))
    expect_true(all(sapply(result, is.numeric)))
})

test_that("plm.density_evaluation_known_c works correctly", {
    RC <- get_model_components('plm',
                               y = y,
                               h = h,
                               c_param = min(h) - 0.1,
                               h_max = max(h),
                               forcepoint = rep(FALSE, length(h)),
                               h_min = min(h))

    theta <- c(log(0.1), 0, runif(5))
    result <- plm.density_evaluation_known_c(theta, RC)

    expect_type(result, "list")
    expect_true(all(c("p", "x", "y_post", "y_post_pred", "sigma_eps", "log_lik") %in% names(result)))
    expect_true(all(sapply(result, is.numeric)))
})

test_that("plm.predict_u_unknown_c works correctly", {
    RC <- get_model_components('plm',
                               y = y,
                               h = h,
                               c_param = NULL,
                               h_max = max(h),
                               forcepoint = rep(FALSE, length(h)),
                               h_min = min(h))

    theta <- c(log(1), log(0.1), 0, runif(5))
    x <- c(1, 2)

    result <- plm.predict_u_unknown_c(theta, x, RC)

    expect_type(result, "list")
    expect_true(all(c("y_post", "y_post_pred", "sigma_eps") %in% names(result)))
    expect_true(all(sapply(result, is.numeric)))
    expect_equal(length(result$y_post), length(RC$h_u))
    expect_equal(length(result$y_post_pred), length(RC$h_u))
})

test_that("plm.predict_u_known_c works correctly", {
    RC <- get_model_components('plm',
                               y = y,
                               h = h,
                               c_param = min(h) - 0.1,
                               h_max = max(h),
                               forcepoint = rep(FALSE, length(h)),
                               h_min = min(h))

    theta <- c(log(0.1), 0, runif(5))
    x <- c(1, 2)

    result <- plm.predict_u_known_c(theta, x, RC)

    expect_type(result, "list")
    expect_true(all(c("y_post", "y_post_pred", "sigma_eps") %in% names(result)))
    expect_true(all(sapply(result, is.numeric)))
    expect_equal(length(result$y_post), length(RC$h_u))
    expect_equal(length(result$y_post_pred), length(RC$h_u))
})

test_that("plm.calc_Dhat works correctly", {
    RC <- get_model_components('plm',
                               y = y,
                               h = h,
                               c_param = NULL,
                               h_max = max(h),
                               forcepoint = rep(FALSE, length(h)),
                               h_min = min(h))

    theta <- matrix(c(log(1), log(0.1), 0, runif(5)), nrow = 8)
    result <- plm.calc_Dhat(theta, RC)

    expect_type(result, "double")
    expect_true(is.finite(result))
})



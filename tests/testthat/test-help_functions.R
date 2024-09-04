context("R and C++ help functions")

# Test C++ functions
test_that("C++ matrix operations work correctly", {
    A <- matrix(1:4, nrow = 2)
    B <- matrix(5:8, nrow = 2)
    c <- matrix(1:2, ncol = 1)

    expect_equal(matMult(A, B), A %*% B)
    expect_equal(choleskyDecomp(A %*% t(A)), t(chol(A %*% t(A))))
    expect_equal(solveArma(A, c), solve(A, c))
    expect_equal(matInverse(A), solve(A))
})

test_that("compute_L and compute_w work correctly", {
    X <- matrix(1:4, nrow = 2)

    L <- compute_L(X, diag(2), diag(2), 1e-6)
    w <- compute_w(L, c(1, 2), X, c(0, 0))

    expect_type(L, "double")
    expect_type(w, "double")
    expect_equal(dim(L), c(2, 2))
    expect_equal(length(w), 2)
})

test_that("get_MCMC_summary_cpp works correctly", {
    X <- matrix(rnorm(100), nrow = 10)

    summary_without_h <- get_MCMC_summary_cpp(X, NULL)
    summary_with_h <- get_MCMC_summary_cpp(X, 1:10)

    expect_s3_class(summary_without_h, "data.frame")
    expect_s3_class(summary_with_h, "data.frame")
    expect_equal(nrow(summary_without_h), 10)
    expect_equal(nrow(summary_with_h), 10)
    expect_true("h" %in% names(summary_with_h))

    expect_error(get_MCMC_summary_cpp(X, 1:5), "Length of h must match the number of rows in X")
})

test_that("variogram_chain works correctly", {
    param_mat1 <- matrix(rnorm(100), nrow = 10)
    param_mat2 <- matrix(rnorm(100), nrow = 10)

    result <- variogram_chain(5, param_mat1, param_mat2, 10, 100)

    expect_type(result, "double")
    expect_equal(dim(result), c(10, 10))
})

test_that("distance_matrix works correctly", {
    x <- 1:5
    dist_mat <- distance_matrix(x)

    expect_type(dist_mat, "double")
    expect_equal(dim(dist_mat), c(5, 5))
    expect_true(all(dist_mat == t(dist_mat)))  # Should be symmetric
})

test_that("create_A_cpp works correctly", {
    A <- create_A_cpp(c(1, 1, 2, 3, 3, 4))

    expect_type(A, "double")
    expect_equal(dim(A), c(6, 4))
    expect_equal(colSums(A), c(2, 1, 2, 1))
})

test_that("pri function works correctly", {
    expect_type(pri("c", c(1, 2)), "double")
    expect_type(pri("sigma_eps2", c(1, 2)), "double")
    expect_type(pri("sigma_b", c(1, 2)), "double")
    expect_type(pri("phi_b", c(1, 2)), "double")
    expect_type(pri("eta_1", c(1, 2)), "double")
    expect_type(pri("eta_minus1", c(1, 2)), "double")
    expect_type(pri("sigma_eta", c(1, 2)), "double")
})

test_that("chain_statistics_cpp works correctly", {
    chains <- matrix(rnorm(100), nrow = 50)
    result <- chain_statistics_cpp(chains)

    expect_type(result, "list")
    expect_true(all(c("W", "var_hat") %in% names(result)))
})

# Test R functions
test_that("priors function works correctly", {
    result_plm0 <- priors("plm0")
    result_plm <- priors("plm")
    result_gplm0 <- priors("gplm0")
    result_gplm <- priors("gplm")

    expect_type(result_plm0, "list")
    expect_type(result_plm, "list")
    expect_type(result_gplm0, "list")
    expect_type(result_gplm, "list")
})

test_that("get_model_components function works correctly", {
    result <- get_model_components("plm0", y, h, NULL, NULL, rep(FALSE, nrow(krokfors)), NULL)

    expect_type(result, "list")
    expect_true(all(c("y", "h", "density_fun", "unobserved_prediction_fun") %in% names(result)))
})

test_that("get_MCMC_summary function works correctly", {
    X <- matrix(rnorm(100), nrow = 10)

    result <- get_MCMC_summary(X, 1:10)

    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 10)
    expect_true(all(c("h", "lower", "median", "upper") %in% names(result)))
})

test_that("get_param_names function works correctly", {
    expect_type(get_param_names("plm0", NULL), "character")
    expect_type(get_param_names("plm", NULL), "character")
    expect_type(get_param_names("gplm0", NULL), "character")
    expect_type(get_param_names("gplm", NULL), "character")
})

test_that("get_param_expression function works correctly", {
    expect_type(get_param_expression("a"), "character")
    expect_type(get_param_expression("b"), "character")
    expect_type(get_param_expression("c"), "character")
    expect_error(get_param_expression("invalid"), "param not found")
})

test_that("get_args_rollout function works correctly", {
    param_vec <- c("a", "b", "c", "sigma_eps")
    result <- get_args_rollout(c("latent_parameters", "hyperparameters"), param_vec)

    expect_type(result, "character")
    expect_equal(length(result), 4)
})

test_that("get_transformed_param function works correctly", {
    expect_type(get_transformed_param(1, "a"), "double")
    expect_type(get_transformed_param(1, "b"), "double")
    expect_type(get_transformed_param(1, "c", h_min = 2), "double")
    expect_error(get_transformed_param(1, "invalid"), "param not found")
})

test_that("h_unobserved function works correctly", {
    RC <- get_model_components("plm0", y, h, NULL, NULL, rep(FALSE, nrow(krokfors)), NULL)
    result <- h_unobserved(RC, 0.5, 11)

    expect_type(result, "double")
    expect_true(all(result >= 0.5 & result <= 11))
})

test_that("B_splines function works correctly", {
    ZZ <- seq(0, 1, length.out = 10)
    result <- B_splines(ZZ)

    expect_type(result, "double")
    expect_equal(dim(result), c(10, 6))
})

test_that("predict_wider function works correctly", {
    p_dat <- data.frame(h = seq(1, 2, by = 0.1), median = runif(11))
    result <- predict_wider(p_dat)

    expect_type(result, "double")
    expect_true(is.matrix(result))
    expect_equal(ncol(result), 10)
})

test_that("chain_statistics function works correctly", {
    chains <- matrix(rnorm(100), nrow = 50)
    result <- chain_statistics(chains)

    expect_type(result, "list")
    expect_true(all(c("W", "var_hat") %in% names(result)))
})

test_that("R_hat function works correctly", {
    chains <- matrix(rnorm(100), nrow = 50)
    result <- R_hat(chains)

    expect_type(result, "double")
    expect_length(result, 1)
})

test_that("LSE and log_mean_LSE functions work correctly", {
    lx <- log(1:10)
    expect_type(LSE(lx), "double")
    expect_type(log_mean_LSE(lx), "double")
})

test_that("various functions exist", {
    expect_true(is.function(get_rhat_dat))
    expect_true(is.function(log_lik_i))
    expect_true(is.function(calc_waic))
    expect_true(is.function(log_ml_harmonic_mean_est))
    expect_true(is.function(post_model_prob_m1))
    expect_true(is.function(SE_Delta_WAIC))
    expect_true(is.function(get_MCMC_output_list))
    expect_true(is.function(run_MCMC))
    expect_true(is.function(get_residuals_dat))
})


test_that("convergence_diagnostics_warnings function works correctly", {
    param_summary <- data.frame(
        r_hat = c(1.05, 1.2),
        eff_n_samples = c(500, 300),
        row.names = c("param1", "param2")
    )

    expect_output(convergence_diagnostics_warnings(param_summary), "Warning:")
})

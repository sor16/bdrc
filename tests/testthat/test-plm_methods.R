context('plm general functions')
test_that("trace plots have known output", {
    skip_on_cran()
    skip_on_ci()
    vdiffr::expect_doppelganger('plm0_trace_plot_hyperparameters', autoplot(plm0.fit, type = 'trace', 'hyperparameters'))
    vdiffr::expect_doppelganger('plm0_trace_plot_hyperparameters_transformed', autoplot(plm0.fit, type = 'trace', 'hyperparameters', transformed = TRUE))
    vdiffr::expect_doppelganger('plm0_trace_plot_single_parameter', autoplot(plm0.fit, type = 'trace', 'sigma_eps'))

    vdiffr::expect_doppelganger('plm_trace_plot_hyperparameters', autoplot(plm.fit, type = 'trace', 'hyperparameters'))
    vdiffr::expect_doppelganger('plm_trace_plot_hyperparameters_transformed', autoplot(plm.fit, type = 'trace', 'hyperparameters', transformed = TRUE))
    vdiffr::expect_doppelganger('plm_trace_plot_single_parameter', autoplot(plm.fit, type = 'trace', 'sigma_eta'))

    vdiffr::expect_doppelganger('gplm0_trace_plot_hyperparameters', autoplot(gplm0.fit, type = 'trace', 'hyperparameters'))
    vdiffr::expect_doppelganger('gplm0_trace_plot_hyperparameters_transformed', autoplot(gplm0.fit, type = 'trace', 'hyperparameters', transformed = TRUE))
    vdiffr::expect_doppelganger('gplm0_trace_plot_single_parameter', autoplot(gplm0.fit, type = 'trace', 'sigma_eps'))

    vdiffr::expect_doppelganger('gplm_trace_plot_hyperparameters', autoplot(gplm.fit, type = 'trace', 'hyperparameters'))
    vdiffr::expect_doppelganger('gplm_trace_plot_hyperparameters_transformed', autoplot(gplm.fit, type = 'trace', 'hyperparameters', transformed = TRUE))
    vdiffr::expect_doppelganger('gplm_trace_plot_single_parameter', autoplot(gplm.fit, type = 'trace', 'sigma_eta'))
})
test_that("histogram have known output", {
    skip_on_cran()
    skip_on_ci()
    vdiffr::expect_doppelganger('plm0_histogram_hyperparameters', autoplot(plm0.fit, type = 'histogram', 'hyperparameters'))
    vdiffr::expect_doppelganger('plm0_histogram_hyperparameters_transformed', autoplot(plm0.fit, type = 'histogram', 'hyperparameters', transformed = TRUE))
    vdiffr::expect_doppelganger('plm0_histogram_single_parameter', autoplot(plm0.fit, type = 'histogram', 'sigma_eps'))

    vdiffr::expect_doppelganger('plm_histogram_plot_hyperparameters', autoplot(plm.fit, type = 'histogram', 'hyperparameters'))
    vdiffr::expect_doppelganger('plm_histogram_plot_hyperparameters_transformed', autoplot(plm.fit, type = 'histogram', 'hyperparameters', transformed = TRUE))
    vdiffr::expect_doppelganger('plm_histogram_plot_single_parameter',autoplot(plm.fit,type='histogram','sigma_eta'))

    vdiffr::expect_doppelganger('gplm0_histogram_plot_hyperparameters', autoplot(gplm0.fit, type = 'histogram', 'hyperparameters'))
    vdiffr::expect_doppelganger('gplm0_histogram_plot_hyperparameters_transformed', autoplot(gplm0.fit, type = 'histogram', 'hyperparameters', transformed = TRUE))
    vdiffr::expect_doppelganger('gplm0_histogram_plot_single_parameter', autoplot(gplm0.fit, type = 'histogram', 'sigma_eps'))

    vdiffr::expect_doppelganger('gplm_histogram_plot_hyperparameters', autoplot(gplm.fit, type = 'histogram', 'hyperparameters'))
    vdiffr::expect_doppelganger('gplm_histogram_plot_hyperparameters_transformed', autoplot(gplm.fit, type = 'histogram', 'hyperparameters', transformed = TRUE))
    vdiffr::expect_doppelganger('gplm_histogram_plot_single_parameter', autoplot(gplm.fit, type = 'histogram', 'sigma_eta'))
})

test_that("rating_curve has known output", {
    skip_on_cran()
    skip_on_ci()
    vdiffr::expect_doppelganger('plm0_rating_curve', autoplot(plm0.fit, type = 'rating_curve'))
    vdiffr::expect_doppelganger('plm0_rating_curve_transformed', autoplot(plm0.fit, type = 'rating_curve', transformed = TRUE))

    vdiffr::expect_doppelganger('plm_rating_curve', autoplot(plm.fit, type = 'rating_curve'))
    vdiffr::expect_doppelganger('plm_rating_curve_transformed', autoplot(plm.fit, type = 'rating_curve', transformed = TRUE))

    vdiffr::expect_doppelganger('gplm0_rating_curve', autoplot(gplm0.fit, type = 'rating_curve'))
    vdiffr::expect_doppelganger('gplm0_rating_curve_transformed', autoplot(gplm0.fit, type = 'rating_curve', transformed = TRUE))

    vdiffr::expect_doppelganger('gplm_rating_curve', autoplot(gplm.fit, type = 'rating_curve'))
    vdiffr::expect_doppelganger('gplm_rating_curve_transformed', autoplot(gplm.fit, type = 'rating_curve', transformed = TRUE))
})

test_that("rating_curve_mean has known output", {
    skip_on_cran()
    skip_on_ci()
    vdiffr::expect_doppelganger('plm0_rating_curve_mean', autoplot(plm0.fit, type = 'rating_curve_mean'))
    vdiffr::expect_doppelganger('plm0_rating_curve_mean_transformed', autoplot(plm0.fit, type = 'rating_curve_mean', transformed = TRUE))

    vdiffr::expect_doppelganger('plm_rating_curve_mean', autoplot(plm.fit, type = 'rating_curve_mean'))
    vdiffr::expect_doppelganger('plm_rating_curve_mean_transformed', autoplot(plm.fit, type = 'rating_curve_mean', transformed = TRUE))

    vdiffr::expect_doppelganger('gplm0_rating_curve_mean', autoplot(gplm0.fit, type = 'rating_curve_mean'))
    vdiffr::expect_doppelganger('gplm0_rating_curve_mean_transformed', autoplot(gplm0.fit, type = 'rating_curve_mean', transformed = TRUE))

    vdiffr::expect_doppelganger('gplm_rating_curve_mean', autoplot(gplm.fit, type = 'rating_curve_mean'))
    vdiffr::expect_doppelganger('gplm_rating_curve_mean_transformed', autoplot(gplm.fit, type = 'rating_curve_mean', transformed = TRUE))
})


test_that("sigma_eps has known output", {
    skip_on_cran()
    skip_on_ci()
    vdiffr::expect_doppelganger('plm0_sigma_eps', autoplot(plm0.fit, type = 'sigma_eps'))
    vdiffr::expect_doppelganger('plm_sigma_eps', autoplot(plm.fit, type = 'sigma_eps'))
    vdiffr::expect_doppelganger('gplm0_sigma_eps', autoplot(gplm0.fit, type = 'sigma_eps'))
    vdiffr::expect_doppelganger('gplm_sigma_eps', autoplot(gplm.fit, type = 'sigma_eps'))
})

test_that("beta has known output", {
    skip_on_cran()
    skip_on_ci()
    vdiffr::expect_doppelganger('gplm0_beta', autoplot(gplm0.fit, type = 'beta'))
    vdiffr::expect_doppelganger('gplm_beta', autoplot(gplm.fit, type = 'beta'))
})

test_that("f has known output", {
    skip_on_cran()
    skip_on_ci()
    vdiffr::expect_doppelganger('plm0_f', autoplot(plm0.fit, type = 'f'))
    vdiffr::expect_doppelganger('plm_f', autoplot(plm.fit, type = 'f'))
    vdiffr::expect_doppelganger('gplm0_f', autoplot(gplm0.fit, type = 'f'))
    vdiffr::expect_doppelganger('gplm_f', autoplot(gplm.fit, type = 'f'))
})

test_that("residuals has known output", {
    skip_on_cran()
    skip_on_ci()
    vdiffr::expect_doppelganger('plm0_residuals', autoplot(plm0.fit, type = 'residuals'))
    vdiffr::expect_doppelganger('plm_residuals', autoplot(plm.fit, type = 'residuals'))
    vdiffr::expect_doppelganger('gplm0_residuals', autoplot(gplm0.fit, type = 'residuals'))
    vdiffr::expect_doppelganger('gplm_residuals', autoplot(gplm.fit, type = 'residuals'))
})

test_that("print methods work for all model types", {
    expect_output(print(plm0.fit), "plm0 - Call:")
    expect_output(print(plm.fit), "plm - Call:")
    expect_output(print(gplm0.fit), "gplm0 - Call:")
    expect_output(print(gplm.fit), "gplm - Call:")
})

test_that("summary methods work for all model types", {
    models <- list(plm0 = plm0.fit, plm = plm.fit, gplm0 = gplm0.fit, gplm = gplm.fit)

    for (model_name in names(models)) {
        summary_output <- capture.output(summary(models[[model_name]]))
        expect_true(any(grepl("Formula:", summary_output)))
        expect_true(any(grepl("Latent parameters:", summary_output)))
        expect_true(any(grepl("Hyperparameters:", summary_output)))
        expect_true(any(grepl("WAIC:", summary_output)))
    }
})

test_that("autoplot methods work for all model types and plot types", {
    models <- list(plm0 = plm0.fit, plm = plm.fit, gplm0 = gplm0.fit, gplm = gplm.fit)
    plot_types <- c("rating_curve", "rating_curve_mean", "f", "sigma_eps", "residuals")

    for (model_name in names(models)) {
        for (plot_type in plot_types) {
            expect_s3_class(autoplot(models[[model_name]], type = plot_type), "ggplot")
        }

        # Test trace and histogram plots
        expect_s3_class(autoplot(models[[model_name]], type = "trace", param = "a"), "ggplot")
        expect_s3_class(autoplot(models[[model_name]], type = "histogram", param = "b"), "ggplot")

        # Test transformed plots
        expect_s3_class(autoplot(models[[model_name]], type = "rating_curve", transformed = TRUE), "ggplot")
    }

    # Test beta plot for gplm0 and gplm
    expect_s3_class(autoplot(gplm0.fit, type = "beta"), "ggplot")
    expect_s3_class(autoplot(gplm.fit, type = "beta"), "ggplot")

    # Test error for beta plot on plm0 and plm
    expect_error(autoplot(plm0.fit, type = "beta"))
    expect_error(autoplot(plm.fit, type = "beta"))
})

test_that("autoplot works with transformed data and different parameters", {
    expect_s3_class(autoplot(plm.fit, type = "trace", param = "hyperparameters"), "ggplot")
    expect_s3_class(autoplot(plm.fit, type = "histogram", param = "latent_parameters"), "ggplot")
    expect_s3_class(autoplot(plm.fit, type = "trace", param = c("a", "b", "sigma_eta")), "ggplot")
})

test_that("predict methods work for all model types", {
    models <- list(plm0 = plm0.fit, plm = plm.fit, gplm0 = gplm0.fit, gplm = gplm.fit)

    for (model_name in names(models)) {
        model <- models[[model_name]]

        pred <- predict(model)
        expect_s3_class(pred, "data.frame")
        expect_equal(colnames(pred), c("h", "lower", "median", "upper"))

        newdata <- seq(min(model$data$W), max(model$data$W), length.out = 10)
        pred_new <- predict(model, newdata = newdata)
        expect_equal(nrow(pred_new), 10)

        pred_wide <- predict(model, wide = TRUE)
        expect_true(is.matrix(pred_wide))
        expect_equal(ncol(pred_wide), 10)  # 10 columns for cm values, 1 for dm

        # Test with a single value
        single_pred <- predict(model, newdata = median(model$data$W))
        expect_equal(nrow(single_pred), 1)

        # Test error handling
        expect_error(predict(model, newdata = c(0, max(model$data$W) + 1)), "newdata must contain values within the range")
        expect_error(predict(model, newdata = c(NA, 1, 2)), "newdata must not include NA")
        expect_error(predict(model, newdata = "invalid"), "newdata must be a vector of type")
    }
})

test_that("plm handles known and unknown c parameter correctly", {
    known_c <- min(krokfors$W) - 0.1
    plm_known_c <- plm(Q ~ W, data = krokfors, c_param = known_c, num_cores = 2)

    expect_null(plm_known_c$c_posterior)
    expect_s3_class(autoplot(plm_known_c, type = "rating_curve"), "ggplot")
    expect_s3_class(autoplot(plm_known_c, type = "f"), "ggplot")

    summary_output <- capture.output(summary(plm_known_c))
    expect_false(any(grepl("c", summary_output)))
})

test_that("error handling works for all model types", {
    models <- list(plm0 = plm0.fit, plm = plm.fit, gplm0 = gplm0.fit, gplm = gplm.fit)

    for (model_name in names(models)) {
        model <- models[[model_name]]

        expect_error(autoplot(model, type = "invalid_type"), "Type argument not recognized")
        expect_error(plot(model, type = "invalid_type"), "Type argument not recognized")
        expect_error(autoplot(model, type = "trace", param = "invalid_param"), "Does not recognize the following input arguments")
        expect_error(plot(model, type = "histogram", param = "invalid_param"), "Does not recognize the following input arguments")
    }
})

test_that("plot methods work silently for all model types", {
    models <- list(plm0 = plm0.fit, plm = plm.fit, gplm0 = gplm0.fit, gplm = gplm.fit)
    plot_types <- c('rating_curve', 'rating_curve_mean', 'f', 'sigma_eps', 'residuals',
                    'trace', 'histogram', 'panel', 'convergence_diagnostics')

    for (model_name in names(models)) {
        model <- models[[model_name]]

        for (plot_type in plot_types) {
            if (plot_type %in% c('trace', 'histogram')) {
                expect_silent(plot(model, type = plot_type, param = 'a'))
            } else {
                expect_silent(plot(model, type = plot_type))
            }

            # Test transformed parameter for applicable plot types
            if (plot_type %in% c('rating_curve', 'rating_curve_mean', 'f', 'sigma_eps', 'residuals')) {
                expect_silent(plot(model, type = plot_type, transformed = TRUE))
            }
        }

        # Test autoplot
        expect_silent(autoplot(model))
        expect_silent(autoplot(model, transformed = TRUE))
    }
})



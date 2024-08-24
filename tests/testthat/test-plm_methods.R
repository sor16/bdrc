context('plm general functions')
test_that("trace plots have known output", {
    skip_on_cran()
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
    vdiffr::expect_doppelganger('plm0_sigma_eps', autoplot(plm0.fit, type = 'sigma_eps'))
    vdiffr::expect_doppelganger('plm_sigma_eps', autoplot(plm.fit, type = 'sigma_eps'))
    vdiffr::expect_doppelganger('gplm0_sigma_eps', autoplot(gplm0.fit, type = 'sigma_eps'))
    vdiffr::expect_doppelganger('gplm_sigma_eps', autoplot(gplm.fit, type = 'sigma_eps'))
})

test_that("beta has known output", {
    skip_on_cran()
    expect_error(autoplot(plm0.fit, type = 'beta'))
    expect_error(autoplot(plm.fit, type = 'beta'))
    vdiffr::expect_doppelganger('gplm0_beta', autoplot(gplm0.fit, type = 'beta'))
    vdiffr::expect_doppelganger('gplm_beta', autoplot(gplm.fit, type = 'beta'))
})

test_that("f has known output", {
    skip_on_cran()
    vdiffr::expect_doppelganger('plm0_f', autoplot(plm0.fit, type = 'f'))
    vdiffr::expect_doppelganger('plm_f', autoplot(plm.fit, type = 'f'))
    vdiffr::expect_doppelganger('gplm0_f', autoplot(gplm0.fit, type = 'f'))
    vdiffr::expect_doppelganger('gplm_f', autoplot(gplm.fit, type = 'f'))
})

test_that("residuals has known output", {
    skip_on_cran()
    vdiffr::expect_doppelganger('plm0_residuals', autoplot(plm0.fit, type = 'residuals'))
    vdiffr::expect_doppelganger('plm_residuals', autoplot(plm.fit, type = 'residuals'))
    vdiffr::expect_doppelganger('gplm0_residuals', autoplot(gplm0.fit, type = 'residuals'))
    vdiffr::expect_doppelganger('gplm_residuals', autoplot(gplm.fit, type = 'residuals'))
})


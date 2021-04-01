context('plm general functions')

test_that("trace plots have known output", {
    expect_doppelganger('bplm0_trace_plot_hyperparameters',autoplot(bplm0.fit,type='trace','hyperparameters'))
    expect_doppelganger('bplm0_trace_plot_hyperparameters_transformed',autoplot(bplm0.fit,type='trace','hyperparameters',transformed=T))
    expect_doppelganger('bplm0_trace_plot_single_parameter',autoplot(bplm0.fit,type='trace','sigma_eps'))

    expect_doppelganger('bplm_trace_plot_hyperparameters',autoplot(bplm.fit,type='trace','hyperparameters'))
    expect_doppelganger('bplm_trace_plot_hyperparameters_transformed',autoplot(bplm.fit,type='trace','hyperparameters',transformed=T))
    expect_doppelganger('bplm_trace_plot_single_parameter',autoplot(bplm.fit,type='trace','sigma_eta'))

    expect_doppelganger('bgplm0_trace_plot_hyperparameters',autoplot(bgplm0.fit,type='trace','hyperparameters'))
    expect_doppelganger('bgplm0_trace_plot_hyperparameters_transformed',autoplot(bgplm0.fit,type='trace','hyperparameters',transformed=T))
    expect_doppelganger('bgplm0_trace_plot_single_parameter',autoplot(bgplm0.fit,type='trace','sigma_eps'))

    expect_doppelganger('bgplm_trace_plot_hyperparameters',autoplot(bgplm.fit,type='trace','hyperparameters'))
    expect_doppelganger('bgplm_trace_plot_hyperparameters_transformed',autoplot(bgplm.fit,type='trace','hyperparameters',transformed=T))
    expect_doppelganger('bgplm_trace_plot_single_parameter',autoplot(bgplm.fit,type='trace','sigma_eta'))
})
test_that("histogram have known output", {
    expect_doppelganger('bplm0_histogram_hyperparameters',autoplot(bplm0.fit,type='histogram','hyperparameters'))
    expect_doppelganger('bplm0_histogram_hyperparameters_transformed',autoplot(bplm0.fit,type='histogram','hyperparameters',transformed=T))
    expect_doppelganger('bplm0_histogram_single_parameter',autoplot(bplm0.fit,type='histogram','sigma_eps'))

    expect_doppelganger('bplm_histogram_plot_hyperparameters',autoplot(bplm.fit,type='histogram','hyperparameters'))
    expect_doppelganger('bplm_histogram_plot_hyperparameters_transformed',autoplot(bplm.fit,type='histogram','hyperparameters',transformed=T))
    expect_doppelganger('bplm_histogram_plot_single_parameter',autoplot(bplm.fit,type='histogram','sigma_eta'))

    expect_doppelganger('bgplm0_histogram_plot_hyperparameters',autoplot(bgplm0.fit,type='histogram','hyperparameters'))
    expect_doppelganger('bgplm0_histogram_plot_hyperparameters_transformed',autoplot(bgplm0.fit,type='histogram','hyperparameters',transformed=T))
    expect_doppelganger('bgplm0_histogram_plot_single_parameter',autoplot(bgplm0.fit,type='histogram','sigma_eps'))

    expect_doppelganger('bgplm_histogram_plot_hyperparameters',autoplot(bgplm.fit,type='histogram','hyperparameters'))
    expect_doppelganger('bgplm_histogram_plot_hyperparameters_transformed',autoplot(bgplm.fit,type='histogram','hyperparameters',transformed=T))
    expect_doppelganger('bgplm_histogram_plot_single_parameter',autoplot(bgplm.fit,type='histogram','sigma_eta'))
})

test_that("rating_curve has known output", {
    expect_doppelganger('bplm0_rating_curve',autoplot(bplm0.fit,type='rating_curve'))
    expect_doppelganger('bplm0_rating_curve_transformed',autoplot(bplm0.fit,type='rating_curve',transformed=T))

    expect_doppelganger('bplm_rating_curve',autoplot(bplm.fit,type='rating_curve'))
    expect_doppelganger('bplm_rating_curve_transformed',autoplot(bplm.fit,type='rating_curve',transformed=T))

    expect_doppelganger('bgplm0_rating_curve',autoplot(bgplm0.fit,type='rating_curve'))
    expect_doppelganger('bgplm0_rating_curve_transformed',autoplot(bgplm0.fit,type='rating_curve',transformed=T))

    expect_doppelganger('bgplm_rating_curve',autoplot(bgplm.fit,type='rating_curve'))
    expect_doppelganger('bgplm_rating_curve_transformed',autoplot(bgplm.fit,type='rating_curve',transformed=T))
})

test_that("rating_curve_mean has known output", {
    expect_doppelganger('bplm0_rating_curve_mean',autoplot(bplm0.fit,type='rating_curve_mean'))
    expect_doppelganger('bplm0_rating_curve_mean_transformed',autoplot(bplm0.fit,type='rating_curve_mean',transformed=T))

    expect_doppelganger('bplm_rating_curve_mean',autoplot(bplm.fit,type='rating_curve_mean'))
    expect_doppelganger('bplm_rating_curve_mean_transformed',autoplot(bplm.fit,type='rating_curve_mean',transformed=T))

    expect_doppelganger('bgplm0_rating_curve_mean',autoplot(bgplm0.fit,type='rating_curve_mean'))
    expect_doppelganger('bgplm0_rating_curve_mean_transformed',autoplot(bgplm0.fit,type='rating_curve_mean',transformed=T))

    expect_doppelganger('bgplm_rating_curve_mean',autoplot(bgplm.fit,type='rating_curve_mean'))
    expect_doppelganger('bgplm_rating_curve_mean_transformed',autoplot(bgplm.fit,type='rating_curve_mean',transformed=T))
})


test_that("sigma_eps has known output", {
    expect_doppelganger('bplm0_sigma_eps',autoplot(bplm0.fit,type='sigma_eps'))
    expect_doppelganger('bplm_sigma_eps',autoplot(bplm.fit,type='sigma_eps'))
    expect_doppelganger('bgplm0_sigma_eps',autoplot(bgplm0.fit,type='sigma_eps'))
    expect_doppelganger('bgplm_sigma_eps',autoplot(bgplm.fit,type='sigma_eps'))
})

test_that("beta has known output", {
    expect_error(autoplot(bplm0.fit,type='beta'))
    expect_error(autoplot(bplm.fit,type='beta'))
    expect_doppelganger('bgplm0_beta',autoplot(bgplm0.fit,type='beta'))
    expect_doppelganger('bgplm_beta',autoplot(bgplm.fit,type='beta'))
})

test_that("f has known output", {
    expect_doppelganger('bplm0_f',autoplot(bplm0.fit,type='f'))
    expect_doppelganger('bplm_f',autoplot(bplm.fit,type='f'))
    expect_doppelganger('bgplm0_f',autoplot(bgplm0.fit,type='f'))
    expect_doppelganger('bgplm_f',autoplot(bgplm.fit,type='f'))
})

test_that("residuals has known output", {
    expect_doppelganger('bplm0_residuals',autoplot(bplm0.fit,type='residuals'))
    expect_doppelganger('bplm_residuals',autoplot(bplm.fit,type='residuals'))
    expect_doppelganger('bgplm0_residuals',autoplot(bgplm0.fit,type='residuals'))
    expect_doppelganger('bgplm_residuals',autoplot(bgplm.fit,type='residuals'))
})





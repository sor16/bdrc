context('Tournament Methods')

test_that("tournament plots can handle different inputs", {
    expect_error(autoplot(t_obj, type = 'tournament_results'))
})

test_that("plot.tournament works correctly", {
    expect_error(plot.tournament(t_obj, type = NULL), "Type not recognized. Possible types are:")
    expect_error(plot.tournament(t_obj, type = "incorrect"), "Type not recognized. Possible types are:")
})

test_that("summary.tournament works correctly", {
    expect_output(summary(t_obj), "Tournament Model Comparison Summary")
    expect_output(summary(t_obj, method = "DIC"), "Method: DIC")
    expect_output(summary(t_obj, winning_criteria = 3), "Winning Criteria: Delta_WAIC > 3")
})

test_that("print.tournament works correctly", {
    expect_output(print(t_obj), "Tournament winner:")
})

test_that("autoplot.tournament works correctly", {
    expect_s3_class(autoplot(t_obj), "ggplot")
    expect_error(autoplot(t_obj, type = 'invalid_type'), "Type argument not recognized")
})


test_that("summary.tournament handles different methods and criteria", {
    expect_output(summary(t_obj, method = "WAIC", winning_criteria = 2), "Method: WAIC")
    expect_output(summary(t_obj, method = "DIC", winning_criteria = 3), "Method: DIC")
    expect_output(summary(t_obj, method = "PMP", winning_criteria = .5), "Method: PMP")

    expect_output(summary(t_obj, winning_criteria = "Delta_WAIC > 5"), "Winning Criteria: Delta_WAIC > 5")
    expect_output(summary(t_obj, method = "PMP"), "The Harmonic Mean Estimator")

})

test_that("tournament_summary_output handles different methods", {
    waic_results <- tournament(t_obj$contestants, method = "WAIC")$summary
    expect_output(tournament_summary_output(waic_results, "WAIC", 2), "Tournament Model Comparison Summary")

    dic_results <- tournament(t_obj$contestants, method = "DIC")$summary
    expect_output(tournament_summary_output(dic_results, "DIC", 2), "Tournament Model Comparison Summary")

    pmp_results <- tournament(t_obj$contestants, method = "PMP")$summary
    expect_output(tournament_summary_output(pmp_results, "PMP", 0.5), "Tournament Model Comparison Summary")
})

test_that("tournament_summary_output handles errors", {
    waic_results <- tournament(t_obj$contestants, method = "WAIC")$summary
    expect_error(tournament_summary_output(waic_results, "Invalid", 2), "Unknown method")
})

test_that("plot_tournament_fun handles boxplot type", {
    p <- plot_tournament_fun(t_obj, type = 'boxplot')
    expect_s3_class(p, "ggplot")
})

test_that("plot_tournament_grob handles different types", {
    types <- c("residuals", "sigma_eps", "f", "rating_curve", "rating_curve_mean",
               "convergence_diagnostics", "panel", "tournament_results")

    for(type in types) {
        p <- plot_tournament_grob(t_obj, type = type)
        expect_true(inherits(p, "grob") | is.list(p))
    }

    expect_error(plot_tournament_grob(t_obj, type = "invalid"), "type is not recognized")
})

test_that("tournament_summary_output formats output correctly", {
    waic_results <- tournament(t_obj$contestants, method = "WAIC")$summary
    output <- capture.output(tournament_summary_output(waic_results, "WAIC", 2))
    expect_true(any(grepl("complexity", output)))
    expect_true(any(grepl("model", output)))
    expect_true(any(grepl("winner", output)))
    expect_true(any(grepl("WAIC", output)))

    dic_results <- tournament(t_obj$contestants, method = "DIC")$summary
    output <- capture.output(tournament_summary_output(dic_results, "DIC", 2))
    expect_true(any(grepl("DIC", output)))

    pmp_results <- tournament(t_obj$contestants, method = "PMP")$summary
    output <- capture.output(tournament_summary_output(pmp_results, "PMP", 0.5))
    expect_true(any(grepl("PMP", output)))
})

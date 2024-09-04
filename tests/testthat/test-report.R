context('Report')


test_that("get_report_components handles different input types correctly", {
    expect_error(get_report_components("invalid_input"), "Please provide a single object of types tournament, gplm, gplm0, plm or plm0.")
    expect_error(get_report_components(plm0.fit, type = 3), "Please input an integer value of 1 or 2 to indicate which type of report is to be produce.")
    expect_error(get_report_components(plm0.fit, type = 2), "It is only possible to produce a type 1 report for a single model object of type gplm, gplm0, plm or plm0.")

    # Test with mock objects
    plm0_report <- get_report_components(plm0.fit, type = 1)
    expect_type(plm0_report, "list")
    expect_equal(plm0_report$obj_class, "plm0")

    tournament_report_1 <- get_report_components(t_obj, type = 1)
    expect_type(tournament_report_1, "list")
    expect_in(tournament_report_1$obj_class, c("plm0","plm","gplm0","gplm"))

    tournament_report_2 <- get_report_components(t_obj, type = 2)
    expect_type(tournament_report_2, "list")
    expect_true("tournament_summary" %in% names(tournament_report_2))
})

test_that("get_report_pages_fun produces correct output", {
    plm0_pages <- get_report_pages_fun(plm0.fit, type = 1)
    expect_type(plm0_pages, "list")
    expect_true(all(sapply(plm0_pages, inherits, "grob")))

    tournament_pages <- get_report_pages_fun(t_obj, type = 2)
    expect_type(tournament_pages, "list")
    expect_true(all(sapply(tournament_pages, inherits, "grob")))
})

test_that("get_report_pages methods work correctly", {
    expect_type(get_report_pages(plm0.fit), "list")
    expect_error(get_report_pages(plm0.fit, type = 2), "It is only possible to produce a type 1 report for a single model object of type gplm, gplm0, plm or plm0.")
    expect_type(get_report_pages(t_obj, type = 2), "list")
})

test_that("report components are correct", {
    report_cmp_1 <- get_report_components(plm0.fit, type = 1)
    t_report_cmp_1 <- get_report_components(t_obj, type = 1)
    t_report_cmp_2 <- get_report_components(t_obj, type = 2)
    expect_equal(names(report_cmp_1), c('main_page_plots', 'main_page_table', 'p_mat_list', 'obj_class'))
    expect_equal(names(t_report_cmp_1), c('main_page_plots', 'main_page_table', 'p_mat_list', 'obj_class'))
    expect_equal(names(t_report_cmp_2), c('main_page_plots', 'main_page_table', 'tournament_summary', 'tournament_plot', 'conv_diag_plots', 'mcmc_hist_list'))
    expect_equal(report_cmp_1$obj_class, 'plm0')
    expect_equal(t_report_cmp_1$obj_class, class(t_obj$winner))
})

test_that("get_report_components handles different input types correctly", {
    expect_error(get_report_components("invalid_input"), "Please provide a single object of types tournament, gplm, gplm0, plm or plm0.")
    expect_error(get_report_components(plm0.fit, type = 3), "Please input an integer value of 1 or 2 to indicate which type of report is to be produce.")
    expect_error(get_report_components(plm0.fit, type = 2), "It is only possible to produce a type 1 report for a single model object of type gplm, gplm0, plm or plm0.")

    plm0_report <- get_report_components(plm0.fit, type = 1)
    expect_type(plm0_report, "list")
    expect_equal(plm0_report$obj_class, "plm0")

    tournament_report_1 <- get_report_components(t_obj, type = 1)
    expect_type(tournament_report_1, "list")
    expect_true(tournament_report_1$obj_class %in% c("plm0","plm","gplm0","gplm"))

    tournament_report_2 <- get_report_components(t_obj, type = 2)
    expect_type(tournament_report_2, "list")
    expect_true("tournament_summary" %in% names(tournament_report_2))
})

test_that("get_report_pages_fun produces correct output", {
    plm0_pages <- get_report_pages_fun(plm0.fit, type = 1)
    expect_type(plm0_pages, "list")
    expect_true(all(sapply(plm0_pages, inherits, "grob")))

    tournament_pages <- get_report_pages_fun(t_obj, type = 2)
    expect_type(tournament_pages, "list")
    expect_true(all(sapply(tournament_pages, inherits, "grob")))
    expect_length(tournament_pages, 10)
})

test_that("get_report_pages methods work correctly", {
    expect_type(get_report_pages(plm0.fit), "list")
    expect_error(get_report_pages(plm0.fit, type = 2), "It is only possible to produce a type 1 report for a single model object of type gplm, gplm0, plm or plm0.")
    expect_type(get_report_pages(t_obj, type = 2), "list")
})

test_that("report components are correct", {
    report_cmp_1 <- get_report_components(plm0.fit, type = 1)
    t_report_cmp_1 <- get_report_components(t_obj, type = 1)
    t_report_cmp_2 <- get_report_components(t_obj, type = 2)

    expect_equal(names(report_cmp_1), c('main_page_plots', 'main_page_table', 'p_mat_list', 'obj_class'))
    expect_equal(names(t_report_cmp_1), c('main_page_plots', 'main_page_table', 'p_mat_list', 'obj_class'))
    expect_equal(names(t_report_cmp_2), c('main_page_plots', 'main_page_table', 'tournament_summary', 'tournament_plot', 'conv_diag_plots', 'mcmc_hist_list'))
    expect_equal(report_cmp_1$obj_class, 'plm0')
    expect_equal(t_report_cmp_1$obj_class, class(t_obj$winner))
})

test_that("get_report methods work correctly", {
    expect_true(is.function(get_report.plm0))
    expect_true(is.function(get_report.plm))
    expect_true(is.function(get_report.gplm0))
    expect_true(is.function(get_report.gplm))
    expect_true(is.function(get_report.tournament))
})

test_that("get_report_components handles all model types correctly", {
    expect_type(get_report_components(plm0.fit, type = 1), "list")
    expect_type(get_report_components(plm.fit, type = 1), "list")
    expect_type(get_report_components(gplm0.fit, type = 1), "list")
    expect_type(get_report_components(gplm.fit, type = 1), "list")
})

test_that("get_report_components handles tournament objects correctly", {
    t_components_1 <- get_report_components(t_obj, type = 1)
    expect_type(t_components_1, "list")
    expect_equal(t_components_1$obj_class, class(t_obj$winner))

    t_components_2 <- get_report_components(t_obj, type = 2)
    expect_type(t_components_2, "list")
    expect_true("tournament_summary" %in% names(t_components_2))
    expect_true("tournament_plot" %in% names(t_components_2))
})


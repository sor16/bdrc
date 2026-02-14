context('Tournament')

test_that("tournament handles errors correctly", {
    mod_list <- list(plm0.fit, plm0.fit, gplm0.fit, gplm.fit)
    expect_error(tournament(data = mod_list))
    expect_error(tournament(model_list = list(plm0.fit, plm.fit, gplm0.fit)))
    expect_error(tournament(mod_list, method = NA ))
    expect_error(tournament(mod_list, method = NULL ))
    expect_error(tournament(mod_list, method = "invalid"))
    expect_error(tournament(mod_list, winning_criteria = "1.5"))
    expect_error(tournament(Q ~ W, krokfors, winning_criteria = "invalid"))
    expect_error(tournament(mod_list, method = "PMP", winning_criteria = 1.5))
    expect_error(tournament(t_obj))

    mod_list_error <- list(plm0.fit, plm0.fit, gplm0.fit, plm0.fit) # two plm0 model objects
    expect_error(tournament(model_list = mod_list_error))
})

test_that("the tournament object is in tact", {
    expect_is(t_obj, "tournament")
    expect_is(t_obj$contestants$plm0, "plm0")
    expect_is(t_obj$contestants$plm, "plm")
    expect_is(t_obj$contestants$gplm0, "gplm0")
    expect_is(t_obj$contestants$gplm, "gplm")
    expect_named(t_obj, c("contestants", "winner", "summary", "info"))
    expect_length(t_obj$contestants, 4)
    expect_true(t_obj$info$method %in% c("WAIC", "DIC", "PMP"))
})

test_that("tournament handles different methods", {
    t_waic <- tournament(t_obj$contestants, method = "WAIC")
    t_dic <- tournament(t_obj$contestants, method = "DIC")
    t_pmp <- tournament(t_obj$contestants, method = "PMP")

    expect_equal(t_waic$info$method, "WAIC")
    expect_equal(t_dic$info$method, "DIC")
    expect_equal(t_pmp$info$method, "PMP")
})

test_that("tournament handles custom winning criteria", {
    t_custom <- tournament(t_obj$contestants, winning_criteria = 3)
    t_expr <- tournament(t_obj$contestants, winning_criteria = "Delta_WAIC > 2 & Delta_WAIC - SE_Delta_WAIC > 0")
    t_dic <- tournament(t_obj$contestants, method="DIC", winning_criteria = 4)
    t_pmp <- tournament(t_obj$contestants, method="PMP", winning_criteria = 0.9)

    expect_equal(t_custom$info$winning_criteria, "Delta_WAIC > 3")
    expect_equal(t_expr$info$winning_criteria, "Delta_WAIC > 2 & Delta_WAIC - SE_Delta_WAIC > 0")
    expect_equal(t_dic$info$winning_criteria, 4)
    expect_equal(t_pmp$info$winning_criteria, 0.9)
})

test_that("evaluate_comparison works correctly", {
    models <- list(
        plm0.fit,
        plm.fit
    )

    waic_result <- evaluate_comparison(models, "WAIC", 2)
    dic_result <- evaluate_comparison(models, "DIC", 2)
    pmp_result <- evaluate_comparison(models, "PMP", 0.75)

    expect_s3_class(waic_result, "data.frame")
    expect_s3_class(dic_result, "data.frame")
    expect_s3_class(pmp_result, "data.frame")

    expect_true("Delta_WAIC" %in% names(waic_result))
    expect_true("Delta_DIC" %in% names(dic_result))
    expect_true("PMP" %in% names(pmp_result))
})

test_that("tournament handles warnings correctly", {
    expect_output(tournament(t_obj$contestants, method = "PMP"), "The Harmonic Mean Estimator")
})

test_that("tournament handles different input types", {
    expect_s3_class(t_obj, "tournament")
    expect_s3_class(tournament(t_obj$contestants), "tournament")
})

test_that("tournament handles verbose input correctly", {
    skip_on_cran()
    t_obj_output <- paste(capture.output(tournament(Q ~ W, krokfors, num_cores = 2)), collapse = "")
    expect_true( grepl("Running tournament", t_obj_output))
    expect_true( grepl("gplm finished", t_obj_output))
    expect_true( grepl("gplm0 finished", t_obj_output))
    expect_true( grepl("plm finished", t_obj_output))
    expect_true( grepl("plm0 finished", t_obj_output))

    t_obj_output <- paste(capture.output(tournament(Q ~ W, krokfors, verbose = FALSE, num_cores = 2)), collapse = "")
    expect_false(grepl("Running tournament", t_obj_output))
})

test_that("tournament summary has correct structure", {
    expect_true(all(c("round", "comparison", "complexity", "model", "winner") %in% names(t_obj$summary)))
    expect_equal(nrow(t_obj$summary), 6)
})

test_that("evaluate_comparison handles string winning criteria for WAIC", {
    models <- list(plm0.fit, plm.fit)
    result <- evaluate_comparison(models, "WAIC", "Delta_WAIC > 2 & Delta_WAIC - SE_Delta_WAIC > 0")
    expect_s3_class(result, "data.frame")
    expect_true("winner" %in% names(result))
})

test_that("tournament handles unexpected arguments", {
    expect_error(tournament(Q ~ W, krokfors, num_cores = 2, unexpected_arg = TRUE),
                 "The following argument\\(s\\) are not recognized: unexpected_arg")
    expect_error(tournament(t_obj, method = "WAIC", winning_criteria = Inf),
                 "For method 'WAIC', when numeric, winning_criteria must be a single finite number.")
    expect_error(tournament(t_obj, method = "WAIC", winning_criteria = NA),
                 "For method 'WAIC', winning_criteria must be either a numeric value or a string expression.")
})

test_that("tournament checks for consistent data in model_list", {
    skip_on_cran()
    # Create a new gplm model with different data
    different_data_plm0 <- plm0(Q ~ W, data = krokfors[1:4, ], num_cores = 2)

    inconsistent_data_models <- list(
        plm0 = different_data_plm0,  # Different data
        plm = plm.fit,
        gplm0 = gplm0.fit,
        gplm = gplm.fit
    )
    expect_error(tournament(model_list = inconsistent_data_models), "The four models added have to be fit on the same data set")
})

test_that("tournament handles invalid winning criteria", {
    skip_on_cran()
    expect_error(tournament(Q ~ W, krokfors, winning_criteria = "invalid = expression"),
                 "For method 'WAIC', when a string, winning_criteria must be a valid R expression")
    expect_error(tournament(Q ~ W, krokfors, winning_criteria = "Delta_WAIC = 2"),  # wrong use of "=="
                 "Use '==' for equality in winning_criteria expressions, not '='")
    expect_error(tournament(Q ~ W, krokfors, method = "DIC", winning_criteria = "invalid"),
                 'For method "DIC", winning_criteria must be a single real number')
    expect_error(tournament(Q ~ W, krokfors, method = "PMP", winning_criteria = 2),
                 "For method 'PMP', winning_criteria must be a single number between 0 and 1")
})


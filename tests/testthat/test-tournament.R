context('Tournament')

test_that("tournament can handle different inputs", {
    expect_error(tournament(plm0.fit,plm0.fit,gplm0.fit,gplm.fit))
    expect_error(tournament(plm0.fit,plm.fit,gplm0.fit))
})

test_that("the tournament object is in tact", {
    expect_is(t_obj,"tournament")
    expect_is(t_obj$contestants$plm0,"plm0")
    expect_is(t_obj$contestants$plm,"plm")
    expect_is(t_obj$contestants$gplm0,"gplm0")
    expect_is(t_obj$contestants$gplm,"gplm")
    #sapply(1:nrow(t_obj$summary),function(i) expect_equal(t_obj$contestants[[t_obj$summary$model[i]]]$DIC,t_obj$summary$DIC[i])) is failing on older R version. Skip for now
})

test_that("Tournament object remains the same", {
    skip_on_cran()
    skip_on_ci()
    skip_on_covr()
    expect_equal_to_reference(t_obj,file='../cached_results/tournament.rds',update=TRUE)
})


context('Tournament')

test_that("tournament can handle different inputs", {
    expect_error(tournament(bplm0.fit,bplm0.fit,bgplm0.fit,bgplm.fit))
    expect_error(tournament(bplm0.fit,bplm.fit,bgplm0.fit))
})

test_that("the tournament object is in tact", {
    expect_is(t_obj,"tournament")
    expect_is(t_obj$contestants$bplm0,"bplm0")
    expect_is(t_obj$contestants$bplm,"bplm")
    expect_is(t_obj$contestants$bgplm0,"bgplm0")
    expect_is(t_obj$contestants$bgplm,"bgplm")
})



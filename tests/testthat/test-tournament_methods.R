context('Tournament Methods')

test_that("tournament plots can handle different inputs", {
    expect_error(autoplot(t_obj,type='tournament_results'))
})

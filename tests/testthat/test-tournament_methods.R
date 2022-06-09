context('Tournament Methods')

test_that("tournament plots can handle different inputs", {
    expect_error(plot(t_obj,type=NULL))
    expect_error(plot(t_obj,type=1))
    expect_error(autoplot(t_obj,type='tournament_results'))
})

context('bplm0')

test_that("bplm0 can handle different inputs", {
    expect_error(bplm0(Q~W,c(1,2,3)))
    expect_error(bplm0('Q~W',V316_river))
    expect_error(bplm0(V~W,V316_river))
    expect_error(bplm0(Q~W+X,V316_river))
    expect_error(bplm0(Q~W,V316_river,c_param=1.5)) # c_param higher than lowest stage measurements
    expect_error(bplm0(Q~W,V316_river,c_param=1L)) # c_param not double
    expect_error(bplm0(Q~W,V316_river,h_max=1.3)) #h_max lower than highest stage measurement
    V316_river_new_names <- V316_river
    names(V316_river_new_names) <- c('t1','t2')
    set.seed(1)
    bplm0.fit_new_names <- bplm0(t2~t1,V316_river_new_names)
    expect_equal(bplm0.fit_new_names$rating_curve,bplm0.fit$rating_curve)
})

test_that("the bplm0 object with unknown c is in tact", {
    expect_is(bplm0.fit,"bplm0")
    #latent parameters
    test_stage_indep_param(bplm0.fit,'a')
    test_stage_indep_param(bplm0.fit,'b')
    #hyperparameters
    test_stage_indep_param(bplm0.fit,'c')
    test_stage_indep_param(bplm0.fit,'sigma_eps')
    #Deviance
    expect_true(is.double(bplm0.fit$Deviance_posterior))
    expect_equal(length(bplm0.fit$Deviance_posterior),bplm0.fit$run_info$num_chains*((bplm0.fit$run_info$nr_iter-bplm0.fit$run_info$burnin)/bplm0.fit$run_info$thin + 1))
    expect_equal(unname(unlist(bplm0.fit$Deviance_summary[1,])),unname(quantile(bplm0.fit$Deviance_posterior,probs=c(0.025,0.5,0.975))))
    #rating curve
    test_stage_dep_component(bplm0.fit,'rating_curve')
    test_stage_dep_component(bplm0.fit,'rating_curve_mean')
    #Other information
    expect_equal(bplm0.fit$formula,Q~W)
    expect_equal(bplm0.fit$data,V316_river[order(V316_river$W),c('Q','W')])
})

test_that("the bplm0 object with known c with a maximum stage value is in tact", {
    expect_is(bplm0.fit_known_c,"bplm0")
    #latent parameters
    test_stage_indep_param(bplm0.fit_known_c,'a')
    test_stage_indep_param(bplm0.fit_known_c,'b')
    #hyperparameters
    expect_true(is.null(bplm0.fit_known_c[['c_posterior']]))
    expect_false('c' %in% row.names(bplm0.fit_known_c))
    test_stage_indep_param(bplm0.fit_known_c,'sigma_eps')
    #Deviance
    expect_true(is.double(bplm0.fit_known_c$Deviance_posterior))
    expect_equal(length(bplm0.fit_known_c$Deviance_posterior),bplm0.fit_known_c$run_info$num_chains*((bplm0.fit_known_c$run_info$nr_iter-bplm0.fit_known_c$run_info$burnin)/bplm0.fit_known_c$run_info$thin + 1))
    expect_equal(unname(unlist(bplm0.fit_known_c$Deviance_summary[1,])),unname(quantile(bplm0.fit_known_c$Deviance_posterior,probs=c(0.025,0.5,0.975))))
    #rating curve
    test_stage_dep_component(bplm0.fit_known_c,'rating_curve')
    test_stage_dep_component(bplm0.fit_known_c,'rating_curve_mean')
    #check if maxmimum stage was in line with output
    expect_equal(max(bplm0.fit_known_c$rating_curve$h),2.5)
    expect_true(max(diff(bplm0.fit_known_c$rating_curve$h))<=(0.05+1e-9)) # added tolerance
})


test_that("bplm0 output remains unchanged", {
    skip_on_cran()
    skip_on_ci()
    skip_on_covr()
    expect_equal_to_reference(bplm0.fit,file='../cached_results/bplm0.fit.rds',update=T)
    expect_equal_to_reference(bplm0.fit_known_c,file='../cached_results/bplm0.fit_known_c.rds',update=T)
})

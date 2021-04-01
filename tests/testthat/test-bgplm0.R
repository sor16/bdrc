context('bgplm0')

test_that("bgplm0 can handle different inputs", {
    expect_error(bgplm0(Q~W,c(1,2,3)))
    expect_error(bgplm0('Q~W',V316_river))
    expect_error(bgplm0(V~W,V316_river))
    expect_error(bgplm0(Q~W+X,V316_river))
    expect_error(bgplm0(Q~W,V316_river,c_param=1.5)) # c_param higher than lowest stage measurements
    expect_error(bgplm0(Q~W,V316_river,c_param=1L)) # c_param not double
    expect_error(bgplm0(Q~W,V316_river,h_max=1.3)) #h_max lower than highest stage measurement
    skip_on_cran()
    V316_river_new_names <- V316_river
    names(V316_river_new_names) <- c('t1','t2')
    set.seed(1)
    bgplm0.fit_new_names <- bgplm0(t2~t1,V316_river_new_names,parallel=F)
    expect_equal(bgplm0.fit_new_names$rating_curve,bgplm0.fit$rating_curve)
})

test_that("the bgplm0 object with unknown c is in tact", {
    expect_is(bgplm0.fit,"bgplm0")
    #latent parameters
    test_stage_indep_param(bgplm0.fit,'a')
    test_stage_indep_param(bgplm0.fit,'b')
    #hyperparameters
    test_stage_indep_param(bgplm0.fit,'c')
    test_stage_indep_param(bgplm0.fit,'sigma_eps')
    test_stage_indep_param(bgplm0.fit,'sigma_beta')
    test_stage_indep_param(bgplm0.fit,'phi_beta')
    #Deviance
    expect_true(is.double(bgplm0.fit$Deviance_posterior))
    expect_equal(length(bgplm0.fit$Deviance_posterior),bgplm0.fit$run_info$num_chains*((bgplm0.fit$run_info$nr_iter-bgplm0.fit$run_info$burnin)/bgplm0.fit$run_info$thin + 1))
    expect_equal(unname(unlist(bgplm0.fit$Deviance_summary[1,])),unname(quantile(bgplm0.fit$Deviance_posterior,probs=c(0.025,0.5,0.975))))
    #rating curve and stage dependent parameters
    test_stage_dep_component(bgplm0.fit,'rating_curve')
    test_stage_dep_component(bgplm0.fit,'rating_curve_mean')
    test_stage_dep_component(bgplm0.fit,'beta')
    test_stage_dep_component(bgplm0.fit,'f')
    #Other information
    expect_equal(bgplm0.fit$formula,Q~W)
    expect_equal(bgplm0.fit$data,V316_river[order(V316_river$W),c('Q','W')])
})

test_that("the bgplm0 object with known c with a maximum stage value is in tact", {
    expect_is(bgplm0.fit_known_c,"bgplm0")
    #latent parameters
    test_stage_indep_param(bgplm0.fit_known_c,'a')
    test_stage_indep_param(bgplm0.fit_known_c,'b')
    #hyperparameters
    expect_true(is.null(bgplm0.fit_known_c[['c_posterior']]))
    expect_false('c' %in% row.names(bgplm0.fit_known_c))
    test_stage_indep_param(bgplm0.fit_known_c,'sigma_eps')
    test_stage_indep_param(bgplm0.fit,'sigma_beta')
    test_stage_indep_param(bgplm0.fit,'phi_beta')
    #Deviance
    expect_true(is.double(bgplm0.fit_known_c$Deviance_posterior))
    expect_equal(length(bgplm0.fit_known_c$Deviance_posterior),bgplm0.fit_known_c$run_info$num_chains*((bgplm0.fit_known_c$run_info$nr_iter-bgplm0.fit_known_c$run_info$burnin)/bgplm0.fit_known_c$run_info$thin + 1))
    expect_equal(unname(unlist(bgplm0.fit_known_c$Deviance_summary[1,])),unname(quantile(bgplm0.fit_known_c$Deviance_posterior,probs=c(0.025,0.5,0.975))))
    #rating curve
    test_stage_dep_component(bgplm0.fit_known_c,'rating_curve')
    test_stage_dep_component(bgplm0.fit_known_c,'rating_curve_mean')
    test_stage_dep_component(bgplm0.fit,'beta')
    test_stage_dep_component(bgplm0.fit,'f')
    #check if maxmimum stage was in line with output
    expect_equal(max(bgplm0.fit_known_c$rating_curve$h),2)
    expect_true(max(diff(bgplm0.fit_known_c$rating_curve$h))<=(0.05+1e-9)) # added tolerance
})


test_that("bgplm0 output remains unchanged", {
    skip_on_cran()
    skip_on_ci()
    skip_on_covr()
    expect_equal_to_reference(bgplm0.fit,file='../cached_results/bgplm0.fit.rds',update=T)
    expect_equal_to_reference(bgplm0.fit_known_c,file='../cached_results/bgplm0.fit_known_c.rds',update=T)
})

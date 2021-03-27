context('bgplm')

test_that("bgplm can handle different inputs", {
    expect_error(bgplm(Q~W,c(1,2,3)))
    expect_error(bgplm('Q~W',V316_river))
    expect_error(bgplm(V~W,V316_river))
    expect_error(bgplm(Q~W+X,V316_river))
    expect_error(bgplm(Q~W,V316_river,c_param=1.5)) # c_param higher than lowest stage measurements
    expect_error(bgplm(Q~W,V316_river,c_param=1L)) # c_param not double
    expect_error(bgplm(Q~W,V316_river,h_max=1.3)) #h_max lower than highest stage measurement
    V316_river_new_names <- V316_river
    names(V316_river_new_names) <- c('t1','t2')
    set.seed(1)
    bgplm.fit_new_names <- bgplm(t2~t1,V316_river_new_names)
    expect_equal(bgplm.fit_new_names$rating_curve,bgplm.fit$rating_curve)
})

test_that("the bgplm object with unknown c is in tact", {
    expect_is(bgplm.fit,"bgplm")
    #latent parameters
    test_stage_indep_param(bgplm.fit,'a')
    test_stage_indep_param(bgplm.fit,'b')
    #hyperparameters
    test_stage_indep_param(bgplm.fit,'c')
    test_stage_indep_param(bgplm.fit,'sigma_beta')
    test_stage_indep_param(bgplm.fit,'phi_beta')
    test_stage_indep_param(bgplm.fit,'sigma_eta')
    test_stage_indep_param(bgplm.fit,'eta_1')
    test_stage_indep_param(bgplm.fit,'eta_2')
    test_stage_indep_param(bgplm.fit,'eta_3')
    test_stage_indep_param(bgplm.fit,'eta_4')
    test_stage_indep_param(bgplm.fit,'eta_5')
    test_stage_indep_param(bgplm.fit,'eta_6')
    #Deviance
    expect_true(is.double(bgplm.fit$Deviance_posterior))
    expect_equal(length(bgplm.fit$Deviance_posterior),bgplm.fit$run_info$num_chains*((bgplm.fit$run_info$nr_iter-bgplm.fit$run_info$burnin)/bgplm.fit$run_info$thin + 1))
    expect_equal(unname(unlist(bgplm.fit$Deviance_summary[1,])),unname(quantile(bgplm.fit$Deviance_posterior,probs=c(0.025,0.5,0.975))))
    #rating curve and stage dependent parameters
    test_stage_dep_component(bgplm.fit,'rating_curve')
    test_stage_dep_component(bgplm.fit,'rating_curve_mean')
    test_stage_dep_component(bgplm.fit,'beta')
    test_stage_dep_component(bgplm.fit,'f')
    test_stage_dep_component(bgplm.fit,'sigma_eps')
    #Other information
    expect_equal(bgplm.fit$formula,Q~W)
    expect_equal(bgplm.fit$data,V316_river[order(V316_river$W),c('Q','W')])
})

test_that("the bgplm object with known c with a maximum stage value is in tact", {
    expect_is(bgplm.fit_known_c,"bgplm")
    #latent parameters
    test_stage_indep_param(bgplm.fit_known_c,'a')
    test_stage_indep_param(bgplm.fit_known_c,'b')
    #hyperparameters
    expect_true(is.null(bgplm.fit_known_c[['c_posterior']]))
    expect_false('c' %in% row.names(bgplm.fit_known_c))
    test_stage_indep_param(bgplm.fit,'sigma_beta')
    test_stage_indep_param(bgplm.fit,'phi_beta')
    test_stage_indep_param(bgplm.fit,'sigma_eta')
    test_stage_indep_param(bgplm.fit,'eta_1')
    test_stage_indep_param(bgplm.fit,'eta_2')
    test_stage_indep_param(bgplm.fit,'eta_3')
    test_stage_indep_param(bgplm.fit,'eta_4')
    test_stage_indep_param(bgplm.fit,'eta_5')
    test_stage_indep_param(bgplm.fit,'eta_6')
    #Deviance
    expect_true(is.double(bgplm.fit_known_c$Deviance_posterior))
    expect_equal(length(bgplm.fit_known_c$Deviance_posterior),bgplm.fit_known_c$run_info$num_chains*((bgplm.fit_known_c$run_info$nr_iter-bgplm.fit_known_c$run_info$burnin)/bgplm.fit_known_c$run_info$thin + 1))
    expect_equal(unname(unlist(bgplm.fit_known_c$Deviance_summary[1,])),unname(quantile(bgplm.fit_known_c$Deviance_posterior,probs=c(0.025,0.5,0.975))))
    #rating curve and stage dependent parameters
    test_stage_dep_component(bgplm.fit,'rating_curve')
    test_stage_dep_component(bgplm.fit,'rating_curve_mean')
    test_stage_dep_component(bgplm.fit,'beta')
    test_stage_dep_component(bgplm.fit,'f')
    test_stage_dep_component(bgplm.fit,'sigma_eps')
    #check if maxmimum stage was in line with output
    expect_equal(max(bgplm.fit_known_c$rating_curve$h),2.5)
    expect_true(max(diff(bgplm.fit_known_c$rating_curve$h))<=(0.05+1e-9)) # added tolerance
})


test_that("bgplm output remains unchanged", {
    skip_on_cran()
    skip_on_ci()
    skip_on_covr()
    expect_equal_to_reference(bgplm.fit,file='../cached_results/bgplm.fit.rds',update=T)
    expect_equal_to_reference(bgplm.fit_known_c,file='../cached_results/bgplm.fit_known_c.rds',update=T)
})

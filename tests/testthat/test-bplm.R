context('bplm')

test_that("bplm can handle different inputs", {
    expect_error(bplm(Q~W,c(1,2,3)))
    expect_error(bplm('Q~W',V316_river))
    expect_error(bplm(V~W,V316_river))
    expect_error(bplm(Q~W+X,V316_river))
    expect_error(bplm(Q~W,V316_river,c_param=1.5)) # c_param higher than lowest stage measurements
    expect_error(bplm(Q~W,V316_river,c_param=1L)) # c_param not double
    expect_error(bplm(Q~W,V316_river,h_max=1.3)) #h_max lower than highest stage measurement
    V316_river_new_names <- V316_river
    names(V316_river_new_names) <- c('t1','t2')
    set.seed(1)
    bplm.fit_new_names <- bplm(t2~t1,V316_river_new_names)
    expect_equal(bplm.fit_new_names$rating_curve,bplm.fit$rating_curve)
})

test_that("the bplm object with unknown c is in tact", {
    expect_is(bplm.fit,"bplm")
    #latent parameters
    test_stage_indep_param(bplm.fit,'a')
    test_stage_indep_param(bplm.fit,'b')
    #hyperparameters
    test_stage_indep_param(bplm.fit,'c')
    test_stage_indep_param(bplm.fit,'sigma_eta')
    test_stage_indep_param(bplm.fit,'eta_1')
    test_stage_indep_param(bplm.fit,'eta_2')
    test_stage_indep_param(bplm.fit,'eta_3')
    test_stage_indep_param(bplm.fit,'eta_4')
    test_stage_indep_param(bplm.fit,'eta_5')
    test_stage_indep_param(bplm.fit,'eta_6')
    #Deviance
    expect_true(is.double(bplm.fit$Deviance_posterior))
    expect_equal(length(bplm.fit$Deviance_posterior),bplm.fit$run_info$num_chains*((bplm.fit$run_info$nr_iter-bplm.fit$run_info$burnin)/bplm.fit$run_info$thin + 1))
    expect_equal(unname(unlist(bplm.fit$Deviance_summary[1,])),unname(quantile(bplm.fit$Deviance_posterior,probs=c(0.025,0.5,0.975))))
    #rating curve and stage dependent parameters
    test_stage_dep_component(bplm.fit,'rating_curve')
    test_stage_dep_component(bplm.fit,'rating_curve_mean')
    test_stage_dep_component(bplm.fit,'sigma_eps')
    #Other information
    expect_equal(bplm.fit$formula,Q~W)
    expect_equal(bplm.fit$data,V316_river[order(V316_river$W),c('Q','W')])
})

test_that("the bplm object with known c with a maximum stage value is in tact", {
    expect_is(bplm.fit_known_c,"bplm")
    #latent parameters
    test_stage_indep_param(bplm.fit_known_c,'a')
    test_stage_indep_param(bplm.fit_known_c,'b')
    #hyperparameters
    expect_true(is.null(bplm.fit_known_c[['c_posterior']]))
    expect_false('c' %in% row.names(bplm.fit_known_c))
    test_stage_indep_param(bplm.fit,'sigma_eta')
    test_stage_indep_param(bplm.fit,'eta_1')
    test_stage_indep_param(bplm.fit,'eta_2')
    test_stage_indep_param(bplm.fit,'eta_3')
    test_stage_indep_param(bplm.fit,'eta_4')
    test_stage_indep_param(bplm.fit,'eta_5')
    test_stage_indep_param(bplm.fit,'eta_6')
    #Deviance
    expect_true(is.double(bplm.fit_known_c$Deviance_posterior))
    expect_equal(length(bplm.fit_known_c$Deviance_posterior),bplm.fit_known_c$run_info$num_chains*((bplm.fit_known_c$run_info$nr_iter-bplm.fit_known_c$run_info$burnin)/bplm.fit_known_c$run_info$thin + 1))
    expect_equal(unname(unlist(bplm.fit_known_c$Deviance_summary[1,])),unname(quantile(bplm.fit_known_c$Deviance_posterior,probs=c(0.025,0.5,0.975))))
    #rating curve and stage dependent parameters
    test_stage_dep_component(bplm.fit,'rating_curve')
    test_stage_dep_component(bplm.fit,'rating_curve_mean')
    test_stage_dep_component(bplm.fit,'sigma_eps')
    #check if maxmimum stage was in line with output
    expect_equal(max(bplm.fit_known_c$rating_curve$h),2.5)
    expect_true(max(diff(bplm.fit_known_c$rating_curve$h))<=(0.05+1e-9)) # added tolerance
})


test_that("bplm output remains unchanged", {
    skip_on_cran()
    skip_on_ci()
    expect_equal_to_reference(bplm.fit,file='../cached_results/bplm.fit.rds',update=T)
    expect_equal_to_reference(bplm.fit_known_c,file='../cached_results/bplm.fit_known_c.rds',update=T)
})

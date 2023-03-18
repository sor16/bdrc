context('plm')

test_that("plm can handle different inputs", {
    expect_error(plm(Q~W,c(1,2,3)))
    expect_error(plm('Q~W',krokfors))
    expect_error(plm(V~W,krokfors))
    expect_error(plm(Q~W+X,krokfors))
    expect_error(plm(Q~W,krokfors,c_param=min(krokfors$W)+0.5)) # c_param higher than lowest stage measurements
    expect_error(plm(Q~W,krokfors,c_param=1L)) # c_param not double
    expect_error(plm(Q~W,krokfors,h_max=max(krokfors$W)-0.5)) #h_max lower than highest stage measurement
    skip_on_cran()
    krokfors_new_names <- krokfors
    names(krokfors_new_names) <- c('t1','t2')
    set.seed(1)
    plm.fit_new_names <- plm(t2~t1,krokfors_new_names,num_cores=2)
    expect_equal(plm.fit_new_names$rating_curve,plm.fit$rating_curve)
})

test_that("the plm object with unknown c is in tact", {
    expect_is(plm.fit,"plm")
    #latent parameters
    test_stage_indep_param(plm.fit,'a')
    test_stage_indep_param(plm.fit,'b')
    #hyperparameters
    test_stage_indep_param(plm.fit,'c')
    test_stage_indep_param(plm.fit,'sigma_eta')
    test_stage_indep_param(plm.fit,'eta_1')
    test_stage_indep_param(plm.fit,'eta_2')
    test_stage_indep_param(plm.fit,'eta_3')
    test_stage_indep_param(plm.fit,'eta_4')
    test_stage_indep_param(plm.fit,'eta_5')
    test_stage_indep_param(plm.fit,'eta_6')
    #Deviance
    expect_true(is.double(plm.fit$Deviance_posterior))
    expect_equal(length(plm.fit$Deviance_posterior),plm.fit$run_info$num_chains*((plm.fit$run_info$nr_iter-plm.fit$run_info$burnin)/plm.fit$run_info$thin + 1))
    expect_equal(unname(unlist(plm.fit$Deviance_summary[1,])),unname(quantile(plm.fit$Deviance_posterior,probs=c(0.025,0.5,0.975))))
    #rating curve and stage dependent parameters
    test_stage_dep_component(plm.fit,'rating_curve')
    test_stage_dep_component(plm.fit,'rating_curve_mean')
    test_stage_dep_component(plm.fit,'sigma_eps')
    #Other information
    expect_equal(plm.fit$formula,Q~W)
    expect_equal(plm.fit$data,krokfors[order(krokfors$W),c('Q','W')])
})

test_that("the plm object with known c with a maximum stage value is in tact", {
    skip_on_cran()
    set.seed(1)
    plm.fit_known_c <- plm(Q~W,krokfors,c_param=known_c,h_max=h_extrap,num_cores=2)
    expect_is(plm.fit_known_c,"plm")
    #latent parameters
    test_stage_indep_param(plm.fit_known_c,'a')
    test_stage_indep_param(plm.fit_known_c,'b')
    #hyperparameters
    expect_true(is.null(plm.fit_known_c[['c_posterior']]))
    expect_false('c' %in% row.names(plm.fit_known_c))
    test_stage_indep_param(plm.fit,'sigma_eta')
    test_stage_indep_param(plm.fit,'eta_1')
    test_stage_indep_param(plm.fit,'eta_2')
    test_stage_indep_param(plm.fit,'eta_3')
    test_stage_indep_param(plm.fit,'eta_4')
    test_stage_indep_param(plm.fit,'eta_5')
    test_stage_indep_param(plm.fit,'eta_6')
    #Deviance
    expect_true(is.double(plm.fit_known_c$Deviance_posterior))
    expect_equal(length(plm.fit_known_c$Deviance_posterior),plm.fit_known_c$run_info$num_chains*((plm.fit_known_c$run_info$nr_iter-plm.fit_known_c$run_info$burnin)/plm.fit_known_c$run_info$thin + 1))
    expect_equal(unname(unlist(plm.fit_known_c$Deviance_summary[1,])),unname(quantile(plm.fit_known_c$Deviance_posterior,probs=c(0.025,0.5,0.975))))
    #rating curve and stage dependent parameters
    test_stage_dep_component(plm.fit,'rating_curve')
    test_stage_dep_component(plm.fit,'rating_curve_mean')
    test_stage_dep_component(plm.fit,'sigma_eps')
    #check if maxmimum stage was in line with output
    expect_equal(max(plm.fit_known_c$rating_curve$h),h_extrap)
    expect_true(max(diff(plm.fit_known_c$rating_curve$h))<=(0.05+1e-9)) # added tolerance
})


# test_that("plm output remains unchanged", {
#     skip_on_cran()
#     skip_on_ci()
#     skip_on_covr()
#     expect_equal_to_reference(plm.fit,file='../cached_results/plm.fit.rds',update=TRUE)
# })

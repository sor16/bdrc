context('gplm0')

test_that("gplm0 can handle different inputs", {
    expect_error(gplm0(Q~W,c(1,2,3)))
    expect_error(gplm0('Q~W',krokfors))
    expect_error(gplm0(V~W,krokfors))
    expect_error(gplm0(Q~W+X,krokfors))
    expect_error(gplm0(Q~W,krokfors,c_param=min(krokfors$W)+0.5)) # c_param higher than lowest stage measurements
    expect_error(gplm0(Q~W,krokfors,c_param=1L)) # c_param not double
    expect_error(gplm0(Q~W,krokfors,h_max=max(krokfors$W)-0.5)) #h_max lower than highest stage measurement
    skip_on_cran()
    krokfors_new_names <- krokfors
    names(krokfors_new_names) <- c('t1','t2')
    set.seed(1)
    gplm0.fit_new_names <- gplm0(t2~t1,krokfors_new_names,num_cores=2)
    expect_equal(gplm0.fit_new_names$rating_curve,gplm0.fit$rating_curve)
})

test_that("the gplm0 object with unknown c is in tact", {
    expect_is(gplm0.fit,"gplm0")
    #latent parameters
    test_stage_indep_param(gplm0.fit,'a')
    test_stage_indep_param(gplm0.fit,'b')
    #hyperparameters
    test_stage_indep_param(gplm0.fit,'c')
    test_stage_indep_param(gplm0.fit,'sigma_eps')
    test_stage_indep_param(gplm0.fit,'sigma_beta')
    test_stage_indep_param(gplm0.fit,'phi_beta')
    #Deviance
    expect_true(is.double(gplm0.fit$Deviance_posterior))
    expect_equal(length(gplm0.fit$Deviance_posterior),gplm0.fit$run_info$num_chains*((gplm0.fit$run_info$nr_iter-gplm0.fit$run_info$burnin)/gplm0.fit$run_info$thin + 1))
    expect_equal(unname(unlist(gplm0.fit$Deviance_summary[1,])),unname(quantile(gplm0.fit$Deviance_posterior,probs=c(0.025,0.5,0.975))))
    #rating curve and stage dependent parameters
    test_stage_dep_component(gplm0.fit,'rating_curve')
    test_stage_dep_component(gplm0.fit,'rating_curve_mean')
    test_stage_dep_component(gplm0.fit,'beta')
    test_stage_dep_component(gplm0.fit,'f')
    #Other information
    expect_equal(gplm0.fit$formula,Q~W)
    expect_equal(gplm0.fit$data,krokfors[order(krokfors$W),c('Q','W')])
})

test_that("the gplm0 object with known c with a maximum stage value is in tact", {
    skip_on_cran()
    set.seed(1)
    gplm0.fit_known_c <- gplm0(Q~W,krokfors,c_param=known_c,h_max=h_extrap,num_cores=2)
    expect_is(gplm0.fit_known_c,"gplm0")
    #latent parameters
    test_stage_indep_param(gplm0.fit_known_c,'a')
    test_stage_indep_param(gplm0.fit_known_c,'b')
    #hyperparameters
    expect_true(is.null(gplm0.fit_known_c[['c_posterior']]))
    expect_false('c' %in% row.names(gplm0.fit_known_c))
    test_stage_indep_param(gplm0.fit_known_c,'sigma_eps')
    test_stage_indep_param(gplm0.fit,'sigma_beta')
    test_stage_indep_param(gplm0.fit,'phi_beta')
    #Deviance
    expect_true(is.double(gplm0.fit_known_c$Deviance_posterior))
    expect_equal(length(gplm0.fit_known_c$Deviance_posterior),gplm0.fit_known_c$run_info$num_chains*((gplm0.fit_known_c$run_info$nr_iter-gplm0.fit_known_c$run_info$burnin)/gplm0.fit_known_c$run_info$thin + 1))
    expect_equal(unname(unlist(gplm0.fit_known_c$Deviance_summary[1,])),unname(quantile(gplm0.fit_known_c$Deviance_posterior,probs=c(0.025,0.5,0.975))))
    #rating curve
    test_stage_dep_component(gplm0.fit_known_c,'rating_curve')
    test_stage_dep_component(gplm0.fit_known_c,'rating_curve_mean')
    test_stage_dep_component(gplm0.fit,'beta')
    test_stage_dep_component(gplm0.fit,'f')
    #check if maxmimum stage was in line with output
    expect_equal(max(gplm0.fit_known_c$rating_curve$h),h_extrap)
    expect_true(max(diff(gplm0.fit_known_c$rating_curve$h))<=(0.05+1e-9)) # added tolerance
})


# test_that("gplm0 output remains unchanged", {
#     skip_on_cran()
#     skip_on_ci()
#     skip_on_covr()
#     expect_equal_to_reference(gplm0.fit,file='../cached_results/gplm0.fit.rds',update=TRUE)
# })

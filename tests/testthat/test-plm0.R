context('plm0')

test_that("plm0 can handle different inputs", {
    expect_error(plm0(Q~W,c(1,2,3)))
    expect_error(plm0('Q~W',krokfors))
    expect_error(plm0(V~W,krokfors))
    expect_error(plm0(Q~W+X,krokfors))
    expect_error(plm0(Q~W,krokfors,c_param=min(krokfors$W)+0.5)) # c_param higher than lowest stage measurements
    expect_error(plm0(Q~W,krokfors,c_param=1L)) # c_param not double
    expect_error(plm0(Q~W,krokfors,h_max=max(krokfors$W)-0.5)) #h_max lower than highest stage measurement
    skip_on_cran()
    krokfors_new_names <- krokfors
    names(krokfors_new_names) <- c('t1','t2')
    set.seed(1)
    plm0.fit_new_names <- plm0(t2~t1,krokfors_new_names,num_cores=2)
    expect_equal(plm0.fit_new_names$rating_curve,plm0.fit$rating_curve)
})

test_that("the plm0 object with unknown c is in tact", {
    expect_is(plm0.fit,"plm0")
    #latent parameters
    test_stage_indep_param(plm0.fit,'a')
    test_stage_indep_param(plm0.fit,'b')
    #hyperparameters
    test_stage_indep_param(plm0.fit,'c')
    test_stage_indep_param(plm0.fit,'sigma_eps')
    #Deviance
    expect_true(is.double(plm0.fit$Deviance_posterior))
    expect_equal(length(plm0.fit$Deviance_posterior),plm0.fit$run_info$num_chains*((plm0.fit$run_info$nr_iter-plm0.fit$run_info$burnin)/plm0.fit$run_info$thin + 1))
    expect_equal(unname(unlist(plm0.fit$Deviance_summary[1,])),unname(quantile(plm0.fit$Deviance_posterior,probs=c(0.025,0.5,0.975))))
    #rating curve
    test_stage_dep_component(plm0.fit,'rating_curve')
    test_stage_dep_component(plm0.fit,'rating_curve_mean')
    #Other information
    expect_equal(plm0.fit$formula,Q~W)
    expect_equal(plm0.fit$data,krokfors[order(krokfors$W),c('Q','W')])
})

test_that("the plm0 object with known c with a maximum stage value is in tact", {
    skip_on_cran()
    set.seed(1)
    plm0.fit_known_c <- plm0(Q~W,krokfors,c_param=known_c,h_max=h_extrap,num_cores=2)
    expect_is(plm0.fit_known_c,"plm0")
    #latent parameters
    test_stage_indep_param(plm0.fit_known_c,'a')
    test_stage_indep_param(plm0.fit_known_c,'b')
    #hyperparameters
    expect_true(is.null(plm0.fit_known_c[['c_posterior']]))
    expect_false('c' %in% row.names(plm0.fit_known_c))
    test_stage_indep_param(plm0.fit_known_c,'sigma_eps')
    #Deviance
    expect_true(is.double(plm0.fit_known_c$Deviance_posterior))
    expect_equal(length(plm0.fit_known_c$Deviance_posterior),plm0.fit_known_c$run_info$num_chains*((plm0.fit_known_c$run_info$nr_iter-plm0.fit_known_c$run_info$burnin)/plm0.fit_known_c$run_info$thin + 1))
    expect_equal(unname(unlist(plm0.fit_known_c$Deviance_summary[1,])),unname(quantile(plm0.fit_known_c$Deviance_posterior,probs=c(0.025,0.5,0.975))))
    #rating curve
    test_stage_dep_component(plm0.fit_known_c,'rating_curve')
    test_stage_dep_component(plm0.fit_known_c,'rating_curve_mean')
    #check if maxmimum stage was in line with output
    expect_equal(max(plm0.fit_known_c$rating_curve$h),h_extrap)
    expect_true(max(diff(plm0.fit_known_c$rating_curve$h))<=(0.05+1e-9)) # added tolerance
})


test_that("plm0 output remains unchanged", {
    skip_on_cran()
    skip_on_ci()
    skip_on_covr()
    expect_equal_to_reference(plm0.fit,file='../cached_results/plm0.fit.rds',update=TRUE)
})

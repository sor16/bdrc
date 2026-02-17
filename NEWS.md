# bdrc 2.0.1

* Fixed phantom observation variance bug in gplm0 and gplm: the phantom observation that centers the GP used zero variance, acting as an infinitely precise constraint that collapsed the posterior of b. Fixed by using the prior variance sig_b^2.
* Fixed eta posterior scaling bug in plm and gplm: the standardized eta coefficients were not scaled by sigma_eta before applying the basis transformation P, resulting in incorrect B-spline coefficient posteriors.
* Fixed sigma_beta back-transform bug in gplm0: stored values were sqrt(sigma_beta) instead of sigma_beta due to an incorrect sqrt(exp()) transform. This affected posterior summaries and plots for the sigma_beta parameter in the gplm0 model.
* Updated reference in DESCRIPTION from arXiv preprint to the published Environmetrics paper (doi:10.1002/env.2711).

# bdrc 2.0.0

* Integrated C++ via Rcpp and RcppArmadillo packages for significant performance improvements.
* Multiple functions rewritten in C++ to speed up the MCMC sampling algorithm and various other tasks.
* The "Deviance" posterior output has been transformed and renamed "Posterior log-likelihood". The Deviance was previously calculated as -2 times the Posterior log-likelihood.
* The plot(tournament_obj, type = "Deviance") figure is now created by evaluating plot(tournament_obj, type = "boxplot").
* The log-likelihood of the models is now computed with the log-transformed discharge observations (normally distributed), rather than with discharge on the real scale (log-normally distributed).
* Pointwise WAIC values (WAIC_i) stored to the model objects.
* Implemented standard error computations for WAIC and Delta_WAIC estimates.
* Applied log-sum-exp trick in WAIC and Bayes factor calculations for numerical stability.
* Added DOI links to references.
* Revised summary() output for tournament objects.
* Updated package vignettes to reflect recent changes and improvements.

# bdrc 1.1.0

* The tournament function now supports information criteria (WAIC & DIC) as the model selection criteria, WAIC being the new default.
* Six new real-world datasets from Iceland and Sweden containing paired observations of stage and discharge have be included in the package.
* The package can now also be used in a user-friendly interactive rating curve builder Shiny application [https://bdrc.shinyapps.io/bdrc/].


# bdrc 1.0.0

* Release of the first version of the package.

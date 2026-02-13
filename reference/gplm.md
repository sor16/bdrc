# Generalized power-law model with variance that varies with stage.

gplm is used to fit a discharge rating curve for paired measurements of
stage and discharge using a generalized power-law model with variance
that varies with stage as described in Hrafnkelsson et al. (2022). See
"Details" for a more elaborate description of the model.

## Usage

``` r
gplm(
  formula,
  data,
  c_param = NULL,
  h_max = NULL,
  parallel = TRUE,
  num_cores = NULL,
  forcepoint = rep(FALSE, nrow(data)),
  verbose = TRUE
)
```

## Arguments

- formula:

  An object of class "formula", with discharge column name as response
  and stage column name as a covariate, i.e. of the form `y`~`x` where
  `y` is discharge in m\\^3/\\s and `x` is stage in m (it is very
  important that the data is in the correct units).

- data:

  A data.frame containing the variables specified in formula.

- c_param:

  The largest stage value for which there is zero discharge. If NULL, it
  is treated as unknown in the model and inferred from the data.

- h_max:

  The maximum stage to which the rating curve should extrapolate to. If
  NULL, the maximum stage value in the data is selected as an upper
  bound.

- parallel:

  A logical value indicating whether to run the MCMC in parallel or not.
  Defaults to TRUE.

- num_cores:

  An integer between 1 and 4 (number of MCMC chains) indicating how many
  cores to use. Only used if parallel=TRUE. If NULL, the number of cores
  available on the device is detected automatically.

- forcepoint:

  A logical vector of the same length as the number of rows in data. If
  an element at index \\i\\ is TRUE it indicates that the rating curve
  should be forced through the \\i\\-th measurement. Use with care, as
  this will strongly influence the resulting rating curve.

- verbose:

  A logical value indicating whether to print progress and diagnostic
  information. If \`TRUE\`, the function will print messages as it runs.
  If \`FALSE\`, the function will run silently. Default is \`TRUE\`.

## Value

gplm returns an object of class "gplm". An object of class "gplm" is a
list containing the following components:

- `rating_curve`:

  A data frame with 2.5%, 50% and 97.5% percentiles of the posterior
  predictive distribution of the rating curve.

- `rating_curve_mean`:

  A data frame with 2.5%, 50% and 97.5% percentiles of the posterior
  distribution of the mean of the rating curve.

- `param_summary`:

  A data frame with 2.5%, 50% and 97.5% percentiles of the posterior
  distribution of latent- and hyperparameters. Additionally contains
  columns with r_hat and the effective number of samples for each
  parameter as defined in Gelman et al. (2013).

- `f_summary`:

  A data frame with 2.5%, 50% and 97.5% percentiles of the posterior
  distribution of \\f(h)\\.

- `beta_summary`:

  A data frame with 2.5%, 50% and 97.5% percentiles of the posterior
  distribution of \\\beta(h)\\.

- `sigma_eps_summary`:

  A data frame with 2.5%, 50% and 97.5% percentiles of the posterior
  distribution of \\\sigma\_\varepsilon(h)\\.

- `posterior_log_likelihood_summary`:

  A data frame with 2.5%, 50% and 97.5% percentiles of the posterior
  log-likelihood values.

- `rating_curve_posterior`:

  A matrix containing the full thinned posterior samples of the
  posterior predictive distribution of the rating curve excluding
  burn-in samples.

- `rating_curve_mean_posterior`:

  A matrix containing the full thinned posterior samples of the
  posterior distribution of the mean of the rating curve excluding
  burn-in samples.

- `a_posterior`:

  A numeric vector containing the full thinned posterior samples of the
  posterior distribution of \\a\\ excluding burn-in samples.

- `b_posterior`:

  A numeric vector containing the full thinned posterior samples of the
  posterior distribution of \\b\\ excluding burn-in samples.

- `c_posterior`:

  A numeric vector containing the full thinned posterior samples of the
  posterior distribution of \\c\\ excluding burn-in samples.

- `sigma_beta_posterior`:

  A numeric vector containing the full thinned posterior samples of the
  posterior distribution of \\\sigma\_\beta\\ excluding burn-in samples.

- `phi_beta_posterior`:

  A numeric vector containing the full thinned posterior samples of the
  posterior distribution of \\\phi\_\beta\\ excluding burn-in samples.

- `sigma_eta_posterior`:

  A numeric vector containing the full thinned posterior samples of the
  posterior distribution of \\\sigma\_\eta\\ excluding burn-in samples.

- `eta_1_posterior`:

  A numeric vector containing the full thinned posterior samples of the
  posterior distribution of \\\eta_1\\ excluding burn-in samples.

- `eta_2_posterior`:

  A numeric vector containing the full thinned posterior samples of the
  posterior distribution of \\\eta_2\\ excluding burn-in samples.

- `eta_3_posterior`:

  A numeric vector containing the full thinned posterior samples of the
  posterior distribution of \\\eta_3\\ excluding burn-in samples.

- `eta_4_posterior`:

  A numeric vector containing the full thinned posterior samples of the
  posterior distribution of \\\eta_4\\ excluding burn-in samples.

- `eta_5_posterior`:

  A numeric vector containing the full thinned posterior samples of the
  posterior distribution of \\\eta_5\\ excluding burn-in samples.

- `eta_6_posterior`:

  A numeric vector containing the full thinned posterior samples of the
  posterior distribution of \\\eta_6\\ excluding burn-in samples.

- `f_posterior`:

  A numeric vector containing the full thinned posterior samples of the
  posterior distribution of \\f(h)\\ excluding burn-in samples.

- `beta_posterior`:

  A numeric vector containing the full thinned posterior samples of the
  posterior distribution of \\\beta(h)\\ excluding burn-in samples.

- `sigma_eps_posterior`:

  A numeric vector containing the full thinned posterior samples of the
  posterior distribution of \\\sigma\_\varepsilon(h)\\ excluding burn-in
  samples.

- `posterior_log_likelihood`:

  A numeric vector containing the full thinned posterior log-likelihood
  values, excluding burn-in samples.

- `D_hat`:

  A statistic defined as -2 times the log-likelihood evaluated at the
  median value of the parameters.

- `effective_num_param_DIC`:

  The effective number of parameters, which is calculated as
  median(-2\*posterior_log_likelihood) minus D_hat.

- `DIC`:

  The Deviance Information Criterion for the model, calculated as D_hat
  plus 2\*effective_num_parameters_DIC.

- `lppd`:

  The log pointwise predictive density of the observed data under the
  model.

- `WAIC`:

  The Widely Applicable Information Criterion for the model, defined as
  -2\*( lppd - effective_num_param_WAIC ).

- `WAIC_i`:

  The pointwise WAIC values, where WAIC := sum(WAIC_i).

- `effective_num_param_WAIC`:

  The effective number of parameters, which is calculated by summing up
  the posterior variance of the log predictive density for each data
  point.

- `autocorrelation`:

  A data frame with the autocorrelation of each parameter for different
  lags.

- `acceptance_rate`:

  The proportion of accepted samples in the thinned MCMC chain
  (excluding burn-in).

- `formula`:

  An object of type "formula" provided by the user.

- `data`:

  The data provided by the user, ordered by stage.

- `run_info`:

  The information about the input arguments and the specific parameters
  used in the MCMC chain.

## Details

The generalized power-law model is of the form \$\$Q=a(h-c)^{f(h)}\$\$
where \\Q\\ is discharge, \\h\\ is stage, \\a\\ and \\c\\ are unknown
constants and \\f\\ is a function of \\h\\, referred to as the
generalized power-law exponent.  
  
The generalized power-law model is here inferred by using a Bayesian
hierarchical model. The function \\f\\ is modeled at the latent level as
a fixed constant \\b\\ plus a continuous stochastic process,
\\\beta(h)\\, which is assumed to be twice differentiable. The model is
on a logarithmic scale \$\$\log(Q_i) = \log(a) + (b + \beta(h_i))
\log(h_i - c) + \varepsilon_i, i = 1,...,n\$\$ where \\\varepsilon_i\\
follows a normal distribution with mean zero and variance
\\\sigma\_\varepsilon(h_i)^2\\ that varies with stage. The stochastic
process \\\beta(h)\\ is assumed a priori to be a Gaussian process
governed by a Matern covariance function with smoothness parameter \\\nu
= 2.5\\. The error variance, \\\sigma\_\varepsilon^2(h)\\, of the
log-discharge data is modeled as an exponential of a B-spline curve,
that is, a linear combination of six B-spline basis functions that are
defined over the range of the stage observations. An efficient posterior
simulation is achieved by sampling from the joint posterior density of
the hyperparameters of the model, and then sampling from the density of
the latent parameters conditional on the hyperparameters.  
  
Bayesian inference is based on the posterior density and summary
statistics such as the posterior mean and 95% posterior intervals are
based on the posterior density. Analytical formulas for these summary
statistics are intractable in most cases and thus they are computed by
generating samples from the posterior density using a Markov chain Monte
Carlo simulation.

## References

Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., and
Rubin, D. B. (2013). Bayesian Data Analysis, Third Edition. Chapman &
Hall/CRC Texts in Statistical Science. Taylor & Francis. doi:
https://doi.org/10.1201/b16018

Hrafnkelsson, B., Sigurdarson, H., Rögnvaldsson, S., Jansson, A. Ö.,
Vias, R. D., and Gardarsson, S. M. (2022). Generalization of the
power-law rating curve using hydrodynamic theory and Bayesian
hierarchical modeling, Environmetrics, 33(2):e2711. doi:
https://doi.org/10.1002/env.2711

Spiegelhalter, D., Best, N., Carlin, B., Van Der Linde, A. (2002).
Bayesian measures of model complexity and fit. Journal of the Royal
Statistical Society: Series B (Statistical Methodology) 64(4), 583–639.
doi: https://doi.org/10.1111/1467-9868.00353

Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation
and widely applicable information criterion in singular learning theory.
J. Mach. Learn. Res. 11, 3571–3594.

## See also

[`summary.gplm`](https://sor16.github.io/bdrc/reference/summary.plm0.md)
for summaries,
[`predict.gplm`](https://sor16.github.io/bdrc/reference/predict.plm0.md)
for prediction and
[`plot.gplm`](https://sor16.github.io/bdrc/reference/plot.plm0.md) for
plots.
[`spread_draws`](https://sor16.github.io/bdrc/reference/spread_draws.md)
and
[`gather_draws`](https://sor16.github.io/bdrc/reference/gather_draws.md)
are also useful to aid further visualization of the full posterior
distributions.

## Examples

``` r
# \donttest{
data(norn)
set.seed(1)
gplm.fit <- gplm(formula=Q~W,data=norn,num_cores=2)
#> Progress:
#> Initializing Metropolis MCMC algorithm...
#> Multiprocess sampling (4 chains in 2 jobs) ...
#> 
#> MCMC sampling completed!
#> 
#> Diagnostics:
#> Acceptance rate: 22.81%.
#> ✔ All chains have mixed well (Rhat < 1.1).
#> ✔ Effective sample sizes sufficient (eff_n_samples > 400).
summary(gplm.fit)
#> 
#> Formula: 
#>  Q ~ W
#> Latent parameters:
#>   lower-2.5% median-50% upper-97.5%
#> a 13.40      24.55      32.51      
#> b  1.83       1.84       1.84      
#> 
#> Hyperparameters:
#>            lower-2.5% median-50% upper-97.5%
#> c          396.421    396.669    396.812    
#> sigma_beta   0.165      0.338      0.806    
#> phi_beta     1.343      2.987      6.894    
#> sigma_eta    0.048      0.436      0.897    
#> eta_1       -6.128     -5.395     -4.399    
#> eta_2       -7.692     -5.974     -4.521    
#> eta_3       -9.093     -7.048     -5.072    
#> eta_4      -10.429     -8.051     -5.135    
#> eta_5      -11.473     -8.717     -5.092    
#> eta_6      -12.436     -9.279     -5.095    
#> 
#> WAIC: -141.1422
# }
```

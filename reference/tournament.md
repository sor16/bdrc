# Tournament - Model comparison

tournament compares four rating curve models of different complexities
and determines the model that provides the best fit of the data at hand.

## Usage

``` r
tournament(
  formula = NULL,
  data = NULL,
  model_list = NULL,
  method = "WAIC",
  winning_criteria = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- formula:

  An object of class "formula", with discharge column name as response
  and stage column name as a covariate.

- data:

  A data.frame containing the variables specified in formula.

- model_list:

  A list of exactly four model objects of types "plm0","plm","gplm0" and
  "gplm" to be used in the tournament. Note that all of the model
  objects are required to be run with the same data and same c_param.

- method:

  A string specifying the method used to estimate the predictive
  performance of the models. The allowed methods are "WAIC", "DIC" and
  "PMP".

- winning_criteria:

  Specifies the criteria for model selection. For "WAIC", it can be a
  numeric value or a string expression. For "DIC", it must be a numeric
  value. For "PMP", it must be a numeric value between 0 and 1. See
  Details section.

- verbose:

  A logical value indicating whether to print progress and diagnostic
  information. If \`TRUE\`, the function will print messages as it runs.
  If \`FALSE\`, the function will run silently. Default is \`TRUE\`.

- ...:

  Optional arguments passed to the model functions.

## Value

An object of type "tournament" with the following elements:

- `contestants`:

  The model objects of types "plm0", "plm", "gplm0" and "gplm" being
  compared.

- `winner`:

  The model object of the tournament winner.

- `info`:

  The specifics about the tournament; the overall winner; the method
  used; and the winning criteria.

- `summary`:

  A data frame with information on results of the different comparisons
  in the power-law tournament. The contents of this data frame depend on
  the method used:

  - For all methods:

    - round: The tournament round

    - comparison: The comparison number

    - complexity: Indicates whether a model is the "more" or "less"
      complex model in a comparison

    - model: The model type

    - winner: Logical value indicating if the model was selected in the
      corresponding comparison

  - Additional columns for method "WAIC":

    - lppd: Log pointwise predictive density

    - eff_num_param: Effective number of parameters (WAIC)

    - WAIC: Widely Applicable Information Criterion

    - SE_WAIC: Standard error of WAIC

    - Delta_WAIC: Difference in WAIC

    - SE_Delta_WAIC: Standard error of the difference in WAIC

  - Additional columns for method "DIC":

    - D_hat: Minus two times the log-likelihood evaluated at the median
      of the posterior samples

    - eff_num_param: Effective number of parameters (DIC)

    - DIC: Deviance Information Criterion

    - Delta_DIC: Difference in DIC

  - Additional columns for method "PMP":

    - log_marg_lik: Logarithm of the marginal likelihood estimated,
      computed with the harmonic-mean estimator

    - PMP: Posterior model probability computed with Bayes factor

## Details

Tournament is a model comparison method that uses WAIC (default method)
to estimate the expected prediction error of the four models and select
the most appropriate model given the data. The first round of model
comparisons sets up model types, "gplm" vs. "gplm0" and "plm" vs.
"plm0". The two comparisons are conducted such that if the WAIC of the
more complex model ("gplm" and "plm", respectively) is smaller than the
WAIC of the simpler models ("gplm0" and "plm0", respectively) by an
input argument called the `winning_criteria` (default value = 2), then
it is chosen as the more appropriate model. If not, the simpler model is
chosen. The more appropriate models move on to the second round and are
compared in the same way. The winner of the second round is chosen as
the overall tournament winner and deemed the most appropriate model
given the data.

The default method "WAIC", or the Widely Applicable Information
Criterion (see Watanabe (2010)), is used to estimate the predictive
performance of the models. This method is a fully Bayesian method that
uses the full set of posterior draws to estimate of the expected log
pointwise predictive density.

Method "DIC", or Deviance Information Criterion (see Spiegelhalter
(2002)), is similar to the "WAIC" but instead of using the full set of
posterior draws to compute the estimate of the expected log pointwise
predictive density, it uses a point estimate of the posterior
distribution.

Method "PMP" uses the posterior model probabilities, calculated with
Bayes factor (see Jeffreys (1961) and Kass and Raftery (1995)), to
compare the models, where all the models are assumed a priori to be
equally likely. This method is not chosen as the default method because
the Bayes factor calculations can be quite unstable.

When method "WAIC" is used, the `winning_criteria` can be either a
numeric value or a string expression. If numeric, it sets the threshold
which the more complex model must exceed to be declared the more
appropriate model. If a string, it must be a valid R expression using
Delta_WAIC and/or SE_Delta_WAIC (e.g., "Delta_WAIC \> 2 & Delta_WAIC -
SE_Delta_WAIC \> 0"). For method "DIC", `winning_criteria` must be a
numeric value. For method "PMP", the winning criteria should be a
numeric value between 0 and 1 (default value = 0.75). This sets the
threshold value for which the posterior probability of the more complex
model, given the data, in each model comparison must exceed to be
declared the more appropriate model. In all cases, the default values
are selected to give the less complex models a slight advantage, which
should give more or less consistent results when applying the tournament
to real world data.

## References

Hrafnkelsson, B., Sigurdarson, H., Rögnvaldsson, S., Jansson, A. Ö.,
Vias, R. D., and Gardarsson, S. M. (2022). Generalization of the
power-law rating curve using hydrodynamic theory and Bayesian
hierarchical modeling, Environmetrics, 33(2):e2711. doi:
https://doi.org/10.1002/env.2711

Jeffreys, H. (1961). Theory of Probability, Third Edition. Oxford
University Press.

Kass, R., and A. Raftery, A. (1995). Bayes Factors. Journal of the
American Statistical Association, 90, 773-795. doi:
https://doi.org/10.1080/01621459.1995.10476572

Spiegelhalter, D., Best, N., Carlin, B., Van Der Linde, A. (2002).
Bayesian measures of model complexity and fit. Journal of the Royal
Statistical Society: Series B (Statistical Methodology) 64(4), 583–639.
doi: https://doi.org/10.1111/1467-9868.00353

Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation
and widely applicable information criterion in singular learning theory.
J. Mach. Learn. Res. 11, 3571–3594.

## See also

[`plm0`](https://sor16.github.io/bdrc/reference/plm0.md)
[`plm`](https://sor16.github.io/bdrc/reference/plm.md),
[`gplm0`](https://sor16.github.io/bdrc/reference/gplm0.md),[`gplm`](https://sor16.github.io/bdrc/reference/gplm.md)
[`summary.tournament`](https://sor16.github.io/bdrc/reference/summary.tournament.md)
and
[`plot.tournament`](https://sor16.github.io/bdrc/reference/plot.tournament.md)

## Examples

``` r
# \donttest{
data(krokfors)
set.seed(1)
t_obj <- tournament(formula = Q ~ W, data = krokfors, num_cores = 2)
#> Running tournament  [                                                ] 0%
#> 
#> Progress:
#> Initializing Metropolis MCMC algorithm...
#> Multiprocess sampling (4 chains in 2 jobs) ...
#> 
#> MCMC sampling completed!
#> 
#> Diagnostics:
#> Acceptance rate: 25.33%.
#> ✔ All chains have mixed well (Rhat < 1.1).
#> ✔ Effective sample sizes sufficient (eff_n_samples > 400).
#> 
#>  ✔  gplm finished   [============                                    ] 25%
#> 
#> Progress:
#> Initializing Metropolis MCMC algorithm...
#> Multiprocess sampling (4 chains in 2 jobs) ...
#> 
#> MCMC sampling completed!
#> 
#> Diagnostics:
#> Acceptance rate: 31.14%.
#> ✔ All chains have mixed well (Rhat < 1.1).
#> ✔ Effective sample sizes sufficient (eff_n_samples > 400).
#> 
#>  ✔  gplm0 finished  [========================                        ] 50%
#> 
#> Progress:
#> Initializing Metropolis MCMC algorithm...
#> Multiprocess sampling (4 chains in 2 jobs) ...
#> 
#> MCMC sampling completed!
#> 
#> Diagnostics:
#> Acceptance rate: 25.66%.
#> ✔ All chains have mixed well (Rhat < 1.1).
#> ✔ Effective sample sizes sufficient (eff_n_samples > 400).
#> 
#>  ✔  plm finished    [====================================            ] 75%
#> 
#> Progress:
#> Initializing Metropolis MCMC algorithm...
#> Multiprocess sampling (4 chains in 2 jobs) ...
#> 
#> MCMC sampling completed!
#> 
#> Diagnostics:
#> Acceptance rate: 36.04%.
#> ✔ All chains have mixed well (Rhat < 1.1).
#> ✔ Effective sample sizes sufficient (eff_n_samples > 400).
#> 
#>  ✔  plm0 finished   [================================================] 100%
t_obj
#> Tournament winner: gplm0
summary(t_obj)
#> 
#> === Tournament Model Comparison Summary ===
#> 
#> Method: WAIC 
#> Winning Criteria: Delta_WAIC > 2 
#> Overall Winner: gplm0 
#> 
#> Comparison 1 Results:
#> -------------------------------------------------------------------------------------------------- 
#> complexity   model  winner lppd     eff_num_param  WAIC       SE_WAIC    Delta_WAIC SE_Delta_WAIC
#> more         gplm          20.7794  6.8706         -27.8176   11.8918    0.5570     0.2416      
#> less         gplm0  <---   20.3710  6.7406         -27.2606   12.0360                           
#> 
#> Comparison 2 Results:
#> -------------------------------------------------------------------------------------------------- 
#> complexity   model  winner lppd     eff_num_param  WAIC       SE_WAIC    Delta_WAIC SE_Delta_WAIC
#> more         plm           5.5842   4.2574         -2.6536    6.6635     -0.4066    0.1904      
#> less         plm0   <---   5.6284   4.0984         -3.0601    6.6931                            
#> 
#> Comparison 3 Results:
#> -------------------------------------------------------------------------------------------------- 
#> complexity   model  winner lppd     eff_num_param  WAIC       SE_WAIC    Delta_WAIC SE_Delta_WAIC
#> more         gplm0  <---   20.3710  6.7406         -27.2606   12.0360    24.2005    9.1834      
#> less         plm0          5.6284   4.0984         -3.0601    6.6931                            
#> 
#> === End of Summary ===

# Using different methods and winning criteria
t_obj_dic <- tournament(Q ~ W,
                        krokfors,
                        num_cores = 2,
                        method = "DIC",
                        winning_criteria = 3)
#> Running tournament  [                                                ] 0%
#> 
#> Progress:
#> Initializing Metropolis MCMC algorithm...
#> Multiprocess sampling (4 chains in 2 jobs) ...
#> 
#> MCMC sampling completed!
#> 
#> Diagnostics:
#> Acceptance rate: 24.81%.
#> ✔ All chains have mixed well (Rhat < 1.1).
#> ✔ Effective sample sizes sufficient (eff_n_samples > 400).
#> 
#>  ✔  gplm finished   [============                                    ] 25%
#> 
#> Progress:
#> Initializing Metropolis MCMC algorithm...
#> Multiprocess sampling (4 chains in 2 jobs) ...
#> 
#> MCMC sampling completed!
#> 
#> Diagnostics:
#> Acceptance rate: 31.42%.
#> ✔ All chains have mixed well (Rhat < 1.1).
#> ✔ Effective sample sizes sufficient (eff_n_samples > 400).
#> 
#>  ✔  gplm0 finished  [========================                        ] 50%
#> 
#> Progress:
#> Initializing Metropolis MCMC algorithm...
#> Multiprocess sampling (4 chains in 2 jobs) ...
#> 
#> MCMC sampling completed!
#> 
#> Diagnostics:
#> Acceptance rate: 26.60%.
#> ✔ All chains have mixed well (Rhat < 1.1).
#> ✔ Effective sample sizes sufficient (eff_n_samples > 400).
#> 
#>  ✔  plm finished    [====================================            ] 75%
#> 
#> Progress:
#> Initializing Metropolis MCMC algorithm...
#> Multiprocess sampling (4 chains in 2 jobs) ...
#> 
#> MCMC sampling completed!
#> 
#> Diagnostics:
#> Acceptance rate: 36.02%.
#> ✔ All chains have mixed well (Rhat < 1.1).
#> ✔ Effective sample sizes sufficient (eff_n_samples > 400).
#> 
#>  ✔  plm0 finished   [================================================] 100%
t_obj_pmp <- tournament(Q ~ W,
                        krokfors,
                        num_cores = 2,
                        method = "PMP",
                        winning_criteria = 0.8)
#> Running tournament  [                                                ] 0%
#> 
#> Progress:
#> Initializing Metropolis MCMC algorithm...
#> Multiprocess sampling (4 chains in 2 jobs) ...
#> 
#> MCMC sampling completed!
#> 
#> Diagnostics:
#> Acceptance rate: 25.10%.
#> ✔ All chains have mixed well (Rhat < 1.1).
#> ✔ Effective sample sizes sufficient (eff_n_samples > 400).
#> 
#>  ✔  gplm finished   [============                                    ] 25%
#> 
#> Progress:
#> Initializing Metropolis MCMC algorithm...
#> Multiprocess sampling (4 chains in 2 jobs) ...
#> 
#> MCMC sampling completed!
#> 
#> Diagnostics:
#> Acceptance rate: 31.05%.
#> ✔ All chains have mixed well (Rhat < 1.1).
#> ✔ Effective sample sizes sufficient (eff_n_samples > 400).
#> 
#>  ✔  gplm0 finished  [========================                        ] 50%
#> 
#> Progress:
#> Initializing Metropolis MCMC algorithm...
#> Multiprocess sampling (4 chains in 2 jobs) ...
#> 
#> MCMC sampling completed!
#> 
#> Diagnostics:
#> Acceptance rate: 25.47%.
#> ✔ All chains have mixed well (Rhat < 1.1).
#> ✔ Effective sample sizes sufficient (eff_n_samples > 400).
#> 
#>  ✔  plm finished    [====================================            ] 75%
#> 
#> Progress:
#> Initializing Metropolis MCMC algorithm...
#> Multiprocess sampling (4 chains in 2 jobs) ...
#> 
#> MCMC sampling completed!
#> 
#> Diagnostics:
#> Acceptance rate: 35.76%.
#> ✔ All chains have mixed well (Rhat < 1.1).
#> ✔ Effective sample sizes sufficient (eff_n_samples > 400).
#> 
#>  ✔  plm0 finished   [================================================] 100%
#> ⚠ Warning: The Harmonic Mean Estimator (HME) is used to estimate the Bayes Factor for the posterior model probability (PMP), which is known to be unstable and potentially unreliable. We recommend using method "WAIC" (Widely Applicable Information Criterion) for model comparison instead.
t_obj_waic_expr <- tournament(Q ~ W,
                              krokfors,
                              num_cores = 2,
                              winning_criteria = "Delta_WAIC > 2 & Delta_WAIC - SE_Delta_WAIC > 0")
#> Running tournament  [                                                ] 0%
#> 
#> Progress:
#> Initializing Metropolis MCMC algorithm...
#> Multiprocess sampling (4 chains in 2 jobs) ...
#> 
#> MCMC sampling completed!
#> 
#> Diagnostics:
#> Acceptance rate: 25.78%.
#> ✔ All chains have mixed well (Rhat < 1.1).
#> ✔ Effective sample sizes sufficient (eff_n_samples > 400).
#> 
#>  ✔  gplm finished   [============                                    ] 25%
#> 
#> Progress:
#> Initializing Metropolis MCMC algorithm...
#> Multiprocess sampling (4 chains in 2 jobs) ...
#> 
#> MCMC sampling completed!
#> 
#> Diagnostics:
#> Acceptance rate: 31.19%.
#> ✔ All chains have mixed well (Rhat < 1.1).
#> ✔ Effective sample sizes sufficient (eff_n_samples > 400).
#> 
#>  ✔  gplm0 finished  [========================                        ] 50%
#> 
#> Progress:
#> Initializing Metropolis MCMC algorithm...
#> Multiprocess sampling (4 chains in 2 jobs) ...
#> 
#> MCMC sampling completed!
#> 
#> Diagnostics:
#> Acceptance rate: 25.85%.
#> ✔ All chains have mixed well (Rhat < 1.1).
#> ✔ Effective sample sizes sufficient (eff_n_samples > 400).
#> 
#>  ✔  plm finished    [====================================            ] 75%
#> 
#> Progress:
#> Initializing Metropolis MCMC algorithm...
#> Multiprocess sampling (4 chains in 2 jobs) ...
#> 
#> MCMC sampling completed!
#> 
#> Diagnostics:
#> Acceptance rate: 35.55%.
#> ✔ All chains have mixed well (Rhat < 1.1).
#> ✔ Effective sample sizes sufficient (eff_n_samples > 400).
#> 
#>  ✔  plm0 finished   [================================================] 100%
# }
```

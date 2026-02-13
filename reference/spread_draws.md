# Spread MCMC chain draws to data.frame on a wide format

Useful to convert MCMC chain draws of particular parameters or output
from the model object to a wide format for further data wrangling

## Usage

``` r
spread_draws(mod, ..., transformed = FALSE)
```

## Arguments

- mod:

  An object of class "plm0", "plm", "gplm0" or "gplm".

- ...:

  Any number of character vectors containing valid names of parameters
  in the model or "rating_curve" and "rating_curve_mean". Also accepts
  "latent_parameters" and "hyperparameters".

- transformed:

  A boolean value determining whether the output is to be represented on
  the transformed scale used for sampling in the MCMC chain or the
  original scale. Defaults to FALSE.

## Value

A data frame with columns:

- `chain`:

  The chain number.

- `iter`:

  The iteration number.

- `param`:

  The parameter name.

- `value`:

  The parameter value.

## References

Hrafnkelsson, B., Sigurdarson, H., Rögnvaldsson, S., Jansson, A. Ö.,
Vias, R. D., and Gardarsson, S. M. (2022). Generalization of the
power-law rating curve using hydrodynamic theory and Bayesian
hierarchical modeling, Environmetrics, 33(2):e2711. doi:
https://doi.org/10.1002/env.2711

## See also

[`plm0`](https://sor16.github.io/bdrc/reference/plm0.md),
[`plm`](https://sor16.github.io/bdrc/reference/plm.md),
[`gplm0`](https://sor16.github.io/bdrc/reference/gplm0.md),
[`gplm`](https://sor16.github.io/bdrc/reference/gplm.md) for further
information on parameters

## Examples

``` r
# \donttest{
 data(krokfors)
 set.seed(1)
 plm0.fit <- plm0(formula=Q~W,data=krokfors,num_cores=2)
#> Progress:
#> Initializing Metropolis MCMC algorithm...
#> Multiprocess sampling (4 chains in 2 jobs) ...
#> 
#> MCMC sampling completed!
#> 
#> Diagnostics:
#> Acceptance rate: 36.09%.
#> ✔ All chains have mixed well (Rhat < 1.1).
#> ✔ Effective sample sizes sufficient (eff_n_samples > 400).
 hyp_samples <- spread_draws(plm0.fit,'hyperparameters')
 head(hyp_samples)
#>   chain iter        c sigma_eps
#> 1     1    1 7.693275 0.2421924
#> 2     1    2 7.708797 0.2256351
#> 3     1    3 7.704167 0.1679654
#> 4     1    4 7.666575 0.1807933
#> 5     1    5 7.658464 0.1570101
#> 6     1    6 7.709366 0.1968574
 rating_curve_samples <- spread_draws(plm0.fit,'rating_curve','rating_curve_mean')
 head(rating_curve_samples)
#>   chain iter        h rating_curve rating_curve_mean
#> 1     1    1 7.673811 0.000000e+00      0.000000e+00
#> 2     1    2 7.673811 0.000000e+00      0.000000e+00
#> 3     1    3 7.673811 0.000000e+00      0.000000e+00
#> 4     1    4 7.673811 1.094841e-06      1.018866e-06
#> 5     1    5 7.673811 6.537614e-06      6.335393e-06
#> 6     1    6 7.673811 0.000000e+00      0.000000e+00
# }
```

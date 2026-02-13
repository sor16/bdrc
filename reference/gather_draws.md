# Gather MCMC chain draws to data.frame on a long format

Useful to convert MCMC chain draws of particular parameters or output
from the model object to a long format for further data wrangling

## Usage

``` r
gather_draws(mod, ..., transformed = F)
```

## Arguments

- mod:

  An object of class "plm0", "plm", "gplm0" or "gplm".

- ...:

  Any number of character vectors containing valid names of parameters
  in the model or "rating_curve" and "rating_curve_mean". Also accepts
  "latent_parameters" and "hyperparameters".

- transformed:

  A boolean value determining whether the parameter is to be represented
  on the transformed scale used for sampling in the MCMC chain or the
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
 hyp_samples <- gather_draws(plm0.fit,'hyperparameters')
 head(hyp_samples)
#>   chain iter name    value
#> 1     1    1    c 7.693275
#> 2     1    2    c 7.708797
#> 3     1    3    c 7.704167
#> 4     1    4    c 7.666575
#> 5     1    5    c 7.658464
#> 6     1    6    c 7.709366
 rating_curve_samples <- gather_draws(plm0.fit,'rating_curve','rating_curve_mean')
 head(rating_curve_samples)
#>   chain iter        h         name        value
#> 1     1    1 7.673811 rating_curve 0.000000e+00
#> 2     1    2 7.673811 rating_curve 0.000000e+00
#> 3     1    3 7.673811 rating_curve 0.000000e+00
#> 4     1    4 7.673811 rating_curve 1.094841e-06
#> 5     1    5 7.673811 rating_curve 6.537614e-06
#> 6     1    6 7.673811 rating_curve 0.000000e+00
# }
```

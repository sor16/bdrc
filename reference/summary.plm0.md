# Summary method for discharge rating curves

Summarize a discharge rating curve model object

## Usage

``` r
# S3 method for class 'plm0'
summary(object, ...)

# S3 method for class 'plm'
summary(object, ...)

# S3 method for class 'gplm0'
summary(object, ...)

# S3 method for class 'gplm'
summary(object, ...)
```

## Arguments

- object:

  an object of class "plm0", "plm", "gplm0" or "gplm".

- ...:

  Not used in this function

## Functions

- `summary(plm0)`: Summary method for plm0

- `summary(plm)`: Summary method for plm

- `summary(gplm0)`: Summary method for gplm0

- `summary(gplm)`: Summary method for gplm

## See also

[`plm0`](https://sor16.github.io/bdrc/reference/plm0.md),
[`plm`](https://sor16.github.io/bdrc/reference/plm.md),
[`gplm0`](https://sor16.github.io/bdrc/reference/gplm0.md) and
[`gplm`](https://sor16.github.io/bdrc/reference/gplm.md) for fitting a
discharge rating curve. It is also useful to look at
[`plot.plm0`](https://sor16.github.io/bdrc/reference/plot.plm0.md),
[`plot.plm`](https://sor16.github.io/bdrc/reference/plot.plm0.md),
[`plot.gplm0`](https://sor16.github.io/bdrc/reference/plot.plm0.md) and
[`plot.gplm`](https://sor16.github.io/bdrc/reference/plot.plm0.md) to
help visualize all aspects of the fitted discharge rating curve.
Additionally,
[`spread_draws`](https://sor16.github.io/bdrc/reference/spread_draws.md)
and
[`spread_draws`](https://sor16.github.io/bdrc/reference/spread_draws.md)
help working directly with the MCMC samples.

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
summary(plm0.fit)
#> 
#> Formula: 
#>  Q ~ W
#> Latent parameters:
#>   lower-2.5% median-50% upper-97.5%
#> a 1.15       1.64       2.12       
#> b 2.46       2.82       3.21       
#> 
#> Hyperparameters:
#>           lower-2.5% median-50% upper-97.5%
#> c         7.580      7.676      7.752      
#> sigma_eps 0.151      0.191      0.253      
#> 
#> WAIC: -2.864096
# }
```

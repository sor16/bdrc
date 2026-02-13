# Predict method for discharge rating curves

Predict the discharge for given stage values based on a discharge rating
curve model object.

## Usage

``` r
# S3 method for class 'plm0'
predict(object, newdata = NULL, wide = FALSE, ...)

# S3 method for class 'plm'
predict(object, newdata = NULL, wide = FALSE, ...)

# S3 method for class 'gplm0'
predict(object, newdata = NULL, wide = FALSE, ...)

# S3 method for class 'gplm'
predict(object, newdata = NULL, wide = FALSE, ...)
```

## Arguments

- object:

  An object of class "plm0", "plm", "gplm0" or "gplm".

- newdata:

  A numeric vector of stage values for which to predict. If omitted, the
  stage values in the data are used.

- wide:

  A logical value denoting whether to produce a wide prediction output.
  If TRUE, then the output is a table with median prediction values for
  an equally spaced grid of stages with 1 cm increments, each row
  containing predictions in a decimeter range of stages.

- ...:

  Not used in this function

## Value

An object of class "data.frame" with four columns:

- `h`:

  The stage.

- `lower`:

  The 2.5% posterior predictive quantile.

- `median`:

  The 50% posterior predictive quantile.

- `upper`:

  The 97.5% posterior predictive quantile.

If wide=TRUE, a matrix as described above (see wide parameter) is
returned.

## Functions

- `predict(plm0)`: Predict method for plm0

- `predict(plm)`: Predict method for plm

- `predict(gplm0)`: Predict method for gplm0

- `predict(gplm)`: Predict method for gplm

## See also

[`plm0`](https://sor16.github.io/bdrc/reference/plm0.md),
[`plm`](https://sor16.github.io/bdrc/reference/plm.md),
[`gplm0`](https://sor16.github.io/bdrc/reference/gplm0.md) and
[`gplm`](https://sor16.github.io/bdrc/reference/gplm.md) for fitting a
discharge rating curve and
[`summary.plm0`](https://sor16.github.io/bdrc/reference/summary.plm0.md),
[`summary.plm`](https://sor16.github.io/bdrc/reference/summary.plm0.md),
[`summary.gplm0`](https://sor16.github.io/bdrc/reference/summary.plm0.md)
and
[`summary.gplm`](https://sor16.github.io/bdrc/reference/summary.plm0.md)
for summaries. It is also useful to look at
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
plm0.fit <- plm0(formula=Q~W,data=krokfors,h_max=10,num_cores=2)
#> Progress:
#> Initializing Metropolis MCMC algorithm...
#> Multiprocess sampling (4 chains in 2 jobs) ...
#> 
#> MCMC sampling completed!
#> 
#> Diagnostics:
#> Acceptance rate: 36.22%.
#> ✔ All chains have mixed well (Rhat < 1.1).
#> ✔ Effective sample sizes sufficient (eff_n_samples > 400).
#predict rating curve on a equally 10 cm spaced grid from 9 to 10 meters
predict(plm0.fit,newdata=seq(9,10,by=0.1))
#>       h     lower    median     upper
#> 1   9.0  2.423761  3.617127  5.355104
#> 2   9.1  2.988084  4.448245  6.528722
#> 3   9.2  3.616549  5.374700  7.940233
#> 4   9.3  4.291666  6.418404  9.515431
#> 5   9.4  5.056743  7.601828 11.323903
#> 6   9.5  5.914344  8.919305 13.169163
#> 7   9.6  6.893257 10.380845 15.570683
#> 8   9.7  7.860588 11.974804 17.826516
#> 9   9.8  9.052186 13.709759 20.592950
#> 10  9.9 10.158058 15.660659 23.351463
#> 11 10.0 11.489777 17.692248 26.874426
# }
```

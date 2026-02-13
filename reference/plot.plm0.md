# Plot method for discharge rating curves

Visualize discharge rating curve model objects

## Usage

``` r
# S3 method for class 'plm0'
plot(
  x,
  type = "rating_curve",
  param = NULL,
  transformed = FALSE,
  title = NULL,
  xlim = NULL,
  ylim = NULL,
  ...
)

# S3 method for class 'plm'
plot(
  x,
  type = "rating_curve",
  param = NULL,
  transformed = FALSE,
  title = NULL,
  xlim = NULL,
  ylim = NULL,
  ...
)

# S3 method for class 'gplm0'
plot(
  x,
  type = "rating_curve",
  param = NULL,
  transformed = FALSE,
  title = NULL,
  xlim = NULL,
  ylim = NULL,
  ...
)

# S3 method for class 'gplm'
plot(
  x,
  type = "rating_curve",
  param = NULL,
  transformed = FALSE,
  title = NULL,
  xlim = NULL,
  ylim = NULL,
  ...
)
```

## Arguments

- x:

  An object of class "plm0", "plm", "gplm0" or "gplm".

- type:

  A character denoting what type of plot should be drawn. Defaults to
  "rating_curve". Possible types are:

  `rating_curve`

  :   Plots the rating curve.

  `rating_curve_mean`

  :   Plots the posterior mean of the rating curve.

  `f`

  :   Plots the power-law exponent.

  `beta`

  :   Plots the random effect in the power-law exponent.

  `sigma_eps`

  :   Plots the standard deviation on the data level.

  `residuals`

  :   Plots the log residuals.

  `trace`

  :   Plots trace plots of parameters given in param.

  `histogram`

  :   Plots histograms of parameters given in param.

  `panel`

  :   Plots a 2x2 panel of plots: "rating curve", "residuals", "f" and
      "sigma_eps".

- param:

  A character vector with the parameters to plot. Defaults to NULL and
  is only used if type is "trace" or "histogram". Allowed values are the
  parameters given in the model summary of x as well as
  "hyperparameters" or "latent_parameters" for specific groups of
  parameters.

- transformed:

  A logical value indicating whether the quantity should be plotted on a
  transformed scale used during the Bayesian inference. Defaults to
  FALSE.

- title:

  A character denoting the title of the plot.

- xlim:

  A numeric vector of length 2, denoting the limits on the x axis of the
  plot. Applicable for types "rating_curve", "rating_curve_mean", "f",
  "beta", "sigma_eps", "residuals".

- ylim:

  A numeric vector of length 2, denoting the limits on the y axis of the
  plot. Applicable for types "rating_curve", "rating_curve_mean", "f",
  "beta", "sigma_eps", "residuals".

- ...:

  Not used in this function

## Value

No return value, called for side effects.

## Functions

- `plot(plm0)`: Plot method for plm0

- `plot(plm)`: Plot method for plm

- `plot(gplm0)`: Plot method for gplm0

- `plot(gplm)`: Plot method for gplm

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
[`spread_draws`](https://sor16.github.io/bdrc/reference/spread_draws.md)
and
[`gather_draws`](https://sor16.github.io/bdrc/reference/gather_draws.md)
to work directly with the MCMC samples.

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

plot(plm0.fit)

plot(plm0.fit,transformed=TRUE)

plot(plm0.fit,type='histogram',param='c')

plot(plm0.fit,type='histogram',param='c',transformed=TRUE)

plot(plm0.fit,type='histogram',param='hyperparameters')

plot(plm0.fit,type='histogram',param='latent_parameters')

plot(plm0.fit,type='residuals')

plot(plm0.fit,type='f')

plot(plm0.fit,type='sigma_eps')

# }
```

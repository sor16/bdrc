# Plot bdrc model objects

Visualize results from model objects in bdrc, plm0, plm, gplm0, gplm

## Usage

``` r
plot_fun(
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

  A character denoting the title of the plot. Defaults to NULL, i.e. no
  title.

- xlim:

  A numeric vector of length 2, denoting the limits on the x axis of the
  plot. Applicable for types "rating_curve", "rating_curve_mean", "f",
  "beta", "sigma_eps", "residuals".

- ylim:

  A numeric vector of length 2, denoting the limits on the y axis of the
  plot. Applicable for types "rating_curve", "rating_curve_mean", "f",
  "beta", "sigma_eps", "residuals".

## Value

Returns an object of class ggplot2.

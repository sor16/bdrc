# Print method for discharge rating curves

Print a discharge rating curve model object

## Usage

``` r
# S3 method for class 'plm0'
print(x, ...)

# S3 method for class 'plm'
print(x, ...)

# S3 method for class 'gplm0'
print(x, ...)

# S3 method for class 'gplm'
print(x, ...)
```

## Arguments

- x:

  an object of class "plm0", "plm", "gplm0" or "gplm".

- ...:

  Not used in this function

## Functions

- `print(plm0)`: Print method for plm0

- `print(plm)`: Print method for plm

- `print(gplm0)`: Print method for gplm0

- `print(gplm)`: Print method for gplm

## See also

[`plm0`](https://sor16.github.io/bdrc/reference/plm0.md),
[`plm`](https://sor16.github.io/bdrc/reference/plm.md),
[`gplm0`](https://sor16.github.io/bdrc/reference/gplm0.md),
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

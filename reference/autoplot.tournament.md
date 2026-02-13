# Autoplot method for discharge rating curve tournament

Compare the four discharge rating curves from the tournament object in
different ways

## Usage

``` r
# S3 method for class 'tournament'
autoplot(object, type = "boxplot", ...)
```

## Arguments

- object:

  An object of class "tournament"

- type:

  A character denoting what type of plot should be drawn. Possible types
  are

- ...:

  Not used in this function

  `boxplot`

  :   Creates a boxplot of the posterior log-likelihood values
      transformed to the deviance scale.

## Value

Returns an object of class "ggplot2".

## See also

[`tournament`](https://sor16.github.io/bdrc/reference/tournament.md) to
run a discharge rating curve tournament and
[`summary.tournament`](https://sor16.github.io/bdrc/reference/summary.tournament.md)
for summaries.

## Examples

``` r
# \donttest{
library(ggplot2)
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
autoplot(t_obj)

# }
```

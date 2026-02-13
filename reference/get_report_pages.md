# Report pages for a discharge rating curve or tournament

Get a list of the pages of a report on a discharge rating curve model or
tournament

## Usage

``` r
get_report_pages(x, type = 1)

# S3 method for class 'plm0'
get_report_pages(x, type = 1)

# S3 method for class 'plm'
get_report_pages(x, type = 1)

# S3 method for class 'gplm0'
get_report_pages(x, type = 1)

# S3 method for class 'gplm'
get_report_pages(x, type = 1)

# S3 method for class 'tournament'
get_report_pages(x, type = 1)
```

## Arguments

- x:

  An object of class "tournament", "plm0", "plm", "gplm0" or "gplm".

- type:

  An integer denoting what type of report is to be produced. Defaults to
  type 1. Possible types are:

  `1`

  :   Produces a report displaying the results of the model (winning
      model if a tournament provided). The first page contains a panel
      of four plots and a summary of the posterior distributions of the
      parameters. On the second page a tabular prediction of discharge
      on an equally spaced grid of stages is displayed. This prediction
      table can span multiple pages.

  `2`

  :   Produces a ten page report and is only permissible for objects of
      class "tournament". The first four pages contain a panel of four
      plots and a summary of the posterior distributions of the
      parameters for each of the four models in the tournament, the
      fifth page shows a summary of the tournament model comparison, the
      sixth page convergence diagnostics plots, and the final four pages
      shows the histograms of the parameters in each of the four models.

## Value

A list of objects of type "grob" that correspond to the pages in a
rating curve report.

## Methods (by class)

- `get_report_pages(plm0)`: Get report pages for plm0 model object

- `get_report_pages(plm)`: Get report pages for plm model object

- `get_report_pages(gplm0)`: Get report pages for gplm0 model object

- `get_report_pages(gplm)`: Get report pages for gplm model object

- `get_report_pages(tournament)`: Get report pages for discharge rating
  curve tournament model object

## See also

[`tournament`](https://sor16.github.io/bdrc/reference/tournament.md) for
running a tournament,
[`summary.tournament`](https://sor16.github.io/bdrc/reference/summary.tournament.md)
for summaries and
[`get_report`](https://sor16.github.io/bdrc/reference/get_report.md) for
generating and saving a report of a tournament object.

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
plm0_pages <- get_report_pages(plm0.fit)
# }
```

# Report for a discharge rating curve or tournament

Save a pdf file with a report of a discharge rating curve object or
tournament.

## Usage

``` r
get_report(x, path = NULL, type = 1, ...)

# S3 method for class 'plm0'
get_report(x, path = NULL, type = 1, ...)

# S3 method for class 'plm'
get_report(x, path = NULL, type = 1, ...)

# S3 method for class 'gplm0'
get_report(x, path = NULL, type = 1, ...)

# S3 method for class 'gplm'
get_report(x, path = NULL, type = 1, ...)

# S3 method for class 'tournament'
get_report(x, path = NULL, type = 1, ...)
```

## Arguments

- x:

  An object of class "tournament", "plm0", "plm", "gplm0" or "gplm".

- path:

  A file path to which the pdf file of the report is saved. If NULL, the
  current working directory is used.

- type:

  An integer denoting what type of report is to be produced. Defaults to
  type 1. Only type 1 is permissible for an object of class "plm0",
  "plm", "gplm0" or "gplm". Possible types are:

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

- ...:

  Further arguments passed to other methods (currently unused).

## Value

No return value, called for side effects.

## Details

This function can only be used in an interactive R session as it asks
permission from the user to write to their file system.

## Methods (by class)

- `get_report(plm0)`: Get report for plm0 model object

- `get_report(plm)`: Get report for plm model object

- `get_report(gplm0)`: Get report for gplm0 model object

- `get_report(gplm)`: Get report for gplm

- `get_report(tournament)`: Get report for discharge rating curve
  tournament

## See also

`get_report` for generating and saving a report.

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
# }
if (FALSE) { # \dontrun{
get_report(plm0.fit)
} # }
```

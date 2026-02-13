# Summary method for a discharge rating curve tournament

Print the summary of a tournament of model comparisons. This function
allows for an efficient and fast re-run of the tournament with different
methods or winning criteria.

## Usage

``` r
# S3 method for class 'tournament'
summary(object, method = NULL, winning_criteria = NULL, ...)
```

## Arguments

- object:

  An object of class "tournament"

- method:

  Optional; a string specifying the method to use for the summary. If
  NULL, uses the method from the original tournament. Options are
  "WAIC", "DIC", or "PMP".

- winning_criteria:

  Optional; specifies new winning criteria for the summary. If NULL,
  uses the criteria from the original tournament. See Details in
  [`tournament`](https://sor16.github.io/bdrc/reference/tournament.md)
  for proper formatting.

- ...:

  Not used in this function

## Value

Prints the summary to the console.

## Details

If either `method` or `winning_criteria` is provided, the function
re-runs the tournament with the new parameters using the fitted models.

## See also

[`tournament`](https://sor16.github.io/bdrc/reference/tournament.md) to
run a discharge rating curve tournament and
[`plot.tournament`](https://sor16.github.io/bdrc/reference/plot.tournament.md)
for visualizing the model comparison

## Examples

``` r
# \donttest{
data(krokfors)
set.seed(1)
t_obj <- tournament(Q ~ W, krokfors, num_cores = 2)
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
summary(t_obj)
#> 
#> === Tournament Model Comparison Summary ===
#> 
#> Method: WAIC 
#> Winning Criteria: Delta_WAIC > 2 
#> Overall Winner: gplm0 
#> 
#> Comparison 1 Results:
#> -------------------------------------------------------------------------------------------------- 
#> complexity   model  winner lppd     eff_num_param  WAIC       SE_WAIC    Delta_WAIC SE_Delta_WAIC
#> more         gplm          20.7794  6.8706         -27.8176   11.8918    0.5570     0.2416      
#> less         gplm0  <---   20.3710  6.7406         -27.2606   12.0360                           
#> 
#> Comparison 2 Results:
#> -------------------------------------------------------------------------------------------------- 
#> complexity   model  winner lppd     eff_num_param  WAIC       SE_WAIC    Delta_WAIC SE_Delta_WAIC
#> more         plm           5.5842   4.2574         -2.6536    6.6635     -0.4066    0.1904      
#> less         plm0   <---   5.6284   4.0984         -3.0601    6.6931                            
#> 
#> Comparison 3 Results:
#> -------------------------------------------------------------------------------------------------- 
#> complexity   model  winner lppd     eff_num_param  WAIC       SE_WAIC    Delta_WAIC SE_Delta_WAIC
#> more         gplm0  <---   20.3710  6.7406         -27.2606   12.0360    24.2005    9.1834      
#> less         plm0          5.6284   4.0984         -3.0601    6.6931                            
#> 
#> === End of Summary ===

# Re-run summary with different method
summary(t_obj, method = "DIC")
#> 
#> === Tournament Model Comparison Summary ===
#> 
#> Method: DIC 
#> Winning Criteria: Delta_DIC > 2 
#> Overall Winner: gplm0 
#> 
#> Comparison 1 Results:
#> -------------------------------------------------------------------------- 
#> complexity   model  winner D_hat      eff_num_param  DIC        Delta_DIC 
#> more         gplm          -42.7732   6.0631         -30.6470   0.8498    
#> less         gplm0  <---   -42.4880   6.3454         -29.7972             
#> 
#> Comparison 2 Results:
#> -------------------------------------------------------------------------- 
#> complexity   model  winner D_hat      eff_num_param  DIC        Delta_DIC 
#> more         plm           -11.3400   2.9739         -5.3923    -0.3058   
#> less         plm0   <---   -11.4712   2.8865         -5.6982              
#> 
#> Comparison 3 Results:
#> -------------------------------------------------------------------------- 
#> complexity   model  winner D_hat      eff_num_param  DIC        Delta_DIC 
#> more         gplm0  <---   -42.4880   6.3454         -29.7972   24.0991   
#> less         plm0          -11.4712   2.8865         -5.6982              
#> 
#> === End of Summary ===

# Re-run summary with different winning criteria
summary(t_obj, winning_criteria = "Delta_WAIC > 3")
#> 
#> === Tournament Model Comparison Summary ===
#> 
#> Method: WAIC 
#> Winning Criteria: Delta_WAIC > 3 
#> Overall Winner: gplm0 
#> 
#> Comparison 1 Results:
#> -------------------------------------------------------------------------------------------------- 
#> complexity   model  winner lppd     eff_num_param  WAIC       SE_WAIC    Delta_WAIC SE_Delta_WAIC
#> more         gplm          20.7794  6.8706         -27.8176   11.8918    0.5570     0.2416      
#> less         gplm0  <---   20.3710  6.7406         -27.2606   12.0360                           
#> 
#> Comparison 2 Results:
#> -------------------------------------------------------------------------------------------------- 
#> complexity   model  winner lppd     eff_num_param  WAIC       SE_WAIC    Delta_WAIC SE_Delta_WAIC
#> more         plm           5.5842   4.2574         -2.6536    6.6635     -0.4066    0.1904      
#> less         plm0   <---   5.6284   4.0984         -3.0601    6.6931                            
#> 
#> Comparison 3 Results:
#> -------------------------------------------------------------------------------------------------- 
#> complexity   model  winner lppd     eff_num_param  WAIC       SE_WAIC    Delta_WAIC SE_Delta_WAIC
#> more         gplm0  <---   20.3710  6.7406         -27.2606   12.0360    24.2005    9.1834      
#> less         plm0          5.6284   4.0984         -3.0601    6.6931                            
#> 
#> === End of Summary ===
# }
```

# Compare two models using a specified model-selection criteria

evaluate_comparison uses the Widely Applicable Information Criterion
(WAIC), the Deviance Information Criterion (DIC), or the posterior model
probabilities (PMP), calculated with Bayes factor, to determine whether
one model is more appropriate than the other given the data at hand.

## Usage

``` r
evaluate_comparison(m, method, winning_criteria)
```

## Arguments

- m:

  A list of two model objects fit on the same dataset. The allowed model
  objects are "gplm", "gplm0", "plm" and "plm0"

- method:

  A string specifying the method used to estimate the predictive
  performance of the models. The allowed methods are "WAIC", "DIC" and
  "PMP".

- winning_criteria:

  For "WAIC", it can be either a numeric value or a string expression.
  For "DIC", it must be a numeric value. For "PMP", it must be a numeric
  value between 0 and 1. This sets the threshold for determining the
  more appropriate model. See Details for more information.

## Value

A data.frame with the summary of the results of each comparison,
including:

- complexity: Indicates whether a model is the "more" or "less" complex
  model in a comparison

- model: The type of model (gplm, gplm0, plm, or plm0)

- Method-specific columns (see Details)

- winner: Logical value indicating if the model was selected

## Details

For "WAIC" method:

- If winning_criteria is numeric, the more complex model wins if
  Delta_WAIC \> winning_criteria

- If winning_criteria is a string, it must be a valid R expression using
  Delta_WAIC and/or SE_Delta_WAIC

- Returns columns: lppd, eff_num_param, WAIC, SE_WAIC, Delta_WAIC,
  SE_Delta_WAIC

For "DIC" method:

- winning_criteria must be numeric

- The more complex model wins if Delta_DIC \> winning_criteria

- Returns columns: D_hat, eff_num_param, DIC, Delta_DIC

For "PMP" method:

- winning_criteria must be a numeric value between 0 and 1

- The more complex model wins if its PMP \> winning_criteria

- Returns columns: log_marg_lik, PMP

## References

Hrafnkelsson, B., Sigurdarson, H., Rögnvaldsson, S., Jansson, A. Ö.,
Vias, R. D., and Gardarsson, S. M. (2022). Generalization of the
power-law rating curve using hydrodynamic theory and Bayesian
hierarchical modeling, Environmetrics, 33(2):e2711. doi:
https://doi.org/10.1002/env.2711

## See also

[`tournament`](https://sor16.github.io/bdrc/reference/tournament.md)

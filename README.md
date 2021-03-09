
<!-- README.md is generated from README.Rmd. Please edit that file -->

## bdrc - Bayesian Discharge Rating Curves <img src="man/figures/logo.png" align="right" alt="" width="140" />

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/sor16/bdrc/branch/master/graph/badge.svg)](https://codecov.io/gh/sor16/bdrc?branch=master)
<!-- badges: end -->

The package implements the following Bayesian hierarchical discharge
rating curve models for paired measurements of stage and discharge in
rivers described in Hrafnkelsson et al.:

`bplm0()` - Power-law model with constant variance

`bplm()` - Power-law model with variance that may vary with stage

`bgplm0()` - Generalized power-law model with constant variance

`bgplm()` - Generalized power-law model with variance that may vary with
stage

## Installation

``` r
# Install release version from CRAN
install.packages("bdrc")
# Install development version from GitHub
devtools::install_github("sor16/bdrc")
```

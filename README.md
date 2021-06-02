
<!-- README.md is generated from README.Rmd. Please edit that file -->

## bdrc - Bayesian Discharge Rating Curves <img src="man/figures/logo.png" align="right" alt="" width="140" />

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/sor16/bdrc/branch/master/graph/badge.svg)](https://codecov.io/gh/sor16/bdrc?branch=master)
[![R build
status](https://github.com/sor16/bdrc/workflows/R-CMD-check/badge.svg)](https://github.com/sor16/bdrc/actions)
<!-- badges: end -->

The package implements the following Bayesian hierarchical discharge
rating curve models for paired measurements of stage and discharge in
rivers described in Hrafnkelsson et al.:

`plm0()` - Power-law model with constant variance

`plm()` - Power-law model with variance that may vary with stage

`gplm0()` - Generalized power-law model with constant variance

`gplm()` - Generalized power-law model with variance that may vary with
stage

## Installation

``` r
# Install release version from CRAN
install.packages("bdrc")
# Install development version from GitHub
devtools::install_github("sor16/bdrc")
```

## Getting started

It is very simple to fit a discharge rating curve with the *bdrc*
package. All you need are two mandatory input arguments, formula and
data. The formula is of the form y~x where y is discharge in \(m^3/s\)
and x is stage in \(m\) (it is very important that the data is in the
correct units). data is a data.frame which must include x and y as
column names. As an example we will use data from the swedish river
*Halla* which is one of the dataset that comes with the package. In this
table, the Q column denotes discharge while W denotes stage:

``` r
gplm.fit <- gplm(Q~W,halla)
```

To dig deeper into the functionality of the package and the different
ways to visualize your discharge rating curve model, we recommend taking
a look at our two vignettes.

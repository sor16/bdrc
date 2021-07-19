
<!-- README.md is generated from README.Rmd. Please edit that file -->

## bdrc - Bayesian Discharge Rating Curves <img src="man/figures/logo.png" align="right" alt="" width="140" />

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/sor16/bdrc/branch/master/graph/badge.svg)](https://codecov.io/gh/sor16/bdrc?branch=master)
[![R build
status](https://github.com/sor16/bdrc/workflows/R-CMD-check/badge.svg)](https://github.com/sor16/bdrc/actions)
<!-- badges: end -->

This software package fits a discharge rating curve based on the
power-law and the generalized power-law from data on paired stage and
discharge measurements in a given river using a Bayesian hierarchical
model as described in Hrafnkelsson et al. (2020). Four models are
implemented:

`plm0()` - Power-law model with a constant variance. This is a Bayesian
hierarchical implementation of the most commonly used discharge rating
curve model in hydrological practice.

`plm()` - Power-law model with variance that varies with stage.

`gplm0()` - Generalized power-law model with a constant variance. The
generalized power-law is introduced in Hrafnkelsson et al. (2020).

`gplm()` - Generalized power-law model with variance that varies with
stage. The generalized power-law is introduced in Hrafnkelsson et
al. (2020).

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
data. The formula is of the form y\~x where y is discharge in
m<sup>3</sup>/s and x is stage in m (it is very important that the data
is in the correct units). data is a data.frame which must include x and
y as column names. As an example, we will use data from the Swedish
gauging station *Krokfors*, which is one of the datasets that come with
the package. In this table, the Q column denotes discharge while W
denotes stage:

``` r
gplm.fit <- gplm(Q~W,krokfors)
```

To dig deeper into the functionality of the package and the different
ways to visualize a discharge rating curve model for your data, we
recommend taking a look at our two vignettes.

## References

Hrafnkelsson, B., Sigurdarson, H., and Gardarsson, S. M. (2020).
*Generalization of the power-law rating curve using hydrodynamic theory
and Bayesian hierarchical modeling*. arXiv preprint 2010.04769.

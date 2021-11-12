
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ssnet

<!-- badges: start -->

[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![Travis build
status](https://travis-ci.org/jmleach-bst/ssnet.svg?branch=master)](https://travis-ci.org/jmleach-bst/ssnet)
<!-- badges: end -->

The goal of `ssnet` is to fit spike-and-slab elastic net GLM’s with or
without spatially structured priors. An expectation maximization (EM)
algorithm is used to fit the models, and spatial structure is
incorporated by using Intrinsic Autoregressions (IAR) as priors for
inclusion probabilities. This allows for variable selection to
incorporate assumptions about spatial clustering of variables that
should (not) be included in the model. Outcome distributions supported
are Binomial, Normal, Poisson, and Multinomial.

## Installation

<!-- You can install the released version of ssnet from [CRAN](https://CRAN.R-project.org) with: -->
<!-- ``` r -->
<!-- install.packages("ssnet") -->
<!-- ``` -->

The `R` package `ssnet` is in development, which version can be
installed from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jmleach-bst/ssnet")
```

## Example

Coming soon - examples and vignette. However, note that the
documentation is thorough (in my opinion), and the most useful
functions’ documentation contain examples.

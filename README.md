
<!-- README.md is generated from README.Rmd. Please edit that file -->

# A Practical Response Adaptive Block Randomization (RABR) Design with Analytic Type I Error Protection

<!-- badges: start -->
<!-- badges: end -->

This package is to evaluate type I error rate, power, and operating
characteristics of RABR via simulations.

## Installation

You can install the released version of RABR from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("RABR")
```

## Example

We provide an example of continuous endpoint. One may refer to the
vignette for more details.

``` r
library(RABR)
library(parallel)
library(doParallel)
#> Loading required package: foreach
#> Loading required package: iterators
RABR.fit = RABRcontinuous(
            MeanVec = c(0.43, 0.48, 0.63, 1.2),
            SdVec = c(1, 1, 1, 1),
            M = 60,
            N = 120,
            R = c(8, 9, 2, 1),
            Nitt = 1000,
            Alpha = 0.025,
            Ncluster = 2,
            Seed = 12345,
            MultiMethod = "dunnett")
##
## Probability of rejecting each elementary null
## hypothesis without multiplicity adjustment
   print(RABR.fit$ProbUnadj)
#> [1] 0.027 0.093 0.877
##
## Probability of rejecting each elementary null
## hypothesis with multiplicity adjustment
   print(RABR.fit$ProbAdj)
#> [1] 0.017 0.062 0.804
##
## Probability of selecting and confirming the
## efficacy of each active treatment group
   print(RABR.fit$ProbAdjSelected)
#> [1] 0.001 0.007 0.802
##
## ProbAdjOverall Probability of rejecting at
## least one elementary null hypothesis
## with multiplicity adjustment
   print(RABR.fit$ProbAdjOverall)
#> [1] 0.81
##
## ASN Average sample size of placebo and active
## treatment groups
   print(RABR.fit$ASN)
#> [1] 39.107 40.746 21.432 18.715
```

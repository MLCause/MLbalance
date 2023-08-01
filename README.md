
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MLbalance

<!-- badges: start -->
<!-- badges: end -->

MLbalance implements a novel machine learning balance test, the balance
permutation test, for experiments with binary, multiarm, and continuous
treatments. The purpose of this test is to detect failures of random
assignment and imbalance across treatment arms. For more detail, see
Rametta and Fuller (2023).

## Installation

You can install the development version of MLbalance from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("CetiAlphaFive/MLbalance")
# OR 
# install.packages("remotes)
remotes::install_github("CetiAlphaFive/MLbalance")
# OR 
# install.packages("pak")
pak::pak("CetiAlphaFive/MLbalance")
```

## Binary Treatment Example

Here is a basic example demonstrating the balance permutation test on a
simulated binary treatment DGP.

``` r
# install.packages("randomizr")
library(MLbalance)
#> 
#> Attaching package: 'MLbalance'
#> The following object is masked from 'package:utils':
#> 
#>     vi
#
# Simple simulation 
n <- 1000
p <- 20
X <- matrix(rnorm(n*p,0,1),n,p)
w_real <- rbinom(n, 1, ifelse(.021 + abs(.4*X[,4] - .5*X[,8]) < 1, .021 + abs(.4*X[,4] - .5*X[,8]), 1))
# install.packages("randomizr")
w_sim <- randomizr::complete_ra(N = n,m = sum(w_real))
e <- rnorm(n,0,1)
y <- 2*w_real*X[,4] + 3*X[,2] -2*X[,8] + e
df <- data.frame(y,w_real,w_sim,X)
#
r.check <- random_check(W_real = df$w_real, #real treatment assignment vector 
                        W_sim  = df$w_sim, #simulated vector, comment out this argument to use permutated real assignment vector instead 
                        X      = subset(df,select = -c(y,w_real,w_sim)) #matrix of pretreatment covariates (or any covariates that SHOULD NOT be related to the assignment process/mechanism
             ); r.check$plot
#> Simulated Assignemnt Vector Provided, Null Distribution Generated Using Simulated Treatment Assignment.
#> 
#> 
#> Simple Count Table(s)
#> 
#> W_real
#>   0   1 
#> 515 485 
#> W_sim
#>   0   1 
#> 515 485 
#> 
#> 
#> Result from difference in variances test (one-sided, greater F-test):
#> 
#>  Statistic p.val Result
#>   61.13291     0   FAIL
#> 
#> 
#> Check diff.var.result in saved output for detailed test result.
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

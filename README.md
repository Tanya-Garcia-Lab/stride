---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->


# stride

<!-- badges: start -->
This R package implements nonparametric estimation of survival models for mixture data.
Mixture data is data collected from multiple populations of interest, but we do not know from which population each observation belongs.  The corresponding references are:

Garcia, T.P. and Parast, L. (2020). Dynamic landmark prediction for mixture data. Biostatistics,  doi:10.1093/biostatistics/kxz052.

Garcia, T.P., Marder, K. and Wang, Y. (2017). Statistical modeling of Huntington disease onset.
In Handbook of Clinical Neurology, vol 144, 3rd Series, editors Andrew Feigin and Karen E. Anderson.

Qing, J., Garcia, T.P., Ma, Y., Tang, M.X., Marder, K., and Wang, Y. (2014).
Combining isotonic regression and EM algorithm to predict genetic risk under monotonicity constraint.
Annals of Applied Statistics, 8(2), 1182-1208.

Wang, Y., Garcia, T.P., and Ma. Y. (2012).  Nonparametric estimation for censored mixture data with
application to the Cooperative Huntington's Observational Research Trial. Journal of the American Statistical Association, 107, 1324-1338.
<!-- badges: end -->



## Installation

You can install the released version of stride from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("stride")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("tpgarcia/stride")
```
## Example

This is a basic example which shows you how to solve a common problem:


```r
library(stride)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so:


```r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" title="plot of chunk pressure" alt="plot of chunk pressure" width="100%" />

In that case, don't forget to commit and push the resulting figure files, so they display on GitHub!

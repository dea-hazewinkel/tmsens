# Summarizing Trimmed Means Linear Model Fits

`summary` method for class "`tm`".

## Usage

``` r
# S3 method for class 'tm'
summary(object, ...)
```

## Arguments

- object:

  an object of class "`tm`"

- ...:

  user specified arguments

## Value

`summary.tm` returns a list of summary statistics of the fitted trimmed
means linear model in `object`, with components

- call:

  the matched call

- n:

  the number of observations per treatment group

- dropout:

  the proportion of dropout per treatment group

- trimfrac:

  the proportion of data that was trimmed away per treatment group

- trimside:

  specifies if lower or higher value trimming was performed

- n_after_trimming:

  the number of observations per treatment group after trimming

- coefficients:

  an array of coefficients with corresponding p-values and 95%
  confidence intervals

- Analysis_details:

  reiterates trimming fraction and side, and, for `adj_est=TRUE`
  specifies if the adjustment was performed on the comparator or
  treatment group.

- SD_outcome:

  an array of the standard deviation per treatment group, for the
  observed outcomes and for the trimmed outcomes

## See also

[`tm`](https://dea-hazewinkel.github.io/tmsens/reference/tm.md). The
function [`coef`](https://rdrr.io/r/stats/coef.html) extracts the array
of regression coefficients with corresponding p-values and 95%
confidence intervals.

## Examples

``` r
set.seed(123456)
test_dat <- as.data.frame(cbind(c(rep(0,500),rep(1,500)),
                          c(sort(rnorm(500,0,1)),sort(rnorm(500,1,1.5))),
                          rbinom(1000,2,0.4), rnorm(1000,0,1)))
colnames(test_dat) <- c("TR", "Y", "U", "U2")
test_dat$Y[1:200] <- NA
tm_obj <- tm(formula= Y ~ TR + U + U2, GR = "TR", trF = 0.5,
             side = "LOW", n_perm = 1000, adj_est = TRUE, data = test_dat)
summary(tm_obj)
#> 
#> Call:
#> tm(formula = Y ~ TR + U + U2, GR = "TR", trF = 0.5, side = "LOW", 
#>     n_perm = 1000, adj_est = TRUE, data = test_dat)
#> 
#> 
#> Analysis details:
#> Treatment (1) group rescaled for adjusted TM estimate
#> 
#> Coefficients:
#>        Estimate      Pval 95% CI LO 95% CI HI
#> TR     1.354308  0.000999  1.177918     1.531
#> U     -0.071352  0.144855 -0.167223     0.025
#> U2     0.055790  0.090909 -0.008575     0.120
#> TRadj  0.969812  0.000999  0.840809     1.099
#> 
#> 
#> Dropout:
#>   1    0  
#>  0%  40%  
#> 
#> Sample size after trimming:
#>   1    0  
#> 250  250  
#> 
#> 
#> Trimming fraction: 
#> 50% lower value trimming
#> 
coef(tm_obj)
#>          Estimate        Pval    95% CI LO  95% CI HI
#> TR     1.35430808 0.000999001  1.177917612 1.53069856
#> U     -0.07135247 0.144855145 -0.167223391 0.02451845
#> U2     0.05579053 0.090909091 -0.008575234 0.12015629
#> TRadj  0.96981237 0.000999001  0.840809363 1.09881537
```

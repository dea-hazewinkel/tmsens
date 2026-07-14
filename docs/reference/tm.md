# Fitting Trimmed Mean Linear Models

`tm` performs a trimmed means analysis for data with a continuous
outcome/response and a binary treatment/exposure variable. Outcomes are
sorted and trimmed per treatment group, and a linear regression is
fitted using [`lm`](https://rdrr.io/r/stats/lm.html).

## Usage

``` r
tm(
  formula,
  GR,
  trF = NULL,
  side = c("LOW", "HIGH"),
  n_perm = 1000,
  adj_est = FALSE,
  data
)
```

## Arguments

- formula:

  an object of class [`formula`](https://rdrr.io/r/stats/formula.html),
  specifying the model, of the form `outcome ~ terms`, where `terms`
  must include the binary treatment variable, with additional variables
  optional.

- GR:

  a string denoting the name of the binary treatment variable. This
  function assumes the lowest value to be the comparator/reference group

- trF:

  a number between 0 and 1, specifying the trimming fraction: the
  proportion of the data that is trimmed away for each treatment group.
  `trF` should be equal to or greater than the largest observed dropout
  proportion. If left unspecified, the largest observed dropout
  proportion is used, or a trimming fraction of 0.5 if there is no
  dropout.

- side:

  specifies if higher value trimming (`"HIGH"`) or lower value trimming
  (`"LOW"`) should be performed. The default is `"LOW"`.

- n_perm:

  the number of permutations performed to obtain the p-value and 95%
  confidence intervals for the estimates. Default is 1000.

- adj_est:

  logical. If `TRUE` the adjusted trimmed means estimate is computed.
  The default is `FALSE`.

- data:

  a data frame containing the variables in the model. `data` should
  contain at least the following: a numeric outcome variable and a
  binary treatment variable (numeric, character or factor).

## Value

`tm` returns an object of class `tm`. The function `summary` is used to
obtain a summary of the results. The generic accessor function
`coefficients` extracts the regression coefficients with corresponding
p-values and 95% confidence intervals.

An object of class "`tm`" is a list containing the following components:

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

## Details

The trimmed means estimate is subject to two assumptions: the strong
MNAR assumption requires that all dropouts (unobserved outcome values)
are located in the fraction of the distribution that is trimmed away;
the location shift assumption requires the group variances of the full
sample to be equal. The adjusted trimmed means estimator relaxes the
latter, but assumes normally distributed outcomes. The adjustment is
performed on the group with the smallest dropout proportion.

The p-value and 95% confidence intervals for the trimmed means estimate
and the adjusted trimmed means estimate are obtained in a permutation
approach.

## Examples

``` r
set.seed(123456)
test_dat <- as.data.frame(cbind(c(rep(0, 500), rep(1, 500)),
                          c(sort(rnorm(500, 0, 1)), sort(rnorm(500, 1, 1.5))),
                          rbinom(1000, 2, 0.4), rnorm(1000, 0, 1)))
colnames(test_dat) <- c("TR", "Y", "U", "U2")
test_dat$Y[1:200] <- NA
# Note that we usually recommend setting n_perm to a larger value, e.g., 1000
tm_obj <- tm(formula= Y ~ TR + U + U2,
             GR = "TR", trF = 0.5, side = "LOW",
             n_perm = 100, adj_est = TRUE, data = test_dat)
print(tm_obj)
#> 
#> Call:
#> tm(formula = Y ~ TR + U + U2, GR = "TR", trF = 0.5, side = "LOW", 
#>     n_perm = 100, adj_est = TRUE, data = test_dat)
#> 
#> Coefficients:
#>        Estimate      Pval 95% CI LO 95% CI HI
#> TR     1.354308  0.009901  1.166162     1.542
#> U     -0.071352  0.108911 -0.163494     0.021
#> U2     0.055791  0.118812 -0.009835     0.121
#> TRadj  0.969812  0.009901  0.850012     1.090
#> 
summary(tm_obj)
#> 
#> Call:
#> tm(formula = Y ~ TR + U + U2, GR = "TR", trF = 0.5, side = "LOW", 
#>     n_perm = 100, adj_est = TRUE, data = test_dat)
#> 
#> 
#> Analysis details:
#> Treatment (1) group rescaled for adjusted TM estimate
#> 
#> Coefficients:
#>        Estimate      Pval 95% CI LO 95% CI HI
#> TR     1.354308  0.009901  1.166162     1.542
#> U     -0.071352  0.108911 -0.163494     0.021
#> U2     0.055791  0.118812 -0.009835     0.121
#> TRadj  0.969812  0.009901  0.850012     1.090
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
```

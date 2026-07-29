# Calculating Bias For Trimmed Mean Linear Models

`tm_bias` calculates the bias and the bias-adjusted estimate for a
trimmed means analysis
([`tm`](https://dea-hazewinkel.github.io/tmsens/reference/tm.md)) of a
given dataset, for a user-specified trimming fraction and dropout
spread. `tm_bias` calculates, under assumption of normally distributed
outcomes, the bias components resulting from violation of the location
shift assumption and violation of the strong MNAR assumption.

## Usage

``` r
tm_bias(
  formula,
  GR,
  trF = NULL,
  side = c("LOW", "HIGH"),
  spread_TG = "max_bias",
  spread_CG = "max_bias",
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

- spread_TG:

  a number between 0 and 1, specifying the dropout spread for the
  treatment group. `spread_TG` should be equal to or greater than the
  observed dropout proportion. If left unspecified, the worst-case
  scenario is assumed, in which dropout is located on the side of the
  distribution opposite from the one that is being trimmed
  (`spread_TG="max_bias"`).

- spread_CG:

  a number between 0 and 1, specifying the dropout spread for the
  comparator group. `spread_CG` should be equal to or greater than the
  observed dropout proportion. If left unspecified, the worst-case
  scenario is assumed, in which dropout is located on the side of the
  distribution opposite from the one that is being trimmed
  (`spread_CG="max_bias"`).

- data:

  a data frame containing the variables in the model. `data` should
  contain at least the following: a numeric outcome variable and a
  binary treatment variable (numeric, character or factor).

## Value

`tm_bias` returns an object of class `tm_bias`.

An object of class "`tm_bias`" is a list containing the following
components:

- call:

  the matched call

- bias_components:

  an array of bias components, including location shift assumption bias
  (LS), Strong MNAR bias in the treatment group (TG) and the comparator
  group (CG)

- total_bias:

  the sum of all bias components

- TM_estimate:

  the trimmed means estimate of the treatment effect

- bias_adj_TM_estimate:

  the bias adjusted trimmed means estimate

- analysis_details:

  the user-specified trimming fraction, trimming side, and dropout
  spread in the treatment (TG) and comparator groups (CG)

- observed_TG_SD:

  observed standard deviation of the treatment group (TG) outcome

- observed_CG_SD:

  observed standard deviation of the comparator group (CG) outcome

- inferred_TG_SD:

  inferred full sample standard deviation of the treatment group (TG)
  outcome

- max_bias_CG:

  an array of bias components, total bias, the bias adjusted estimate,
  and inferred full sample group standard deviations, calculated under
  the assumption of worst-case scenario dropout, with dropout in the
  comparator group (CG) on the opposite side of the distribution from
  the one that is being trimmed

- max_bias_TG:

  an array of bias components, total bias, the bias adjusted estimate,
  and inferred full sample group standard deviations, calculated under
  the assumption of worst-case scenario dropout, with dropout in the
  treatment group (TG) on the opposite side of the distribution from the
  one that is being trimmed

- groups:

  a named vector giving the values of the treatment variable that define
  the treatment (TG) and comparator (CG) groups

## Details

The trimmed means estimate is subject to two assumptions: the strong
MNAR assumption requires that all dropouts (unobserved outcome values)
are located in the fraction of the distribution that is trimmed away;
the location shift assumption requires the group variances of the full
sample to be equal. The bias resulting from the violation of either
assumption can be calculated under assumption of normally distributed
outcomes.

Obtaining the strong MNAR assumption bias requires an additional
assumption about the distribution of the dropouts: it is assumed that
the dropouts are spread homogeneously across the specified dropout
spread. For example, under lower value trimming (`side="LOW"`), and a
treatment group dropout spread of 0.6 (`spread_TG=0.6`), any value in
the bottom 60% of the treatment group outcome distribution is equally
likely to be missing.

The specified dropout spread for a given treatment group has
implications for the unobserved full sample variance that is inferred
from the observed data. For example, for an observed dropout of 0.4 and
an assumed dropout spread of 0.5, the inferred full sample variance will
be larger than for an assumed dropout spread of e.g., 0.8.

In addition to calculating the bias for a user-specified dropout spread,
`tm_bias` also calculates the maximal bias. For example, for lower value
trimming (`side="LOW"`), the worst-case scenario would involve lower
value dropout in the treatment group (TG) and higher value dropout in
the comparator group (CG), and vice versa. Bias components are
calculated for both scenarios. If the dropout spread (`spread_TG`,
`spread_CG`) is left unspecified for either treatment group, the
function will return only these quantities.

## Examples

``` r
test_dat <- as.data.frame(cbind(c(rep(0, 500), rep(1, 500)),
  c(sort(rnorm(500, 0, 1)), sort(rnorm(500, 1, 1.5)))))
colnames(test_dat) <- c("TR", "Y")
test_dat$Y[which(test_dat$TR == 0)[1:150]] <- NA
test_dat$Y[which(test_dat$TR == 1)[sample(seq(1, 400), 200, replace = FALSE)]] <- NA
tm_bias_obj <- tm_bias(formula = Y ~ TR, "TR", trF = 0.5,
                       side = "LOW", spread_TG = 0.4,
                       spread_CG = 0.6, data = test_dat)
print(tm_bias_obj)
#> 
#> Call:
#> tm_bias(formula = Y ~ TR, GR = "TR", trF = 0.5, side = "LOW", 
#>     spread_TG = 0.4, spread_CG = 0.6, data = test_dat)
#> 
#> Analysis details:
#>                                              
#>  50% trimming,                               
#>  under assumption of lower value dropout,    
#>  with 40% dropout spread in the TG group (1),
#>  and 60% dropout spread in the CG group (0)  
#> 
#> Bias components:
#>   Bias type                 Bias   
#>   LS                        1.41339
#>   Strong MNAR TG group (1)  0.00000
#>   Strong MNAR CG group (0)  0.01761
#> 
#> Total bias:
#>       
#>  1.431
#> 
#> TM estimate:
#>        
#>  0.9402
#> 
#> Bias adjusted TM estimate:
#>         
#>  -0.4908
#> 
#> Observed SDs:
#>                                           
#>  Observed TG group ( 1) SD for 40% dropout
#>  1.605                                    
#>                                           
#>  Observed CG group ( 0) SD for 30% dropout
#>  0.6966                                   
#> 
#> Inferred SDs:
#>                                                                                            
#>  Inferred TG group ( 1) full sample SD for 40% dropout in the lower 40% of the distribution
#>  2.47                                                                                      
#>                                                                                            
#>  Inferred CG group ( 0) full sample SD for 30% dropout in the lower 60% of the distribution
#>  0.6988                                                                                    
#> 
#> Bias under maximal violation of the strong MNAR assumption in the CG group (0):
#>                                 
#> Strong MNAR bias          0.6615
#> LS bias                   1.1792
#> Total bias                1.8407
#> Bias adjusted estimate   -0.9005
#> Inferred SD TG group (1)  2.4702
#> Inferred SD CG group (0)  0.9923
#> 
#> Bias under maximal violation of the strong MNAR assumption in the TG group (1):
#>                                 
#> Strong MNAR bias         -2.5105
#> LS bias                   1.1792
#> Total bias               -1.3313
#> Bias adjusted estimate    2.2715
#> Inferred SD TG group (1)  2.4702
#> Inferred SD CG group (0)  0.9923
```

# Changelog

## tmsens 0.4.0

- Fixed the adjusted trimmed means estimator for higher value trimming
  (`side = "HIGH"`), which produced incorrect estimates for outcomes not
  centred near zero
- The permutation p-value is now two-sided, so negative treatment
  effects are detected, and uses the `(b+1)/(n_perm+1)` estimator so it
  is never exactly zero
- Fixed an off-by-one in higher value trimming which removed one more
  observation per group than the specified trimming fraction
- Fixed `print.tm_bias()` erroring when the dropout spread was left
  unspecified
- [`tm()`](https://dea-hazewinkel.github.io/tmsens/reference/tm.md) and
  [`tm_bias()`](https://dea-hazewinkel.github.io/tmsens/reference/tm_bias.md)
  now work with any treatment variable name, not just `TR`
- `side` now defaults to `"LOW"` and invalid values give a clear error
- [`tm_bias()`](https://dea-hazewinkel.github.io/tmsens/reference/tm_bias.md)
  gains the documented `trF` default
- [`tm_bias()`](https://dea-hazewinkel.github.io/tmsens/reference/tm_bias.md)
  no longer errors when a treatment group has no dropout
- [`summary.tm()`](https://dea-hazewinkel.github.io/tmsens/reference/summary.tm.md)
  output now shows dropout and post-trimming sample sizes for both
  treatment groups
- [`tm()`](https://dea-hazewinkel.github.io/tmsens/reference/tm.md) no
  longer perturbs the random number generator state beyond the draws
  needed for the permutation test
- The trimming fraction and treatment variable are now validated up
  front with clear error messages
- Treatment groups are defined by factor level order, fixing numeric
  group codes such as 2 and 10 and respecting user-set reference levels
- The permutation test is around 9 times faster
- `tm_bias` objects gain a `groups` component giving the treatment and
  comparator group values
- Corrected the author name spelling and publication year in the package
  citation

## tmsens 0.3.1

CRAN release: 2024-08-29

- Minor fixes to improve the code, the documentation, and the package
  website

## tmsens 0.3.0

CRAN release: 2023-05-10

- Simplified the NAMESPACE so functions from other packages are not
  imported

- Modified the
  [`tm()`](https://dea-hazewinkel.github.io/tmsens/reference/tm.md)
  helpfile example to run faster

## tmsens 0.2.0

CRAN release: 2023-03-30

- Bump roxygen2 version

- Update maintainer email address

## tmsens 0.1.0

CRAN release: 2022-02-03

- Initial submission to CRAN

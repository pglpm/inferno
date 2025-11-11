# Calculate quantile width through batches

Modified from from
https://github.com/LaplacesDemonR/LaplacesDemon/blob/master/R/ESS.R

## Usage

``` r
funMCEI(x, fn, p = c(0.055, 0.945), ...)
```

## Arguments

- x:

  A matrix, rows being MC samples and columns being quantities whose
  MCSE is to be estimated.

## Value

Estimates of the MC standard error for each trace. Division by sqrt(N)
is already performed.

# Calculate MC effective sample size using LaplacesDemon's algorithm

Modified from from
https://github.com/LaplacesDemonR/LaplacesDemon/blob/master/R/ESS.R

## Usage

``` r
funESSLD(x)
```

## Arguments

- x:

  A matrix, rows being MC samples and columns being quantities whose
  MCSE is to be estimated.

## Value

Estimates of the effective sample size for each trace.

# Calculate MC standard error using LaplacesDemon's batch means

Modified from from
https://github.com/LaplacesDemonR/LaplacesDemon/blob/master/R/ESS.R

## Usage

``` r
funMCSELD(x)
```

## Arguments

- x:

  A matrix, rows being MC samples and columns being quantities whose
  MCSE is to be estimated.

## Value

Estimates of the MC standard error for each trace. Division by sqrt(N)
is already performed.

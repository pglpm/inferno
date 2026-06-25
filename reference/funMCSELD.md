# Calculate MC standard error using LaplacesDemon's batch means

Modified from from
https://github.com/LaplacesDemonR/LaplacesDemon/blob/master/R/ESS.R

## Usage

``` r
funMCSELD(x)
```

## Arguments

- x:

  matrix, each row being a "trace", that is a set of MC samples, whose
  MCSE is to be estimated.

## Value

MCSE estimates, one for each trace. Division by sqrt(N) is already
performed.

## Details

Used in 'util_combineYX()' in 'Pr()'.

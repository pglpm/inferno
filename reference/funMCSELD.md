# Calculate MC standard error using LaplacesDemon's batch means

This function gives a good approximation of the "true" standard
deviation in the case of independent samples. Multiply by `qnorm(x)` to
obtain the `x`-quantile.

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

Modified from
<https://github.com/LaplacesDemonR/LaplacesDemon/blob/master/R/ESS.R>.

Tested also on t-distributions with df=1.1 and Pareto with a=1.5 (mean
exists, variance infinite).

`sd() / sqrt(funESS3()` gives essentially identical results to
`funMCSELD()`, but it's 20 times slower.

Used in 'util_combineYX()' in 'Pr()'.

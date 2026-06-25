# Compute ESS

Modified from 'rstan'
<https://github.com/stan-dev/rstan/blob/develop/rstan/rstan/R/monitor.R>

## Usage

``` r
funESS3(x)
```

## Arguments

- x:

  Vector of MC samples.

## Value

Effective Sample Size.

## Details

Used in 'workerfun()' in 'learn()', and in 'funMCEQ()'.

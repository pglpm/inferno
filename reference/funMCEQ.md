# Calculate credibility quantiles on estimated quantile

Modified from Vehtari et al.

## Usage

``` r
funMCEQ(x, prob = c(0.055, 0.945), Qpair = pnorm(c(-1, 1)))
```

## Arguments

- x:

  A vector of MC samples

- prob::

  Quantile whose error intervalis being estimated

- Qpair::

  Lower and higher credibility-quantiles requested

## Value

Estimates lower and higher credibility-quantiles on estimated quantile

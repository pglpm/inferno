# Calculate credibility quantiles on estimated quantile

Calculates the lower and upper bound of a credibility interval, for
various quantiles of the empirical distribution of a vector of MC
samples.

## Usage

``` r
funMCEQ(x, prob = c(0.055, 0.945), Qpair = pnorm(c(-1, 1)))
```

## Arguments

- x:

  A vector of MC samples

- prob:

  numeric vector of probabilities: quantiles whose error interval is
  being estimated.

- Qpair:

  vector of length two (further elements are ignored): lower and higher
  credibility-quantiles requested. Default yields a credibility interval
  of 68%, or one nominal normal standard deviation.

## Value

A matrix with two rows and as many columns as elements in 'prob'. Forr
each column, the first and second row determine the lower and upper
bound of the credibility interval of width `Qpair[2] - Qpair[2]`.

## Details

Used in 'workerfun()' in 'learn()'

# Generate datapoints

This function generate datapoints according to the posterior probability
`Pr(Y | X, data)` calculated with
[`learn()`](https://pglpm.github.io/inferno/reference/learn.md), for the
variates specified in the argument `Y`, and conditional on the variate
values specified in the argument `X`. If `X` is omitted or `NULL`, then
the posterior probability `Pr(Y | data)` is used. Each variate in the
argument `X` can be specified either as a point-value `X = x` or as a
left-open interval `X ≤ x` or as a right-open interval `X ≥ x`, through
the argument `tails`.

## Usage

``` r
rPr(n, Ynames, X = NULL, learnt, tails = NULL, mcsamples = NULL)
```

## Arguments

- n:

  Positive integer: number of samples to draw.

- Ynames:

  Character vector: names of variates to draw jointly

- X:

  List or data.table or `NULL`: set of values of variates on which we
  want to condition the joint probability for `Y`. If `NULL` (default),
  no conditioning is made. Any rows beyond the first are discarded

- learnt:

  Either a character with the name of a directory or full path for a
  'learnt.rds' object, produced by the
  [`learn()`](https://pglpm.github.io/inferno/reference/learn.md)
  function, or such an object itself.

- tails:

  Named vector or list, or `NULL` (default). The names must match some
  or all of the variates in arguments `X`. For variates in this list,
  the probability conditional is understood in an semi-open interval
  sense: `X ≤ x` or `X ≥ x`, an so on. See analogous argument in
  [`Pr()`](https://pglpm.github.io/inferno/reference/Pr.md).

- mcsamples:

  Vector of integers, or `'all'`, or `NULL` (default): which Monte Carlo
  samples calculated by the
  [`learn()`](https://pglpm.github.io/inferno/reference/learn.md)
  function should be used to draw the variate values. The default is to
  choose a random subset if `n` is smaller than their number, otherwise
  to recycle them as necessary.

## Value

A data frame of joint draws of the variates `Ynames` from the posterior
distribution, conditional on `X`. The row names of the data frame report
the Monte Carlo sample (from
[`learn()`](https://pglpm.github.io/inferno/reference/learn.md)) used
for that draw, and the total number of draws from that sample so far.

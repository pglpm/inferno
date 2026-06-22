# Generate datapoints

This function generate datapoints according to the posterior probability
`Pr(Y | X, data)`, for the variates specified in the argument `Y`, and
conditional on the variate values specified in the argument `X`. It is
somewhat analogous to the `r`-variants of R distribution functions, such
as [`stats::rnorm()`](https://rdrr.io/r/stats/Normal.html). If `X` is
omitted or `NULL`, then the posterior probability `Pr(Y | data)` is
used. Each variate in the argument `X` can be specified either as a
point-value `X = x` or as a left-open interval `X ≤ x` or as a
right-open interval `X ≥ x`, through the argument `tails`.

## Usage

``` r
rPr(
  n,
  Ynames,
  X = NULL,
  learnt,
  tails = NULL,
  mcsamples = NULL,
  parallel = NULL
)
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
  [`learn()`](https://pglpm.github.io/prova/reference/learn.md)
  function, or such an object itself.

- tails:

  Named vector or list, or `NULL` (default). The names must match some
  or all of the variates in arguments `X`. For variates in this list,
  the probability conditional is understood in an semi-open interval
  sense: `X ≤ x` or `X ≥ x`, an so on. See analogous argument in
  [`Pr()`](https://pglpm.github.io/prova/reference/Pr.md).

- mcsamples:

  Vector of integers, or `'all'`, or `NULL` (default): which Monte Carlo
  samples calculated by the
  [`learn()`](https://pglpm.github.io/prova/reference/learn.md) function
  should be used to draw the variate values. The default is to choose a
  random subset if `n` is smaller than their number, otherwise to
  recycle them as necessary.

## Value

A data frame of joint draws of the variates `Ynames` from the posterior
distribution, conditional on `X`. The row names of the data frame report
the Monte Carlo sample (from
[`learn()`](https://pglpm.github.io/prova/reference/learn.md)) used for
that draw, and the total number of draws from that sample so far.

## See also

[`learn()`](https://pglpm.github.io/prova/reference/learn.md), which
generates the `learnt` objects required by
[`qPr()`](https://pglpm.github.io/prova/reference/qPr.md).

[`Pr()`](https://pglpm.github.io/prova/reference/Pr.md) to calculate
joint and conditional probabilities.

[`qPr()`](https://pglpm.github.io/prova/reference/qPr.md) to calculate
quantiles.

## Examples

``` r
## Load the example `learnt` object included in the package
learnt <- learntExample

## ## Example 1:
## Generate 10 values of the 'species' variate,
## according to the frequency distribution estimated from the data

datapoints <- rPr(
  n = 10,
  Ynames = 'species',
  learnt = learnt
)

c(datapoints)
#> $species
#>  [1] "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"    "Gentoo"    "Chinstrap"
#>  [7] "Gentoo"    "Adelie"    "Adelie"    "Adelie"   
#> 


## ## Example 2:
## Generate 5 joint values of the 'species' and 'bill_len' variates.

datapoints <- rPr(
  n = 5,
  Ynames = c('species', 'bill_len'),
  learnt = learnt
)

print(datapoints, row.names = FALSE) ## row names give MCMC information
#>    species bill_len
#>     Adelie     36.0
#>     Adelie     38.4
#>     Adelie     35.8
#>  Chinstrap     42.1
#>  Chinstrap     48.6


## ## Example 3:
## Generate 5 values of the 'species' variate,
## for the subpopulation of penguins having bill length shorter than 40 mm

datapoints <- rPr(
  n = 5,
  Ynames = 'species',
  X = data.frame(bill_len = 40),
  tails = list(bill_len = -1),
  learnt = learnt
)

c(datapoints)
#> $species
#> [1] "Adelie" "Adelie" "Adelie" "Adelie" "Adelie"
#> 
```

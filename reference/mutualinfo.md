# Calculate mutual information between groups of joint variates

This function calculates various entropic information measures of two
variates (each variate may consist of joint variates): the mutual
information, the conditional entropies, and the entropies.

## Usage

``` r
mutualinfo(
  Y1names,
  Y2names = NULL,
  X = NULL,
  learnt,
  tails = NULL,
  n = NULL,
  unit = "Sh",
  parallel = TRUE,
  silent = FALSE
)
```

## Arguments

- Y1names:

  Character vector: first group of joint variates

- Y2names:

  Character vector or `NULL`: second group of joint variates

- X:

  Matrix or data.frame or `NULL`: values of some variates conditional on
  which we want the probabilities.

- learnt:

  Either a character with the name of a directory or full path for an
  'learnt.rds' object, or such an object itself.

- tails:

  Named vector or list, or `NULL` (default). The names must match some
  or all of the variates in arguments `X`. For variates in this list,
  the probability conditional is understood in an semi-open interval
  sense: `X ≤ x` or `X ≥ x`, an so on. See analogous argument in
  [`Pr()`](https://pglpm.github.io/prova/reference/Pr.md).

- n:

  Integer or `NULL` (default): number of samples from which to
  approximately calculate the mutual information. Default as many as
  Monte Carlo samples in `learnt`.

- unit:

  Either one of 'Sh' for *shannon* (default), 'Hart' for *hartley*,
  'nat' for *natural unit*, or a positive real indicating the base of
  the logarithms to be used.

- parallel:

  Logical or positive integer or cluster object. `TRUE` (default): use
  roughly half of available cores; `FALSE`: use serial computation;
  integer: use this many cores. It can also be a cluster object
  previously created with
  [`parallel::makeCluster()`](https://rdrr.io/r/parallel/makeCluster.html);
  in this case the parallel computation will use this object.

- silent:

  Logical: give warnings or updates in the computation?

## Value

A list consisting of the following elements:

- `MI`, a vector of `value` and `accuracy`: the mutual information
  between (joint) variates `Y1names` and (joint) variates `Y2names`.

- `CondEn12`, `CondEn21`, vectors of `value` and `accuracy`: the
  conditional entropy of the first variate given the second, and vice
  versa.

- `En1`, `En2`, vectors of `value` and `accuracy`: the (differential)
  entropies of the first and second variates.

- `MI.rGauss`, a vector of `value` and `accuracy`: the absolute value of
  the Pearson correlation coefficient of a *multivariate Gaussian
  distribution* having mutual information `MI` (the two are related by
  `MI = -log(1 - MI.rGauss^2)/2`); it may provide a vague intuition for
  the `MI` value for people more familiar with Pearson's correlation,
  but should be taken with a grain of salt.

- `unit`, `Y1names`, `Y1names`: same as the input arguments, included
  for the user's convenience.

## See also

[`Pr()`](https://pglpm.github.io/prova/reference/Pr.md) to calculate
probabilities and their variability.

[`learn()`](https://pglpm.github.io/prova/reference/learn.md), which
generates the `learnt` objects required by `mutualinfo()`.

## Examples

``` r
## Load the example `learnt` object calculated from the "penguins" dataset;
## variates: 'species' and 'bill_len'
learnt <- learntExample

## mutual information between variates 'species' and 'bill_len'
MI <- mutualinfo(Y1names = 'species', Y2names = 'bill_len',
  learnt = learnt, parallel = 1)
#> Registered socket cluster with 1 nodes on host ‘localhost’.
#> Closing connections to cores.

paste0(MI$MI, ' ', MI$unit, collapse = ' +/- ')
#> [1] "0.699139790987302 Sh +/- 0.053 Sh"

## Shannon entropy of variate 'species'
paste0(MI$En1, ' ', MI$unit, collapse = ' +/- ')
#> [1] "1.5404796693525 Sh +/- 0.029 Sh"


## Shannon entropy of variate 'species',
## conditional on a bill length of 30 mm:
entr <- mutualinfo(
  Y1names = 'species',
  X = data.frame(bill_len = 30),
  learnt = learnt, parallel = 1
)
#> Registered socket cluster with 1 nodes on host ‘localhost’.
#> Closing connections to cores.

paste0(entr$En1, ' ', entr$unit, collapse = ' +/- ')
#> [1] "0.440800870784225 Sh +/- 0.081 Sh"

## the entropy is now lower; indeed a penguin with a short bill length
## is more probably of the 'Adelie' species:
probs <- Pr(
  Y = data.frame(species = c('Adelie', 'Gentoo', 'Chinstrap')),
  X = data.frame(bill_len = 30),
  learnt = learnt, parallel = 1
)
#> Registered socket cluster with 1 nodes on host ‘localhost’.
#> Closing connections to cores.

probs$values
#>            X
#> Y                   30
#>   Adelie    0.92985096
#>   Gentoo    0.03626731
#>   Chinstrap 0.03388172

```

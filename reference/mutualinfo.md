# Calculate mutual information between groups of joint variates

This function calculates various entropic information measures of two
variates (each variate may consist of joint variates): the mutual
information, the conditional entropies, and the entropies.

## Usage

``` r
mutualinfo(
  Y1names,
  Y2names,
  X = NULL,
  learnt,
  tails = NULL,
  n = NULL,
  unit = "Sh",
  parallel = NULL,
  silent = FALSE
)
```

## Arguments

- Y1names:

  Character vector: first group of joint variates

- Y2names:

  Character vector or NULL: second group of joint variates

- X:

  Matrix or data.frame or NULL: values of some variates conditional on
  which we want the probabilities.

- learnt:

  Either a character with the name of a directory or full path for an
  'learnt.rds' object, or such an object itself.

- tails:

  Named vector or list, or `NULL` (default). The names must match some
  or all of the variates in arguments `X`. For variates in this list,
  the probability conditional is understood in an semi-open interval
  sense: `X ≤ x` or `X ≥ x`, an so on. See analogous argument in
  [`Pr()`](https://pglpm.github.io/inferno/reference/Pr.md).

- n:

  Integer or `NULL` (default): number of samples from which to
  approximately calculate the mutual information. Default as many as
  Monte Carlo samples in `learnt`.

- unit:

  Either one of 'Sh' for *shannon* (default), 'Hart' for *hartley*,
  'nat' for *natural unit*, or a positive real indicating the base of
  the logarithms to be used.

- parallel:

  Logical or positive integer or cluster object. `TRUE`: use roughly
  half of available cores; `FALSE`: use serial computation; integer: use
  this many cores. It can also be a cluster object previously created
  with
  [`parallel::makeCluster()`](https://rdrr.io/r/parallel/makeCluster.html);
  in this case the parallel computation will use this object.

- silent:

  Logical: give warnings or updates in the computation?

## Value

A list consisting of the elements `MI`, `CondEn12`, `CondEn21`, `En1`,
`En2`, `MI.rGauss`, `unit`, `Y1names`, `Y1names`. All elements except
`unit`, `Y1names`, `Y2names` are a vector of `value` and `accuracy`.
Element `MI` is the mutual information between (joint) variates
`Y1names` and (joint) variates `Y2names`. Element`CondEn12` is the
conditional entropy of the first variate given the second, and vice
versa for `CondEn21`. Elements `En1` and `En1` are the (differential)
entropies of the first and second variates. Elements `unit`, `Y1names`,
`Y2names` are identical to the same inputs. Element `MI.rGauss` is the
absolute value of the Pearson correlation coefficient of a *multivariate
Gaussian distribution* having mutual information `MI` (the two are
related by `MI = -log(1 - MI.rGauss^2)/2`); it may provide a vague
intuition for the `MI` value for people more familiar with Pearson's
correlation, but should be taken with a grain of salt.

# Calculate posterior probabilities

This function calculates the posterior probability `Pr(Y | X, data)`,
where `Y` and `X` are two (non overlapping) sets of joint variate
values. If `X` is omitted or `NULL`, then the posterior probability
`Pr(Y | data)` is calculated. The function also gives quantiles about
the possible variability of the probability `Pr(Y | X, newdata, data)`
that we could have if more learning data were provided, as well as a
number of samples of the possible values of such probabilities. If
several joint values are given for `Y` or `X`, the function will create
a 2D grid of results for all possible compbinations of the given `Y` and
`X` values. This function also allows for base-rate or other
prior-probability corrections: If a prior (for instance, a base rate)
for `Y` is given, the function will calculate the
`Pr(Y | X, data, prior)` from `Pr(X | Y, data)` and the prior by means
of Bayes's theorem. Each variate in each argument `Y`, `X` can be
specified either as a point-value `Y = y` or as a left-open interval
`Y ≤ y` or as a right-open interval `Y ≥ y`, through the argument
`tails`.

## Usage

``` r
Pr(
  Y,
  X = NULL,
  learnt,
  tails = NULL,
  priorY = NULL,
  nsamples = "all",
  quantiles = c(0.055, 0.25, 0.75, 0.945),
  parallel = NULL,
  silent = FALSE,
  keepYX = TRUE
)
```

## Arguments

- Y:

  Matrix or data.table: set of values of variates of which we want the
  joint probability of. One variate per column, one set of values per
  row.

- X:

  Matrix or data.table or `NULL` (default): set of values of variates on
  which we want to condition the joint probability of `Y`. If `NULL`, no
  conditioning is made (except for conditioning on the learning dataset
  and prior assumptions). One variate per column, one set of values per
  row.

- learnt:

  Either a character with the name of a directory or full path for a
  'learnt.rds' object, produced by the
  [`learn()`](https://pglpm.github.io/inferno/reference/learn.md)
  function, or such an object itself.

- tails:

  Named vector or list, or `NULL` (default). The names must match some
  or all of the variates in arguments `Y` and `X`. For variates in this
  list, the probability arguments are understood in an semi-open
  interval sense: `Y ≤ y` or `Y ≥ y`, an so on. This is true for
  variates on the left and on the right of the conditional sign `\|`. A
  left-open interval `Y ≤ y` is indicated by the values `'<='` or
  `'left'` or `-1`; a right-open interval `Y ≥ y` is indicated by the
  values `'>='` or `'right'` or `+1`. Values `NULL`, `'=='`, `0`
  indicate that a point value `Y = y` (not an interval) should be
  calculated. **NB**: the semi-open intervals *always* include the given
  value; this is important for ordinal or rounded variates. For
  instance, if `Y` is an integer variate, then to calculate `P(Y < 3)`
  you should require `P(Y <= 2)`; for this reason we also have that
  `P(Y <= 2)` and `P(Y >= 2)` generally add up to *more* than 1.

- priorY:

  Numeric vector with the same length as the rows of `Y`, or `TRUE`, or
  `NULL` (default): prior probabilities or base rates for the `Y`
  values. If `TRUE`, the prior probabilities are assumed to be all
  equal.

- nsamples:

  Integer or `NULL` or `'all'` (default): desired number of samples of
  the variability of the probability for `Y`. If `NULL`, no samples are
  reported. If `'all'` (or `Inf`), all samples obtained by the
  [`learn()`](https://pglpm.github.io/inferno/reference/learn.md)
  function are used.

- quantiles:

  Numeric vector, between 0 and 1, or `NULL`: desired quantiles of the
  variability of the probability for `Y`. Default
  `c(0.055, 0.25, 0.75, 0.945)`, that is, the 5.5%, 25%, 75%, 94.5%
  quantiles (these are typical quantile values in the Bayesian
  literature: they give 50% and 89% credibility intervals, which
  correspond to 1 shannons and 0.5 shannons of uncertainty). If `NULL`,
  no quantiles are calculated.

- parallel:

  Logical or positive integer or cluster object. `TRUE`: use roughly
  half of available cores; `FALSE`: use serial computation; integer: use
  this many cores. It can also be a cluster object previously created
  with
  [`parallel::makeCluster()`](https://rdrr.io/r/parallel/makeCluster.html);
  in this case the parallel computation will use this object.

- silent:

  Logical, default `FALSE`: give warnings or updates in the computation?

- keepYX:

  Logical, default `TRUE`: keep a copy of the `Y` and `X` arguments in
  the output? This is used for the plot method.

## Value

A list of class `probability`, consisting of the elements `values`,
`quantiles` (possibly `NULL`), `samples` (possibly `NULL`),
`values.MCaccuracy`, `quantiles.MCaccuracy` (possibly `NULL`), `Y`, `X`.
Element `values`: a matrix with the probabilities
P(Y\|X,data,assumptions), for all combinations of values of `Y` (rows)
and `X` (columns). Element `quantiles`: an array with the variability
quantiles (3rd dimension of the array) for such probabilities. Element
`samples`: an array with the variability samples (3rd dimension of the
array) for such probabilities. Elements `values.MCaccuracy` and
`quantiles.MCaccuracy`: arrays with the numerical accuracies (roughly
speaking a standard deviation) of the Monte Carlo calculations for the
`values` and `quantiles` elements. Elements `Y`, `X`: copies of the `Y`
and `X` arguments.

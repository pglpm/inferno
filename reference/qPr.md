# Calculate quantiles

This function calculates the quantiles of `Pr(Y | X, data)` at specified
cumulative-probability levels (that is, the values of `Y` having
specified cumulative probabilities), as well as the variability of those
quantiles if more learning data were provided. It is somewhat analogous
to the `q`-variants of R distribution functions, such as
[`stats::qnorm()`](https://rdrr.io/r/stats/Normal.html). The variability
can be expressed in the form of quantiles, samples, or both, as in the
[`Pr()`](https://pglpm.github.io/prova/reference/Pr.md) function. If
several joint values are given for the probability levels and for `X`,
the function creates a 2D grid of results for all possible combinations
of the given probability levels and `X` values. Each variate in the
argument `X` can be specified either as a point-value `X = x` or as a
left-open interval `X ≤ x` or as a right-open interval `X ≥ x`, through
the argument `tails`.

## Usage

``` r
qPr(
  p = c(0.055, 0.5, 0.945),
  Yname,
  X = NULL,
  learnt,
  tails = NULL,
  priorY = NULL,
  nsamples = "all",
  quantiles = c(0.055, 0.5, 0.945),
  parallel = NULL,
  silent = FALSE,
  keepYX = TRUE,
  tol = .Machine$double.eps * 10
)
```

## Arguments

- p:

  Numeric vector of probability levels. Default: `c(0.055, 0.5, 0.945)`.

- Yname:

  Character vector: name of variate whose quantiles will be computed.

- X:

  Matrix or data.table or `NULL` (default): set of values of variates on
  which we want to condition. If `NULL`, no conditioning is made (except
  for conditioning on the learning dataset and prior assumptions). One
  variate per column, one set of values per row.

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

- priorY:

  Numeric vector with the same length as the rows of `Y`, or `TRUE`, or
  `NULL` (default): prior probabilities or base rates for the `Y`
  values. If `TRUE`, the prior probabilities are assumed to be all
  equal. For the moment only the value `NULL` is accepted.

- nsamples:

  Integer or `NULL` or `'all'` (default): desired number of samples of
  the variability of the quantile for `Y`. If `NULL`, no samples are
  reported. If `'all'` (or `Inf`), all samples obtained by the
  [`learn()`](https://pglpm.github.io/prova/reference/learn.md) function
  are used.

- quantiles:

  Numeric vector, between 0 and 1, or `NULL`: desired quantiles of the
  variability of the quantile for `Y`. Default
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

- tol:

  numeric positive: tolerance in the calculation of quantiles. Default:
  `.Machine$double.eps * 10` (typically `2.22045e-15`).

## Value

A list of the elements `values`, `quantiles` (possibly `NULL`),
`samples` (possibly `NULL`), `Y`, `X`. Element `values`: a matrix with
the requested `Y`-quantiles conditional on the requested `X`-values, for
all combinations of `p` (rows) and `X` (columns). Element `quantiles`:
an array with the variability quantiles (3rd dimension of the array).
Element `samples`: an array with the variability samples (3rd dimension
of the array). Elements `Y`, `X`: copies of the `p` and `X` arguments.

## References

P.G.L. Porta Mana: *What's special about 89% credibility intervals?*.
2025 <doi:10.5281/zenodo.17072199>.

## See also

[`learn()`](https://pglpm.github.io/prova/reference/learn.md), which
generates the `learnt` objects required by `qPr()`.

[`Pr()`](https://pglpm.github.io/prova/reference/Pr.md) to calculate
joint and conditional probabilities.

[`rPr()`](https://pglpm.github.io/prova/reference/rPr.md) to generate
datapoints.

## Examples

``` r

## Load the example `learnt` object included in the package
learnt <- learntExample

## ## Example 1:
## Calculate the 5.5%-, 50%-, and 94.5%-quantiles for the variate "bill lengt",
## that is, the bill lengths having such cumulative probabilities

quants <- qPr(
  Yname = 'bill_len',
  learnt = learnt, parallel = 1
)
#> Registered socket cluster with 1 nodes on host ‘localhost’.
#> Closing connections to cores.

## display the quantile values
quants$values
#>         X
#> bill_len [,1]
#>    0.055 35.7
#>    0.5   44.3
#>    0.945 52.1

## verify these values using Pr():
probs <- Pr(
  Y = data.frame(bill_len = c(quants$values)),
  tails = list(bill_len = -1),
  learnt = learnt, parallel = 1)
#> Registered socket cluster with 1 nodes on host ‘localhost’.
#> Closing connections to cores.

## the cumulative probabilities are indeed 0.055, 0.5, 0.945 within numerical error:
probs$values
#>       X
#> Y            [,1]
#>   35.7 0.05541489
#>   44.3 0.50269977
#>   52.1 0.94784616

## display the uncertainty about the quantiles
quants$quantiles
#> , ,  = 5.5%
#> 
#>         X
#> bill_len [,1]
#>    0.055 35.3
#>    0.5   43.3
#>    0.945 51.5
#> 
#> , ,  = 50%
#> 
#>         X
#> bill_len [,1]
#>    0.055 35.7
#>    0.5   44.3
#>    0.945 52.1
#> 
#> , ,  = 94.5%
#> 
#>         X
#> bill_len [,1]
#>    0.055 36.2
#>    0.5   45.1
#>    0.945 52.8
#> 
```

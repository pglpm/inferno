# Calculate mutual information between groups of joint variates

This function calculates various entropic information measures between
two grops of joint variates: the mutual information, the conditional
entropies, and the entropies.

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
  sense: \\X \le x\\ or \\X \ge x\\, an so on. See analogous argument in
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
  the Pearson correlation coefficient \\r\\ of a *multivariate Gaussian
  distribution* having mutual information `MI`; the two are related by
  \\\mathrm{MI} = -\ln(1 - r^2)/2\\. It may provide a vague intuition
  for the `MI` value for people more familiar with Pearson's
  correlation, but should be taken with a grain of salt.

- `unit`, `Y1names`, `Y1names`: same as the input arguments, included
  for the user's convenience.

## Details

If \\Y_1\\ and \\Y_2\\ are two variates, each of which can be a joint
variate such as \\Y_1 = (Y\_{1,1}, Y\_{1,2}, \dotsc)\\, and \\X\\ a
third, also possibly join, variate, then the mutual information
\\\mathit{MI}\\ between \\Y_1\\ and \\Y_2\\, conditional on \\X = x\\,
is given by \$\$\mathit{MI}(Y_1, Y_2 \vert X = x) \coloneqq \sum\_{y_1,
y_2} \mathrm{Pr}(Y_1 = y_1, Y_2 = y_2 \vert X = x, \text{data})
\log_2\frac{ \mathrm{Pr}(Y_1 = y_1, Y_2 = y_2 \vert X = x, \text{data})
}{ \mathrm{Pr}(Y_1 = y_1 \vert X = x, \text{data}) \cdot \mathrm{Pr}(Y_2
= y_2 \vert X = x, \text{data}) } \\ \mathrm{Sh} \$\$ an expression
which can also be written in several other equivalent ways. It is an
information-theoretic measure of association that is model-free, that
is, does not depend on assumptions such as linearity, gaussianity, and
similar. See
[`vignette('mutualinfo')`](https://pglpm.github.io/prova/articles/mutualinfo.md)
for discussion and example uses, and also the "References" section. If
\\Y_1, Y_2\\ are *jointly gaussian variates*, then there is a
mathematical correspondence between their mutual information and their
Pearson correlation coefficient; see output `MI.rGauss` in the "Value"
section.

The conditional entropy of \\Y_1\\ with respect to \\Y_2\\, conditional
on \\X = x\\, is given by \$\$\mathit{CondEn12}(Y_1, Y_2 \vert X = x)
\coloneqq -\sum\_{y_1, y_2} \mathrm{Pr}(Y_1 = y_1 \vert Y_2 = y_2, X =
x, \text{data}) \log_2 \mathrm{Pr}(Y_1 = y_1 \vert Y_2 = y_2, X = x,
\text{data}) \cdot \mathrm{Pr}(Y_2 = y_2 \vert X = x, \text{data}) \\
\mathrm{Sh} \$\$

The (differential) entropy of \\Y_1\\, conditional on \\X = x\\, is
given by \$\$\mathit{En1}(Y_1 \vert X = x) \coloneqq -\sum\_{y_1}
\mathrm{Pr}(Y_1 = y_1 \vert X = x, \text{data}) \log_2 \mathrm{Pr}(Y_1 =
y_1 \vert X = x, \text{data}) \\ \mathrm{Sh} \$\$

see "References" section for discussions about entropy and conditional
entropy.

The function `mutualinfo()` calculates the quantities above for the
joint variates specified in the arguments `Y1names` and `Y2names`,
conditional on the values of the variates specified in the data frame
`X`. If `X` is omitted or `NULL`, then the posterior probabilities
\\\mathrm{Pr}(Y_1 \| \text{data})\\ etc. are used. Each variate in the
argument `X` can be specified either as a point-value \\X = x\\ or as a
left-open interval \\X \le x\\ or as a right-open interval \\X \ge x\\,
through the argument `tails`.

The computation of these quantities is done via Monte Carlo integration,
using the samples produced by the
[`learn()`](https://pglpm.github.io/prova/reference/learn.md) function.
The present function also output the numerical error associated with
this computation.

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

# \donttest{
## the entropy is now lower; indeed a penguin with a short bill length
## is most probably of the 'Adelie' species:
probs <- Pr(
  Y = data.frame(species = c('Adelie', 'Gentoo', 'Chinstrap')),
  X = data.frame(bill_len = 30),
  learnt = learnt, parallel = 1
)
#> Registered socket cluster with 1 nodes on host ‘localhost’.
#> Closing connections to cores.

print(probs)
#> , , |bill_len = 30
#> 
#>            prob. & vrb.
#> species     value   Q5.5%   Q25%  Q75% Q94.5%
#>   Adelie    0.930 0.62600 0.9430 0.992  0.998
#>   Gentoo    0.036 0.00021 0.0025 0.022  0.150
#>   Chinstrap 0.034 0.00019 0.0016 0.022  0.170
#> 
# }
```

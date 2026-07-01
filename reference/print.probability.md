# Print an object of class "probability"

This [`base::print()`](https://rdrr.io/r/base/print.html) method is a
utility to display selected elements of a "probability" object obtained
with [`Pr()`](https://pglpm.github.io/prova/reference/Pr.md); typically
its posterior probabilies (element `$values`) and their variabilities
(element `$quantiles`). If the `Y` or `X` variates are joint variates,
this method also allow to display only selected values of them

## Usage

``` r
# S3 method for class 'probability'
print(x, elements = NULL, subset = NULL, digits = TRUE, ...)
```

## Arguments

- x:

  Object of class "probability", obtained with
  [`Pr()`](https://pglpm.github.io/prova/reference/Pr.md).

- elements:

  character or integer vector, or `NULL` (default): elements of the
  "probability" object to display. The syntax is the same as with
  [`[`](https://rdrr.io/r/base/Extract.html). If `NULL`, the elements
  `$values` and `$quantiles` are displayed together in a special way.

- subset:

  Named list or named vector: which variate values to display. For the
  variates corresponding to the names in this list, only the vector of
  values corresponding to that variate is displayed.

- digits:

  positive number or `NULL` or `TRUE` (default): minimal number of
  significant digits, see
  [`base::print.default()`](https://rdrr.io/r/base/print.default.html).
  If value is `TRUE`, then the significant digits for elements `$values`
  and `$quantiles` are determined from their respective
  `$values.MCaccuracy` and `$quantiles.MCaccuracy` elements of the
  `probability` object, see
  [`Pr()`](https://pglpm.github.io/prova/reference/Pr.md); whereas
  `$samples` elements use 2 significant digits.

- ...:

  Other parameters to be passed to
  [`base::print()`](https://rdrr.io/r/base/print.html).

## Value

Its `x` argument, [invisibly](https://rdrr.io/r/base/invisible.html);
see [`base::print()`](https://rdrr.io/r/base/print.html).

## See also

[`Pr()`](https://pglpm.github.io/prova/reference/Pr.md) to calculate
posterior probabilities and quantiles.

[`plot.probability()`](https://pglpm.github.io/prova/reference/plot.probability.md)
to plot probabilities and quantiles calculated by \`Pr()'.
[`hist.probability()`](https://pglpm.github.io/prova/reference/hist.probability.md)
to plot the variability of the probabilities as a distribution.

## Examples

``` r
## Load the example `learnt` object calculated from the "penguins" dataset;
## variates: 'species' and 'bill_len'
learnt <- learntExample

## Calculate the 3 x 2 probabilities for the 3 species
## given bill-lengths of 43 mm and 44 mm

Y <- data.frame(species = c('Adelie', 'Chinstrap', 'Gentoo'))
X <- data.frame(bill_len = c(43, 44))

probs <- Pr(Y = Y, X = X, learnt = learnt, parallel = 1)

## display the values and variabilities of these probabilities
print(probs)
#> , , |bill_len = 43
#> 
#>            prob. & vrb.
#> species     value Q5.5%  Q25%  Q75% Q94.5%
#>   Adelie    0.465 0.367 0.421 0.513  0.568
#>   Chinstrap 0.146 0.081 0.117 0.172  0.219
#>   Gentoo    0.389 0.299 0.348 0.430  0.481
#> 
#> , , |bill_len = 44
#> 
#>            prob. & vrb.
#> species     value Q5.5%  Q25%  Q75% Q94.5%
#>   Adelie    0.222 0.140 0.187 0.255  0.307
#>   Chinstrap 0.205 0.120 0.172 0.243  0.297
#>   Gentoo    0.572 0.467 0.527 0.620  0.672
#> 

## diplay 'values' only, and only for the species value 'Gentoo'
print(probs, elements = 'values', subset = list(species = 'Gentoo'))
#> $values
#>         |bill_len
#> species     43    44
#>   Gentoo 0.389 0.572
#> 
```

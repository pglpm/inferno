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
print(x, elements = NULL, subset = NULL, digits = 2, ...)
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

  positive number, default 2: minimal number of significant digits, see
  [`base::print.default()`](https://rdrr.io/r/base/print.default.html).

- ...:

  Other parameters to be passed to
  [`base::print()`](https://rdrr.io/r/base/print.html).

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
#> Registered socket cluster with 1 nodes on host ‘localhost’.
#> Closing connections to cores.

## display the values and variabilities of these probabilities
print(probs)
#> , , |bill_len = 43
#> 
#>            prob. & vrb.
#> species     value Q5.5% Q25% Q75% Q94.5%
#>   Adelie     0.46 0.367 0.42 0.51   0.57
#>   Chinstrap  0.15 0.081 0.12 0.17   0.22
#>   Gentoo     0.39 0.299 0.35 0.43   0.48
#> 
#> , , |bill_len = 44
#> 
#>            prob. & vrb.
#> species     value Q5.5% Q25% Q75% Q94.5%
#>   Adelie     0.22  0.14 0.19 0.26   0.31
#>   Chinstrap  0.21  0.12 0.17 0.24   0.30
#>   Gentoo     0.57  0.47 0.53 0.62   0.67
#> 

## diplay 'values' only, and only for the species value 'Gentoo'
print(probs, elements = 'values', subset = list(species = 'Gentoo'))
#> $values
#>         |bill_len
#> species    43   44
#>   Gentoo 0.39 0.57
#> 
```

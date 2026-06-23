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
print(x, elements = NULL, subset = NULL, ...)
```

## Arguments

- x:

  Object of class "probability", obtained with
  [`Pr()`](https://pglpm.github.io/prova/reference/Pr.md).

- elements:

  character or integer vector, or `NULL` (default): elements of the
  "probability" object to display. The syntax is the same as
  \[base::\[\]. If `NULL`, display elements `$values` and `$quantiles`.

- subset:

  Named list or named vector: which variate values to display. For the
  variates corresponding to the names in this list, only the vector of
  values corresponding to that variate is displayed.

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
#> $values
#>            bill_len
#> species            43        44
#>   Adelie    0.4647433 0.2223224
#>   Chinstrap 0.1458345 0.2054491
#>   Gentoo    0.3894222 0.5722285
#> 
#> $quantiles
#> , , Q = 5.5%
#> 
#>            bill_len
#> species             43        44
#>   Adelie    0.36692485 0.1447719
#>   Chinstrap 0.08096214 0.1197566
#>   Gentoo    0.29852092 0.4671262
#> 
#> , , Q = 25%
#> 
#>            bill_len
#> species            43        44
#>   Adelie    0.4211860 0.1865739
#>   Chinstrap 0.1167195 0.1716849
#>   Gentoo    0.3476311 0.5269590
#> 
#> , , Q = 75%
#> 
#>            bill_len
#> species            43        44
#>   Adelie    0.5131288 0.2551308
#>   Chinstrap 0.1722602 0.2429360
#>   Gentoo    0.4302261 0.6203536
#> 
#> , , Q = 94.5%
#> 
#>            bill_len
#> species            43        44
#>   Adelie    0.5678666 0.3069388
#>   Chinstrap 0.2190872 0.2965559
#>   Gentoo    0.4811995 0.6718201
#> 
#> 

## diplay 'values' only, and only for the species value 'Gentoo'
print(probs, elements = 'values', subset = list(species = 'Gentoo'))
#> $values
#>         bill_len
#> species         43        44
#>   Gentoo 0.3894222 0.5722285
#> 
```

# Plot pairs of quantiles

Utility function to plot pairs of quantiles obtained with
[`Pr()`](https://pglpm.github.io/prova/reference/Pr.md).

## Usage

``` r
plotquantiles(
  x,
  y,
  xdomain = NULL,
  alpha.f = 0.25,
  col = palette(),
  border = NA,
  type = "n",
  ...
)
```

## Arguments

- x:

  Numeric or character: vector of x-coordinates. See
  [`flexiplot()`](https://pglpm.github.io/prova/reference/flexiplot.md).

- y:

  Numeric: a matrix having as many rows as `x` and an even number of
  columns, with one column per quantile. Typically these quantiles have
  been obtained with
  [`Pr()`](https://pglpm.github.io/prova/reference/Pr.md), as their
  `$quantiles` value. This value is a three-dimensional array, and one
  of its columns (corresponding to the possible values of the `X`
  argument of [`Pr()`](https://pglpm.github.io/prova/reference/Pr.md))
  or one of its rows (corresponding to the possible values of the `Y`
  argument of [`Pr()`](https://pglpm.github.io/prova/reference/Pr.md))
  should be selected before being used as `y` input.

- xdomain:

  Character or numeric or `NULL` (default): vector of possible values of
  the variable represented in the x-axis, if the `x` argument is a
  character vector. The ordering of the values is respected. If `NULL`,
  then `unique(x)` is used.

- alpha.f:

  Numeric, default 0.25: opacity of the quantile bands, `0` being
  completely invisible and `1` completely opaque.

- col:

  Fill colour of the quantile bands. Can be specified in any of the
  usual ways, see for instance
  [`grDevices::col2rgb()`](https://rdrr.io/r/grDevices/col2rgb.html).
  Default `#4477AA`.

- border:

  Fill colour of the quantile bands. Can be specified in any of the
  usual ways, see for instance
  [`grDevices::col2rgb()`](https://rdrr.io/r/grDevices/col2rgb.html). If
  `NA` (default), no border is drawn.

- type:

  see analogous argument in
  [`flexiplot()`](https://pglpm.github.io/prova/reference/flexiplot.md).

- ...:

  Other parameters to be passed to
  [`flexiplot()`](https://pglpm.github.io/prova/reference/flexiplot.md).

## See also

[`Pr()`](https://pglpm.github.io/prova/reference/Pr.md) to calculate
posterior probabilities and quantiles.

[`plot.probability()`](https://pglpm.github.io/prova/reference/plot.probability.md)
to directly plot posterior probabilities and quantiles contained in a
probability object.

[`flexiplot()`](https://pglpm.github.io/prova/reference/flexiplot.md)
for more general plots.

## Examples

``` r
## Load the example `learnt` object calculated from the "penguins" dataset;
## variates: 'species' and 'bill_len'
learnt <- learntExample

## create a grid of values for variate "bill length",
## based on the information in the dataset and metadata:
values <- vrtgrid(vrt = 'bill_len', learnt = learnt)

## calculate the probabilities and quantiles
probs <- Pr(Y = data.frame(bill_len = values), learnt = learnt, parallel = 1)
#> Registered socket cluster with 1 nodes on host ‘localhost’.
#> Closing connections to cores.

## plot the quantiles, setting lower plot range to zero
plotquantiles(x = values, y = probs$quantiles[, 1, ], ylim = c(0, NA),
  xlab = 'bill length', ylab = 'probability')

## add a plot of the probabilities in thick dashed red
flexiplot(x = values, y = probs$values, lwd = 5, lty = 2, col = 2, add = TRUE)

```

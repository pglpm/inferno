# Plot pairs of quantiles

Utility function to plot pair of quantiles obtained with
[`Pr()`](https://pglpm.github.io/inferno/reference/Pr.md).

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
  [`flexiplot()`](https://pglpm.github.io/inferno/reference/flexiplot.md).

- y:

  Numeric: a matrix having as many rows as `x` and an even number of
  columns, with one column per quantile. Typically these quantiles have
  been obtained with
  [`Pr()`](https://pglpm.github.io/inferno/reference/Pr.md), as their
  `$quantiles` value. This value is a three-dimensional array, and one
  of its columns (corresponding to the possible values of the `X`
  argument of [`Pr()`](https://pglpm.github.io/inferno/reference/Pr.md))
  or one of its rows (corresponding to the possible values of the `Y`
  argument of [`Pr()`](https://pglpm.github.io/inferno/reference/Pr.md))
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

- ...:

  Other parameters to be passed to
  [`flexiplot()`](https://pglpm.github.io/inferno/reference/flexiplot.md).

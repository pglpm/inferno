# Plot numeric or character values

Plot function that modifies and expands the **graphics** package's
[`graphics::matplot()`](https://rdrr.io/r/graphics/matplot.html)
function in several ways.

## Usage

``` r
flexiplot(
  x,
  y,
  type = NULL,
  lty = c(1, 2, 4, 3, 6, 5),
  lwd = 2,
  pch = c(1, 2, 0, 5, 6, 3),
  col = palette(),
  xlab = NULL,
  ylab = NULL,
  xlim = NULL,
  ylim = NULL,
  add = FALSE,
  xdomain = NULL,
  ydomain = NULL,
  alpha.f = 1,
  xjitter = NULL,
  yjitter = NULL,
  grid = TRUE,
  cex.main = 1,
  ...
)
```

## Arguments

- x:

  Numeric or character: vector of x-coordinates. If missing, a numeric
  vector `1:...` is created having as many values as the rows of `y`.

- y:

  Numeric or character: vector of y coordinates. If missing, a numeric
  vector `1:...` is created having as many values as the rows of `x`.

- xlim, ylim:

  `NULL` (default) or a vector of two values. In the latter case, if any
  of the two values is not finite (including `NA` or `NULL`), then the
  `min` or `max` `x`- or `y`-coordinates of the plotted points are used.

- xdomain, ydomain:

  Character or numeric or `NULL` (default): vector of possible values of
  the variables represented in the `x`- and `y`-axes, in case the `x` or
  `y` argument is a character vector. The ordering of the values is
  respected. If `NULL`, then `unique(x)` or `unique(y)` is used.

- alpha.f:

  Numeric, default 1: opacity of the colours, `0` being completely
  invisible and `1` completely opaque.

- xjitter, yjitter:

  Logical or `NULL` (default): add
  [`base::jitter()`](https://rdrr.io/r/base/jitter.html) to `x`- or
  `y`-values? Useful when plotting discrete variates. If `NULL`, jitter
  is added if the values are of character class.

- grid:

  Logical: whether to plot a light grid. Default `TRUE`.

- ...:

  Other parameters to be passed to
  [`graphics::matplot()`](https://rdrr.io/r/graphics/matplot.html).

## Details

Some of the additional features provided by `flexiplot` are the
following. First, either or both `x` and `y` arguments can be of class
[`base::character`](https://rdrr.io/r/base/character.html). In this
case, axes labels corresponding to the unique values are used (see
arguments `xdomain` and `ydomain`). A jitter can also be added to the
generated points, via the `xjitter` and `yjitter` switches. Second, it
allows for the specification of only a lower or upper limit in the
`xlim` and `ylim` arguments. Third, it uses a cleaner plotting style and
a default argument `type = 'l'` (line plot) rather than `type = 'p'`
(point plot).

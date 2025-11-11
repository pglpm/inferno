# Plot an object of class "probability"

This [`base::plot()`](https://rdrr.io/r/base/plot.html) method is a
utility to plot probabilities obtained with
[`Pr()`](https://pglpm.github.io/inferno/reference/Pr.md), as well as
their variabilities. The probabilities are plotted either against `Y`,
with one curve for each value of `X`, or vice versa.

## Usage

``` r
# S3 method for class 'probability'
plot(
  p,
  variability = NULL,
  PvsY = NULL,
  legend = "top",
  lty = c(1, 2, 4, 3, 6, 5),
  lwd = 2,
  col = palette(),
  type = NULL,
  alpha.f = 1,
  var.alpha.f = NULL,
  xlab = NULL,
  ylab = NULL,
  main = NULL,
  ylim = c(0, NA),
  grid = TRUE,
  add = FALSE,
  ...
)
```

## Arguments

- p:

  Object of class "probability", obtained with
  [`Pr()`](https://pglpm.github.io/inferno/reference/Pr.md).

- variability:

  One of the values `'quantiles'`, `'samples'`, `'none'` (equivalent to
  `NA` or `FALSE`), or `NULL` (default), in which case the variability
  available in `p` is used. This argument chooses how to represent the
  variability of the probability; see
  [`Pr()`](https://pglpm.github.io/inferno/reference/Pr.md). If the
  requested variability is not available in the object `p`, then a
  warning is issued and no variability is plotted.

- PvsY:

  Logical or `NULL`: should probabilities be plotted against their `Y`
  argument? If `NULL`, the argument between `Y` and `X` having larger
  number of values is chosen. As many probability curves will be plotted
  as the number of values of the other argument.

- legend:

  One of the values `'bottomright'`, `'bottom'`, `'bottomleft'`,
  `'left'`, `'topleft'`, `'top'`, `'topright'`, `'right'`, `'center'`
  (see [`graphics::legend()`](https://rdrr.io/r/graphics/legend.html)):
  plot a legend at that position. A value `FALSE` or any other does not
  plot any legend. Default `'top'`.

- alpha.f:

  Numeric, default 0.25: opacity of the colours, `0` being completely
  invisible and `1` completely opaque.

- var.alpha.f:

  Numeric: opacity of the quantile bands or of the samples, `0` being
  completely invisible and `1` completely opaque.

- ...:

  Other parameters to be passed to
  [`flexiplot()`](https://pglpm.github.io/inferno/reference/flexiplot.md).

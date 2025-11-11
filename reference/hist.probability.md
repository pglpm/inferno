# Plot the variability of an object of class "probability" as a histogram

This [`graphics::hist()`](https://rdrr.io/r/graphics/hist.html)ogram
method is a utility to visualize the variability of the probabilities
obtained with [`Pr()`](https://pglpm.github.io/inferno/reference/Pr.md),
which can also be interpreted as the probability density for the
whole-population frequencies.

## Usage

``` r
# S3 method for class 'probability'
hist(
  p,
  breaks = NULL,
  legend = "top",
  lty = c(1, 2, 4, 3, 6, 5),
  lwd = 2,
  col = palette(),
  alpha.f = 1,
  fill.alpha.f = 0.125,
  showmean = TRUE,
  xlab = NULL,
  ylab = NULL,
  xlim = NULL,
  ylim = c(0, NA),
  main = NULL,
  grid = TRUE,
  add = FALSE,
  ...
)
```

## Arguments

- p:

  Object of class "probability", obtained with
  [`Pr()`](https://pglpm.github.io/inferno/reference/Pr.md).

- breaks:

  `NULL` or as in function
  [`graphics::hist()`](https://rdrr.io/r/graphics/hist.html). If `NULL`
  (default), an optimal number of breaks for each probability
  distribution is computed.

- legend:

  One of the values `"bottomright"`, `"bottom"`, `"bottomleft"`,
  `"left"`, `"topleft"`, `"top"`, `"topright"`, `"right"`, `"center"`
  (see [`graphics::legend()`](https://rdrr.io/r/graphics/legend.html)):
  plot a legend at that position. A value `FALSE` or any other does not
  plot any legend. Default `"top"`.

- fill.alpha.f:

  Numeric, default 0.125: opacity of the histogram filling. `0` means no
  filling.

- showmean:

  Logical, default `TRUE`: show the means of the probability
  distributions? The means correspond to the probabilities about the
  next observed unit.

- ...:

  Other parameters to be passed to
  [`flexiplot()`](https://pglpm.github.io/inferno/reference/flexiplot.md).

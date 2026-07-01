# Plot the variability of an object of class "probability" as a histogram

The posterior probabilities calculated with the
[`Pr()`](https://pglpm.github.io/prova/reference/Pr.md) function, and
outputted as a `probability` object, have an associated variability that
comes from the finite size of the data sample. This variability can be
interpreted in two ways:

- How the probabilities would change, if we could collect a very large
  (infinite) amount of additional data, and how likely would such change
  be;

- The relative frequency of a particular variate value in the full
  (sampled and unsampled) population is unknown; we can quantify our
  uncertainty about this relative frequency with a probability
  distribution.

The [`hist()`](https://rdrr.io/r/graphics/hist.html) method for a
`probability` object is a utility to visualize this kind of variability,
in the form of a distribution.

## Usage

``` r
# S3 method for class 'probability'
hist(
  x,
  subset = NULL,
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

- x:

  Object of class "probability", obtained with
  [`Pr()`](https://pglpm.github.io/prova/reference/Pr.md).

- subset:

  Named list or named vector: which variate values to display. For the
  variates corresponding to the names in this list, only the vector of
  values corresponding to that variate is displayed.

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

- lty, lwd, col, alpha.f, xlab, ylab, xlim, ylim, main, grid, add:

  see analogous arguments in
  [`graphics::matplot()`](https://rdrr.io/r/graphics/matplot.html)

- fill.alpha.f:

  Numeric, default 0.125: opacity of the histogram filling. `0` means no
  filling.

- showmean:

  Logical, default `TRUE`: show the means of the probability
  distributions? The means correspond to the probabilities about the
  next observed unit.

- ...:

  Other parameters to be passed to
  [`flexiplot()`](https://pglpm.github.io/prova/reference/flexiplot.md).

## Value

[Invisibly](https://rdrr.io/r/base/invisible.html), an object of class
["histogram"](https://rdrr.io/r/graphics/hist.html).

## See also

[`Pr()`](https://pglpm.github.io/prova/reference/Pr.md) to calculate
posterior probabilities and quantiles.

[`plot.probability()`](https://pglpm.github.io/prova/reference/plot.probability.md)
to plot the posterior probabilities.

[`flexiplot()`](https://pglpm.github.io/prova/reference/flexiplot.md)
(on which `hist.probability()` is based) for more general plots.

[`plotquantiles()`](https://pglpm.github.io/prova/reference/plotquantiles.md)
to plot quantile ranges.

## Examples

``` r
## Load the example `learnt` object calculated from the "penguins" dataset;
## variates: 'species' and 'bill_len'
learnt <- learntExample

## calculate the probability, and its variability,
## for the value 'Adelie' of the "species" variate
probs <- Pr(Y = data.frame(species = 'Adelie'), learnt = learnt, parallel = 1)
probs$values
#>         
#> species      [,1]
#>   Adelie 0.440685

## show the variability of this probability; equivalently show
## the probability distribution for the relative frequency of
## 'Adelie' penguins in the full population
hist(probs, legend = 'topright')

```

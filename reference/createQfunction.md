# Calculate and save transformation function for ordinal variates

NB: the functional form of this function does not depend on the number
of components, minalpha, and maxalpha parameters

## Usage

``` r
createQfunction(
  nint = 3600,
  nsamples = 2^24L,
  mean = 0,
  sd = 3,
  shapelo = 0.5,
  shapehi = 0.5,
  rate = 1,
  file = paste0("__Qfunction", nint, "_", sd),
  save = TRUE,
  plot = FALSE
)
```

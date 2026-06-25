# Calculate and save transformation function for ordinal variates

It creates the interpolation functions 'util_Q', 'util_invQ',
'util_invDQ' and saves them into 'sysdata.rda'.

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
  plot = FALSE
)
```

## Details

Those three functions are used to transform variates having bounded
domains into variates with unbounded domains. See
<https://github.com/pglpm/prova/raw/main/development/manual/pglpm2024-bayes_nonparam.pdf>.

NB: the functional form of this function does not depend on the number
of components, minalpha, and maxalpha parameters

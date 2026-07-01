# Join '\_\_\_\_tempPtraces-' files

Join '\_\_\_\_tempPtraces-' files

## Usage

``` r
util_joinPtraces(path)
```

## Value

A [data frame](https://rdrr.io/r/base/data.frame.html) of MCMC traces.

## Details

For deeper monitoring of the MCMC, the user can require the 'learn()'
function not to clean intermediate MCMC-related files generated during
the computation. The files with prefix '\_\_\_\_tempPtraces-' contain
chunks of MCMC traces.

The present function can be used to join them into a single trace.

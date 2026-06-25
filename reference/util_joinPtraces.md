# Join '\_\_\_\_tempPtraces-' files

For deeper monitoring of the MCMC, the user can require the 'learn()'
function not to clean intermediate MCMC-related files generated during
the computation. The files with prefix '\_\_\_\_tempPtraces-' contain
chunks of MCMC traces.

## Usage

``` r
util_joinPtraces(path)
```

## Details

The present function can be used to join them into a single trace.

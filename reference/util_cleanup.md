# Cleanup a learn()-output directory

For deeper monitoring of the MCMC, the user can require the 'learn()'
function not to clean intermediate MCMC-related files generated during
the computation.

## Usage

``` r
util_cleanup(path)
```

## Details

The present function can be used to remove these intermediate files from
the output directory created by 'learn()'.

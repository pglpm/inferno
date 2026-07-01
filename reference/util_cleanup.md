# Cleanup a learn()-output directory

Cleanup a learn()-output directory

## Usage

``` r
util_cleanup(path)
```

## Value

No return value; called for side effects.

## Details

For deeper monitoring of the MCMC, the user can require the 'learn()'
function not to clean intermediate MCMC-related files generated during
the computation.

The present function can be used to remove these intermediate files from
the output directory created by 'learn()'.

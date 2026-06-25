# Format datapoints used for MCMC monitoring

Used in 'util_Pcheckpoints()' within 'learn()'.

## Usage

``` r
util_prepPcheckpoints(x, auxmetadata, pointsid = NULL)
```

## Arguments

- x:

  Datapoints to be used for checking MCMC progress

- auxmetadata:

  auxmetadata object

- pointsid:

  Id of datapoints

## Value

some arguments to be repeatedly used in util_Pcheckpoints

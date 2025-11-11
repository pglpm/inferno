# Create a grid of values for a variate

This function create a set of values for a variate `vrt`, based on the
metadata stored in a `learnt` object created by the
[`learn()`](https://pglpm.github.io/inferno/reference/learn.md)
function). The set of values depends on the type of variate (nominal or
continuous, rounded, and so on, see
[`metadata`](https://pglpm.github.io/inferno/reference/metadatatemplate.md)).
The range of values is chosen to include, and extend slightly beyond,
the range observed in the data used in the
[`learn()`](https://pglpm.github.io/inferno/reference/learn.md)
function. Variate domains are always respected.

## Usage

``` r
vrtgrid(vrt, learnt, length.out = 129)
```

## Arguments

- vrt:

  Character: name of the variate, must match one of the names in the
  `metadata` file provided to the
  [`learn()`](https://pglpm.github.io/inferno/reference/learn.md)
  function.

- learnt:

  Either a character with the name of a directory or full path for a
  'learnt.rds' object, produced by the
  [`learn()`](https://pglpm.github.io/inferno/reference/learn.md)
  function, or such an object itself.

- length.out:

  Numeric, positive (default 129): number of values to be created; used
  only for continuous, non-rounded variates (see
  [`metadata`](https://pglpm.github.io/inferno/reference/metadatatemplate.md)).

## Value

A numeric or character vector of values.

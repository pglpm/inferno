# Subset variates of an object of class `probability`

An object of class `probability`, obtained with the
[`Pr()`](https://pglpm.github.io/inferno/reference/Pr.md) function,
holds the probabilities for all possible combinations of values of a set
of joint variates `Y` conditional on a set of joint variates `X`. In
some cases one may wish to exclude some of the values of the `Y` or `X`
variates. For instance `Y` in the probability-class object could include
the variate "age" with values from 18 to 100, and one may want to retain
the values from 60 to 80.

## Usage

``` r
# S3 method for class 'probability'
subset(x, subset)
```

## Arguments

- x:

  Object of class "probability", obtained with
  [`Pr()`](https://pglpm.github.io/inferno/reference/Pr.md).

- subset:

  Named list or named vector: variates to subset, given as list names,
  and corresponding values to subset.

## Value

A list of class `probability`, identical to the original object `x`
except for a reduced range of values in some if its variates.

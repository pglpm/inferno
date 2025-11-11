# Transforms variates to different representations

Transforms variates to different representations

## Usage

``` r
vtransform(
  x,
  auxmetadata,
  Rout = NULL,
  Cout = NULL,
  Dout = NULL,
  Bout = NULL,
  Oout = NULL,
  Nout = NULL,
  variates = NULL,
  logjacobianOr = NULL
)
```

## Arguments

- x:

  data.table object containing data to be transformed

- auxmetadata:

  auxmetadata object

- Rout:

  Character, output of R-type variate, with possible values:
  'normalized': for internal MCMC use 'mi': for use in mutualinfo()
  'original': original representation

- Cout:

  Character, output of C-type variate, with possible values: 'init': for
  internal MCMC use (init input) 'left', 'right': for internal MCMC use
  'aux', 'lat': for internal MCMC use 'boundnormalized': for sampling
  functions 'boundisinf': for sampling functions 'mi': for use in
  mutualinfo() 'original': original representation

- Dout:

  Character, output of D-type variate, with possible values: 'init': for
  internal MCMC use (init input) 'left', 'right': for internal MCMC use
  'aux': for internal MCMC use 'boundisinf': for sampling functions
  'normalized': for sampling functions 'mi': for use in mutualinfo()
  'original': original representation

- Bout:

  Character, output of B-type variate, with possible values: 'numeric':
  for internal MCMC use, values 0,1 'original': original representation

- Oout:

  Character, output of O-type variate, with possible values: 'numeric':
  for internal MCMC use, values 1,2,... 'original': original
  representation

- Nout:

  Character, output of N-type variate, with possible values: 'numeric':
  for internal MCMC use, values 1,2,... 'original': original
  representation

- variates:

  Character vector, names of variates corresponding to columns of x (in
  case x misses column names)

- logjacobianOr:

  Logical or `NULL`: output is the log-Jacobian in orginal or
  transformed domain? `NULL` (default) means do not calculate the
  log-Jacobians

## Value

data frame of transformed variates, or their log-Jacobians

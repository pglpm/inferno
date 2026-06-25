# Calculate quantiles for continuous Y by bisection

Used in 'qPr()'.

## Usage

``` r
util_qYXcont(
  iyx,
  params1,
  params2,
  auxmetadata,
  temporarydir,
  usememory = TRUE,
  doquantiles,
  quantiles,
  dosamples,
  nsamples,
  Qerror,
  tol = .Machine$double.eps * 3
)
```

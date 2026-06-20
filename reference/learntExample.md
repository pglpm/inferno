# Example `learnt` object produced by learn()

An example `learnt` object obtained by means of the
[`learn()`](https://pglpm.github.io/prova/reference/learn.md) function,
using the [datasets::penguins](https://rdrr.io/r/datasets/penguins.html)
dataset and the metadata in
[metadataExample](https://pglpm.github.io/prova/reference/metadataExample.md),
according to the call

    learn(data = penguins, metadata = metadataExample,
      nsamples = 225, nchains = 15)

It is a list that essentially contains posterior hyperparameters for
drawing statistical inferences about the variates `species` and
`bill_len`.

**Note** that the
[`learn()`](https://pglpm.github.io/prova/reference/learn.md) function
that produced `learntExample` was called with the option to create only
a limited number (225) of Monte Carlo samples, in order to reduce its
memory size. Thus the numerical error associated with the Monte Carlo
approximation is relatively in inferences drawn from the posterior
hyperparameters saved in `learntExample`. It is only meant to be used
for illustration purposes of the package's capabilities.

## Usage

``` r
learntExample
```

## Format

### `learntExample`

A list containing results from Markov-chain Monte Carlo computation,
including diagnostics and variate metadata.

## See also

[`learn()`](https://pglpm.github.io/prova/reference/learn.md), which
produces this kind of object.

[`Pr()`](https://pglpm.github.io/prova/reference/Pr.md),
[`qPr()`](https://pglpm.github.io/prova/reference/qPr.md),
[`rPr()`](https://pglpm.github.io/prova/reference/rPr.md),
[`mutualinfo()`](https://pglpm.github.io/prova/reference/mutualinfo.md):
functions that require this kind of object in order to calculate
probabilities and quantiles, generate data points, and calculate mutual
information.

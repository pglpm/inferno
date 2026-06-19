# Example `learnt` object produced by learn()

An example `learnt` object obtained by means of the
[`learn()`](https://pglpm.github.io/prova/reference/learn.md) function,
using the [datasets::penguins](https://rdrr.io/r/datasets/penguins.html)
dataset and the metadata in
[metadataExample](https://pglpm.github.io/prova/reference/metadataExample.md).
It is a list that essentially contains posterior hyperparameters for
drawing statistical inferences about the variates `species` and
`bill_len`.

**Note** that the
[`learn()`](https://pglpm.github.io/prova/reference/learn.md) function
that produced `learntExample` was called with the option to create only
a limited number (360) of Monte Carlo samples, in order to reduce its
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

[`learn()`](https://pglpm.github.io/prova/reference/learn.md) which
produces this kind of object,
[`Pr()`](https://pglpm.github.io/prova/reference/Pr.md) which calculates
variousposterior probabilities based on this kind of object.

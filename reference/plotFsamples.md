# Plot one-dimensional posterior probabilities

Plot one-dimensional posterior probabilities

## Usage

``` r
plotFsamples(
  filename,
  learnt,
  data,
  plotprobability = TRUE,
  plotvariability = "samples",
  nFsamples = NULL,
  datahistogram = !(missing(data) || is.null(data)),
  datascatter = !(missing(data) || is.null(data)),
  parallel = TRUE,
  silent = FALSE
)
```

## Arguments

- filename:

  Character: name of plot output file

- data:

  data.table object or filepath: datapoints

- plotprobability:

  Logical: plot the resulting probability curve

- plotvariability:

  Character, either 'samples' or 'quantiles': how to plot the
  variability of the probability distribution with new samples

- nFsamples:

  Positive number: if plotvariability='samples', then number of samples
  of representative frequency distributions to display as variability;
  if plotvariability='quantiles', then the quantiles (in range 0 to 0.5)
  to show

- datahistogram:

  Logical: plot the data as histogram?

- datascatter:

  Logical: plot the data as scatterplot along the x-axis?

- parallel:

  Logical or numeric: whether to use pre-existing parallel workers, or
  how many to create and use

- silent:

  Logical: give warnings or updates in the computation

- learned:

  Either a character with the name of a directory or full path for an
  'learnt.rds' object, or such an object itself

## Value

A list with the mutual information, its error, and its unit

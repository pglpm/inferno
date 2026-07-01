# Build augmented metadata file

Builds an augmented metadata [data
frame](https://rdrr.io/r/base/data.frame.html) from the metadata and
data given to 'learn()'. This augmented metadata object is saved in the
'learnt' object produced by 'learn()'.

## Usage

``` r
buildauxmetadata(data, metadata, Dthreshold = 1, tscalefactor = 4.266)
```

## Arguments

- data:

  data.frame object

- metadata:

  data.frame object

- Dthreshold:

  Positive number: threshold of fraction of unique datapoints to total
  datapoints, to decide whether to treat a rounded variate as continuous

- tscalefactor:

  Positive number: scaling factor for variate conversion

## Value

A [data frame](https://rdrr.io/r/base/data.frame.html) with auxmetadata.

## Details

In addition to the original metadata it contains info about transformed
variates and their domains, estimated location- and scale-parameters,
and similar metadata.

Used in 'learn()'.

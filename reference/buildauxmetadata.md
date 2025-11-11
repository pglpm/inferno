# Build preliminary metadata flie

Build preliminary metadata flie

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

an auxmetadata data.frame object

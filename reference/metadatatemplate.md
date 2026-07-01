# Metadata and helper function for metadata

Metadata and helper function to create a template metadata file or
object.

## Usage

``` r
metadatatemplate(
  data,
  file = NULL,
  includevrt = NULL,
  excludevrt = NULL,
  addsummary2metadata = FALSE,
  backupfiles = FALSE,
  verbose = TRUE
)
```

## Arguments

- data:

  A dataset, given as a [data
  frame](https://rdrr.io/r/base/data.frame.html) or as a file path to a
  csv file.

- file:

  Character or `NULL` (default): name of csv file where the metadata
  should be saved; if `NULL`: output metadata as `VALUE`.

- includevrt:

  Character or `NULL`: name of variates in dataset to be included.

- excludevrt:

  Character or `NULL`: name of variates in dataset to be excluded.

- addsummary2metadata:

  Logical: also output some diagnostic statistics in the metadata?
  Default `FALSE`.

- backupfiles:

  Logical: rename previous metadata file if it exists? Default `TRUE`.

- verbose:

  Logical: output heuristics for each variate? Default `TRUE`.

## Value

A preliminary [data frame](https://rdrr.io/r/base/data.frame.html)
containing the metadata,
[invisibly](https://rdrr.io/r/base/invisible.html) if `file = NULL`. If
argument `file` is a character, a preliminary metadata file is also
created with that name or path.

## Details

The [`learn()`](https://pglpm.github.io/prova/reference/learn.md)
function needs metadata about the variates present in the data. Such
metadata can be provided either as a `csv` file or as a
[`base::data.frame()`](https://rdrr.io/r/base/data.frame.html). The
function `buildmetadata` creates a template metadata csv-file, or
outputs a metadata data.frame, by trying to *guess* metadata information
from the dataset.The guesses may be very incorrect (as already said,
metadata is information not contained in the data, so no algorithm can
exist that extracts it from the data). **The user *must* modify and
correct this template, using it as a starting point to prepare the
correct metadata information.**

## Metadata information and format

In order to correctly learn from a dataset, the
[`learn()`](https://pglpm.github.io/prova/reference/learn.md) function
needs information that is not contained in the data themeselves; that
is, it needs *meta*data. Metadata are provided either as a `csv` file or
as a [`base::data.frame()`](https://rdrr.io/r/base/data.frame.html).

A metadata file or data.frame must contain one row for each simple
variate in the given inference problem, and the following fields
(columns), even if some of them may be empty:

`name`, `type`, `domainmin`, `domainmax`, `datastep`, `minincluded`,
`maxincluded`, `V1`, `V2`, (possibly additional `V`-fields, sequentially
numbered)

The `type` field has three possible values: `nominal`, `ordinal`,
`continuous`. The remaining fields that must be filled in depend on the
`type` field. Here is a list of requirements:

- **`nominal`** and **`ordinal`**: require *either* `V1`, `V2`, ...
  fields *or* `domainmin`, `domainmax`, `datastep` (all three) fields.
  No other fields are required.

- **`continuous`**: requires `domainmin`, `domainmax`, `datastep`,
  `minincluded`, `maxincluded`.

Here are the meanings and possible values of the fields:

**`name`**: The name of the variate. This must be the same character
string as it appears in the dataset (be careful about upper- and
lower-case).

**`type`**: The data type of variate `name`. Possible values are
`nominal`, `ordinal`, `continuous`.

- A *nominal* (also called *categorical*) variate has a discrete, finite
  number of possible values which have no intrinsic ordering. Examples
  could be a variate related to colour, with values "red", "green",
  "blue", and so on; or a variate related to cat breeds, with values
  "Siamese", "Abyssinian", "Persian", and so on. The possible values of
  the variate must be given in the fields `V1`, `V2`, and so on. It is
  important to include values that are possible but are *not* present in
  the dataset. A variate having only two possible values (binary
  variate), for example "yes" and "no", can be specified as nominal.

- An *ordinal* variate has a discrete, finite number of possible values
  which do have an intrinsic ordering. Examples could be a Likert-scaled
  variate for the results of a survey, with values "very dissatisfied",
  "dissatisfied", "satisfied", "very satisfied"; or a variate related to
  the levels of some quantities, with values "low", "medium", "high"; or
  a variate having a numeric scale with values from 1 to 10. Whether a
  variate is nominal or ordinal often depends on the context. The
  possible values of the variate but be given in either one (but not
  both) or two ways: (1) in the fields `V1`, `V2`, ..., as for nominal
  variates; (2) as the fields `domainmin`, `domainmax`, `datastep`.
  Option (2) only works with numeric, equally spaced values: it assumes
  that the first value is `domainmin`, the second is
  `domainmin`+`datastep`, the third is `domainmin`+2\*`datastep`, and so
  on up to the last value, `domainmax`.

- A *continuous* variate has a continuum of values with an intrinsic
  ordering. Examples could be a variate related to the width of an
  object; or to the age of a person; or one coordinate of an object in a
  particular reference system. A continuous variate requires
  specification of the fields `domainmin`, `domainmax`, `datastep`,
  `minincluded`, `maxincluded`. Some naturally continuous variates are
  often rounded to a given precision; for instance, the age of a person
  might be reported as rounded to the nearest year (25 years, 26 years,
  and so on); or the length of an object might be reported to the
  nearest centimetre (1 m, 1.01 m, 1.02 m, and so on). The minimum
  distance between such rounded values **must** be reported in the
  `datastep` field; this would be `1` in the age example and `0.01` in
  the length example above. See below for further explanation of why
  reporting such rounding is important.

**`domainmin`**: The minimum value that the variate (ordinal or
continuous) can take on. Possible values are a real number or an empty
value, which is then interpreted as `-Inf` (explicit values like `-Inf`,
`-inf`, `-infinity` should also work). Some continuous variates, like
age or distance or temperature, are naturally positive, and therefore
have `domainmin` equal `0`. But in other contexts the minimum value
could be different. For instance, if a given inference problem only
involves people of age 18 or more, then `domainmin` would be set to
`18`.

**`domainmax`**: The maximum value that the variate (ordinal or
continuous) can take on. Possible values are a real number, or an empty
value, which is then interpreted as `+Inf` (explicit values like `Inf`,
`inf`, `infinity` should also work). As with `domainmin`, the maximum
value depends on the context. An age-related variate could theoretically
have `domainmax` equal to infinity (empty value in the metadata file);
but if a given study categorizes some people as "90 years old or older",
then `domainmax` should be set to `90`.

**`datastep`**: The minimum distance between the values of a variate
(ordinal or continuous). Possible values are a positive real number or
an empty value, which is then interpreted as 0 (the explicit value `0`
is also accepted). For a numeric ordinal variate, `datastep` is the step
between consecutive values. For a continuous *rounded* variate,
`datastep` is the minimum distance between different values that occurs
because of rounding; see the examples given above. The function
`buildmetadata` has some heuristics to determine whether the variate is
rounded or not. See further details under the section Rounding below.

**`minincluded`**, **`maxincluded`**: Whether the minimum (`domainmin`)
and maximum(`domainmax`) values of a *continuous* variate can really
appear in the data or not. Possible values are `true` (or `t` or `yes`)
or `false` (or `f`, `no`, or an empty field); upper- or lower-case is
irrelevant. Here are some examples about the meaning of these fields.
(a) A continuous *unrounded* variate such as temperature has 0 as a
minimum possible value `domainmin`, but this value itself is physically
impossible and can never appear in data; in this case `minincluded` is
empty (or set to `false` or `no`). (b) A variate related to the
*unrounded* length, in metres, of some objects may take on any positive
real value; but suppose that all objects of length 5 or less are grouped
together under the value `5`. It is then possible for several datapoints
to have value `5`: one such datapoint could originally have the value
3.782341...; another the value 4.929673..., and so on. In this case
`domainmin` is set to `5`, and `minincluded` is set to `true` (or
`yes`). Similarly for the maximum value of a variate and `maxincluded`.
Note that if `domainmin` is minus-infinity (empty value in the metadata
file), then `minincluded` is automatically empty (that is, `false`), and
similarly for `maxincluded` if `domainmax` is infinity.

## See also

[`learn()`](https://pglpm.github.io/prova/reference/learn.md), which
generates the information necessary to calculate posterior
probabilities, based on data and metadata.

## Examples

``` r
## Create a preliminary data frame of metadata for the `penguins` dataset
metadata <- metadatatemplate(data = datasets::penguins, file = NULL)
#> Converting factors to characters
#> Analyzing8variates for344datapoints.
#> * "species" variate:
#>   - 3 different values detected:
#> "Adelie", "Chinstrap", "Gentoo"
#>   which do not seem to refer to an ordered scale.
#>   Assuming variate to be NOMINAL.
#> * "island" variate:
#>   - 3 different values detected:
#> "Biscoe", "Dream", "Torgersen"
#>   which do not seem to refer to an ordered scale.
#>   Assuming variate to be NOMINAL.
#> * "bill_len" variate:
#>   - Numeric values between 32.1 and 59.6
#>   Assuming variate to be CONTINUOUS.
#>   - Distance between datapoints is a multiple of 0.1
#>   Assuming variate to be ROUNDED.
#>   - All values are positive
#>   Assuming "domainmin" to be 0
#> * "bill_dep" variate:
#>   - Numeric values between 13.1 and 21.5
#>   Assuming variate to be CONTINUOUS.
#>   - Distance between datapoints is a multiple of 0.1
#>   Assuming variate to be ROUNDED.
#>   - All values are positive
#>   Assuming "domainmin" to be 0
#> * "flipper_len" variate:
#>   - Numeric values between 172 and 231
#>   Assuming variate to be CONTINUOUS.
#>   - Distance between datapoints is a multiple of 1
#>   Assuming variate to be ROUNDED.
#>   - All values are positive
#>   Assuming "domainmin" to be 0
#> * "body_mass" variate:
#>   - Numeric values between 2700 and 6300
#>   Assuming variate to be CONTINUOUS.
#>   - Distance between datapoints is a multiple of 25
#>   Assuming variate to be ROUNDED.
#>   - All values are positive
#>   Assuming "domainmin" to be 0
#> * "sex" variate:
#>   - 2 different values detected:
#> "female", "male"
#>   which do not seem to refer to an ordered scale.
#>   Assuming variate to be NOMINAL.
#> * "year" variate:
#>   - Only3 different numeric values detected:
#> from 2007 to 2009 in steps of 1
#>   Assuming variate to be ORDINAL.
#> =========
#> WARNINGS - please make sure to check these variates in the metadata file:
#> 
#> * "flipper_len" variate appears to be continuous and rounded,
#> but it could also be an ordinal variate
#> 
#> * "body_mass" variate appears to be continuous and rounded,
#> but it could also be an ordinal variate
#> 
#> * "year" variate appears to have been rounded
#> and then transformed to logarithmic scale.
#> This may lead to problems in the inference.
#> Preferably, transform it back to non-logarithmic scale.
#> =========

## Note how the preliminary data frame includes additional spots
## for values of nominal and ordinal variates
## which could be missing from the data
print(metadata)
#>          name       type domainmin domainmax datastep minincluded maxincluded
#> 1     species    nominal        NA        NA       NA          NA          NA
#> 2      island    nominal        NA        NA       NA          NA          NA
#> 3    bill_len continuous         0        NA      0.1          NA          NA
#> 4    bill_dep continuous         0        NA      0.1          NA          NA
#> 5 flipper_len continuous         0        NA      1.0          NA          NA
#> 6   body_mass continuous         0        NA     25.0          NA          NA
#> 7         sex    nominal        NA        NA       NA          NA          NA
#> 8        year    ordinal      2007      2009      1.0          NA          NA
#>       V1        V2        V3 V4 V5 V6 V7 V8 V9 V10 V11
#> 1 Adelie Chinstrap    Gentoo NA NA NA NA NA NA  NA  NA
#> 2 Biscoe     Dream Torgersen NA NA NA NA NA NA  NA  NA
#> 3   <NA>      <NA>      <NA> NA NA NA NA NA NA  NA  NA
#> 4   <NA>      <NA>      <NA> NA NA NA NA NA NA  NA  NA
#> 5   <NA>      <NA>      <NA> NA NA NA NA NA NA  NA  NA
#> 6   <NA>      <NA>      <NA> NA NA NA NA NA NA  NA  NA
#> 7 female      male      <NA> NA NA NA NA NA NA  NA  NA
#> 8   <NA>      <NA>      <NA> NA NA NA NA NA NA  NA  NA


## Create a preliminary data frame of metadata for the `penguins` dataset,
## including only the 'species' and 'bill_len' variates:
metadata2 <- metadatatemplate(
  data = datasets::penguins, file = NULL,
  includevrt = c('species', 'bill_len')
)
#> Converting factors to characters
#> Analyzing2variates for344datapoints.
#> * "species" variate:
#>   - 3 different values detected:
#> "Adelie", "Chinstrap", "Gentoo"
#>   which do not seem to refer to an ordered scale.
#>   Assuming variate to be NOMINAL.
#> * "bill_len" variate:
#>   - Numeric values between 32.1 and 59.6
#>   Assuming variate to be CONTINUOUS.
#>   - Distance between datapoints is a multiple of 0.1
#>   Assuming variate to be ROUNDED.
#>   - All values are positive
#>   Assuming "domainmin" to be 0

print(metadata2)
#>       name       type domainmin domainmax datastep minincluded maxincluded
#> 1  species    nominal        NA        NA       NA          NA          NA
#> 2 bill_len continuous         0        NA      0.1          NA          NA
#>       V1        V2     V3 V4 V5 V6 V7 V8 V9 V10 V11
#> 1 Adelie Chinstrap Gentoo NA NA NA NA NA NA  NA  NA
#> 2   <NA>      <NA>   <NA> NA NA NA NA NA NA  NA  NA


## Create a preliminary data frame of metadata for the `penguins` dataset,
## excluding the 'year' variate:
metadata3 <- metadatatemplate(
  data = datasets::penguins, file = NULL,
  excludevrt = 'year'
)
#> Converting factors to characters
#> Analyzing7variates for344datapoints.
#> * "species" variate:
#>   - 3 different values detected:
#> "Adelie", "Chinstrap", "Gentoo"
#>   which do not seem to refer to an ordered scale.
#>   Assuming variate to be NOMINAL.
#> * "island" variate:
#>   - 3 different values detected:
#> "Biscoe", "Dream", "Torgersen"
#>   which do not seem to refer to an ordered scale.
#>   Assuming variate to be NOMINAL.
#> * "bill_len" variate:
#>   - Numeric values between 32.1 and 59.6
#>   Assuming variate to be CONTINUOUS.
#>   - Distance between datapoints is a multiple of 0.1
#>   Assuming variate to be ROUNDED.
#>   - All values are positive
#>   Assuming "domainmin" to be 0
#> * "bill_dep" variate:
#>   - Numeric values between 13.1 and 21.5
#>   Assuming variate to be CONTINUOUS.
#>   - Distance between datapoints is a multiple of 0.1
#>   Assuming variate to be ROUNDED.
#>   - All values are positive
#>   Assuming "domainmin" to be 0
#> * "flipper_len" variate:
#>   - Numeric values between 172 and 231
#>   Assuming variate to be CONTINUOUS.
#>   - Distance between datapoints is a multiple of 1
#>   Assuming variate to be ROUNDED.
#>   - All values are positive
#>   Assuming "domainmin" to be 0
#> * "body_mass" variate:
#>   - Numeric values between 2700 and 6300
#>   Assuming variate to be CONTINUOUS.
#>   - Distance between datapoints is a multiple of 25
#>   Assuming variate to be ROUNDED.
#>   - All values are positive
#>   Assuming "domainmin" to be 0
#> * "sex" variate:
#>   - 2 different values detected:
#> "female", "male"
#>   which do not seem to refer to an ordered scale.
#>   Assuming variate to be NOMINAL.
#> =========
#> WARNINGS - please make sure to check these variates in the metadata file:
#> 
#> * "flipper_len" variate appears to be continuous and rounded,
#> but it could also be an ordinal variate
#> 
#> * "body_mass" variate appears to be continuous and rounded,
#> but it could also be an ordinal variate
#> =========

print(metadata3)
#>          name       type domainmin domainmax datastep minincluded maxincluded
#> 1     species    nominal        NA        NA       NA          NA          NA
#> 2      island    nominal        NA        NA       NA          NA          NA
#> 3    bill_len continuous         0        NA      0.1          NA          NA
#> 4    bill_dep continuous         0        NA      0.1          NA          NA
#> 5 flipper_len continuous         0        NA      1.0          NA          NA
#> 6   body_mass continuous         0        NA     25.0          NA          NA
#> 7         sex    nominal        NA        NA       NA          NA          NA
#>       V1        V2        V3 V4 V5 V6 V7 V8 V9 V10 V11
#> 1 Adelie Chinstrap    Gentoo NA NA NA NA NA NA  NA  NA
#> 2 Biscoe     Dream Torgersen NA NA NA NA NA NA  NA  NA
#> 3   <NA>      <NA>      <NA> NA NA NA NA NA NA  NA  NA
#> 4   <NA>      <NA>      <NA> NA NA NA NA NA NA  NA  NA
#> 5   <NA>      <NA>      <NA> NA NA NA NA NA NA  NA  NA
#> 6   <NA>      <NA>      <NA> NA NA NA NA NA NA  NA  NA
#> 7 female      male      <NA> NA NA NA NA NA NA  NA  NA

## Generate 10 points for a continuous variate in (0, 1)
dataset <- runif(10)

## `metadatatemplate` correctly guesses the variate minimum,
## but not the maximum (`NA` is equivalent to `+Inf`)
metadata <- metadatatemplate(data = dataset, file = NULL)
#> Analyzing1variates for10datapoints.
#> * "data" variate:
#>   - Numeric values between 0.00249242829158902 and 0.996412934735417
#>   Assuming variate to be CONTINUOUS.
#>   - All values are positive
#>   Assuming "domainmin" to be 0
#>   with 0 excluded from domain.
print(metadata)
#>   name       type domainmin domainmax datastep minincluded maxincluded V1 V2 V3
#> 1 data continuous         0        NA       NA          NA          NA NA NA NA
#>   V4 V5 V6 V7 V8 V9 V10 V11
#> 1 NA NA NA NA NA NA  NA  NA
```

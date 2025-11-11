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

  A dataset, given as a
  [`base::data.frame()`](https://rdrr.io/r/base/data.frame.html) or as a
  file path to a csv file.

- file:

  Character: name of csv file where the metadata should be saved; if
  `NULL`: output metadata as `VALUE`.

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

If `file = NULL`, a preliminary metadata file is created and `VALUE` is
`NULL`; otherwise `VALUE` is a
[`base::data.frame()`](https://rdrr.io/r/base/data.frame.html)
containing the metadata.

## Details

The [`learn()`](https://pglpm.github.io/inferno/reference/learn.md)
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
[`learn()`](https://pglpm.github.io/inferno/reference/learn.md) function
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

## Rounded continuous variates

To be written.

## Necessity of metadata

To be written.

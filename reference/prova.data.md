# Write and read CSV files in **Prova**

Utility functions to read and write CSV files in the format required by
**Prova**

## Usage

``` r
pwrite.csv(x, file, ...)

pread.csv(file, ...)
```

## Arguments

- x:

  The object to be written. Preferably a matrix or data frame; if not,
  it is attempted to coerce `x` to a data frame. See
  [`utils::write.table()`](https://rdrr.io/r/utils/write.table.html).

- file:

  Either a character naming a file or a connection open for writing or
  reading. See
  [`utils::write.table()`](https://rdrr.io/r/utils/write.table.html) and
  [`utils::read.table()`](https://rdrr.io/r/utils/read.table.html).

- ...:

  Other arguments to be passed to
  [`utils::write.table()`](https://rdrr.io/r/utils/write.table.html) or
  [`utils::read.table()`](https://rdrr.io/r/utils/read.table.html).
  Arguments 'row.names', 'quote', 'na', 'na.strings', 'tryLogical',
  'sep', 'dec' are not allowed.

## Details

The functions
[`learn()`](https://pglpm.github.io/prova/reference/learn.md) and
[`metadatatemplate()`](https://pglpm.github.io/prova/reference/metadatatemplate.md)
accept CSV files formatted as follows:

- Decimal values should be separated by a *dot*; no comma should be used
  to separate thousands etc. Example: `86342.75` .

- Character and names should be quoted in single or double quotes.
  Example: `"female"`.

- Values should be separated by *commas*, not by tabs or semicolons.

- Missing values should be simply *empty*, not denoted by "NA",
  "missing", "-", or similar.

- Preferably there should not be factors (see
  [base::factor](https://rdrr.io/r/base/factor.html)); use character
  names instead.

The utility functions `pwrite.csv()` and `pread.csv()` are wrappers to
[`utils::write.table()`](https://rdrr.io/r/utils/write.table.html) and
[`utils::read.table()`](https://rdrr.io/r/utils/read.table.html) that
set appropriate default parameters according to the formatting rules
above.

## See also

[`metadatatemplate()`](https://pglpm.github.io/prova/reference/metadatatemplate.md)
to help writing metadata files.

[`learn()`](https://pglpm.github.io/prova/reference/learn.md), which
needs a metadata data-frame or CSV file.

## Examples

``` r
## Save the 'penguins' dataset in a (temporary) file
filename <- tempfile(fileext = '.csv')

pwrite.csv(penguins, file = filename)

## check first few lines of the raw file
writeLines(readLines(filename, n = 10))
#> "species","island","bill_len","bill_dep","flipper_len","body_mass","sex","year"
#> "Adelie","Torgersen",39.1,18.7,181,3750,"male",2007
#> "Adelie","Torgersen",39.5,17.4,186,3800,"female",2007
#> "Adelie","Torgersen",40.3,18,195,3250,"female",2007
#> "Adelie","Torgersen",,,,,,2007
#> "Adelie","Torgersen",36.7,19.3,193,3450,"female",2007
#> "Adelie","Torgersen",39.3,20.6,190,3650,"male",2007
#> "Adelie","Torgersen",38.9,17.8,181,3625,"female",2007
#> "Adelie","Torgersen",39.2,19.6,195,4675,"male",2007
#> "Adelie","Torgersen",34.1,18.1,193,3475,,2007
```

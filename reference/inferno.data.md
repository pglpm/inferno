# Write and read data in **inferno**

Utility functions to read and write CSV files in the format required by
**inferno**

## Usage

``` r
write.csvi(x, file, ...)

read.csvi(file, ...)
```

## Arguments

- x:

  The object to be written, preferably a matrix or data frame. If not,
  it is attempted to coerce `x` to a data frame. See
  [`utils::write.table()`](https://rdrr.io/r/utils/write.table.html).

- file:

  Either a character naming a file or a connection open for writing or
  reading. See
  [`utils::write.table()`](https://rdrr.io/r/utils/write.table.html) and
  [`utils::read.table()`](https://rdrr.io/r/utils/read.table.html).

- ...:

  Other parameters to be passed to
  [`utils::write.table()`](https://rdrr.io/r/utils/write.table.html) or
  [`utils::read.table()`](https://rdrr.io/r/utils/read.table.html).

## Details

The functions of the **inferno** package accept CSV files formatted as
follows:

- Decimal values should be separated by a *dot*; no comma should be used
  to separate thousands etc. Example: `86342.75` .

- Character and names should be quoted in single or double quotes.
  Example: `"female"`.

- Values should be separated by *commas*, not by tabs or semicolons.

- Missing values should be simply *empty*, not denoted by "NA",
  "missing", "-", or similar.

- Preferably there should not be
  [`base::factor`](https://rdrr.io/r/base/factor.html)s; use character
  names instead.

The utility functions `write.csvi()` and `read.csvi()` are wrappers to
[`utils::write.table()`](https://rdrr.io/r/utils/write.table.html) or
[`utils::read.table()`](https://rdrr.io/r/utils/read.table.html) that
sets appropriate default parameters

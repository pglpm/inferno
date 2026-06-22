#' @name prova.data
#' @rdname prova.data
#'
#' @title Write and read CSV files in **Prova**
#'
#' @description Utility functions to read and write CSV files in the format required by **Prova**
#'
#' @param x The object to be written. Preferably a matrix or data frame; if not, it is attempted to coerce `x` to a data frame. See [utils::write.table()].
#'
#' @param file Either a character naming a file or a connection open for writing or reading. See [utils::write.table()] and [utils::read.table()].
#' @param ... Other arguments to be passed to [utils::write.table()] or [utils::read.table()]. Arguments 'row.names', 'quote', 'na', 'na.strings', 'tryLogical', 'sep', 'dec' are not allowed.
#'
#' @details
#' The functions [learn()] and [metadatatemplate()] accept CSV files formatted as follows:
#'
#' - Decimal values should be separated by a *dot*; no comma should be used to separate thousands etc. Example: `86342.75`.
#' - Character and names should be quoted in single or double quotes. Example: `"female"`.
#' - Values should be separated by *commas*, not by tabs or semicolons.
#' - Missing values should be simply *empty*, not denoted by "NA", "missing", "-", or similar.
#' - Preferably there should not be factors (see [base::factor]); use character names instead.
#'
#' The utility functions `pwrite.csv()` and `pread.csv()` are wrappers to [utils::write.table()] and [utils::read.table()] that set appropriate default parameters according to the formatting rules above.
#'
#' @seealso
#' [metadatatemplate()] to help writing metadata files.
#'
#' [learn()], which needs a metadata data-frame or CSV file.
#'
#' @examples
#' ## Save the 'penguins' dataset in a (temporary) file
#' filename <- tempfile(fileext = '.csv')
#'
#' pwrite.csv(penguins, file = filename)
#'
#' ## check first few lines of the raw file
#' writeLines(readLines(filename, n = 10))
#'
NULL

#' @rdname prova.data
#'
#' @import utils
#'
#' @export
pwrite.csv <- function(x, file, ...){
    write.csv(x = x, file = file,
        row.names = FALSE,
        quote = TRUE,
        # sep = ", ",
        # dec = ".",
        na = '',
        ...)
}

#' @rdname prova.data
#'
#' @import utils
#'
#' @export
pread.csv <- function(file, ...){
    read.csv(file = file,
        na.strings = '',
        stringsAsFactors = FALSE,
        tryLogical = FALSE,
        sep = ",",
        dec = ".",
        ...)
}



## TO BE WRITTEN
## This function checks whether data and metadata are mutually consistent.
## Below are some snippets of checks that were included in other scripts

    ## if (xinfo$type == 'binary') { # seems binary variate
    ##   if (length(unique(x)) != 2) {
    ##     cat('Warning: inconsistencies with variate', xn, '\n')
##   }
## }





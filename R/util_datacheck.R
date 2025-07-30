#' @name inferno.data
#' @rdname inferno.data
#'
#' @title Write and read data in **inferno**
#'
#' @description Utility functions to read and write CSV files in the format required by **inferno**
#'
#' @param x The object to be written, preferably a matrix or data frame. If not, it is attempted to coerce `x` to a data frame. See [utils::write.table()].
#'
#' @param file Either a character naming a file or a connection open for writing or reading. See [utils::write.table()] and [utils::read.table()].
#' @param ... Other parameters to be passed to [utils::write.table()] or [utils::read.table()].
#'
#' @details
#' The functions of the **inferno** package accept CSV files formatted as follows:
#'
#' - Decimal values should be separated by a *dot*; no comma should be used to separate thousands etc. Example: `86342.75` .
#' - Character and names should be quoted in single or double quotes. Example: `"female"`.
#' - Values should be separated by *commas*, not by tabs or semicolons.
#' - Missing values should be simply *empty*, not denoted by "NA", "missing", "-", or similar.
#' - Preferably there should not be [`base::factor`]s; use character names instead.
#'
#' The utility functions [write.csvi()] and [read.csvi()] are wrappers to [utils::write.table()] or [utils::read.table()] that sets appropriate default parameters
NULL

#' @rdname inferno.data
#' @export
write.csvi <- function(x, file, ...){
    write.csv(x = x, file = file,
        row.names = FALSE,
        quote = TRUE,
        ## sep = ",",
        ## dec = ".",
        na = '',
        ...)
}

#' @rdname inferno.data
#' @export
read.csvi <- function(file, ...){
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





#' Subset variates of an object of class `probability`
#'
#' An object of class `probability`, calculated by the \code{\link{Pr}} or \code{\link{tailPr}} functions, holds the probabilities for all possible combinations of values of a set of joint variates `Y` conditional on a set of joint variates `X`. In some cases one may wish to exclude some of the values of the `Y` or `X` variates. For instance `Y` in the probability-class object could include the variate "age" with values from 18 to 100, and one may want to retain the values from 60 to 80.
#'
#' @param x object of class "probability", obtained with \code{\link{Pr}} or \code{\link{tailPr}}.
#' @param subset named list or named vector: variates to subset, given as list names, and corresponding values to subset.
#'
#' @return A list of class `probability`, identical to the original object `x` except for a reduced range of values in some if its variates.
#'
#' @export
subset.probability <- function(
    x,
    subset
){
    Ynames <- names(x$Y)
    Xnames <- names(x$X)
    subset <- as.list(subset)
    vrtnames <- names(subset)
    
    if(!all(vrtnames %in% c(Ynames, Xnames))){
        stop("probability object does not contain some of the given variates")
    }
    
    ## subset Y
    for(vrt in vrtnames[vrtnames %in% Ynames]){
        selvals <- x$Y[[vrt]] %in% subset[[vrt]]
        x$values <- x$values[selvals, , drop = FALSE]
        if(!is.null(x$quantiles)){
            x$quantiles <- x$quantiles[selvals, , , drop = FALSE]
        }
        if(!is.null(x$samples)){
            x$samples <- x$samples[selvals, , , drop = FALSE]
        }
        x$Y <- x$Y[selvals, , drop = FALSE]
    }

    ## subset X
    for(vrt in vrtnames[vrtnames %in% Xnames]){
        selvals <- x$X[[vrt]] %in% subset[[vrt]]
        x$values <- x$values[, selvals, drop = FALSE]
        if(!is.null(x$quantiles)){
            x$quantiles <- x$quantiles[, selvals, , drop = FALSE]
        }
        if(!is.null(x$samples)){
            x$samples <- x$samples[, selvals, , drop = FALSE]
        }
        x$X <- x$X[selvals, , drop = FALSE]
    }

    x
}

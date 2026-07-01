#' Subset variates of an object of class "probability"
#'
#' An object of class "probability", obtained with the [Pr()] function, holds the probabilities for all possible combinations of values of a set of joint variates `Y` conditional on a set of joint variates `X`, together with the variabilities of these probabilities and some other information. In some cases one may wish to exclude some of the values of the `Y` or `X` variates. For instance `Y` in the probability-class object could include the variate "age" with values from 18 to 100, and one may want to retain the values from 60 to 80.
#'
#' @param x Object of class "probability", obtained with [Pr()].
#' @param subset Named list or named vector: variates to subset, given as list names, and corresponding values to subset.
#'
#' @return An object of class "probability", identical to the original object `x` except for a reduced range of values in some if its variates.
#'
## #' @seealso
## #' [Pr()], which generates probability objects.
## #'
## #' [plot.probability()] to plot probabilities and quantiles calculated by `Pr()'.
## #'
## #' [hist.probability()] to plot histograms of the probability distributions calculated by `Pr()`.
## #'
## #' @examples
## #' ## Load the example `learnt` object calculated from the "penguins" dataset;
## #' ## variates: 'species' and 'bill_len'
## #' learnt <- learntExample
## #'
## #' ## Calculate the probability object for the three values of variate 'species',
## #' ## given values 43 and 44 of variate 'bill_len';
## #' ## this object contains probabilities, quantiles, and other information
## #' probs <- Pr(
## #'   Y = data.frame(species = c('Adelie', 'Chinstrap', 'Gentoo')),
## #'   X = data.frame(bill_len = c(43, 44)),
## #'   learnt = learnt, parallel = 1
## #' )
## #'
## #' probs$values
## #'
## #' ## Subset by retaining the values 'Adelie' and 'Gentoo' for species,
## #' ## and 44 for bill length
## #' newprobs <- prsubset(
## #'   probs,
## #'   subset = list(species = c('Adelie', 'Gentoo'), bill_len = 43)
## #' )
## #'
## #' newprobs$values
## #'
## #' ## Plot these conditional probabilities and their variabilities
## #' plot(newprobs)
## #'
## #' hist(newprobs)
## #'
#' @keywords internal
prsubset <- function(
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
        if(!is.null(x$values.MCaccuracy)){
            x$values.MCaccuracy <- x$values.MCaccuracy[selvals, , drop = FALSE]
        }
        if(!is.null(x$quantiles)){
            x$quantiles <- x$quantiles[selvals, , , drop = FALSE]
        }
        if(!is.null(x$quantiles.MCaccuracy)){
            x$quantiles.MCaccuracy <- x$quantiles.MCaccuracy[selvals, , , drop = FALSE]
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
        if(!is.null(x$values.MCaccuracy)){
            x$values.MCaccuracy <- x$values.MCaccuracy[, selvals, drop = FALSE]
        }
        if(!is.null(x$quantiles)){
            x$quantiles <- x$quantiles[, selvals, , drop = FALSE]
        }
        if(!is.null(x$quantiles.MCaccuracy)){
            x$quantiles.MCaccuracy <- x$quantiles.MCaccuracy[, selvals, , drop = FALSE]
        }
        if(!is.null(x$samples)){
            x$samples <- x$samples[, selvals, , drop = FALSE]
        }
        x$X <- x$X[selvals, , drop = FALSE]
    }

    x
}

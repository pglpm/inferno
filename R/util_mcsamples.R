#' Eliminate samples from mcsamples object
#'
#' Used in 'learn()'.
#'
#' @keywords internal
mcsubset <- function(learnt, subsamples) {
    lapply(learnt, function(xx) {
        do.call('[', c(
            list(xx),
            rep(TRUE, length(dim(xx)) - 1),
            list(subsamples),
            list(drop = FALSE))
        )
    })
}


#' Concatenate mcsample objects
#'
#' Used in 'learn()'.
#'
#' @keywords internal
mcjoin <- function(x, y){
    if(is.null(x)){
        y
    } else {
        mapply(
            function(xx, yy) {
                temp <- c(xx, yy)
                dx <- dim(yy)[-length(dim(yy))]
                dim(temp) <- c(dx, length(temp) / prod(dx))
                temp
            },
            x, y,
            SIMPLIFY = FALSE
        )
    }
}


#' Bind 3D arrays by first dimension
#'
#' Used in 'util_checkpoints()' within 'learn()', and in various functions in 'util_lprobs.R'.
#'
#' NB: the following variant is slower:
#'
#' ```
#' function(x, y) {
#'     out <- c(aperm(x), aperm(y))
#'     dim(out) <- c(rev(dim(x)[-1]), dim(x)[1] + dim(y)[1])
#'     aperm(out)
#' }
#' ```
#'
#' @keywords internal
learnbind <- function(x, y) {
    if(is.null(x)) {
        y
    } else {
        nrx <- dim(x)[1]
        nry <- dim(y)[1]
        out <- array(data = NA, dim = c(nrx + nry, dim(x)[-1]),
            dimnames = NULL)
        out[seq_len(nrx), ,] <- x
        out[nrx + seq_len(nry), ,] <- y
        out
    }
}

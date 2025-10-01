#' Eliminate samples from a 'learnt' object
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
#' @keywords internal
mcjoin <- function(x, y){
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

#' Bind 3D arrays by first dimension
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
## ## old, slower variant
## learnbind2 <- function(x, y) {
##     out <- c(aperm(x), aperm(y))
##     dim(out) <- c(rev(dim(x)[-1]), dim(x)[1] + dim(y)[1])
##     aperm(out)
## }


#' Cumulative sum along first dimension
#'
#' @keywords internal
rowcumsum <- function(x){
    for(i in 2:(dim(x)[1])){
        x[i,,] <- x[i,,] + x[i-1,,]
    }
    x
}

#' Inverse cumulative sum along first dimension
#'
#' @keywords internal
rowinvcumsum <- function(x){
    for(i in (dim(x)[1] - 1):1){
        x[i,,] <- x[i,,] + x[i+1,,]
    }
    x
}

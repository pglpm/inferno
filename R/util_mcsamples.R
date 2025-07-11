#### Eliminate samples from a 'learnt' object
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

#### Concatenate mcsample objects
#' @keywords internal
mcjoin <- function(x, y){
    mapply(
        function(xx, yy) {
            temp <- c(xx, yy)
            dx <- dim(xx)[-length(dim(xx))]
            dim(temp) <- c(dx, length(temp) / prod(dx))
            temp
        },
        x, y,
        SIMPLIFY = FALSE
    )
}

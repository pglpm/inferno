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

#### Bind 3D arrays by first dimension
#' @keywords internal
learnbind <- function(x, y) {
    out <- c(aperm(x), aperm(y))
    dim(out) <- c(dim(x)[-1], dim(x)[1] + dim(y)[1])
    aperm(out)
}

#### Cumulative sum along first dimension
#' @keywords internal
rowcumsum <- function(x){
    for(i in 2:(dim(x)[1])){
        x[i,,] <- x[i,,] + x[i-1,,]
    }
    x
}

#### Inverseumulative sum along first dimension
#' @keywords internal
rowinvcumsum <- function(x){
    for(i in (dim(x)[1] - 1):1){
        x[i,,] <- x[i,,] + x[i+1,,]
    }
    x
}

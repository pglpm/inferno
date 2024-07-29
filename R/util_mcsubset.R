#### Eliminate samples from an mcoutput object
mcsubset <- function(mcoutput, subsamples) {
    lapply(mcoutput, function(xx) {
        do.call('[', c(
            list(xx),
            rep(TRUE, length(dim(xx)) - 1),
            list(subsamples),
            list(drop = FALSE))
        )
    })
}

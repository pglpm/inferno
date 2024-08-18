#### Eliminate samples from an learned object
mcsubset <- function(learned, subsamples) {
    lapply(learned, function(xx) {
        do.call('[', c(
            list(xx),
            rep(TRUE, length(dim(xx)) - 1),
            list(subsamples),
            list(drop = FALSE))
        )
    })
}

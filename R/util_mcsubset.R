#### Eliminate samples from an learnt object
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

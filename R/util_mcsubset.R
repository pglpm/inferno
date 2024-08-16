#### Eliminate samples from an agent object
mcsubset <- function(agent, subsamples) {
    lapply(agent, function(xx) {
        do.call('[', c(
            list(xx),
            rep(TRUE, length(dim(xx)) - 1),
            list(subsamples),
            list(drop = FALSE))
        )
    })
}

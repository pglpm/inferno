#' Eliminate samples from a learnt object
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

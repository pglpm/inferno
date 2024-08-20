#### Function for calculating MCMC standard error
#' @keywords internal
funMCSE <- function(x) {
    N <- nrow(x)
    b <- floor(sqrt(N))
    a <- floor(N/b)
    Ys <- rbind(sapply(seq_len(a), function(k) {
        colMeans(x[((k - 1) * b + 1):(k * b), , drop = FALSE])
    }))
    ##
    sqrt(b * rowSums((Ys - colMeans(x))^2) / ((a - 1) * N))
}

#### Function for calculating number of needed MCMC iterations
#' @keywords internal
mcmcstop <- function(
    traces,
    nsamples,
    availiter,
    relerror,
    ## diagnESS,
    ## diagnIAT,
    ## diagnBMK,
    ## diagnMCSE,
    ## diagnStat,
    ## diagnBurn,
    ## diagnBurn2,
    ## diagnThin,
    thinning
) {
    ## Based on doi.org/10.1080/10618600.2015.1044092

    relmcse <- funMCSE(traces) / apply(traces, 2, sd)
    ess <- 1/relmcse^2
    maxrelmcse <- max(relmcse)
    autothinning <- ceiling(1.5 * nrow(traces)/ess)

    if(is.null(thinning)) {
        thinning <- max(autothinning)
    }
    missingsamples <- thinning * (nsamples - 1) - availiter

    if(max(relmcse) < relerror) {
        ## sampling could be stopped
        reqiter <- 0
    } else {
        ## sampling should continue
        reqiter <- ceiling(thinning * sqrt(nsamples))
    }

    if(missingsamples > 0) {
        ## not enough samples
        reqiter <- max(reqiter, missingsamples)
    }

    list(
        reqiter = reqiter,
        proposed.thinning = thinning,
        toprint = list(
            'rel. MC standard error' = relmcse,
            'eff. sample size' = ess,
            'needed thinning' = autothinning
        )
    )
}

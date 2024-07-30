#### Function for calculating MCMC standard error
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
    if(is.null(thinning)) {
        autothinning <- ceiling(1.5 * nrow(traces)/min(ess))
    } else {
        autothinning <- thinning
    }
    missingsamples <- autothinning * (nsamples - 1) - availiter

    if(max(relmcse) < relerror) {
        ## sampling could be stopped
        reqiter <- 0
    } else {
        ## sampling should continue
        reqiter <- ceiling(autothinning * sqrt(nsamples))
    }

    if(missingsamples > 0) {
        ## not enough samples
        reqiter <- max(reqiter, missingsamples)
    }

    list(
        reqiter = reqiter,
        proposed.thinning = autothinning,
        toprint = list(
            'rel. MC standard error' = relmcse,
            'eff. sample size' = ess
        )
        )
}

## oldmcmclength <- function(nsamplesperchain, nitertot, thinning,
##                         diagnESS, diagnIAT, diagnBMK, diagnMCSE,
##                         diagnStat, diagnBurn, diagnBurn2, diagnThin) {
##   ## This function uses (or can potentially use) various diagnostics from
##   ## the LaplacesDemon package in order to determine
##   ## the total number of needed MCMC iterations as well as
##   ## the needed thinning
##   ## Current method: no. of iterations is given by:
##   ## 2 * diagnosed burn-in +
##   ## 2 * (maximum between needed thinning and integrated autocorrelation time) *
##   ## number of independent samples needed (minus last one)
##   list(reqiter = ceiling(3 * max(diagnBurn2) +
##           ( 2 * ceiling(max(diagnIAT, diagnThin)) *
##            (nsamplesperchain - 1L) ) ),
##        thinning = (if (is.null(thinning)) {
##          2 * ceiling(max(diagnIAT, diagnThin))
##        } else {
##          thinning
##        })
##        )
## }

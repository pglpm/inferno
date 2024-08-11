#' Calculate joint frequencies for checkpoints in inferpopulation()
#'
#' @param Y matrix or data.table: values of some already-transformed
#'   variates of which we want the joint probability; one variate per column
#' @param mcsamples object internal to `inferpopulation()`,
#'   containing partial Monte Carlo draws
#' @param auxmetadata object internal to `inferpopulation()`,
#'   containing processed metadata information
#'
#' @return The joint frequencies of Y correspoinding to the Monte Carlo samples
#'
#' @import foreach
Pcheckpoints <- function(
    Y,
    mcsamples,
    auxmetadata
) {

    nsamples <- ncol(mcsamples$W)
    ncomponents <- nrow(mcsamples$W)

### Guide to indices:
    ## .i. = order in Y corresponding to appearance in vnames
    ## .t. = vnames present in Y, kept in their vnames-order

    Yv <- colnames(Y)

#### Type R
    vnames <- auxmetadata[auxmetadata$mcmctype == 'R', 'name']
    YiR <- match(vnames, Yv)
    YtR <- which(!is.na(YiR))
    YiR <- YiR[YtR]
    YnR <- length(YiR)
    if (YnR > 0) {
        mcsamples$Rvar <- sqrt(mcsamples$Rvar)
    }

#### Type C
    vnames <- auxmetadata[auxmetadata$mcmctype == 'C', 'name']
    YiC <- match(vnames, Yv)
    YtC <- which(!is.na(YiC))
    YiC <- YiC[YtC]
    YnC <- length(YiC)
    if (YnC > 0) {
        mcsamples$Cvar <- sqrt(mcsamples$Cvar)
        Clefts <- auxmetadata[match(vnames, auxmetadata$name), 'tleftbound']
        Crights <- auxmetadata[match(vnames, auxmetadata$name), 'trightbound']
    }

#### Type D
    vnames <- auxmetadata[auxmetadata$mcmctype == 'D', 'name']
    YiD <- match(vnames, Yv)
    YtD <- which(!is.na(YiD))
    YiD <- YiD[YtD]
    YnD <- length(YiD)
    if (YnD > 0) {
        mcsamples$Dvar <- sqrt(mcsamples$Dvar)
        Dsteps <- auxmetadata[match(vnames, auxmetadata$name), 'halfstep'] /
            auxmetadata[match(vnames, auxmetadata$name), 'tscale']
        Dlefts <- auxmetadata[match(vnames, auxmetadata$name), 'tleftbound']
        Drights <- auxmetadata[match(vnames, auxmetadata$name), 'trightbound']
    }

#### Type O
    vnames <- auxmetadata[auxmetadata$mcmctype == 'O', 'name']
    YiO <- match(vnames, Yv)
    YtO <- which(!is.na(YiO))
    YiO <- YiO[YtO]
    YnO <- length(YiO)

#### Type N
    vnames <- auxmetadata[auxmetadata$mcmctype == 'N', 'name']
    YiN <- match(vnames, Yv)
    YtN <- which(!is.na(YiN))
    YiN <- YiN[YtN]
    YnN <- length(YiN)

#### Type B
    vnames <- auxmetadata[auxmetadata$mcmctype == 'B', 'name']
    YiB <- match(vnames, Yv)
    YtB <- which(!is.na(YiB))
    YiB <- YiB[YtB]
    YnB <- length(YiB)


    foreach(y = t(Y), .combine = `cbind`) %do% {

#### the loop is over the columns of y
#### each instance is a 1-column vector
        if (all(is.na(y))) {
            lprobY <- array(NA, dim = c(ncomponents, nsamples))
        } else {
            lprobY <- (if (YnR > 0) { # continuous
                colSums(
                    dnorm(
                        x = y[YiR, ],
                        mean = mcsamples$Rmean[YtR, , , drop = FALSE],
                        sd = mcsamples$Rvar[YtR, , , drop = FALSE],
                        log = TRUE
                    ), na.rm = TRUE)
            } else {
                0
            }) +
                (if (YnC > 0) { # censored
                    isfin <- is.finite(y[YiC, ])
                    indf <- which(isfin)
                    indi <- which(!isfin)
                    (if (length(indf) > 0) {
                        colSums(
                            dnorm(
                                x = y[YiC[indf], ],
                                mean = mcsamples$Cmean[YtC[indf], , , drop = FALSE],
                                sd = mcsamples$Cvar[YtC[indf], , , drop = FALSE],
                                log = TRUE
                            ), na.rm = TRUE)
                    } else {
                        0
                    }) +
                        (if (length(indi) > 0) {
                            vt <- YtC[indi]
                            vx <- pmin(
                                pmax(Clefts[vt], y[YiC[indi], ]),
                                Crights[vt])
                            ## for upper tail, take opposite mean and value
                            colSums(
                                pnorm(
                                    q = -abs(vx),
                                    mean = -sign(vx) *
                                        mcsamples$Cmean[vt, , , drop = FALSE],
                                    sd = mcsamples$Cvar[vt, , , drop = FALSE],
                                    log.p = TRUE
                                ), na.rm = TRUE)
                        } else {
                            0
                        })
                } else {
                    0
                }) +
                (if (YnD > 0) { # discretized
                    vrights <- y[YiD, ] + Dsteps[YtD]
                    vrights[vrights >= Drights[YtD]] <- +Inf
                    vlefts <- y[YiD, ] - Dsteps[YtD]
                    vlefts[vlefts <= Dlefts[YtD]] <- -Inf
                    colSums(log(
                        pnorm(
                            q = vrights,
                            mean = mcsamples$Dmean[YtD, , , drop = FALSE],
                            sd = mcsamples$Dvar[YtD, , , drop = FALSE]
                        ) -
                            pnorm(
                                q = vlefts,
                                mean = mcsamples$Dmean[YtD, , , drop = FALSE],
                                sd = mcsamples$Dvar[YtD, , , drop = FALSE]
                            )
                    ), na.rm = TRUE)
                } else {
                    0
                }) +
                (if (YnO > 0) { # nominal
                    colSums(log(
                        aperm(
                            vapply(seq_len(YnO), function(v) {
                                mcsamples$Oprob[YtO[v], , y[YiO[v], ], ]
                            }, mcsamples$W),
                            c(3, 1, 2))
                    ), na.rm = TRUE)
                } else {
                    0
                }) +
                (if (YnN > 0) { # nominal
                    colSums(log(
                        aperm(
                            vapply(seq_len(YnN), function(v) {
                                mcsamples$Nprob[YtN[v], , y[YiN[v], ], ]
                            }, mcsamples$W),
                            c(3, 1, 2))
                    ), na.rm = TRUE)
                } else {
                    0
                }) +
                (if (YnB > 0) { # binary
                    colSums(log(
                    (y[YiB, ] * mcsamples$Bprob[YtB, , , drop = FALSE]) +
                        ((1 - y[YiB, ]) *
                             (1 - mcsamples$Bprob[YtB, , , drop = FALSE]))
                    ), na.rm = TRUE)
                } else {
                    0
                })
        }
#### Output: rows=components, columns=samples
        lprobX <- apply(log(mcsamples$W), 2, function(xx) {
            xx - max(xx[is.finite(xx)])
        })
        colSums(exp(lprobX + lprobY)) / colSums(exp(lprobX))
    }
}

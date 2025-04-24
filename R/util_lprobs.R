#' Calculate collection of log-probabilities for different components and samples
#'
#' @param X numerical matrix: transformed variates
#' @param learnt: Monte-Carlo-output object
#' @param nR etc: Parameters containing appropriate indices
#'
#' @return Matrix with as many rows as components and as many cols as samples
#' @keywords internal
util_lprobs <- function(
    x,
    learnt,
    ##
    nR, iR, tR,
    nLR, iLR, tLR,
    nUR, iUR, tUR,
    ##
    nC, iC, tC,
    nLC, iLC, tLC,
    nUC, iUC, tUC,
    Clefts, Crights,
    ##
    nD, iD, tD,
    nLD, iLD, tLD,
    nUD, iUD, tUD,
    Dsteps, Dlefts, Drights,
    ##
    nO, iO, tO,
    nLO, iLO, tLO,
    nUO, iUO, tUO,
    ##
    nN, iN, tN,
    nB, iB, tB
) {
    (
        (if (nR > 0) { # continuous
            colSums(
                dnorm(
                    x = x[iR, ],
                    mean = learnt$Rmean[tR, , , drop = FALSE],
                    sd = learnt$Rvar[tR, , , drop = FALSE],
                    log = TRUE
                ), na.rm = TRUE)
        } else {
            0
        }) +
            (if (nC > 0) { # censored
                isfin <- is.finite(x[iC, ])
                indf <- which(isfin)
                indi <- which(!isfin)
                (if (length(indf) > 0) {
                    colSums(
                        dnorm(
                            x = x[iC[indf], ],
                            mean = learnt$Cmean[tC[indf], , , drop = FALSE],
                            sd = learnt$Cvar[tC[indf], , , drop = FALSE],
                            log = TRUE
                        ), na.rm = TRUE)
                } else {
                    0
                }) +
                    (if (length(indi) > 0) {
                        vt <- tC[indi]
                        vx <- pmin(
                            pmax(Clefts[vt], x[iC[indi], ]),
                            Crights[vt])
                        ## for upper tail, take opposite mean and value
                        colSums(
                            pnorm(
                                q = -abs(vx),
                                mean = -sign(vx) *
                                    learnt$Cmean[vt, , , drop = FALSE],
                                sd = learnt$Cvar[vt, , , drop = FALSE],
                                log.p = TRUE
                            ), na.rm = TRUE)
                    } else {
                        0
                    })
            } else {
                0
            }) +
            (if (nD > 0) { # discretized
                vrights <- pmax(x[iD, ] + Dsteps[tD], Dlefts[tD])
                vrights[vrights - Dsteps[tD] >= Drights[tD]] <- +Inf
                vlefts <- pmin(x[iD, ] - Dsteps[tD], Drights[tD])
                vlefts[vlefts + Dsteps[tD] <= Dlefts[tD]] <- -Inf
                colSums(log(
                    pnorm(
                        q = vrights,
                        mean = learnt$Dmean[tD, , , drop = FALSE],
                        sd = learnt$Dvar[tD, , , drop = FALSE]
                    ) -
                        pnorm(
                            q = vlefts,
                            mean = learnt$Dmean[tD, , , drop = FALSE],
                            sd = learnt$Dvar[tD, , , drop = FALSE]
                        )
                ), na.rm = TRUE)
            } else {
                0
            }) +
            (if (nO > 0) { # nominal
                colSums(log(
                    aperm(
                        vapply(X = seq_len(nO), FUN = function(v) {
                            learnt$Oprob[tO[v], , x[iO[v], ], ]
                        }, FUN.VALUE = learnt$W),
                        c(3, 1, 2))
                ), na.rm = TRUE)
            } else {
                0
            }) +
            (if (nN > 0) { # nominal
                colSums(log(
                    aperm(
                        vapply(X = seq_len(nN), FUN = function(v) {
                            learnt$Nprob[tN[v], , x[iN[v], ], ]
                        }, FUN.VALUE = learnt$W),
                        c(3, 1, 2))
                ), na.rm = TRUE)
            } else {
                0
            }) +
            (if (nB > 0) { # binary
                ## Bprob is the probability that x=1
                colSums(log(
                    x[iB, ] * learnt$Bprob[tB, , , drop = FALSE] +
                        (1 - x[iB, ]) *
                        (1 - learnt$Bprob[tB, , , drop = FALSE])
                ), na.rm = TRUE)
            } else {
                0
            }) +
    ## lower-tails
        (if (nLR > 0) { # continuous
            colSums(
                pnorm(
                    q = x[iLR, ],
                    mean = learnt$Rmean[tLR, , , drop = FALSE],
                    sd = learnt$Rvar[tLR, , , drop = FALSE],
                    log.p = TRUE,
                    lower.tail= TRUE
                ), na.rm = TRUE)
        } else {
            0
        }) +
            (if (nLC > 0) { # censored
                if(lower.tail) {
                    vx <- pmax(Clefts[tLC], x[iLC, ])
                } else {
                    vx <- pmin(Crights[tLC], x[iLC, ])
                }
                colSums(
                    pnorm(
                        q = vx,
                        mean = learnt$Cmean[tLC, , , drop = FALSE],
                        sd = learnt$Cvar[tLC, , , drop = FALSE],
                        log.p = TRUE,
                        lower.tail= TRUE
                    ), na.rm = TRUE)
            } else {
                0
            }) +
            (if (nLD > 0) { # discretized
                ## if(lower.tail) {
                ##     vx <- pmax(x[iLD, ] + Dsteps[tLD], Dlefts[tLD])
                ##     vx[vx - Dsteps[tLD] >= Drights[tLD]] <- +Inf
                ## } else {
                ##     vx <- pmin(x[iLD, ] + Dsteps[tLD], Drights[tLD])
                ##     vx[vx + Dsteps[tLD] <= Dlefts[tLD]] <- -Inf
                ## }
                ## ## Old, before 'eq' argument
                ## vx <- pmax(x[iLD, ] + Dsteps[tLD], Dlefts[tLD])
                ## vx[vx - Dsteps[tLD] >= Drights[tLD]] <- +Inf
                vx <- x[iLD, ]
                vx[vx <= Dlefts[tLD]] <- -Inf
                vx[vx >= Drights[tLD]] <- +Inf
                colSums(
                    pnorm(
                        q = vx,
                        mean = learnt$Dmean[tLD, , , drop = FALSE],
                        sd = learnt$Dvar[tLD, , , drop = FALSE],
                        log.p = TRUE,
                        lower.tail= TRUE
                    ), na.rm = TRUE)
            } else {
                0
            }) +
            (if (nLO > 0) { # nominal
                seqO <- seq_len(dim(learnt$Oprob)[3])
                colSums(log(
                    aperm(
                        vapply(seq_len(nLO), function(v) {
                            if(lower.tail) {
                                ## vx <- seq_len(x[iLO[v], ])
                                vx <- (seqO <= x[iLO[v], ])
                            } else {
                                ## vx <- seq(1+x[iLO[v], ], 20)
                                vx <- (seqO > x[iLO[v], ])
                            }
                            apply(X = learnt$Oprob[tLO[v], , vx, , drop = F],
                                MARGIN = -3, FUN = sum, na.rm = TRUE)
                        }, learnt$W),
                        c(3, 1, 2))
                ), na.rm = TRUE)
            } else {
                0
            }) +
    ## upper-tails
        (if (nUR > 0) { # continuous
            colSums(
                pnorm(
                    q = x[iUR, ],
                    mean = learnt$Rmean[tUR, , , drop = FALSE],
                    sd = learnt$Rvar[tUR, , , drop = FALSE],
                    log.p = TRUE,
                    lower.tail= FALSE
                ), na.rm = TRUE)
        } else {
            0
        }) +
            (if (nUC > 0) { # censored
                if(lower.tail) {
                    vx <- pmax(Clefts[tUC], x[iUC, ])
                } else {
                    vx <- pmin(Crights[tUC], x[iUC, ])
                }
                colSums(
                    pnorm(
                        q = vx,
                        mean = learnt$Cmean[tUC, , , drop = FALSE],
                        sd = learnt$Cvar[tUC, , , drop = FALSE],
                        log.p = TRUE,
                        lower.tail= FALSE
                    ), na.rm = TRUE)
            } else {
                0
            }) +
            (if (nUD > 0) { # discretized
                ## if(lower.tail) {
                ##     vx <- pmax(x[iUD, ] + Dsteps[tUD], Dlefts[tUD])
                ##     vx[vx - Dsteps[tUD] >= Drights[tUD]] <- +Inf
                ## } else {
                ##     vx <- pmin(x[iUD, ] + Dsteps[tUD], Drights[tUD])
                ##     vx[vx + Dsteps[tUD] <= Dlefts[tUD]] <- -Inf
                ## }
                ## ## Old, before 'eq' argument
                ## vx <- pmax(x[iUD, ] + Dsteps[tUD], Dlefts[tUD])
                ## vx[vx - Dsteps[tUD] >= Drights[tUD]] <- +Inf
                vx <- x[iUD, ]
                vx[vx <= Dlefts[tUD]] <- -Inf
                vx[vx >= Drights[tUD]] <- +Inf
                colSums(
                    pnorm(
                        q = vx,
                        mean = learnt$Dmean[tUD, , , drop = FALSE],
                        sd = learnt$Dvar[tUD, , , drop = FALSE],
                        log.p = TRUE,
                        lower.tail= FALSE
                    ), na.rm = TRUE)
            } else {
                0
            }) +
            (if (nUO > 0) { # nominal
                seqO <- seq_len(dim(learnt$Oprob)[3])
                colSums(log(
                    aperm(
                        vapply(seq_len(nUO), function(v) {
                            if(lower.tail) {
                                ## vx <- seq_len(x[iUO[v], ])
                                vx <- (seqO <= x[iUO[v], ])
                            } else {
                                ## vx <- seq(1+x[iUO[v], ], 20)
                                vx <- (seqO > x[iUO[v], ])
                            }
                            apply(X = learnt$Oprob[tUO[v], , vx, , drop = F],
                                MARGIN = -3, FUN = sum, na.rm = TRUE)
                        }, learnt$W),
                        c(3, 1, 2))
                ), na.rm = TRUE)
            } else {
                0
            })
    )
}


#' Calculate collection of log-probabilities for different components and samples
#'
#' @param X numerical matrix: transformed variates
#' @param agent: Monte-Carlo-output object
#' @param nR etc: Parameters containing appropriate indices
#'
#' @return Matrix with as many rows as components and as many cols as samples
util_lprob <- function(
    x,
    agent,
    nR, iR, tR,
    nC, iC, tC, Clefts, Crights,
    nD, iD, tD, Dsteps, Dlefts, Drights,
    nO, iO, tO,
    nN, iN, tN,
    nB, iB, tB
) {
    (
        (if (nR > 0) { # continuous
            colSums(
                dnorm(
                    x = x[iR, ],
                    mean = agent$Rmean[tR, , , drop = FALSE],
                    sd = agent$Rvar[tR, , , drop = FALSE],
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
                            mean = agent$Cmean[tC[indf], , , drop = FALSE],
                            sd = agent$Cvar[tC[indf], , , drop = FALSE],
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
                                    agent$Cmean[vt, , , drop = FALSE],
                                sd = agent$Cvar[vt, , , drop = FALSE],
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
                vrights[vrights >= Drights[tD]] <- +Inf
                vlefts <- pmin(x[iD, ] - Dsteps[tD], Drights[tD])
                vlefts[vlefts <= Dlefts[tD]] <- -Inf
                colSums(log(
                    pnorm(
                        q = vrights,
                        mean = agent$Dmean[tD, , , drop = FALSE],
                        sd = agent$Dvar[tD, , , drop = FALSE]
                    ) -
                        pnorm(
                            q = vlefts,
                            mean = agent$Dmean[tD, , , drop = FALSE],
                            sd = agent$Dvar[tD, , , drop = FALSE]
                        )
                ), na.rm = TRUE)
            } else {
                0
            }) +
            (if (nO > 0) { # nominal
                colSums(log(
                    aperm(
                        vapply(seq_len(nO), function(v) {
                            agent$Oprob[tO[v], , x[iO[v], ], ]
                        }, agent$W),
                        c(3, 1, 2))
                ), na.rm = TRUE)
            } else {
                0
            }) +
            (if (nN > 0) { # nominal
                colSums(log(
                    aperm(
                        vapply(seq_len(nN), function(v) {
                            agent$Nprob[tN[v], , x[iN[v], ], ]
                        }, agent$W),
                        c(3, 1, 2))
                ), na.rm = TRUE)
            } else {
                0
            }) +
            (if (nB > 0) { # binary
                ## Bprob is the probability that x=1
                colSums(log(
                    x[iB, ] * agent$Bprob[tB, , , drop = FALSE] +
                        (1 - x[iB, ]) *
                        (1 - agent$Bprob[tB, , , drop = FALSE])
                ), na.rm = TRUE)
            } else {
                0
            })
    )
}

#' Calculate collection of log-probabilities for different components and samples
#'
#' @param X Numerical matrix: transformed variates
#' @param learnt: Monte-Carlo-output object
#' @param nR etc: Parameters containing appropriate indices
#'
#' @return Matrix with as many rows as components and as many cols as samples
#' @keywords internal
util_lprob <- function(
    x,
    learnt,
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
                    learnt$Oprob[x[iO], , , drop = FALSE]
                ), na.rm = TRUE)
            } else {
                0
            }) +
            (if (nN > 0) { # nominal
                colSums(log(
                    learnt$Nprob[x[iN], , , drop = FALSE]
                ), na.rm = TRUE)
            } else {
                0
            }) +
            (if (nB > 0) { # binary
                ## Bprob is the probability that x=1
                pv <- learnt$Bprob[tB, , , drop = FALSE]
                colSums(log(
                    x[iB, ] * pv + (1 - x[iB, ]) * (1 - pv)
                ), na.rm = TRUE)
            } else {
                0
            })
    )
}

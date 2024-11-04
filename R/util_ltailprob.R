#' Calculate collection of log-probabilities for different components and samples
#'
#' @param X numerical matrix: transformed variates
#' @param learnt: Monte-Carlo-output object
#' @param nR etc: Parameters containing appropriate indices
#'
#' @return Matrix with as many rows as components and as many cols as samples
#' @keywords internal
util_ltailprob <- function(
    x,
    learnt,
    nR, iR, tR,
    nC, iC, tC, Clefts, Crights,
    nD, iD, tD, Dsteps, Dlefts, Drights,
    nO, iO, tO,
    ## nN, iN, tN,
    ## nB, iB, tB,
    lower.tail = TRUE
) {
    (
        (if (nR > 0) { # continuous
            colSums(
                pnorm(
                    q = x[iR, ],
                    mean = learnt$Rmean[tR, , , drop = FALSE],
                    sd = learnt$Rvar[tR, , , drop = FALSE],
                    log.p = TRUE,
                    lower.tail = lower.tail
                ), na.rm = TRUE)
        } else {
            0
        }) +
            (if (nC > 0) { # censored
                if(lower.tail) {
                    vx <- pmax(Clefts[tC], x[iC, ])
                } else {
                    vx <- pmin(Crights[tC], x[iC, ])
                }
                colSums(
                    pnorm(
                        q = vx,
                        mean = learnt$Cmean[tC, , , drop = FALSE],
                        sd = learnt$Cvar[tC, , , drop = FALSE],
                        log.p = TRUE,
                        lower.tail = lower.tail
                    ), na.rm = TRUE)
            } else {
                0
            }) +
            (if (nD > 0) { # discretized
                ## if(lower.tail) {
                ##     vx <- pmax(x[iD, ] + Dsteps[tD], Dlefts[tD])
                ##     vx[vx - Dsteps[tD] >= Drights[tD]] <- +Inf
                ## } else {
                ##     vx <- pmin(x[iD, ] + Dsteps[tD], Drights[tD])
                ##     vx[vx + Dsteps[tD] <= Dlefts[tD]] <- -Inf
                ## }
                ## ## Old, before 'eq' argument
                ## vx <- pmax(x[iD, ] + Dsteps[tD], Dlefts[tD])
                ## vx[vx - Dsteps[tD] >= Drights[tD]] <- +Inf
                vx <- x[iD, ]
                vx[vx <= Dlefts[tD]] <- -Inf
                vx[vx >= Drights[tD]] <- +Inf
                colSums(
                    pnorm(
                        q = vx,
                        mean = learnt$Dmean[tD, , , drop = FALSE],
                        sd = learnt$Dvar[tD, , , drop = FALSE],
                        log.p = TRUE,
                        lower.tail = lower.tail
                    ), na.rm = TRUE)
            } else {
                0
            }) +
            (if (nO > 0) { # nominal
                seqO <- seq_len(dim(learnt$Oprob)[3])
                colSums(log(
                    aperm(
                        vapply(seq_len(nO), function(v) {
                            if(lower.tail) {
                                ## vx <- seq_len(x[iO[v], ])
                                vx <- (seqO <= x[iO[v], ])
                            } else {
                                ## vx <- seq(1+x[iO[v], ], 20)
                                vx <- (seqO > x[iO[v], ])
                            }
                            apply(X = learnt$Oprob[tO[v], , vx, , drop = F],
                                MARGIN = -3, FUN = sum, na.rm = TRUE)
                        }, learnt$W),
                        c(3, 1, 2))
                ), na.rm = TRUE)
            } else {
                0
            })
    )
}


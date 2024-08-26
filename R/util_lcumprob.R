#' Calculate collection of log-probabilities for different components and samples
#'
#' @param X numerical matrix: transformed variates
#' @param learnt: Monte-Carlo-output object
#' @param nR etc: Parameters containing appropriate indices
#'
#' @return Matrix with as many rows as components and as many cols as samples
#' @keywords internal
util_lcumprob <- function(
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
                        mean = learnt$Cmean[vt, , , drop = FALSE],
                        sd = learnt$Cvar[vt, , , drop = FALSE],
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
                vx <- pmax(x[iD, ] + Dsteps[tD], Dlefts[tD])
                vx[vx - Dsteps[tD] >= Drights[tD]] <- +Inf
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
                colSums(log(
                    aperm(
                        vapply(seq_len(nO), function(v) {
                            if(lower.tail) {
                                vx <- seq_len(x[iO[v], ])
                            } else {
                                vx <- -seq_len(x[iO[v], ])
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






## lower.tail <- T
## testv <- vapply(seq_len(3), function(v) {
##     if(lower.tail) {
##         vx <- seq_len(v+1)
##     } else {
##         vx <- -seq_len(v - 1L)
##     }
##     str(vx)
##     out <- apply(testl2$Nprob[v, , vx, , drop = F],
##         c(1,2,4), sum, na.rm = TRUE)
##     str(out)
##     out
## }, testl2$W)
##
## lower.tail <- F
## testv2 <- vapply(seq_len(3), function(v) {
##     if(lower.tail) {
##         vx <- seq_len(v)
##     } else {
##         vx <- -seq_len(v)
##     }
##     str(vx)
##     out <- apply(testl2$Nprob[v, , vx, , drop = F],
##         -3, sum, na.rm = TRUE)
##     str(out)
##     out
## }, testl2$W)
##
## testv <- vapply(seq_len(3), function(v) {
##     if(lower.tail) {
##         vx <- seq_len(v)
##     } else {
##         vx <- -seq_len(v - 1L)
##     }
##     str(vx)
##     out <- testl2$Nprob[v, , v, ]
##     str(out)
##     out
## }, testl2$W)

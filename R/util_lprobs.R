#' Calculate collection of log-probabilities for different components and samples
#'
#' @param X numerical matrix: transformed variates
#' @param learnt: Monte-Carlo-output object
#' @param nR etc: Parameters containing appropriate indices
#'
#' @return Matrix with as many rows as components and as many cols as samples
#' @keywords internal
util_lprobs <- function(
    nV0, xV0, V0mean, V0sd,
    nV1, xV1, V1mean, V1sd,
    nV2, xV2, V2mean, V2sd, V2steps,
    nVN, VNprobs, xVN,
    nVB, VBprobs, xVB
) {
    out <- 0
    ## point probability density
    if(nV0) {
        out <- out + colSums(
            x = dnorm(x = xV0, mean = V0mean, sd = V0sd, log = TRUE),
            na.rm = TRUE, dims = 1)
    },
    ## tail probability
    if(nV1) {
        out <- out + colSums(
            x = pnorm(q = xV1, mean = V1mean, sd = V1sd,
                log.p = TRUE, lower.tail = TRUE),
            na.rm = TRUE, dims = 1)
    },
    ## interval probability
    if(nV2) {
        out <- pnorm(q = xV2 - V2steps, mean = V2mean, sd = V2sd,
            log.p = TRUE, lower.tail = TRUE)
        ##
        out <- out + colSums(
            x = out + log(expm1(
                pnorm(q = xV2 + V2steps, mean = V2mean, sd = V2sd,
                    log.p = TRUE, lower.tail = TRUE) - out
            )),
            na.rm = TRUE, dims = 1)
    },
    ##
    if(nVN) {
        out <- out + colSums(
            x = log(VNprobs[xVN, , , drop = FALSE]),
            na.rm = TRUE, dims = 1)
    }
    ##
    if(nVB) {
        ## VBprob is the probability that x = 1
        out <- out + colSums(
            x = log(1 - xVB - VBprobs + 2 * xVB * VBprobs),
            na.rm = TRUE, dims = 1)
    }
    out
}

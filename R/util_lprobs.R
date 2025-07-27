#' Calculate collection of log-probabilities for different components and samples
#'
#' @param X numerical matrix: transformed variates
#' @param learnt: Monte-Carlo-output object
#' @param nR etc: Parameters containing appropriate indices
#'
#' @return Matrix with as many rows as components and as many cols as samples
#' @keywords internal
util_lprobs <- function(
    nV0, V0mean, V0sd, xV0,
    nV1, V1mean, V1sd, xV1,
    nV2, V2mean, V2sd, V2steps, xV2,
    nVN, VNprobs, xVN,
    nVB, VBprobs, xVB
) {
    out <- 0
    ## point probability density
    if(nV0) {
        out <- out + colSums(
            x = dnorm(x = xV0, mean = V0mean, sd = V0sd, log = TRUE),
            na.rm = TRUE, dims = 1)
    }
    ## tail probability
    if(nV1) {
        out <- out + colSums(
            x = pnorm(q = xV1, mean = V1mean, sd = V1sd,
                log.p = TRUE, lower.tail = TRUE),
            na.rm = TRUE, dims = 1)
    }
    ## interval probability
    if(nV2) {
        pright <- pnorm(q = xV2 + V2steps, mean = V2mean, sd = V2sd,
            log.p = TRUE, lower.tail = TRUE)
        ##
        out <- out + colSums(
            x = pright + log(-expm1(
                pnorm(q = xV2 - V2steps, mean = V2mean, sd = V2sd,
                    log.p = TRUE, lower.tail = TRUE) - pright
            )),
            na.rm = TRUE, dims = 1)
        ## ## this alternate form leads to infinities in some cases
        ## pleft <- pnorm(q = xV2 - V2steps, mean = V2mean, sd = V2sd,
        ##     log.p = TRUE, lower.tail = TRUE)
        ## ##
        ## out <- out + colSums(
        ##     x = pleft + log(expm1(
        ##         pnorm(q = xV2 + V2steps, mean = V2mean, sd = V2sd,
        ##             log.p = TRUE, lower.tail = TRUE) - pleft
        ##     )),
        ##     na.rm = TRUE, dims = 1)
    }
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

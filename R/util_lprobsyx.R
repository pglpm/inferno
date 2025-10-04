#' Calculate collection of log-probabilities for different components and samples
#' @return Matrix with as many rows as components and as many cols as samples
#' @keywords internal
util_lprobsave <- function(xVs, params, logW = 0, temporarydir, lab) {
    with(c(xVs, params), {
    out <- logW
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

    saveRDS(out,
        file.path(temporarydir,
            paste0(lab, ii, '__.rds'))
    )
    })
}

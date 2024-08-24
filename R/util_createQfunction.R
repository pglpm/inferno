## NB: the functional form of this function does not depend on
##   the number of components, minalpha, and maxalpha parameters

#' Calculate and save transformation function for ordinal variates
#' @keywords internal
createQfunction <- function(
    nint = 3600,
    nsamples = 2^24L,
    mean = 0,
    sd = 3,
    shapelo = 0.5,
    shapehi = 0.5,
    rate = 1,
    file = paste0('__Qfunction', nint, '_', sd),
    save = TRUE, # save using usethis::use_data()
    plot = F
) {
    ##
    seqnint <- (1:(nint - 1)) / nint
    means <- rnorm(nsamples, mean = mean, sd = sd)
    sds <- sqrt(
        nimble::rinvgamma(nsamples, shape = shapelo, rate = nimble::rinvgamma(nsamples, shape = shapehi, rate = rate))
        ## extraDistr::rbetapr(nsamples, shape1=shapehi, shape2=shapelo, scale=scale)
    )
    xsamples <- rnorm(nsamples,
        mean = means,
        sd = sds
    )
    ##
    thismad <- mad(xsamples, constant = 1)
    oquants <- c(
        NULL,
        quantile(x = xsamples, probs = seqnint, na.rm = TRUE,
            type = 6, names = FALSE),
        NULL
    )
    rm(xsamples)
    oquants <- (oquants - rev(oquants)) / 2
    util_Q <- approxfun(x = seqnint, y = oquants, yleft = -Inf, yright = +Inf)
    util_invQ <- approxfun(y = seqnint, x = oquants, yleft = 0, yright = 1)

    if (is.character(file)) {
        saveRDS(util_Q, paste0(file, '.rds'))
        saveRDS(util_invQ, paste0(file,'_inv', '.rds'))
    }
    ##
    xsamples <- util_Q(seqnint)
    oquants <- foreach(x = xsamples, .combine = c) %dopar% {
        mean(dnorm(x, mean = means, sd = sds))
    }
    oquants <- (oquants + rev(oquants)) / 2
    ##
    ##
    util_invDQ <- approxfun(x = xsamples, y = oquants, yleft = 0, yright = 0)
    if (is.character(file)) {
        saveRDS(util_invDQ, paste0('D', file, '.rds'))
    }
    if (is.character(plot)) {
        ##
        nint <- 256
        xgrid <- seq(1 / nint, (nint - 1) / nint, length.out = nint - 1)
        pdff(plot)
        tplot(
            x = xgrid, y = list(
                util_Q(xgrid), qnorm(xgrid, sd = thismad / qnorm(3 / 4)), qcauchy(xgrid, scale = thismad) # ,qlogis(xgrid,scale=1/qlogis(3/4))
            ),
            lwd = c(3, 2, 2, 5), lty = c(1, 2, 4, 3), alpha = c(0, rep(0.25, 3)),
            ylim = range(util_Q(xgrid)),
            ## xticks=c(0,0.25,0.5,0.75,1),xlabels=c(0,expression(italic(m)/4),expression(italic(m)/2),expression(3*italic(m)/4),expression(italic(m))),
            xlab = expression(italic(x)), ylab = expression(italic(Q)(italic(x))),
            mar = c(NA, 5, 1, 1)
        )
        legend('topleft', c(expression(italic(Q)), 'Gauss', 'Cauchy'), lwd = c(3, 2, 2, 5), lty = c(1, 2, 4, 3), col = c(1, 2, 3, 4), bty = 'n')
        dev.off()
    }

    if(save){
        usethis::use_data(util_Q, util_invQ, util_invDQ,
            internal = TRUE, overwrite = TRUE)
    }

    list(util_Q = util_Q, util_invQ = util_invQ, util_invDQ = util_invDQ)
}

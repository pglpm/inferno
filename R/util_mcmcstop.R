#### Function for calculating MCMC standard error
#' @keywords internal
funMCSE <- function(x) {
    x <- as.matrix(x)
    N <- nrow(x)
    b <- floor(sqrt(N))
    a <- floor(N/b)
    i0 <- N - a * b
    Ys <- rbind(sapply(seq_len(a), function(k) {
        colMeans(x[((k - 1) * b + 1 + i0):(k * b + i0), , drop = FALSE])
    }))
    ## ## Note: this is 'w" in doi.org/10.1080/10618600.2015.1044092
    sqrt(b * rowSums((Ys - rowMeans(Ys))^2) / ((a - 1) * N))
}


#### Function for calculating effective sample size
#### from https://github.com/LaplacesDemonR/LaplacesDemon/blob/master/R/ESS.R
#' @keywords internal
funESS <- function(x){
    x <- as.matrix(x)
    N <- nrow(x)
    M <- ncol(x)
    v0 <- order <- rep(0, M)
    names(v0) <- names(order) <- colnames(x)
    z <- 1:N
    for (i in 1:M) {
        lm.out <- lm(x[, i] ~ z)
        ## if(!isTRUE(all.equal(sd(residuals(lm.out)), 0))) {
            ar.out <- try(ar(x[,i], aic=TRUE), silent=TRUE)
            if(!inherits(ar.out, "try-error")) {
                v0[i] <- ar.out$var.pred / {1 - sum(ar.out$ar)}^2
                ## order[i] <- ar.out$order
            }
        ## }
    }
    ## spec <- list(spec=v0, order=order)
    ## spec <- spec$spec
    Y <- x - matrix(colMeans(x), N, M, byrow = TRUE)
    temp <- N * (N * colMeans(Y * Y) / (N - 1)) / v0
    v0[which(v0 != 0)] <- temp[which(v0 != 0)]
    v0[which(v0 < .Machine$double.eps)] <- .Machine$double.eps
    v0[which(v0 > N)] <- N
    v0
}

## #### Function for calculating number of needed MCMC iterations
## #' @keywords internal
## mcmcstop <- function(
##     traces,
##     nsamples,
##     availiter,
##     relerror,
##     ## diagnESS,
##     ## diagnIAT,
##     ## diagnBMK,
##     ## diagnMCSE,
##     ## diagnStat,
##     ## diagnBurn,
##     ## diagnBurn2,
##     ## diagnThin,
##     thinning
## ) {
##     ## Based on doi.org/10.1080/10618600.2015.1044092
## 
##     ## ## 'mcse' is 'w' or 'sigma/sqrt(n)' in doi.org/10.1080/10618600.2015.1044092
##     ## mcse <- funMCSE(traces)
##     ## N <- nrow(traces)
##     ## ## 'sds' is 'lambda' in doi.org/10.1080/10618600.2015.1044092
##     ## sds <- apply(traces, 2, sd)
##     ## avg <- apply(traces, 2, mean)
## 
##     relmcse <- funMCSE(traces) / apply(traces, 2, sd)
##     ## relmcse2 <- (mcse + 1/N) / sds
## 
##     ess <- funESS(traces)
## 
##     ## autothinning <- ceiling(1.5 * nrow(traces)/ess)
##     autothinning <- ceiling(nrow(traces)/ess)
## 
##     if(is.null(thinning)) {
##         thinning <- max(autothinning)
##     }
## 
##     missingsamples <- thinning * (nsamples - 1) - availiter
## 
##     if(max(relmcse) <= relerror) {
##         ## sampling could be stopped,
##         ## unless we still lack the required number of samples
##         reqiter <- max(0, missingsamples)
##     } else {
##         ## sampling should continue
##         reqiter <- max(ceiling(thinning * sqrt(nsamples)),
##             missingsamples)
##     }
## 
##     list(
##         reqiter = reqiter,
##         proposed.thinning = thinning,
##         toprint = list(
##             'rel. MC standard error' = relmcse,
##             'eff. sample size' = ess,
##             'needed thinning' = autothinning,
##             'average' = apply(traces, 2, mean)
##         )
##     )
## }

## #### Function for calculating number of needed MCMC iterations
## #' @keywords internal
## mcmcstopess <- function(
##     traces,
##     nsamples,
##     availiter,
##     reqess,
##     ## diagnESS,
##     ## diagnIAT,
##     ## diagnBMK,
##     ## diagnMCSE,
##     ## diagnStat,
##     ## diagnBurn,
##     ## diagnBurn2,
##     ## diagnThin,
##     thinning
## ) {
##     ## Based on doi.org/10.1080/10618600.2015.1044092
## 
##     ## ## 'mcse' is 'w' or 'sigma/sqrt(n)' in doi.org/10.1080/10618600.2015.1044092
##     ## mcse <- funMCSE(traces)
##     ## N <- nrow(traces)
##     ## ## 'sds' is 'lambda' in doi.org/10.1080/10618600.2015.1044092
##     ## sds <- apply(traces, 2, sd)
##     ## avg <- apply(traces, 2, mean)
## 
##     relmcse <- funMCSE(traces) / apply(traces, 2, sd)
##     ## relmcse2 <- (mcse + 1/N) / sds
## 
##     ess <- funESS(traces)
## 
##     ## autothinning <- ceiling(1.5 * nrow(traces)/ess)
##     autothinning <- ceiling(nrow(traces)/ess)
## 
##     if(is.null(thinning)) {
##         thinning <- max(autothinning)
##     }
## 
##     missingsamples <- thinning * (nsamples - 1) - availiter
##     reqsamples <- thinning * (reqess + 2) - availiter
## 
##     if(min(ess) >= reqess + 2) {
##         ## sampling could be stopped,
##         ## unless we still lack the required number of samples
##         reqiter <- max(0, missingsamples)
##     } else {
##         ## sampling should continue
##         reqiter <- max(reqsamples, missingsamples)
##     }
## 
##     list(
##         reqiter = reqiter,
##         proposed.thinning = thinning,
##         toprint = list(
##             'rel. MC standard error' = relmcse,
##             'eff. sample size' = ess,
##             'needed thinning' = autothinning,
##             'average' = apply(traces, 2, mean)
##         )
##     )
## }

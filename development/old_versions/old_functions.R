## flexiplot <- function(
##     x, y,
##     type = 'l',
##     lty = c(1, 2, 4, 3, 6, 5),
##     lwd = 2,
##     pch = c(1, 2, 0, 5, 6, 3), #, 4,
##     col = palette(),
##     xlab = NULL, ylab = NULL,
##     xlim = NULL, ylim = NULL,
##     add = FALSE,
##     xdomain = NULL, ydomain = NULL,
##     alpha.f = 1,
##     xjitter = NULL,
##     yjitter = NULL,
##     ## c( ## Tol's colour-blind-safe scheme
##     ##     '#4477AA',
##     ##     '#EE6677',
##     ##     '#228833',
##     ##     '#CCBB44',
##     ##     '#66CCEE',
##     ##     '#AA3377' #, '#BBBBBB'
##     ## ),
##     grid = TRUE,
##     cex.main = 1,
##     ...
## ){
##     xat <- yat <- NULL
## 
##     if(missing('x') && !missing('y')){
##         x <- numeric(NROW(y))
##         if(is.null(xdomain) && is.null(xlim)){
##             xat <- 0
##             xdomain <- NA
##             if(!is.null(xjitter)){
##                 xlim <- c(-0.04, 0.04)
##             }
##             if(is.null(xlab)){ xlab <- NA }
##             if(is.null(ylab)){ ylab <- deparse1(substitute(y)) }
##         }
##     } else if(!missing('x') && missing('y')){
##         y <- numeric(NROW(x))
##         if(is.null(ydomain) && is.null(ylim)){
##             yat <- 0
##             ydomain <- NA
##             if(!is.null(yjitter)){
##                 ylim <- c(-0.04, 0.04)
##             }
##             if(is.null(ylab)){ ylab <- NA }
##             if(is.null(xlab)){ xlab <- deparse1(substitute(x)) }
##         }
##     } else if(!missing('x') && !missing('y')){
##         if(is.null(xlab)){ xlab <- deparse1(substitute(x)) }
##         if(is.null(ylab)){ ylab <- deparse1(substitute(y)) }
##     } else {
##         stop('Arguments "x" and "y" cannot both be missing')
##     }
## 
##     ## if x is character, convert to numeric
##     if(is.character(x)){
##         if(is.null(xdomain)){ xdomain <- unique(x) }
##         ## we assume the user has sorted the values in a meaningful order
##         ## because the lexical order may not be correct
##         ## (think of values like 'low', 'medium', 'high')
##         x <- as.numeric(factor(x, levels = xdomain))
##         if(is.null(xjitter)){xjitter <- TRUE}
##         xat <- seq_along(xdomain)
##     }
##     if(isTRUE(xjitter)){x <- jitter(x)}
## 
##     ## if y is character, convert to numeric
##     if(is.character(y)){
##         if(is.null(ydomain)){ ydomain <- unique(y) }
##         ## we assume the user has sorted the values in a meaningful order
##         ## because the lexical order may not be correct
##         ## (think of values like 'low', 'medium', 'high')
##         y <- as.numeric(factor(y, levels = ydomain))
##         if(is.null(yjitter)){yjitter <- TRUE}
##         yat <- seq_along(ydomain)
##     }
##     if(isTRUE(yjitter)){y <- jitter(y)}
## 
##     ## Syntax of xlim and ylim that allows
##     ## for the specification of only upper- or lower-bound
##     if(length(xlim) == 2){
##         if(is.null(xlim[1]) || !is.finite(xlim[1])){ xlim[1] <- min(x[is.finite(x)]) }
##         if(is.null(xlim[2]) || !is.finite(xlim[2])){ xlim[2] <- max(x[is.finite(x)]) }
##     }
##     if(length(ylim) == 2){
##         if(is.null(ylim[1]) || !is.finite(ylim[1])){ ylim[1] <- min(y[is.finite(y)]) }
##         if(is.null(ylim[2]) || !is.finite(ylim[2])){ ylim[2] <- max(y[is.finite(y)]) }
##     }
## 
##     if(is.null(xlab) && !missing(x)) {
##         xlab <- deparse1(substitute(x))
##     }
##     if(is.na(alpha.f)){alpha.f <- 1}
##     col <- adjustcolor(col, alpha.f = alpha.f)
## 
##     graphics::matplot(x, y, xlim = xlim, ylim = ylim, type = type, axes = F,
##         col = col, lty = lty, lwd = lwd, pch = pch, cex.main = cex.main, add = add, xlab = xlab, ylab = ylab, ...)
##     if(!add){
##         graphics::axis(1, at = xat, labels = xdomain, tick = !grid,
##             col = 'black', lwd = 1, lty = 1, ...)
##         graphics::axis(2, at = yat, labels = ydomain, tick = !grid,
##             col = 'black', lwd = 1, lty = 1, ...)
##         if(grid){
##             graphics::grid(nx = NULL, ny = NULL, lty = 1, col = '#BBBBBB80')
##         }
##     }
## }





#### Possibly for future versions
## #' Summary for an object of class 'probability'
## #'
## #' Should this be 'print'?
## #'
## #' @export
## summary.probability <- function(x, ...){print.default(x, ...)}





#### The following functions are not used for the moment,
#### but may be useful in future versions.

## #' Calculate quantile width through batches
## #'
## #' Modified from
## #' from https://github.com/LaplacesDemonR/LaplacesDemon/blob/master/R/ESS.R
## #'
## #' @param x A matrix, rows being MC samples and columns being quantities whose MCSE is to be estimated.
## #'
## #' @return Estimates of the MC standard error for each trace. Division by sqrt(N) is already performed.
## #'
## #' @import stats
## #'
## #' @keywords internal
## funMCEI <- function(x, fn, p = c(0.055, 0.945), ...) {
##     N <- length(x)
##     a <- floor(sqrt(N))
##     b <- N %/% a
##     y <- x[rep(x = seq_len(a), each = b) +
##                round(seq(from = 0, to = N-a, length.out = b))]
##     dim(y) <- c(b, a)
##     quantile(x = apply(y, 2, FUN = fn, ...),
##         probs = p, na.rm = FALSE, names = FALSE, type = 6)
## }


## #' Calculate MC effective sample size using LaplacesDemon's algorithm
## #'
## #' Modified from
## #' from https://github.com/LaplacesDemonR/LaplacesDemon/blob/master/R/ESS.R
## #'
## #' @param x A matrix, rows being MC samples and columns being quantities whose MCSE is to be estimated.
## #'
## #' @return Estimates of the effective sample size for each trace.
## #'
## #' @import stats
## #'
## #' @keywords internal
## funESSLD <- function(x){
##     x <- as.matrix(x)
##     N <- nrow(x)
##     M <- ncol(x)
##     v0 <- order <- rep(0, M)
##     names(v0) <- names(order) <- colnames(x)
##     z <- 1:N
##     for (i in 1:M) {
##         lm.out <- lm(x[, i] ~ z)
##         ## if(!isTRUE(all.equal(sd(residuals(lm.out)), 0))) {
##             ar.out <- try(ar(x[,i], aic=TRUE), silent=TRUE)
##             if(!inherits(ar.out, "try-error")) {
##                 v0[i] <- ar.out$var.pred / {1 - sum(ar.out$ar)}^2
##                 ## order[i] <- ar.out$order
##             }
##         ## }
##     }
##     ## spec <- list(spec=v0, order=order)
##     ## spec <- spec$spec
##     Y <- x - matrix(colMeans(x), N, M, byrow = TRUE)
##     temp <- N * (N * colMeans(Y * Y) / (N - 1)) / v0
##     v0[which(v0 != 0)] <- temp[which(v0 != 0)]
##     v0[which(v0 < .Machine$double.eps)] <- .Machine$double.eps
##     v0[which(v0 > N)] <- N
##     v0
## }


## #' Calculate MC standard error, from Geyer's mcmc package
## #'
## #' @param x A matrix, rows being MC samples and columns being quantities whose MCSE is to be estimated.
## #'
## #' @return Estimates of the MC standard error for each trace. Division by sqrt(N) is already performed.
## #'
## #' @keywords internal
## funMCSEGeyer <- function(x){
##     x <- as.matrix(x)
##     N <- nrow(x)
##     apply(x, 2, function(atrace){
##         sqrt(mcmc::initseq(atrace)$var.con / N)
##     })
## }

## #' Function for calculating MC standard error, from Geyer's mcmc package
## #'
## #' @param x A matrix, rows being MC samples and columns being quantities whose ESS is to be estimated.
## #'
## #' @return Estimates of ESS for each trace.
## #'
## #' @keywords internal
## funESSGeyer <- function(x){
##     x <- as.matrix(x)
##     (apply(x, 2, sd) / funMCSE2(x))^2
## }


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








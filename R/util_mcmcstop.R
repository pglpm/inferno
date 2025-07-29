## #' Notes on the MCSE and ESS functions below
## #'
## #' `funMCSELD()` gives a good approximation of the "true" standard deviation in the case of independent samples. Multiply by `qnorm(x)` to obtain the `x`-quantile.
## #'
## #' `sd() / sqrt(funESS3()` gives essentially identical results to `funMCSELD()`, but it's 20 times slower.
## #'
## #' `funMCEQ()` gives a very good approximation of the "true" credibility quantiles in the case of independent samples.
## #'
## #' All above tested on t-distributions with df=1.1 and Pareto with a=1.5 (mean exists, variance infinite).


#' Calculate quantile width through batches
#'
#' Modified from
#' from https://github.com/LaplacesDemonR/LaplacesDemon/blob/master/R/ESS.R
#'
#' @param x: a matrix, rows being MC samples and columns being quantities whose MCSE is to be estimated.
#'
#' @return Estimates of the MC standard error for each trace. Division by sqrt(N) is already performed.
#'
#' @keywords internal
funMCEI <- function(x, fn, p = c(0.055, 0.945), ...) {
    N <- length(x)
    a <- floor(sqrt(N))
    b <- N %/% a
    y <- x[rep(x = seq_len(a), each = b) +
               round(seq(from = 0, to = N-a, length.out = b))]
    dim(y) <- c(b, a)
    quantile(apply(y, 2, FUN = fn, ...),
        probs = p, na.rm = FALSE, names = FALSE, type = 6)
}

#' Calculate credibility quantiles on estimated quantile
#'
#' Modified from Vehtari et al.
#'
#' @param x: a vector of MC samples
#' @param prob: quantile whose error intervalis being estimated
#' @param Qpair: lower and higher credibility-quantiles requested
#' 
#' @return Estimates lower and higher credibility-quantiles on estimated quantile
#'
#' @keywords internal
funMCEQ <- function(x, prob = c(0.055, 0.945), Qpair = pnorm(c(-1, 1))){
    N <- length(x)
    straces <- sort(x)
    sapply(prob, function(aprob) {
        Xlo <- quantile(x, aprob, na.rm = FALSE, names = FALSE, type = 6)
        ##
        essXlo <- funESS3(x <= Xlo)
        ##
        a <- qbeta(Qpair, essXlo * aprob + 1, essXlo * (1 - aprob) + 1)
        c(straces[max(round(a[1] * N), 1)], straces[min(round(a[2] * N), N)])
    })
}
## funMCEQ0 <- function(x, prob, Qpair){
##     N <- length(x)
##     straces <- sort(x)
##         Xlo <- quantile(x, prob, na.rm = FALSE, names = FALSE, type = 6)
##         ##
##         essXlo <- funESS3(x <= Xlo)
##         ##
##         a <- qbeta(Qpair, essXlo * prob + 1, essXlo * (1 - prob) + 1)
##         c(straces[max(round(a[1] * N), 1)], straces[min(round(a[2] * N), N)])
## }


#' Calculate MC standard error using LaplacesDemon's batch means
#'
#' Modified from
#' from https://github.com/LaplacesDemonR/LaplacesDemon/blob/master/R/ESS.R
#'
#' @param x: a matrix, rows being MC samples and columns being quantities whose MCSE is to be estimated.
#'
#' @return Estimates of the MC standard error for each trace. Division by sqrt(N) is already performed.
#'
#' @keywords internal
funMCSELD <- function(x) {
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


#' Calculate MC effective sample size using LaplacesDemon's algorithm
#'
#' Modified from
#' from https://github.com/LaplacesDemonR/LaplacesDemon/blob/master/R/ESS.R
#'
#' @param x: a matrix, rows being MC samples and columns being quantities whose MCSE is to be estimated.
#'
#' @return Estimates of the effective sample size for each trace.
#'
#' @keywords internal
funESSLD <- function(x){
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

## #' Calculate MC standard error, from Geyer's mcmc package
## #'
## #' @param x: a matrix, rows being MC samples and columns being quantities whose MCSE is to be estimated.
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
## #' @param x: a matrix, rows being MC samples and columns being quantities whose ESS is to be estimated.
## #'
## #' @return Estimates of ESS for each trace.
## #'
## #' @keywords internal
## funESSGeyer <- function(x){
##     x <- as.matrix(x)
##     (apply(x, 2, sd) / funMCSE2(x))^2
## }


#' Find optimal FFT size
#'
#' From **rstan** <https://github.com/stan-dev/rstan/blob/develop/rstan/rstan/R/monitor.R>
#'
#' @param N: integer.
#'
#' @return Optimal FFT size
#'
#' @keywords internal
fftNGS <- function(N) {
  # Find the optimal next size for the FFT so that
  # a minimum number of zeros are padded.
  if (N <= 2)
    return(2)
  while (TRUE) {
    m <- N
    while ((m %% 2) == 0) m <- m / 2
    while ((m %% 3) == 0) m <- m / 3
    while ((m %% 5) == 0) m <- m / 5
    if (m <= 1)
      return(N)
    N <- N + 1
  }
}

#' Compute autocovariance
#'
#' From **rstan** <https://github.com/stan-dev/rstan/blob/develop/rstan/rstan/R/monitor.R>
#'
#' @param y: time series
#'
#' @return Autocovariances at different lags
#'
#' @keywords internal
funAC <- function(y) {
  N <- length(y)
  Mt2 <- 2 * fftNGS(N)
  yc <- y - mean(y)
  yc <- c(yc, rep.int(0, Mt2 - N))
  transform <- fft(yc)
  ac <- fft(Conj(transform) * transform, inverse = TRUE)
  # use "biased" estimate as recommended by Geyer (1992)
  ac <- Re(ac)[1:N] / (N^2 * 2)
  ac
}


#' Compute ESS
#'
#' From **rstan** <https://github.com/stan-dev/rstan/blob/develop/rstan/rstan/R/monitor.R>
#'
#' @param x: vector of Monte Carlo samples
#'
#' @return Effective Sample Size
#'
#' @keywords internal
funESS3 <- function(x){
  N <- length(x)
  if (N < 3L || any(!is.finite(x))) {
    return(NA_real_)
  }
  acov <- funAC(x)
  chain_mean <- mean(x)
  var_plus <- acov[1]
  mean_var <- var_plus * N / (N - 1)
  ##
  ## Geyer's initial positive sequence
  rho_hat_t <- rep.int(0, N)
  t <- 0
  rho_hat_even <- 1
  rho_hat_t[t + 1] <- rho_hat_even
  rho_hat_odd <- 1 - (mean_var - acov[t + 2]) / var_plus
  rho_hat_t[t + 2] <- rho_hat_odd
  while (t < length(acov) - 5 &&
             !is.nan(rho_hat_even + rho_hat_odd) &&
             (rho_hat_even + rho_hat_odd > 0)) {
                 t <- t + 2
                 rho_hat_even = 1 - (mean_var - acov[t + 1]) / var_plus
                 rho_hat_odd = 1 - (mean_var - acov[t + 2]) / var_plus
                 if ((rho_hat_even + rho_hat_odd) >= 0) {
                     rho_hat_t[t + 1] <- rho_hat_even
                     rho_hat_t[t + 2] <- rho_hat_odd
                 }
             }
  max_t <- t
  ## this is used in the improved estimate
  if (rho_hat_even>0){rho_hat_t[max_t + 1] <- rho_hat_even}
  ##
  ## Geyer's initial monotone sequence
  t <- 0
  while (t <= max_t - 4) {
    t <- t + 2
    if (rho_hat_t[t + 1] + rho_hat_t[t + 2] > rho_hat_t[t - 1] + rho_hat_t[t]) {
      rho_hat_t[t + 1] = (rho_hat_t[t - 1] + rho_hat_t[t]) / 2
      rho_hat_t[t + 2] = rho_hat_t[t + 1]
    }
  }
  ess <- N
  # Geyer's truncated estimate
  # tau_hat <- -1 + 2 * sum(rho_hat_t[1:max_t])
  # Improved estimate reduces variance in antithetic case
  tau_hat <- -1 + 2 * sum(rho_hat_t[1:max_t]) + rho_hat_t[max_t+1]
  # Safety check for negative values and with max ess equal to ess*log10(ess)
  tau_hat <- max(tau_hat, 1/log10(ess))
  ess <- ess / tau_hat
  ess
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

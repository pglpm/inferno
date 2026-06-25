#' Calculate credibility quantiles on estimated quantile
#'
#' Calculates the lower and upper bound of a credibility interval, for various quantiles of the empirical distribution of a vector of MC samples.
#'
#' Tests show that it gives a very good approximation of the "true" credibility quantiles in the case of independent samples.
#'
#' Tested also on t-distributions with df=1.1 and Pareto with a=1.5 (mean exists, variance infinite).
#'
#' Used in 'workerfun()' in 'learn()'
#'
#' @param x A vector of MC samples
#' @param prob numeric vector of probabilities: quantiles whose error interval is being estimated.
#' @param Qpair vector of length two (further elements are ignored): lower and higher credibility-quantiles requested. Default yields a credibility interval of 68%, or one nominal normal standard deviation.
#'
#' @return A matrix with two rows and as many columns as elements in 'prob'. Forr each column, the first and second row determine the lower and upper bound of the credibility interval of width `Qpair[2] - Qpair[2]`.
#'
#' @import stats
#'
#' @keywords internal
funMCEQ <- function(x, prob = c(0.055, 0.945), Qpair = pnorm(c(-1, 1))){
    N <- length(x)
    straces <- sort(x)
    sapply(prob, function(aprob) {
        Xlo <- quantile(x = x, probs = aprob,
            na.rm = FALSE, names = FALSE, type = 6)
        ##
        essXlo <- funESS3(x <= Xlo)
        ##
        a <- qbeta(Qpair, essXlo * aprob + 1, essXlo * (1 - aprob) + 1)
        c(straces[max(round(a[1] * N), 1)], straces[min(round(a[2] * N), N)])
    })
}

#' Compute ESS
#'
#' Modified from 'rstan' <https://github.com/stan-dev/rstan/blob/develop/rstan/rstan/R/monitor.R>
#'
#' Used in 'workerfun()' in 'learn()', and in 'funMCEQ()'.
#'
#' @param x Vector of MC samples.
#'
#' @return Effective Sample Size.
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


#' Calculate MC standard error using LaplacesDemon's batch means
#'
#' This function gives a good approximation of the "true" standard deviation in the case of independent samples. Multiply by `qnorm(x)` to obtain the `x`-quantile.
#'
#' Modified from <https://github.com/LaplacesDemonR/LaplacesDemon/blob/master/R/ESS.R>.
#'
#' Tested also on t-distributions with df=1.1 and Pareto with a=1.5 (mean exists, variance infinite).
#'
#' `sd() / sqrt(funESS3()` gives essentially identical results to `funMCSELD()`, but it's 20 times slower.
#'
#' Used in 'util_combineYX()' in 'Pr()'.
#'
#' @param x matrix, each row being a "trace", that is a set of MC samples, whose MCSE is to be estimated.
#'
#' @return MCSE estimates, one for each trace. Division by sqrt(N) is already performed.
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


#' Compute autocovariance
#'
#' Modified from rstan <https://github.com/stan-dev/rstan/blob/develop/rstan/rstan/R/monitor.R>.
#'
#' Used in 'funESS3()'.
#'
#' @param y Time series
#'
#' @return Autocovariances at different lags
#'
#' @import stats
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


#' Find optimal FFT size
#'
#' Modified from rstan <https://github.com/stan-dev/rstan/blob/develop/rstan/rstan/R/monitor.R>.
#'
#' Used in 'funAC()'.
#'
#' @param N Integer.
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

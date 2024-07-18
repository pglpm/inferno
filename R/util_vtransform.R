## Transformation from variate to internal variable
vtransform <- function(x, auxmetadata,
                       Rout = '',
                       Cout = 'init',
                       Dout = 'data',
                       Lout = 'data',
                       Bout = 'numeric',
                       Oout = 'numeric',
                       Nout = 'numeric',
                       variates = NULL,
                       invjacobian = FALSE,
                       useLquantiles = FALSE) {
  ## Qf <- readRDS('Qfunction3600_3.rds')
  ## DQf <- readRDS('DQfunction3600_3.rds')
  ## invQf <- readRDS('invQfunction3600_3.rds')
  x <- as.data.frame(cbind(x))
  if (!is.null(variates)) {
    colnames(x) <- variates
  }
  matrix(sapply(colnames(x), function(v) {
    ##
    ## datum <- unlist(x[, v, with = F]) # old version, remove
    datum <- x[[v]]
    info <- as.list(auxmetadata[auxmetadata$name == v, ])
    ##
    if (invjacobian) {
#### Calculation of reciprocal Jacobian factors
      if (info$mcmctype %in% c('B', 'N', 'O', 'L')) {
        datum <- rep(1L, length(datum))
      } else {
        if (info$transform == 'log') {
          datum <- (datum - info$domainmin) * info$tscale
        } else if (info$transform == 'logminus') {
          datum <- (info$domainmax - datum) * info$tscale
        } else if (info$transform == 'Q') {
          datum <- Qf((datum - info$domainmin) / (info$domainmax - info$domainmin))
          datum <- DQf(datum) * info$tscale * (info$domainmax - info$domainmin)
        } else {
          datum <- rep(info$tscale, length(datum))
        }
        xv <- data.matrix(x[, v, drop = FALSE])
        if (info$mcmctype %in% c('C', 'D')) {
          datum[(xv >= info$censormax) | (xv <= info$censormin)] <- 1L
        }
        datum[is.na(xv)] <- 1L
      }
    } else {
#### Transformation to internal value for MCMC

#### Continuous, open domain

      if (info$mcmctype == 'R') {
        if(Rout == 'id') { # used in mutualinfo()
          datum
        } else { # transformation for internal MCMC use
          if (info$transform == 'log') {
            datum <- log(datum - info$domainmin)
          } else if (info$transform == 'logminus') {
            datum <- log(info$domainmax - datum)
          } else if (info$transform == 'Q') {
            datum <- Qf((datum - info$domainmin) / (info$domainmax - info$domainmin))
          }
          datum <- (datum - info$tlocation) / info$tscale
        }

#### Continuous, closed domain
      } else if (info$mcmctype == 'C') {
        if (info$transform == 'log') {
          datum <- log(datum - info$domainmin)
          censormin <- log(info$censormin - info$domainmin)
          censormax <- log(info$censormin - info$domainmin)
        } else if (info$transform == 'logminus') {
          datum <- log(info$domainmax - datum)
          censormin <- log(info$domainmax - censormin)
          censormax <- log(info$domainmax - censormax)
        } else if (info$transform == 'Q') {
          datum <- Qf(0.5 +
                      (datum - (info$domainmin + info$domainmax)/2) /
                      (info$domainmax - info$domainmin))
          censormin <- Qf(0.5 +
                      (info$censormin - (info$domainmin + info$domainmax)/2) /
                      (info$domainmax - info$domainmin))
          censormax <- Qf(0.5 +
                      (info$censormax - (info$domainmin + info$domainmax)/2) /
                      (info$domainmax - info$domainmin))
        } else {
          censormin <- info$censormin
          censormax <- info$censormax

        }
        ## datum <- (datum-info$tlocation)/info$tscale
        ## censormax <- (censormax-info$tlocation)/info$tscale
        ## censormin <- (censormin-info$tlocation)/info$tscale
        if (Cout == 'left') { # in MCMC
          xv <- data.matrix(x[, v, drop = FALSE])
          sel <- is.na(xv) | (xv < info$censormax)
          datum[sel] <- -Inf
          datum[!sel] <- censormax
        } else if (Cout == 'right') { # in MCMC
          xv <- data.matrix(x[, v, drop = FALSE])
          sel <- is.na(xv) | (xv > info$censormin)
          datum[sel] <- +Inf
          datum[!sel] <- censormin
        } else if (Cout == 'lat') { # latent variable in MCMC
          xv <- data.matrix(x[, v, drop = FALSE])
          sel <- is.na(xv) | (xv >= info$censormax) | (xv <= info$censormin)
          datum[sel] <- NA
        } else if (Cout == 'init') { # init in MCMC
          xv <- data.matrix(x[, v, drop = FALSE])
          datum[is.na(xv)] <- 0L
          datum[!is.na(xv) & (xv <= info$censormin)] <- censormin - 0.125 * info$tscale
          datum[!is.na(xv) & (xv >= info$censormax)] <- censormax + 0.125 * info$tscale
          datum[!is.na(xv) & (xv < info$censormax) & (xv > info$censormin)] <- NA
        } else if (Cout == 'aux') { # aux variable in MCMC
          xv <- data.matrix(x[, v, drop = FALSE])
          sel <- is.na(xv)
          datum[sel] <- NA
          datum[!sel] <- 1L
        } else if (Cout == 'boundisinf') { # in sampling functions
          xv <- data.matrix(x[, v, drop = FALSE])
          datum[xv >= info$censormax] <- +Inf
          datum[xv <= info$censormin] <- -Inf
        } else if (Cout == 'sleft') { # in sampling functions
          datum <- rep(censormin, length(datum))
        } else if (Cout == 'sright') { # in sampling functions
          datum <- rep(censormax, length(datum))
        } else if (Cout == 'idboundinf') { # in mutualinfo
          datum <- xv <- data.matrix(x[, v, drop = FALSE])
          datum[xv >= info$censormax] <- +Inf
          datum[xv <= info$censormin] <- -Inf
        }
        ##
        if (Cout != 'aux' && Cout != 'idboundinf') {
          datum <- (datum - info$tlocation) / info$tscale
        }

#### Continuous rounded
      } else if (info$mcmctype == 'D') {
        censormin <- info$censormin
        censormax <- info$censormax
        leftbound <- pmax(datum - info$step, info$domainmin, na.rm = T)
        leftbound[leftbound <= censormin] <- info$domainmin
        leftbound[leftbound >= censormax] <- censormax
        rightbound <- pmin(datum + info$step, info$domainmax, na.rm = T)
        rightbound[rightbound >= censormax] <- info$domainmax
        rightbound[rightbound <= censormin] <- censormin
        if (info$transform == 'log') {
          datum <- log(datum - info$domainmin)
          leftbound <- log(leftbound - info$domainmin)
          rightbound <- log(rightbound - info$domainmin)
          censormin <- log(censormin - info$domainmin)
          censormax <- log(censormax - info$domainmin)
        } else if (info$transform == 'logminus') {
          datum <- log(info$domainmax - datum)
          leftbound <- log(info$domainmax - leftbound)
          rightbound <- log(info$domainmax - rightbound)
          censormin <- log(info$domainmax - censormin)
          censormax <- log(info$domainmax - censormax)
        } else if (info$transform == 'Q') {
          datum <- Qf((datum - info$domainmin) / (info$domainmax - info$domainmin))
          leftbound <- Qf((leftbound - info$domainmin) / (info$domainmax - info$domainmin))
          rightbound <- Qf((rightbound - info$domainmin) / (info$domainmax - info$domainmin))
          censormin <- Qf((censormin - info$domainmin) / (info$domainmax - info$domainmin))
          censormax <- Qf((censormax - info$domainmin) / (info$domainmax - info$domainmin))
        }
        ## datum <- (datum-info$tlocation)/info$tscale
        ## rightbound <- (rightbound-info$tlocation)/info$tscale
        ## leftbound <- (leftbound-info$tlocation)/info$tscale
        ## censormax <- (censormax-info$tlocation)/info$tscale
        ## censormin <- (censormin-info$tlocation)/info$tscale
        if (Dout == 'left') {
          datum <- leftbound
        } else if (Dout == 'right') {
          datum <- rightbound
        } else if (Dout == 'init') { # init in MCMC
          xv <- data.matrix(x[, v, drop = FALSE])
          datum[is.na(xv)] <- 0L
          datum[!is.na(xv) & (xv <= info$censormin)] <- censormin - 0.125 * info$tscale
          datum[!is.na(xv) & (xv >= info$censormax)] <- censormax + 0.125 * info$tscale
        } else if (Dout == 'aux') { # aux variable in MCMC
          xv <- data.matrix(x[, v, drop = FALSE])
          sel <- is.na(xv)
          datum[sel] <- NA
          datum[!sel] <- 1L
        } else if (Dout == 'boundisinf') { # in sampling functions
          xv <- data.matrix(x[, v, drop = FALSE])
          datum[xv >= info$censormax] <- +Inf
          datum[xv <= info$censormin] <- -Inf
        } else if (Dout == 'sleft') { # in sampling functions
          datum <- rep(censormin, length(datum))
        } else if (Dout == 'sright') { # in sampling functions
          datum <- rep(censormax, length(datum))
        } else if (Dout == 'idboundinf') { # in mutualinfo
          datum <- xv <- data.matrix(x[, v, drop = FALSE])
          datum[xv >= info$censormax] <- +Inf
          datum[xv <= info$censormin] <- -Inf
        }
        ##
        if (Dout != 'aux' && Dout != 'idboundinf') {
          datum <- (datum - info$tlocation) / info$tscale
        }

#### Latent
      } else if (info$mcmctype == 'L') {
        olocation <- (info$Nvalues * info$domainmin - info$domainmax) / (info$Nvalues - 1)
        oscale <- (info$domainmax - info$domainmin) / (info$Nvalues - 1)
        ##
        ## output is in range 1 to Nvalues
        datum <- round((datum - olocation) / oscale)
        if (Lout == 'init') { # in sampling functions or init MCMC
          datum[is.na(datum)] <- info$Nvalues / 2 + 0.5
          datum <- Qf((datum - 0.5) / info$Nvalues)
          datum[datum == +Inf] <- 1e6
          datum[datum == -Inf] <- -1e6
          if (useLquantiles) {
            datum <- (datum - info$tlocation) / info$tscale
          }
        } else if (Lout == 'left') { # as left for MCMC
          datum <- Qf(pmax(0, datum - 1L) / info$Nvalues)
          if (useLquantiles) {
            datum <- (datum - info$tlocation) / info$tscale
          }
          datum[is.na(datum)] <- -Inf
        } else if (Lout == 'right') { # as right for MCMC
          datum <- Qf(pmin(info$Nvalues, datum) / info$Nvalues)
          if (useLquantiles) {
            datum <- (datum - info$tlocation) / info$tscale
          }
          datum[is.na(datum)] <- +Inf
        } else if (Lout == 'aux') { # aux variable in MCMC
          sel <- is.na(datum)
          datum[sel] <- NA
          datum[!sel] <- 1L
        } else if (Lout == 'boundisinf') { # in output functions
          datum <- datum - 1L
        } else if (Lout == 'integer') { # in mutualinfo
          datum <- data.matrix(x[, v, drop = FALSE])
          if (useLquantiles) {
            datum <- datum * info$tscale + info$tlocation
          }
          datum <- round(invQf(datum) * info$Nvalues + 0.5)
          datum[datum < 1] <- 1L
          datum[datum > info$Nvalues] <- info$Nvalues
        }

      ## Ordinal
      } else if (info$mcmctype == 'O') {
        bvalues <- 1:info$Nvalues
        names(bvalues) <- unlist(info[paste0('V', bvalues)])
        if (Oout == 'numeric') {
          datum <- bvalues[as.character(datum)]
        } else if (Oout == 'original') {
          datum <- names(bvalues[datum])
        }

      ## Nominal
      } else if (info$mcmctype == 'N') {
        bvalues <- 1:info$Nvalues
        names(bvalues) <- unlist(info[paste0('V', bvalues)])
        if (Nout == 'numeric') {
          datum <- bvalues[as.character(datum)]
        } else if (Nout == 'original') {
          datum <- names(bvalues[datum])
        }

      ## Binary
      } else if (info$mcmctype == 'B') {
        bvalues <- 0:1
        names(bvalues) <- unlist(info[c('V1', 'V2')])
        if (Bout == 'numeric') {
          datum <- bvalues[as.character(datum)]
        } else if (Bout == 'original') {
          datum <- names(bvalues[datum + 1])
        }
      }
    }
    ##
    datum
  }), ncol = ncol(x), dimnames = list(NULL, colnames(x)))
}

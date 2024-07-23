#' Transforms variates to different representations
#'
#' @param x data.table object containing data to be transformed
#' @param auxmetadata auxmetadata object
#' @param Rout string, output of R-type variate, with possible values:
#'  'normalized': for internal MCMC use
#'  'mi': for use in mutualinfo()
#'  'original': original representation
#' @param Cout string, output of C-type variate, with possible values:
#'  'init': for internal MCMC use (init input)
#'  'left', 'right': for internal MCMC use
#'  'aux', 'lat': for internal MCMC use
#'  'boundisinf': for sampling functions
#'  'mi': for use in mutualinfo()
#'  'original': original representation
#' @param Dout string, output of D-type variate, with possible values:
#'  'init': for internal MCMC use (init input)
#'  'left', 'right': for internal MCMC use
#'  'aux': for internal MCMC use
#'  'boundisinf': for sampling functions
#'  'mi': for use in mutualinfo()
#'  'original': original representation
#' @param Lout string, output of L-type variate, with possible values:
#'  'init': for internal MCMC use (init input)
#'  'left', 'right': for internal MCMC use
#'  'aux': for internal MCMC use
#'  'boundisinf': for sampling functions
#'  'mi': for use in mutualinfo()
#'  'original': original representation
#' @param Oout string, output of O-type variate, with possible values:
#'  'numeric': for internal MCMC use, values 1,2,...
#'  'original': original representation
#' @param Nout string, output of N-type variate, with possible values:
#'  'numeric': for internal MCMC use, values 1,2,...
#'  'original': original representation
#' @param Bout string, output of B-type variate, with possible values:
#'  'numeric': for internal MCMC use, values 0,1
#'  'original': original representation
#' @param variates string vector, names of variates
#'   corresponding to columns of x (in case x misses column names)
#' @param invjacobian logical: calculate Jacobian factor?
#' @param useLquantiles logical: use L quantiles in transformation?
vtransform <- function(x, auxmetadata,
                       Rout = NULL,
                       Cout = NULL,
                       Dout = NULL,
                       Lout = NULL,
                       Bout = NULL,
                       Oout = NULL,
                       Nout = NULL,
                       variates = NULL,
                       invjacobian = FALSE,
                       useLquantiles = FALSE) {
  ## Qf <- readRDS('Qfunction3600_3.rds')
  ## DQf <- readRDS('DQfunction3600_3.rds')
  ## invQf <- readRDS('invQfunction3600_3.rds')
  x <- as.data.frame(cbind(x))
  if (!missing(variates)) {
    colnames(x) <- variates
  }
  matrix(sapply(colnames(x), function(v) {
    datum <- x[[v]]
    info <- as.list(auxmetadata[auxmetadata$name == v, ])

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
          datum <- Qf(0.5 +
                      (datum - (info$domainmin + info$domainmax)/2) /
                      (info$domainmax - info$domainmin)
                      )
          datum <- DQf(datum) * info$tscale * (info$domainmax - info$domainmin)
        } else if (info$transform == 'identity') {
          datum[] <- info$tscale
        } else {
          stop('Unknown transformation for variate', v)
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

        if (Rout == 'normalized') {
          ## used for MCMC
          if (info$transform == 'log') {
            datum <- log(datum - info$domainmin)
          } else if (info$transform == 'logminus') {
            datum <- log(info$domainmax - datum)
          } else if (info$transform == 'Q') {
            datum <- Qf(0.5 +
                        (datum - (info$domainmin + info$domainmax)/2) /
                        (info$domainmax - info$domainmin)
                        )
          }
          datum <- (datum - info$tlocation) / info$tscale

        } else if (Rout == 'original') {
          ## transformation from MCMC representation
          ## to original domain
          datum <- datum * info$tscale + info$tlocation
          if (info$transform == 'log') {
            datum <- exp(datum) + info$domainmin
          } else if (info$transform == 'logminus') {
            datum <- info$domainmax - exp(datum)
          } else if (info$transform == 'Q') {
            datum <- (invQf(datum) - 0.5) * (info$domainmax - info$domainmin) +
              (info$domainmin + info$domainmax)/2
          }

        } else if (Rout != 'mi'){
          ## used in mutualinfo()
          stop('Unknown transformation for variate', v)
        }

#### Continuous, closed domain
      } else if (info$mcmctype == 'C') {

        if (Cout == 'init') {
          ## init value used for MCMC
          selna <- is.na(datum)
          selmin <- !is.na(datum) & (datum <= info$censormin)
          selmax <- !is.na(datum) & (datum >= info$censormax)
          selmid <- !is.na(datum) &
            (datum < info$censormax) & (datum > info$censormin)
          datum[selna] <- 0L
          datum[selmin] <- info$tcensormin - 0.125
          datum[selmax] <- info$tcensormax + 0.125
          datum[selmid] <- NA # so that it's discarded by Nimble

        } else if (Cout == 'left') {
          ## used for MCMC
          sel <- is.na(datum) | (datum < info$censormax)
          datum[sel] <- -Inf
          datum[!sel] <- info$tcensormax

        } else if (Cout == 'right') {
          ## used for MCMC
          sel <- is.na(datum) | (datum > info$censormin)
          datum[sel] <- +Inf
          datum[!sel] <- info$tcensormin

        } else if (Cout == 'lat') {
          ## latent variable used for MCMC
          datum[(datum >= info$censormax) | (datum <= info$censormin)] <- NA

          if (info$transform == 'log') {
            datum <- log(datum - info$domainmin)
          } else if (info$transform == 'logminus') {
            datum <- log(info$domainmax - datum)
          } else if (info$transform == 'Q') {
            datum <- Qf(0.5 +
                        (datum -
                         (info$domainmin + info$domainmax)/2) /
                        (info$domainmax - info$domainmin)
                        )
          }
          ##
          datum <- (datum - info$tlocation) / info$tscale

        } else if (Cout == 'aux') {
          ## aux variable used for MCMC
          sel <- is.na(datum)
          datum[sel] <- NA
          datum[!sel] <- 1L

        } else if (Cout == 'boundisinf') {
          ## used in sampling functions
          selmax <- (datum >= info$censormax)
          selmin <- (datum <= info$censormin)
          ##
          if (info$transform == 'log') {
            datum <- log(datum - info$domainmin)
          } else if (info$transform == 'logminus') {
            datum <- log(info$domainmax - datum)
          } else if (info$transform == 'Q') {
            datum <- Qf(0.5 +
                        (datum -
                         (info$domainmin + info$domainmax)/2) /
                        (info$domainmax - info$domainmin)
                        )
          }
          ##
          datum[selmax] <- +Inf
          datum[selmin] <- -Inf
          datum <- (datum - info$tlocation) / info$tscale

        ## } else if (Cout == 'sleft') {
        ##   ## used in sampling functions
        ##   datum[] <- info$tcensormin
        ## 
        ## } else if (Cout == 'sright') {
        ##   ## used in sampling functions
        ##   datum[] <- info$tcensormax

        } else if (Cout == 'mi') {
          ## used in mutualinfo
          datum[datum >= info$tcensormax] <- +Inf
          datum[datum <= info$tcensormin] <- -Inf

        } else if (Cout == 'original') {
          ## transformation from MCMC representation
          ## to original domain
          datum <- datum * info$tscale + info$tlocation
          if (info$transform == 'log') {
            datum <- exp(datum) + info$domainmin
          } else if (info$transform == 'logminus') {
            datum <- info$domainmax - exp(datum)
          } else if (info$transform == 'Q') {
            datum <- (invQf(datum) - 0.5) * (info$domainmax - info$domainmin) +
              (info$domainmin + info$domainmax)/2
          }
          datum[datum <= info$censormin] <- censormin
          datum[datum >= info$censormax] <- censormax

        } else {
          stop('Unknown transformation for variate', v)
        }



#### Continuous rounded
      } else if (info$mcmctype == 'D') {
        ## tcensormin <- info$censormin
        ## tcensormax <- info$censormax
        ## leftbound <- pmax(datum - info$step, info$domainmin, na.rm = T)
        ## leftbound[leftbound <= info$censormin] <- info$domainmin
        ## leftbound[leftbound >= info$censormax] <- info$censormax
        ## rightbound <- pmin(datum + info$step, info$domainmax, na.rm = T)
        ## rightbound[rightbound >= info$censormax] <- info$domainmax
        ## rightbound[rightbound <= info$censormin] <- info$censormin
        ## if (info$transform == 'log') {
        ##   datum <- log(datum - info$domainmin)
        ##   leftbound <- log(leftbound - info$domainmin)
        ##   rightbound <- log(rightbound - info$domainmin)
        ##   tcensormin <- log(tcensormin - info$domainmin)
        ##   tcensormax <- log(tcensormax - info$domainmin)
        ## } else if (info$transform == 'logminus') {
        ##   datum <- log(info$domainmax - datum)
        ##   leftbound <- log(info$domainmax - leftbound)
        ##   rightbound <- log(info$domainmax - rightbound)
        ##   tcensormin <- log(info$domainmax - tcensormin)
        ##   tcensormax <- log(info$domainmax - tcensormax)
        ## } else if (info$transform == 'Q') {
        ##   datum <- Qf((datum - info$domainmin) / (info$domainmax - info$domainmin))
        ##   leftbound <- Qf((leftbound - info$domainmin) / (info$domainmax - info$domainmin))
        ##   rightbound <- Qf((rightbound - info$domainmin) / (info$domainmax - info$domainmin))
        ##   tcensormin <- Qf((tcensormin - info$domainmin) / (info$domainmax - info$domainmin))
        ##   tcensormax <- Qf((tcensormax - info$domainmin) / (info$domainmax - info$domainmin))
        ## }
        ## datum <- (datum-info$tlocation)/info$tscale
        ## rightbound <- (rightbound-info$tlocation)/info$tscale
        ## leftbound <- (leftbound-info$tlocation)/info$tscale
        ## tcensormax <- (tcensormax-info$tlocation)/info$tscale
        ## tcensormin <- (tcensormin-info$tlocation)/info$tscale

        if (Dout == 'init') {
          ## init value used for MCMC
          selna <- is.na(datum)
          selmin <- !is.na(datum) & (datum <= info$censormin)
          selmax <- !is.na(datum) & (datum >= info$censormax)
          datum[selna] <- 0L
          datum[selmin] <- info$tcensormin - 0.125
          datum[selmax] <- info$tcensormax + 0.125

        } else if (Dout == 'left') {
          ## used for MCMC
          datum <- pmax(datum - info$step, info$domainmin, na.rm = T)
          datum[datum <= info$censormin] <- info$domainmin
          datum[datum >= info$censormax] <- info$censormax
          ##
          if (info$transform == 'log') {
            datum <- log(datum - info$domainmin)
          } else if (info$transform == 'logminus') {
            datum <- log(info$domainmax - datum)
          } else if (info$transform == 'Q') {
            datum <- Qf(0.5 +
                        (datum -
                         (info$domainmin + info$domainmax)/2) /
                        (info$domainmax - info$domainmin)
                        )
          }
          datum <- (datum - info$tlocation) / info$tscale

        } else if (Dout == 'right') {
          ## used for MCMC
          datum <- pmin(datum + info$step, info$domainmax, na.rm = T)
          datum[datum >= info$censormax] <- info$domainmax
          datum[datum <= info$censormin] <- info$censormin
          ##
          if (info$transform == 'log') {
            datum <- log(datum - info$domainmin)
          } else if (info$transform == 'logminus') {
            datum <- log(info$domainmax - datum)
          } else if (info$transform == 'Q') {
            datum <- Qf(0.5 +
                        (datum -
                         (info$domainmin + info$domainmax)/2) /
                        (info$domainmax - info$domainmin)
                        )
          }
          datum <- (datum - info$tlocation) / info$tscale

        } else if (Dout == 'aux') {
          ## aux variable used for MCMC
          sel <- is.na(datum)
          datum[sel] <- NA
          datum[!sel] <- 1L

        } else if (Dout == 'boundisinf') {
          ## used in sampling functions
          selmax <- (datum >= info$censormax)
          selmin <- (datum <= info$censormin)
          ##
          if (info$transform == 'log') {
            datum <- log(datum - info$domainmin)
          } else if (info$transform == 'logminus') {
            datum <- log(info$domainmax - datum)
          } else if (info$transform == 'Q') {
            datum <- Qf(0.5 +
                        (datum -
                         (info$domainmin + info$domainmax)/2) /
                        (info$domainmax - info$domainmin)
                        )
          }
          ##
          datum[selmax] <- +Inf
          datum[selmin] <- -Inf
          datum <- (datum - info$tlocation) / info$tscale

        ## } else if (Dout == 'sleft') {
        ##   ## in sampling functions
        ##   datum[] <- info$tcensormin
        ## 
        ## } else if (Dout == 'sright') {
        ##   ## in sampling functions
        ##   datum[] <- info$tcensormax

        } else if (Dout == 'mi') {
          ## used in mutualinfo
          datum[datum >= info$tcensormax] <- +Inf
          datum[datum <= info$tcensormin] <- -Inf

        } else if (Dout == 'original') {
          ## transformation from MCMC representation
          ## to original domain
          datum <- datum * info$tscale + info$tlocation
          if (info$transform == 'log') {
            datum <- exp(datum) + info$domainmin
          } else if (info$transform == 'logminus') {
            datum <- info$domainmax - exp(datum)
          } else if (info$transform == 'Q') {
            datum <- (invQf(datum) - 0.5) * (info$domainmax - info$domainmin) +
              (info$domainmin + info$domainmax)/2
          }
          datum[datum <= info$censormin] <- censormin
          datum[datum >= info$censormax] <- censormax

        } else {
          stop('Unknown transformation for variate', v)
        }

#### Latent *** NOT USED AT THE MOMENT***
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

        } else if (Lout == 'mi') { # in mutualinfo
          datum <- data.matrix(x[, v, drop = FALSE])
          if (useLquantiles) {
            datum <- datum * info$tscale + info$tlocation
          }
          datum <- round(invQf(datum) * info$Nvalues + 0.5)
          datum[datum < 1] <- 1L
          datum[datum > info$Nvalues] <- info$Nvalues

        } else if (Lout != 'normalized') { # in sampling functions
          stop('Unknown transformation for variate', v)
        }

        ## Ordinal
      } else if (info$mcmctype == 'O') {
        bvalues <- 1:info$Nvalues
        names(bvalues) <- unlist(info[paste0('V', bvalues)])

        if (Oout == 'numeric') {
          datum <- bvalues[as.character(datum)]

        } else if (Oout == 'original') {
          datum <- names(bvalues[datum])

        } else if (Oout != 'mi') {
          stop('Unknown transformation for variate', v)
        }

        ## Nominal
      } else if (info$mcmctype == 'N') {
        bvalues <- 1:info$Nvalues
        names(bvalues) <- unlist(info[paste0('V', bvalues)])

        if (Nout == 'numeric') {
          datum <- bvalues[as.character(datum)]

        } else if (Nout == 'original') {
          datum <- names(bvalues[datum])

        } else if (Nout != 'mi') {
          stop('Unknown transformation for variate', v)
        }

        ## Binary
      } else if (info$mcmctype == 'B') {
        bvalues <- 0:1
        names(bvalues) <- unlist(info[c('V1', 'V2')])

        if (Bout == 'numeric') {
          datum <- bvalues[as.character(datum)]

        } else if (Bout == 'original') {
          datum <- names(bvalues[datum + 1])

        } else if (Bout != 'mi') {
          stop('Unknown transformation for variate', v)
        }
      }

    }
    datum
  }), ncol = ncol(x), dimnames = list(NULL, colnames(x)))
}

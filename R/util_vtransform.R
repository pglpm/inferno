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
as.data.frame(lapply(colnames(x), function(v) {
    datum <- x[[v]]
    xinfo <- as.list(auxmetadata[auxmetadata$name == v, ])

    if (invjacobian) {
#### Calculation of reciprocal Jacobian factors
      if (xinfo$mcmctype %in% c('B', 'N', 'O', 'L')) {
        datum <- rep(1L, length(datum))
      } else {
        if (xinfo$transform == 'log') {
          datum <- (datum - xinfo$domainmin) * xinfo$tscale
        } else if (xinfo$transform == 'logminus') {
          datum <- (xinfo$domainmax - datum) * xinfo$tscale
        } else if (xinfo$transform == 'Q') {
          datum <- Qf(0.5 +
                      (datum - (xinfo$domainmin + xinfo$domainmax)/2) /
                      (xinfo$domainmax - xinfo$domainmin)
                      )
          datum <- DQf(datum) * xinfo$tscale * (xinfo$domainmax - xinfo$domainmin)
        } else if (xinfo$transform == 'identity') {
          datum[] <- xinfo$tscale
        } else {
          stop('Unknown transformation for variate', v)
        }
        xv <- data.matrix(x[, v, drop = FALSE])
        if (xinfo$mcmctype %in% c('C', 'D')) {
          datum[(xv >= xinfo$censormax) | (xv <= xinfo$censormin)] <- 1L
        }
        datum[is.na(xv)] <- 1L
      }
    } else {
#### Transformation to internal value for MCMC

#### Continuous, open domain

      if (xinfo$mcmctype == 'R') {

        if (Rout == 'normalized') {
          ## used for MCMC
          if (xinfo$transform == 'log') {
            datum <- log(datum - xinfo$domainmin)
          } else if (xinfo$transform == 'logminus') {
            datum <- log(xinfo$domainmax - datum)
          } else if (xinfo$transform == 'Q') {
            datum <- Qf(0.5 +
                        (datum - (xinfo$domainmin + xinfo$domainmax)/2) /
                        (xinfo$domainmax - xinfo$domainmin)
                        )
          }
          datum <- (datum - xinfo$tlocation) / xinfo$tscale

        } else if (Rout == 'original') {
          ## transformation from MCMC representation
          ## to original domain
          datum <- datum * xinfo$tscale + xinfo$tlocation
          if (xinfo$transform == 'log') {
            datum <- exp(datum) + xinfo$domainmin
          } else if (xinfo$transform == 'logminus') {
            datum <- xinfo$domainmax - exp(datum)
          } else if (xinfo$transform == 'Q') {
            datum <- (invQf(datum) - 0.5) * (xinfo$domainmax - xinfo$domainmin) +
              (xinfo$domainmin + xinfo$domainmax)/2
          }

        } else if (Rout != 'mi'){
          ## used in mutualinfo()
          stop('Unknown transformation for variate', v)
        }

#### Continuous, closed domain
      } else if (xinfo$mcmctype == 'C') {

        if (Cout == 'init') {
          ## init value used for MCMC
          selna <- is.na(datum)
          selmin <- !is.na(datum) & (datum <= xinfo$censormin)
          selmax <- !is.na(datum) & (datum >= xinfo$censormax)
          selmid <- !is.na(datum) &
            (datum < xinfo$censormax) & (datum > xinfo$censormin)
          datum[selna] <- 0L
          datum[selmin] <- xinfo$tcensormin - 0.125
          datum[selmax] <- xinfo$tcensormax + 0.125
          datum[selmid] <- NA # so that it's discarded by Nimble

        } else if (Cout == 'left') {
          ## used for MCMC
          sel <- is.na(datum) | (datum < xinfo$censormax)
          datum[sel] <- -Inf
          datum[!sel] <- xinfo$tcensormax

        } else if (Cout == 'right') {
          ## used for MCMC
          sel <- is.na(datum) | (datum > xinfo$censormin)
          datum[sel] <- +Inf
          datum[!sel] <- xinfo$tcensormin

        } else if (Cout == 'lat') {
          ## latent variable used for MCMC
          datum[(datum >= xinfo$censormax) | (datum <= xinfo$censormin)] <- NA

          if (xinfo$transform == 'log') {
            datum <- log(datum - xinfo$domainmin)
          } else if (xinfo$transform == 'logminus') {
            datum <- log(xinfo$domainmax - datum)
          } else if (xinfo$transform == 'Q') {
            datum <- Qf(0.5 +
                        (datum -
                         (xinfo$domainmin + xinfo$domainmax)/2) /
                        (xinfo$domainmax - xinfo$domainmin)
                        )
          }
          ##
          datum <- (datum - xinfo$tlocation) / xinfo$tscale

        } else if (Cout == 'aux') {
          ## aux variable used for MCMC
          sel <- is.na(datum)
          datum[sel] <- NA
          datum[!sel] <- 1L

        } else if (Cout == 'boundisinf') {
          ## used in sampling functions
          selmax <- (datum >= xinfo$censormax)
          selmin <- (datum <= xinfo$censormin)
          ##
          if (xinfo$transform == 'log') {
            datum <- log(datum - xinfo$domainmin)
          } else if (xinfo$transform == 'logminus') {
            datum <- log(xinfo$domainmax - datum)
          } else if (xinfo$transform == 'Q') {
            datum <- Qf(0.5 +
                        (datum -
                         (xinfo$domainmin + xinfo$domainmax)/2) /
                        (xinfo$domainmax - xinfo$domainmin)
                        )
          }
          ##
          datum[selmax] <- +Inf
          datum[selmin] <- -Inf
          datum <- (datum - xinfo$tlocation) / xinfo$tscale

        ## } else if (Cout == 'sleft') {
        ##   ## used in sampling functions
        ##   datum[] <- xinfo$tcensormin
        ## 
        ## } else if (Cout == 'sright') {
        ##   ## used in sampling functions
        ##   datum[] <- xinfo$tcensormax

        } else if (Cout == 'mi') {
          ## used in mutualinfo
          datum[datum >= xinfo$tcensormax] <- +Inf
          datum[datum <= xinfo$tcensormin] <- -Inf

        } else if (Cout == 'original') {
          ## transformation from MCMC representation
          ## to original domain
          datum <- datum * xinfo$tscale + xinfo$tlocation
          if (xinfo$transform == 'log') {
            datum <- exp(datum) + xinfo$domainmin
          } else if (xinfo$transform == 'logminus') {
            datum <- xinfo$domainmax - exp(datum)
          } else if (xinfo$transform == 'Q') {
            datum <- (invQf(datum) - 0.5) * (xinfo$domainmax - xinfo$domainmin) +
              (xinfo$domainmin + xinfo$domainmax)/2
          }
          datum[datum <= xinfo$censormin] <- xinfo$censormin
          datum[datum >= xinfo$censormax] <- xinfo$censormax

        } else {
          stop('Unknown transformation for variate', v)
        }



#### Continuous rounded
      } else if (xinfo$mcmctype == 'D') {
        ## tcensormin <- xinfo$censormin
        ## tcensormax <- xinfo$censormax
        ## leftbound <- pmax(datum - xinfo$step, xinfo$domainmin, na.rm = T)
        ## leftbound[leftbound <= xinfo$censormin] <- xinfo$domainmin
        ## leftbound[leftbound >= xinfo$censormax] <- xinfo$censormax
        ## rightbound <- pmin(datum + xinfo$step, xinfo$domainmax, na.rm = T)
        ## rightbound[rightbound >= xinfo$censormax] <- xinfo$domainmax
        ## rightbound[rightbound <= xinfo$censormin] <- xinfo$censormin
        ## if (xinfo$transform == 'log') {
        ##   datum <- log(datum - xinfo$domainmin)
        ##   leftbound <- log(leftbound - xinfo$domainmin)
        ##   rightbound <- log(rightbound - xinfo$domainmin)
        ##   tcensormin <- log(tcensormin - xinfo$domainmin)
        ##   tcensormax <- log(tcensormax - xinfo$domainmin)
        ## } else if (xinfo$transform == 'logminus') {
        ##   datum <- log(xinfo$domainmax - datum)
        ##   leftbound <- log(xinfo$domainmax - leftbound)
        ##   rightbound <- log(xinfo$domainmax - rightbound)
        ##   tcensormin <- log(xinfo$domainmax - tcensormin)
        ##   tcensormax <- log(xinfo$domainmax - tcensormax)
        ## } else if (xinfo$transform == 'Q') {
        ##   datum <- Qf((datum - xinfo$domainmin) / (xinfo$domainmax - xinfo$domainmin))
        ##   leftbound <- Qf((leftbound - xinfo$domainmin) / (xinfo$domainmax - xinfo$domainmin))
        ##   rightbound <- Qf((rightbound - xinfo$domainmin) / (xinfo$domainmax - xinfo$domainmin))
        ##   tcensormin <- Qf((tcensormin - xinfo$domainmin) / (xinfo$domainmax - xinfo$domainmin))
        ##   tcensormax <- Qf((tcensormax - xinfo$domainmin) / (xinfo$domainmax - xinfo$domainmin))
        ## }
        ## datum <- (datum-xinfo$tlocation)/xinfo$tscale
        ## rightbound <- (rightbound-xinfo$tlocation)/xinfo$tscale
        ## leftbound <- (leftbound-xinfo$tlocation)/xinfo$tscale
        ## tcensormax <- (tcensormax-xinfo$tlocation)/xinfo$tscale
        ## tcensormin <- (tcensormin-xinfo$tlocation)/xinfo$tscale

        if (Dout == 'init') {
          ## init value used for MCMC
          selna <- is.na(datum)
          selmin <- !is.na(datum) & (datum <= xinfo$censormin)
          selmax <- !is.na(datum) & (datum >= xinfo$censormax)
          datum[selna] <- 0L
          datum[selmin] <- xinfo$tcensormin - 0.125
          datum[selmax] <- xinfo$tcensormax + 0.125

        } else if (Dout == 'left') {
          ## used for MCMC
          datum <- pmax(datum - xinfo$step, xinfo$domainmin, na.rm = T)
          datum[datum <= xinfo$censormin] <- xinfo$domainmin
          datum[datum >= xinfo$censormax] <- xinfo$censormax
          ##
          if (xinfo$transform == 'log') {
            datum <- log(datum - xinfo$domainmin)
          } else if (xinfo$transform == 'logminus') {
            datum <- log(xinfo$domainmax - datum)
          } else if (xinfo$transform == 'Q') {
            datum <- Qf(0.5 +
                        (datum -
                         (xinfo$domainmin + xinfo$domainmax)/2) /
                        (xinfo$domainmax - xinfo$domainmin)
                        )
          }
          datum <- (datum - xinfo$tlocation) / xinfo$tscale

        } else if (Dout == 'right') {
          ## used for MCMC
          datum <- pmin(datum + xinfo$step, xinfo$domainmax, na.rm = T)
          datum[datum >= xinfo$censormax] <- xinfo$domainmax
          datum[datum <= xinfo$censormin] <- xinfo$censormin
          ##
          if (xinfo$transform == 'log') {
            datum <- log(datum - xinfo$domainmin)
          } else if (xinfo$transform == 'logminus') {
            datum <- log(xinfo$domainmax - datum)
          } else if (xinfo$transform == 'Q') {
            datum <- Qf(0.5 +
                        (datum -
                         (xinfo$domainmin + xinfo$domainmax)/2) /
                        (xinfo$domainmax - xinfo$domainmin)
                        )
          }
          datum <- (datum - xinfo$tlocation) / xinfo$tscale

        } else if (Dout == 'aux') {
          ## aux variable used for MCMC
          sel <- is.na(datum)
          datum[sel] <- NA
          datum[!sel] <- 1L

        } else if (Dout == 'boundisinf') {
          ## used in sampling functions
          selmax <- (datum >= xinfo$censormax)
          selmin <- (datum <= xinfo$censormin)
          ##
          if (xinfo$transform == 'log') {
            datum <- log(datum - xinfo$domainmin)
          } else if (xinfo$transform == 'logminus') {
            datum <- log(xinfo$domainmax - datum)
          } else if (xinfo$transform == 'Q') {
            datum <- Qf(0.5 +
                        (datum -
                         (xinfo$domainmin + xinfo$domainmax)/2) /
                        (xinfo$domainmax - xinfo$domainmin)
                        )
          }
          ##
          datum[selmax] <- +Inf
          datum[selmin] <- -Inf
          datum <- (datum - xinfo$tlocation) / xinfo$tscale

        ## } else if (Dout == 'sleft') {
        ##   ## in sampling functions
        ##   datum[] <- xinfo$tcensormin
        ## 
        ## } else if (Dout == 'sright') {
        ##   ## in sampling functions
        ##   datum[] <- xinfo$tcensormax

        } else if (Dout == 'mi') {
          ## used in mutualinfo
          datum[datum >= xinfo$tcensormax] <- +Inf
          datum[datum <= xinfo$tcensormin] <- -Inf

        } else if (Dout == 'original') {
          ## transformation from MCMC representation
          ## to original domain
          datum <- datum * xinfo$tscale + xinfo$tlocation
          if (xinfo$transform == 'log') {
            datum <- exp(datum) + xinfo$domainmin
          } else if (xinfo$transform == 'logminus') {
            datum <- xinfo$domainmax - exp(datum)
          } else if (xinfo$transform == 'Q') {
            datum <- (invQf(datum) - 0.5) * (xinfo$domainmax - xinfo$domainmin) +
              (xinfo$domainmin + xinfo$domainmax)/2
          }
          datum[datum <= xinfo$censormin] <- xinfo$censormin
          datum[datum >= xinfo$censormax] <- xinfo$censormax

        } else {
          stop('Unknown transformation for variate', v)
        }

#### Latent *** NOT USED AT THE MOMENT***
      } else if (xinfo$mcmctype == 'L') {
        olocation <- (xinfo$Nvalues * xinfo$domainmin - xinfo$domainmax) / (xinfo$Nvalues - 1)
        oscale <- (xinfo$domainmax - xinfo$domainmin) / (xinfo$Nvalues - 1)
        ##
        ## output is in range 1 to Nvalues
        datum <- round((datum - olocation) / oscale)

        if (Lout == 'init') { # in sampling functions or init MCMC
          datum[is.na(datum)] <- xinfo$Nvalues / 2 + 0.5
          datum <- Qf((datum - 0.5) / xinfo$Nvalues)
          datum[datum == +Inf] <- 1e6
          datum[datum == -Inf] <- -1e6
          if (useLquantiles) {
            datum <- (datum - xinfo$tlocation) / xinfo$tscale
          }

        } else if (Lout == 'left') { # as left for MCMC
          datum <- Qf(pmax(0, datum - 1L) / xinfo$Nvalues)
          if (useLquantiles) {
            datum <- (datum - xinfo$tlocation) / xinfo$tscale
          }
          datum[is.na(datum)] <- -Inf

        } else if (Lout == 'right') { # as right for MCMC
          datum <- Qf(pmin(xinfo$Nvalues, datum) / xinfo$Nvalues)
          if (useLquantiles) {
            datum <- (datum - xinfo$tlocation) / xinfo$tscale
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
            datum <- datum * xinfo$tscale + xinfo$tlocation
          }
          datum <- round(invQf(datum) * xinfo$Nvalues + 0.5)
          datum[datum < 1] <- 1L
          datum[datum > xinfo$Nvalues] <- xinfo$Nvalues

        } else if (Lout != 'normalized') { # in sampling functions
          stop('Unknown transformation for variate', v)
        }

        ## Ordinal
      } else if (xinfo$mcmctype == 'O') {
        bvalues <- 1:xinfo$Nvalues
        names(bvalues) <- unlist(xinfo[paste0('V', bvalues)])

        if (Oout == 'numeric') {
          datum <- bvalues[as.character(datum)]

        } else if (Oout == 'original') {
          datum <- names(bvalues[datum])

        } else if (Oout != 'mi') {
          stop('Unknown transformation for variate', v)
        }

        ## Nominal
      } else if (xinfo$mcmctype == 'N') {
        bvalues <- 1:xinfo$Nvalues
        names(bvalues) <- unlist(xinfo[paste0('V', bvalues)])

        if (Nout == 'numeric') {
          datum <- bvalues[as.character(datum)]

        } else if (Nout == 'original') {
          datum <- names(bvalues[datum])

        } else if (Nout != 'mi') {
          stop('Unknown transformation for variate', v)
        }

        ## Binary
      } else if (xinfo$mcmctype == 'B') {
        bvalues <- 0:1
        names(bvalues) <- unlist(xinfo[c('V1', 'V2')])

        if (Bout == 'numeric') {
          datum <- bvalues[as.character(datum)]

        } else if (Bout == 'original') {
          datum <- names(bvalues[datum + 1L])

        } else if (Bout != 'mi') {
          stop('Unknown transformation for variate', v)
        }
      }

    }
    datum
  }), row.names = NULL,#ncol = ncol(x), dimnames = list(NULL, colnames(x)))
  col.names = colnames(x))
}

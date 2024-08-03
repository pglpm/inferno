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
#'
#' @return data frame of transformed variates
vtransform <- function(
    x,
    auxmetadata,
    Rout = NULL,
    Cout = NULL,
    Dout = NULL,
    Lout = NULL,
    Bout = NULL,
    Oout = NULL,
    Nout = NULL,
    variates = NULL,
    invjacobian = FALSE
) {
    useLquantiles <- TRUE

    x <- as.data.frame(cbind(x))
    if (!(missing(variates) || is.null(variates))) {
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
                    datum <- util_Q(0.5 +
                                (datum - (xinfo$domainmin + xinfo$domainmax)/2) /
                                (xinfo$domainmax - xinfo$domainmin)
                    )
                    datum <- util_DQ(datum) * xinfo$tscale * (xinfo$domainmax - xinfo$domainmin)
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
                        datum <- util_Q(0.5 +
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
                        datum <- (util_invQ(datum) - 0.5) * (xinfo$domainmax - xinfo$domainmin) +
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
                    ## output:
                    ## for datapoints with known value within domain
                    ## it must be set to NA (will be ignored by Nimble)
                    ## because it's already given as data;
                    ## for datapoints with unknown value
                    ## we set it to 0;
                    ## for datapoints at boundaries
                    ## we set it slightly outside the boundaries.
                    selna <- is.na(datum)
                    selmin <- !is.na(datum) & (datum <= xinfo$censormin)
                    selmax <- !is.na(datum) & (datum >= xinfo$censormax)
                    selmid <- !is.na(datum) &
                        (datum < xinfo$censormax) & (datum > xinfo$censormin)
                    datum[selna] <- 0L
                    datum[selmin] <- xinfo$tcensormin - 0.125
                    datum[selmax] <- xinfo$tcensormax + 0.125
                    datum[selmid] <- NA

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
                    datum[(datum >= xinfo$censormax) |
                          (datum <= xinfo$censormin)] <- NA

                    if (xinfo$transform == 'log') {
                        datum <- log(datum - xinfo$domainmin)
                    } else if (xinfo$transform == 'logminus') {
                        datum <- log(xinfo$domainmax - datum)
                    } else if (xinfo$transform == 'Q') {
                        datum <- util_Q(0.5 +
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
                        datum <- util_Q(0.5 +
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
                        datum <- (util_invQ(datum) - 0.5) * (xinfo$domainmax - xinfo$domainmin) +
                            (xinfo$domainmin + xinfo$domainmax)/2
                    }
                    datum[datum <= xinfo$censormin] <- xinfo$censormin
                    datum[datum >= xinfo$censormax] <- xinfo$censormax

                } else {
                    stop('Unknown transformation for variate', v)
                }



#### Continuous rounded
            } else if (xinfo$mcmctype == 'D') {

                if (Dout == 'init') {
                    ## init value used for MCMC
                    ## output:
                    ## for datapoints with known value within domain
                    ## we set it to their transformed value;
                    ## for datapoints with unknown value
                    ## we set it to 0;
                    ## for datapoints at boundaries
                    ## we set it slightly outside the boundaries.
                    selna <- is.na(datum)
                    selmin <- !is.na(datum) & (datum <= xinfo$censormin)
                    selmax <- !is.na(datum) & (datum >= xinfo$censormax)
                    if (xinfo$transform == 'log') {
                        datum <- log(datum - xinfo$domainmin)
                    } else if (xinfo$transform == 'logminus') {
                        datum <- log(xinfo$domainmax - datum)
                    } else if (xinfo$transform == 'Q') {
                        datum <- util_Q(0.5 +
                                    (datum -
                                     (xinfo$domainmin + xinfo$domainmax)/2) /
                                    (xinfo$domainmax - xinfo$domainmin)
                        )
                    }
                    ##
                    datum <- (datum - xinfo$tlocation) / xinfo$tscale
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
                        datum <- util_Q(0.5 +
                                    (datum -
                                     (xinfo$domainmin + xinfo$domainmax)/2) /
                                    (xinfo$domainmax - xinfo$domainmin)
                        )
                    }
                    datum <- (datum - xinfo$tlocation) / xinfo$tscale

                } else if (Dout == 'right') {
                    ## used for MCMC
                    datum <- pmin(datum + xinfo$step, xinfo$domainmax, na.rm = T)
                    datum[datum <= xinfo$censormin] <- xinfo$censormin
                    datum[datum >= xinfo$censormax] <- xinfo$domainmax
                    ##
                    if (xinfo$transform == 'log') {
                        datum <- log(datum - xinfo$domainmin)
                    } else if (xinfo$transform == 'logminus') {
                        datum <- log(xinfo$domainmax - datum)
                    } else if (xinfo$transform == 'Q') {
                        datum <- util_Q(0.5 +
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
                        datum <- util_Q(0.5 +
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
                        datum <- (util_invQ(datum) - 0.5) * (xinfo$domainmax - xinfo$domainmin) +
                            (xinfo$domainmin + xinfo$domainmax)/2
                    }
                    datum[datum <= xinfo$censormin] <- xinfo$censormin
                    datum[datum >= xinfo$censormax] <- xinfo$censormax

                } else {
                    stop('Unknown transformation for variate', v)
                }

#### Latent *** NOT USED AT THE MOMENT***
            } else if (xinfo$mcmctype == 'L') {

                if (Lout == 'init') {
                    ## init value used for MCMC
                    olocation <- (xinfo$Nvalues * xinfo$domainmin -
                                  xinfo$domainmax) / (xinfo$Nvalues - 1)
                    oscale <- (xinfo$domainmax - xinfo$domainmin) /
                        (xinfo$Nvalues - 1)
                    ## output is in range 1 to Nvalues
                    datum <- round((datum - olocation) / oscale)

                    datum <- util_Q((datum - 0.5) / xinfo$Nvalues)
                    datum[is.na(datum)] <- 0
                    datum[datum == +Inf] <- 1e6
                    datum[datum == -Inf] <- -1e6
                    if (useLquantiles) {
                        datum <- (datum - xinfo$tlocation) / xinfo$tscale
                    }

                } else if (Lout == 'left') { # as left for MCMC
                    olocation <- (xinfo$Nvalues * xinfo$domainmin -
                                  xinfo$domainmax) / (xinfo$Nvalues - 1)
                    oscale <- (xinfo$domainmax - xinfo$domainmin) /
                        (xinfo$Nvalues - 1)
                    ## output is in range 1 to Nvalues
                    datum <- round((datum - olocation) / oscale)

                    datum <- util_Q(pmax(0, datum - 1L) / xinfo$Nvalues)
                    if (useLquantiles) {
                        datum <- (datum - xinfo$tlocation) / xinfo$tscale
                    }
                    datum[is.na(datum)] <- -Inf

                } else if (Lout == 'right') { # as right for MCMC
                    olocation <- (xinfo$Nvalues * xinfo$domainmin -
                                  xinfo$domainmax) / (xinfo$Nvalues - 1)
                    oscale <- (xinfo$domainmax - xinfo$domainmin) /
                        (xinfo$Nvalues - 1)
                    ## output is in range 1 to Nvalues
                    datum <- round((datum - olocation) / oscale)

                    datum <- util_Q(pmin(xinfo$Nvalues, datum) / xinfo$Nvalues)
                    if (useLquantiles) {
                        datum <- (datum - xinfo$tlocation) / xinfo$tscale
                    }
                    datum[is.na(datum)] <- +Inf

                } else if (Lout == 'aux') { # aux variable in MCMC
                    sel <- is.na(datum)
                    datum[sel] <- NA
                    datum[!sel] <- 1L

                ## } else if (Lout == 'boundisinf') { # in output functions
                ##     olocation <- (xinfo$Nvalues * xinfo$domainmin -
                ##                   xinfo$domainmax) / (xinfo$Nvalues - 1)
                ##     oscale <- (xinfo$domainmax - xinfo$domainmin) /
                ##         (xinfo$Nvalues - 1)
                ##     ## output is in range 1 to Nvalues
                ##     datum <- round((datum - olocation) / oscale) - 1L

                } else if (Lout == 'mi') { # in mutualinfo
                    if (useLquantiles) {
                        datum <- datum * xinfo$tscale + xinfo$tlocation
                    }
                    datum <- round(util_invQ(datum) * xinfo$Nvalues + 0.5)
                    datum[datum < 1] <- 1L
                    datum[datum > xinfo$Nvalues] <- xinfo$Nvalues

                } else if (Lout == 'normalized') { # in sampling functions
                    olocation <- (xinfo$Nvalues * xinfo$domainmin -
                                  xinfo$domainmax) / (xinfo$Nvalues - 1)
                    oscale <- (xinfo$domainmax - xinfo$domainmin) /
                        (xinfo$Nvalues - 1)
                    ## output is in range 1 to Nvalues
                    datum <- round((datum - olocation) / oscale)

                } else {
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
    }), row.names = NULL,  col.names = colnames(x))
}

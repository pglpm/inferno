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
#'  'boundnormalized': for sampling functions
#'  'boundisinf': for sampling functions
#'  'mi': for use in mutualinfo()
#'  'original': original representation
#' @param Dout string, output of D-type variate, with possible values:
#'  'init': for internal MCMC use (init input)
#'  'left', 'right': for internal MCMC use
#'  'aux': for internal MCMC use
#'  'boundisinf': for sampling functions
#'  'normalized': for sampling functions
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
#' @param logjacobianOr logical or `NULL`: output is the log-Jacobian in orginal or transformed domain? `NULL` (default) means do not calculate the log-Jacobians
#'
#' @return data frame of transformed variates, or their log-Jacobians
#' @keywords internal
vtransform <- function(
    x,
    auxmetadata,
    Rout = NULL,
    Cout = NULL,
    Dout = NULL,
    Bout = NULL,
    Oout = NULL,
    Nout = NULL,
    variates = NULL,
    logjacobianOr = NULL
) {
    useLquantiles <- TRUE

    if(!is.data.frame(x)){
        if(is.vector(x)){
            dim(x) <- c(length(x), 1)
        }
        x <- as.data.frame(x)
    }
    if (!is.null(variates)) {
        colnames(x) <- variates
    }
    as.data.frame(lapply(colnames(x), function(v) {
        datum <- x[[v]]
        with(data = as.list(auxmetadata[auxmetadata$name == v, ]),
        {

            if (is.null(logjacobianOr)) {
#### Transformation to internal value for MCMC

#### Continuous, open domain
                if (mcmctype == 'R') {

                    if (Rout == 'normalized') {
                        ## used for MCMC
                        if (transform == 'log') {
                            datum <- log(datum - domainmin)
                        } else if (transform == 'logminus') {
                            datum <- -log(domainmax - datum)
                        } else if (transform == 'Q') {
                            datum <- util_Q(0.5 +
                                                (datum - (domainmin + domainmax)/2) /
                                                (domainmax - domainmin)
                            )
                        }
                        datum <- (datum - tlocation) / tscale

                    } else if (Rout == 'original') {
                        ## transformation from MCMC representation
                        ## to original domain
                        datum <- datum * tscale + tlocation
                        if (transform == 'log') {
                            datum <- exp(datum) + domainmin
                        } else if (transform == 'logminus') {
                            datum <- domainmax - exp(-datum)
                        } else if (transform == 'Q') {
                            datum <- (util_invQ(datum) - 0.5) * (domainmax - domainmin) +
                                (domainmin + domainmax)/2
                        }

                    } else if (Rout != 'mi'){
                        ## used in mutualinfo()
                        stop('Unknown transformation for variate', v)
                    }

#### Continuous, closed domain
                } else if (mcmctype == 'C') {

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
                        selmin <- !is.na(datum) & (datum <= domainmin)
                        selmax <- !is.na(datum) & (datum >= domainmax)
                        selmid <- !is.na(datum) &
                            (datum > domainmin) & (datum < domainmax)
                        datum[selna] <- 0L
                        datum[selmin] <- tdomainmin - 0.125
                        datum[selmax] <- tdomainmax + 0.125
                        datum[selmid] <- NA

                    } else if (Cout == 'left') {
                        ## used for MCMC
                        sel <- is.na(datum) | (datum < domainmax)
                        datum[sel] <- -Inf
                        datum[!sel] <- tdomainmax

                    } else if (Cout == 'right') {
                        ## used for MCMC
                        sel <- is.na(datum) | (datum > domainmin)
                        datum[sel] <- +Inf
                        datum[!sel] <- tdomainmin

                    } else if (Cout == 'lat') {
                        ## latent variable used for MCMC
                        ## Boundary values are not fixed, free to roam
                        datum[(datum >= domainmax) |
                                  (datum <= domainmin)] <- NA

                        ## Non-boundary values are fixed data
                        if (transform == 'log') {
                            datum <- log(datum - domainmin)
                        } else if (transform == 'logminus') {
                            datum <- -log(domainmax - datum)
                        } else if (transform == 'Q') {
                            datum <- util_Q(0.5 +
                                                (datum -
                                                     (domainmin + domainmax)/2) /
                                                (domainmax - domainmin))
                        }
                        ##
                        datum <- (datum - tlocation) / tscale

                    } else if (Cout == 'aux') {
                        ## aux, indicator variable used for MCMC
                        sel <- is.na(datum)
                        datum[sel] <- NA
                        datum[!sel] <- 1L

                    } else if (Cout == 'boundisna') {
                        ## used in sampling functions
                        ## points outside or on boundaries are moved to infinity
                        datum[datum <= domainmin] <- NA
                        datum[datum >= domainmax] <- NA
                        ##
                        if (transform == 'log') {
                            datum <- log(datum - domainmin)
                        } else if (transform == 'logminus') {
                            datum <- -log(domainmax - datum)
                        } else if (transform == 'Q') {
                            datum <- util_Q(0.5 +
                                                (datum -
                                                     (domainmin + domainmax)/2) /
                                                (domainmax - domainmin))
                        }
                        ##
                        datum <- (datum - tlocation) / tscale

                    } else if (Cout == 'leftbound') {
                        ## used in Pr()
                        ## non-boundary points are st to NA
                        ## boundaries are enforced
                        datum[datum > domainmin] <- NA
                        datum[datum < domainmin] <- domainmin
                        ##
                        if (transform == 'log') {
                            datum <- log(datum - domainmin)
                        } else if (transform == 'logminus') {
                            datum <- -log(domainmax - datum)
                        } else if (transform == 'Q') {
                            datum <- util_Q(0.5 +
                                                (datum -
                                                     (domainmin + domainmax)/2) /
                                                (domainmax - domainmin))
                        }
                        ##
                        datum <- (datum - tlocation) / tscale

                    } else if (Cout == 'rightbound') {
                        ## used in Pr()
                        ## non-boundary points are st to NA
                        ## boundaries are enforced
                        datum[datum < domainmax] <- NA
                        datum[datum > domainmax] <- domainmax
                        ##
                        if (transform == 'log') {
                            datum <- log(datum - domainmin)
                        } else if (transform == 'logminus') {
                            datum <- -log(domainmax - datum)
                        } else if (transform == 'Q') {
                            datum <- util_Q(0.5 +
                                                (datum -
                                                     (domainmin + domainmax)/2) /
                                                (domainmax - domainmin))
                        }
                        ##
                        datum <- (datum - tlocation) / tscale

                    } else if (Cout == 'boundisinf') {
                        ## used in sampling functions
                        ## points outside or on boundaries are moved to infinity
                        selmax <- (datum >= domainmax)
                        selmin <- (datum <= domainmin)
                        ##
                        if (transform == 'log') {
                            datum <- log(datum - domainmin)
                        } else if (transform == 'logminus') {
                            datum <- -log(domainmax - datum)
                        } else if (transform == 'Q') {
                            datum <- util_Q(0.5 +
                                                (datum -
                                                     (domainmin + domainmax)/2) /
                                                (domainmax - domainmin))
                        }
                        ##
                        datum <- (datum - tlocation) / tscale
                        datum[selmax] <- +Inf
                        datum[selmin] <- -Inf

                    } else if (Cout == 'boundnormalized') {
                        ## used in Pr()
                        ## points outside boundaries are moved to boundaries
                        datum[datum < domainmin] <- domainmin
                        datum[datum > domainmax] <- domainmax
                        ##
                        if (transform == 'log') {
                            datum <- log(datum - domainmin)
                        } else if (transform == 'logminus') {
                            datum <- -log(domainmax - datum)
                        } else if (transform == 'Q') {
                            datum <- util_Q(0.5 +
                                                (datum -
                                                     (domainmin + domainmax)/2) /
                                                (domainmax - domainmin))
                        }
                        ##
                        datum <- (datum - tlocation) / tscale

                    } else if (Cout == 'miboundisna') {
                        ## used in mutualinfo
                        datum[datum >= tdomainmax] <- NA
                        datum[datum <= tdomainmin] <- NA

                    } else if (Cout == 'mileftbound') {
                        ## used in mutualinfo
                        datum[datum > tdomainmin] <- NA
                        datum[datum < tdomainmin] <- tdomainmin

                    } else if (Cout == 'mirightbound') {
                        ## used in mutualinfo
                        datum[datum < tdomainmax] <- NA
                        datum[datum > tdomainmax] <- tdomainmax

                    } else if (Cout == 'original') {
                        ## transformation from MCMC representation
                        ## to original domain
                        datum <- datum * tscale + tlocation
                        if (transform == 'log') {
                            datum <- exp(datum) + domainmin
                        } else if (transform == 'logminus') {
                            datum <- domainmax - exp(-datum)
                        } else if (transform == 'Q') {
                            datum <- (util_invQ(datum) - 0.5) *
                                (domainmax - domainmin) +
                                (domainmin + domainmax) / 2
                        }
                        datum[datum < domainmin] <- domainmin
                        datum[datum > domainmax] <- domainmax

                    } else {
                        stop('Unknown transformation for variate ', v)
                    }

#### Continuous rounded
                } else if (mcmctype == 'D') {
                    ## NOTA BENE:
                    ## except for the 'original' (inverse-transformation) case,
                    ## we expect D-data to be given with correctly rounded values

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
                        datum <- (datum - tlocation) / tscale
                        datum[selna] <- 0L

                    } else if (Dout == 'left') {
                        ## used for MCMC
                        datum <- datum - halfstep
                        selmin <- is.na(datum) | (datum <= domainmin)
                        selmax <- (datum >= domainmaxminushs)
                        datum <- (datum - tlocation) / tscale
                        datum[selmin] <- -Inf
                        datum[selmax] <- tdomainmaxminushs

                    } else if (Dout == 'right') {
                        ## used for MCMC
                        datum <- datum + halfstep
                        selmin <- (datum <= domainminplushs)
                        selmax <- is.na(datum) | (datum >= domainmax)
                        datum <- (datum - tlocation) / tscale
                        datum[selmin] <- tdomainminplushs
                        datum[selmax] <- +Inf

                    } else if (Dout == 'aux') {
                        ## aux variable used for MCMC
                        sel <- is.na(datum)
                        datum[sel] <- NA
                        datum[!sel] <- 1L

                        ## ## not used anymore
                        ## } else if (Dout == 'boundisinf') {
                        ##     ## used in sampling functions
                        ##     selmin <- (datum - halfstep <= domainmin)
                        ##     selmax <- (datum + halfstep >= domainmax)
                        ##     ##
                        ##     datum <- (datum - tlocation) / tscale
                        ##     datum[selmin] <- -Inf
                        ##     datum[selmax] <- +Inf

                    } else if (Dout == 'boundnormalized') {
                        ## used in sampling functions
                        datum[datum < domainmin] <- domainmin
                        datum[datum > domainmax] <- domainmax
                        ##
                        datum <- (datum - tlocation) / tscale

                    } else if (Dout == 'boundisna') {
                        ## used in Pr()
                        datum[datum <= domainminplushs] <- NA
                        datum[datum >= domainmaxminushs] <- NA
                        datum <- (datum - tlocation) / tscale

                    } else if (Dout == 'leftbound') {
                        ## used in Pr()
                        datum[datum > domainminplushs] <- NA
                        datum[datum < domainmin] <- domainmin
                        datum <- (datum - tlocation) / tscale

                    } else if (Dout == 'rightbound') {
                        ## used in Pr()
                        datum[datum < domainmaxminushs] <- NA
                        datum[datum > domainmax] <- domainmax
                        datum <- (datum - tlocation) / tscale

                    } else if (Dout == 'miboundisna') {
                        ## used in mutualinfo
                        datum[datum <= tdomainminplushs] <- NA
                        datum[datum >= tdomainmaxminushs] <- NA
                        datum <- round((datum * tscale) / (2 * halfstep)) *
                            2 * halfstep / tscale

                    } else if (Dout == 'mileftbound') {
                        ## used in mutualinfo
                        datum[datum > tdomainminplushs] <- NA
                        datum[datum < tdomainmin] <- tdomainmin
                        datum <- round((datum * tscale) / (2 * halfstep)) *
                            2 * halfstep / tscale

                    } else if (Dout == 'mirightbound') {
                        ## used in mutualinfo
                        datum[datum < tdomainmaxminushs] <- NA
                        datum[datum > tdomainmax] <- tdomainmax
                        datum <- round((datum * tscale) / (2 * halfstep)) *
                            2 * halfstep / tscale

                    } else if (Dout == 'original') {
                        ## transformation from MCMC representation
                        ## to original domain
                        datum <- round((datum * tscale) / (2 * halfstep)) *
                            2 * halfstep + tlocation
                        datum[datum < domainmin] <- domainmin
                        datum[datum > domainmax] <- domainmax
                        ## datum[datum <= domainmin] <- domainmin
                        ## datum[datum >= domainmax] <- domainmax

                    } else {
                        stop('Unknown transformation for variate ', v)
                    }

                    ## Ordinal
                } else if (mcmctype == 'O') {

                    if (Oout == 'numeric') {
                        bvalues <- seq_len(Nvalues)
                        names(bvalues) <- sapply(bvalues, function(x)get(paste0('V', x)))
                        ## names(bvalues) <- unlist(xinfo[paste0('V', bvalues)])
                        datum <- bvalues[as.character(datum)]

                    } else if (Oout == 'index') {
                        bvalues <- seq_len(Nvalues)
                        names(bvalues) <- sapply(bvalues, function(x)get(paste0('V', x)))
                        ## names(bvalues) <- unlist(xinfo[paste0('V', bvalues)])
                        datum <- bvalues[as.character(datum)] + indexpos

                    } else if (Oout == 'original') {
                        bvalues <- seq_len(Nvalues)
                        names(bvalues) <- sapply(bvalues, function(x)get(paste0('V', x)))
                        ## names(bvalues) <- unlist(xinfo[paste0('V', bvalues)])
                        datum <- names(bvalues[datum])

                    } else if (Oout != 'mi') {
                        stop('Unknown transformation for variate ', v)
                    }

                    ## Nominal
                } else if (mcmctype == 'N') {

                    if (Nout == 'numeric') {
                        bvalues <- seq_len(Nvalues)
                        names(bvalues) <- sapply(bvalues, function(x)get(paste0('V', x)))
                        ## names(bvalues) <- unlist(xinfo[paste0('V', bvalues)])
                        datum <- bvalues[as.character(datum)]

                    } else if (Nout == 'index') {
                        bvalues <- seq_len(Nvalues)
                        names(bvalues) <- sapply(bvalues, function(x)get(paste0('V', x)))
                        ## names(bvalues) <- unlist(xinfo[paste0('V', bvalues)])
                        datum <- bvalues[as.character(datum)] + indexpos

                    } else if (Nout == 'original') {
                        bvalues <- seq_len(Nvalues)
                        names(bvalues) <- sapply(bvalues, function(x)get(paste0('V', x)))
                        ## names(bvalues) <- unlist(xinfo[paste0('V', bvalues)])
                        datum <- names(bvalues[datum])

                    } else if (Nout != 'mi') {
                        stop('Unknown transformation for variate', v)
                    }

                    ## Binary
                } else if (mcmctype == 'B') {
                    ## V1 corresponds to 0
                    ## V2 corresponds to 1
                    ## in the MC sampler, Bprob is the probability of V2
                    if (Bout == 'numeric') {
                        bvalues <- 0:1
                        names(bvalues) <- sapply(bvalues + 1L,
                            function(x)get(paste0('V', x)))
                        datum <- bvalues[as.character(datum)]

                    } else if (Bout == 'original') {
                        bvalues <- 0:1
                        names(bvalues) <- sapply(bvalues + 1L,
                            function(x)get(paste0('V', x)))
                        datum <- names(bvalues[datum + 1L])

                    } else if (Bout != 'mi') {
                        stop('Unknown transformation for variate ', v)
                    }
                }

#### Calculation of reciprocal log-Jacobian factors
            } else if (logjacobianOr) {
                ## 'datum' is given in the original domain

                if (mcmctype %in% c('R', 'C')) {
                    if (mcmctype %in% c('C')) {
                        selma <- datum <= domainminplushs | datum >= domainmaxminushs
                    }
                    if (transform == 'log') {
                        datum <- -log(datum - domainmin) - log(tscale)
                    } else if (transform == 'logminus') {
                        datum <- -log(domainmax - datum) - log(tscale)
                    } else if (transform == 'Q') {
                        datum <- util_Q(0.5 +
                                            (datum - (domainmin + domainmax)/2) /
                                            (domainmax - domainmin)
                        )
                        datum <- -log(util_invDQ(datum)) -
                            log(tscale) - log(domainmax - domainmin)
                    } else if (transform == 'identity') {
                        datum[!is.na(datum)] <- -log(tscale)
                    } else {
                        stop('Unknown transformation for variate ', v)
                    }
                    if (mcmctype %in% c('C')) {
                        datum[selma] <- 0
                    }
                    datum[is.na(datum)] <- 0
                } else {
                    datum <- numeric(length(datum))
                }

            } else if (!logjacobianOr) {

                ## 'datum' is given in the transformed domain

                if (mcmctype %in% c('R', 'C')) {
                    if (mcmctype %in% c('C')) {
                        selma <- datum <= tdomainminplushs | datum >= tdomainmaxminushs
                    }
                    if (transform == 'log') {
                        datum <- tscale * datum + tlocation - log(tscale)
                    } else if (transform == 'logminus') {
                        datum <- -tscale * datum - tlocation - log(tscale)
                    } else if (transform == 'Q') {
                        datum <- -log(util_invDQ(tscale * datum + tlocation)) -
                            log(tscale) - log(domainmax - domainmin)
                    } else if (transform == 'identity') {
                        datum[!is.na(datum)] <- -log(tscale)
                    } else {
                        stop('Unknown transformation for variate ', v)
                    }
                    if (mcmctype %in% c('C')) {
                        datum[selma] <- 0
                    }
                    datum[is.na(datum)] <- 0
                } else {
                    datum <- numeric(length(datum))
                }

            } else {
                stop('Unknown "logjacobianOr" value')
            }

            datum
        })
    }), row.names = NULL,  col.names = colnames(x))
}

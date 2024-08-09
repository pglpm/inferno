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
        with(data = as.list(auxmetadata[auxmetadata$name == v, ]),
        {

        if (invjacobian) {
#### Calculation of reciprocal Jacobian factors
            if (mcmctype %in% c('R', 'C')) {
                if (mcmctype %in% c('C')) {
                    selma <- datum <= leftbound | datum >= rightbound
                }
                if (transform == 'log') {
                    datum <- (datum - domainmin) * tscale
                } else if (transform == 'logminus') {
                    datum <- (domainmax - datum) * tscale
                } else if (transform == 'Q') {
                    datum <- util_Q(0.5 +
                                    (datum - (domainmin + domainmax)/2) /
                                    (domainmax - domainmin)
                    )
                    datum <- util_DQ(datum) * tscale * (domainmax - domainmin)
                } else if (transform == 'identity') {
                    datum[] <- tscale
                } else {
                    stop('Unknown transformation for variate ', v)
                }
                if (mcmctype %in% c('C')) {
                    datum[selma] <- 1
                }
                datum[is.na(datum)] <- 1
            } else {
                datum <- rep(1, length(datum))
            }
        } else {
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
                    selmin <- !is.na(datum) & (datum <= leftbound)
                    selmax <- !is.na(datum) & (datum >= rightbound)
                    selmid <- !is.na(datum) &
                        (datum > leftbound) & (datum < rightbound)
                    datum[selna] <- 0L
                    datum[selmin] <- tleftbound - 0.125
                    datum[selmax] <- trightbound + 0.125
                    datum[selmid] <- NA

                } else if (Cout == 'left') {
                    ## used for MCMC
                    sel <- is.na(datum) | (datum < rightbound)
                    datum[sel] <- -Inf
                    datum[!sel] <- trightbound

                } else if (Cout == 'right') {
                    ## used for MCMC
                    sel <- is.na(datum) | (datum > leftbound)
                    datum[sel] <- +Inf
                    datum[!sel] <- tleftbound

                } else if (Cout == 'lat') {
                    ## latent variable used for MCMC
                    ## Boundary values are not fixed, free to roam
                    datum[(datum >= rightbound) |
                          (datum <= leftbound)] <- NA

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

                } else if (Cout == 'boundisinf') {
                    ## used in sampling functions
                    selmax <- (datum >= rightbound)
                    selmin <- (datum <= leftbound)
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

                } else if (Cout == 'normalized') {
                    ## used only for debugging
                    datum[datum <= leftbound] <- domainmin
                    datum[datum >= rightbound] <- domainmax
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

                } else if (Cout == 'mi') {
                    ## used in mutualinfo
                    datum[datum >= trightbound] <- +Inf
                    datum[datum <= tleftbound] <- -Inf

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
                    datum[datum <= leftbound] <- leftbound
                    datum[datum >= rightbound] <- rightbound

                } else {
                    stop('Unknown transformation for variate ', v)
                }



#### Continuous rounded
            } else if (mcmctype == 'D') {

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
                    selmin <- is.na(datum) | (datum <= leftbound)
                    selmax <- (datum >= rightbound)
                    datum <- (datum - halfstep - tlocation) / tscale
                    datum[selmin] <- -Inf
                    datum[selmax] <- (rightbound - tlocation) / tscale

                } else if (Dout == 'right') {
                    ## used for MCMC
                    selmin <- (datum <= leftbound)
                    selmax <- is.na(datum) | (datum >= rightbound)
                    datum <- (datum + halfstep - tlocation) / tscale
                    datum[selmin] <- (leftbound - tlocation) / tscale
                    datum[selmax] <- +Inf

                } else if (Dout == 'aux') {
                    ## aux variable used for MCMC
                    sel <- is.na(datum)
                    datum[sel] <- NA
                    datum[!sel] <- 1L

                } else if (Dout == 'boundisinf') {
                    ## used in sampling functions
                    selmin <- (datum - halfstep <= domainmin)
                    selmax <- (datum + halfstep >= domainmax)
                    ##
                    datum <- (datum - tlocation) / tscale
                    datum[selmin] <- -Inf
                    datum[selmax] <- +Inf

                } else if (Dout == 'normalized') {
                    ## used in sampling functions
                    datum[datum <= leftbound] <- domainmin
                    datum[datum >= rightbound] <- domainmax
                    datum <- (datum - tlocation) / tscale

                } else if (Dout == 'mi') {
                    ## used in mutualinfo
                    datum[datum <= tleftbound] <- -Inf
                    datum[datum >= trightbound] <- +Inf

                } else if (Dout == 'original') {
                    ## transformation from MCMC representation
                    ## to original domain
                    datum <- tlocation + 2 * halfstep *
                        round((datum * tscale) / (2 * halfstep))
                    datum[datum <= leftbound] <- domainmin
                    datum[datum >= rightbound] <- domainmax

                } else {
                    stop('Unknown transformation for variate ', v)
                }

                ## Ordinal
            } else if (mcmctype == 'O') {

                if (Oout == 'numeric') {
                bvalues <- 1:Nvalues
                names(bvalues) <- sapply(bvalues, function(x)get(paste0('V',x)))
                ## names(bvalues) <- unlist(xinfo[paste0('V', bvalues)])
                    datum <- bvalues[as.character(datum)]

                } else if (Oout == 'original') {
                bvalues <- 1:Nvalues
                names(bvalues) <- sapply(bvalues, function(x)get(paste0('V',x)))
                ## names(bvalues) <- unlist(xinfo[paste0('V', bvalues)])
                    datum <- names(bvalues[datum])

                } else if (Oout != 'mi') {
                    stop('Unknown transformation for variate ', v)
                }

                ## Nominal
            } else if (mcmctype == 'N') {

                if (Nout == 'numeric') {
                bvalues <- 1:Nvalues
                names(bvalues) <- sapply(bvalues, function(x)get(paste0('V',x)))
                ## names(bvalues) <- unlist(xinfo[paste0('V', bvalues)])
                    datum <- bvalues[as.character(datum)]

                } else if (Nout == 'original') {
                bvalues <- 1:Nvalues
                names(bvalues) <- sapply(bvalues, function(x)get(paste0('V',x)))
                ## names(bvalues) <- unlist(xinfo[paste0('V', bvalues)])
                    datum <- names(bvalues[datum])

                } else if (Nout != 'mi') {
                    stop('Unknown transformation for variate', v)
                }

                ## Binary
            } else if (mcmctype == 'B') {

                if (Bout == 'numeric') {
                bvalues <- 0:1
                names(bvalues) <- sapply(bvalues+1, function(x)get(paste0('V',x)))
                ## names(bvalues) <- unlist(xinfo[c('V1', 'V2')])
                    datum <- bvalues[as.character(datum)]

                } else if (Bout == 'original') {
                bvalues <- 0:1
                names(bvalues) <- sapply(bvalues+1, function(x)get(paste0('V',x)))
                ## names(bvalues) <- unlist(xinfo[c('V1', 'V2')])
                    datum <- names(bvalues[datum + 1L])

                } else if (Bout != 'mi') {
                    stop('Unknown transformation for variate ', v)
                }
            }

        }
            datum
        })
    }), row.names = NULL,  col.names = colnames(x))
}

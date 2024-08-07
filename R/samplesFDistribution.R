#' Calculate posterior probabilities
#'
#' @param Y matrix or data.table: values of some variates of which we want
#'   the joint probability one variate per column
#' @param X matrix or data.table: values of some variates conditional on
#'   which we want the joint probability one variate per column
#' @param mcoutput Either a string with the name of a directory or full
#'   path for a 'FDistribution.rds' object, or such an object itself
#' @param subsamples numeric: number of Monte Carlo samples to use
#' @param jacobian include the Jacobian in the output probability
#' @param fn NULL or function to apply to the group of MCsamples,
#'   example 'identity' or 'mean'
#' @param combine how to combine the output for the different variate values
#' @param parallel logical or numeric: whether to use pre-existing parallel
#'   workers, or how many to create and use
#' @param silent logical: give warnings or updates in the computation
#'
#' @return A list with the mutual information, its error, and its unit
#'
#' @import parallel foreach doParallel
#'
#' @export
samplesFDistribution <- function(
    Y,
    X,
    mcoutput,
    subsamples,
    jacobian = TRUE,
    fn = NULL,
    combine = `rbind`,
    parallel = TRUE,
    silent = FALSE
) {
    if (!silent) {
        cat('\n')
    }

#### Determine the status of parallel processing
    if (is.logical(parallel) && parallel) {
        if (foreach::getDoParRegistered()) {
            if (!silent) {
                cat('Using already registered', foreach::getDoParName(),
                    'with', foreach::getDoParWorkers(), 'workers\n')
            }
            ncores <- foreach::getDoParWorkers()
        } else {
            if (!silent) {
                cat('No parallel backend registered.\n')
            }
            ncores <- 1
        }
    } else if (is.numeric(parallel) && parallel >= 2) {
        if (foreach::getDoParRegistered()) {
            ncores <- min(foreach::getDoParWorkers(), parallel)
            if (!silent) {
                cat('Using already registered', foreach::getDoParName(),
                    'with', foreach::getDoParWorkers(), 'workers\n')
                if(parallel > ncores) {
                    cat('NOTE: fewer pre-registered cores',
                        'than requested in the "parallel" argument.\n')
                }
            }
        } else {
            ## ##
            ## ## Alternative way to register cores;
            ## ## might need to be used for portability to Windows?
            ## registerDoSEQ()
            ## cl <- makePSOCKcluster(ncores)
            ## ##
            cl <- parallel::makeCluster(parallel)
            doParallel::registerDoParallel(cl)
            if (!silent) {
                cat('Registered', foreach::getDoParName(),
                    'with', foreach::getDoParWorkers(), 'workers\n')
            }
            ncores <- parallel
            closecoresonexit <- function(){
                cat('\nClosing connections to cores.\n')
                foreach::registerDoSEQ()
                parallel::stopCluster(cl)
                env <- foreach:::.foreachGlobals
                rm(list=ls(name=env), pos=env)
            }
            on.exit(closecoresonexit())
        }
    } else {
        if (!silent) {
            cat('No parallel backend registered.\n')
        }
        ncores <- 1
    }

#### Determine the status of parallel processing
    ## if (!missing(parallel) && is.logical(parallel) && parallel) {
    ##   if (foreach::getDoParRegistered()) {
    ##     if (!silent) {
    ##       cat('Using already registered', foreach::getDoParName(),
    ##           'with', foreach::getDoParWorkers(), 'workers\n')
    ##     }
    ##     ncores <- foreach::getDoParWorkers()
    ##   } else {
    ##     if (!silent) {
    ##       cat('No parallel backend registered.\n')
    ##     }
    ##     ncores <- 1
    ##   }
    ## } else if (!missing(parallel) && is.numeric(parallel) && parallel >= 2) {
    ##   if (foreach::getDoParRegistered()) {
    ##     ncores <- min(foreach::getDoParWorkers(), parallel)
    ##     if (!silent) {
    ##       cat('Using already registered', foreach::getDoParName(),
    ##           'with', ncores, 'workers\n')
    ##     }
    ##   } else {
    ##     cl <- parallel::makeCluster(parallel)
    ##     doParallel::registerDoParallel(cl)
    ##     if (!silent) {
    ##       cat('Registered', foreach::getDoParName(), 'with',
    ##           foreach::getDoParWorkers(), 'workers\n')
    ##     }
    ##   }
    ## } else {
    ##   if (!silent) {
    ##     cat('No parallel backend registered.\n')
    ##   }
    ##   ncores <- 1
    ## }

    if (ncores < 2) {
        `%dochains%` <- `%do%`
    } else {
        `%dochains%` <- `%dopar%`
    }

    ## Extract Monte Carlo output & auxmetadata
    ## If mcoutput is a string, check if it's a folder name or file name
    if (is.character(mcoutput)) {
                                        # Check if 'mcoutput' is a folder containing Fdistribution.rds
        if (file_test('-d', mcoutput) &&
            file.exists(file.path(mcoutput, 'Fdistribution.rds'))) {
            mcoutput <- readRDS(file.path(mcoutput, 'Fdistribution.rds'))
        } else {
            ## Assume 'mcoutput' the full path of Fdistributions.rds
            ## possibly without the file extension '.rds'
            mcoutput <- paste0(sub('.rds$', '', mcoutput), '.rds')
            if (file.exists(mcoutput)) {
                mcoutput <- readRDS(mcoutput)
            } else {
                stop("The argument 'mcoutput' must be a folder containing Fdistribution.rds, or the path to an rds-file containing the output from 'inferpopulation'.")
            }
        }
    }
    ## Add check to see that mcoutput is correct type of object?
    auxmetadata <- mcoutput$auxmetadata
    mcoutput$auxmetadata <- NULL
    mcoutput$auxinfo <- NULL

    ## Consistency checks
    if (length(dim(Y)) != 2) {
        stop('Y must have two dimensions')
    }
    if (missing(X)) {
        X <- NULL
    }
    if (!is.null(X) && length(dim(X)) != 2) {
        stop('X must be NULL or have two dimensions')
    }
    ##
    if (!is.null(X) && ncol(X) == 0) {
        X <- NULL
    }

    ## More consistency checks
    Yv <- colnames(Y)
    if (!all(Yv %in% auxmetadata$name)) {
        stop('unknown Y variates\n')
    }
    if (length(unique(Yv)) != length(Yv)) {
        stop('duplicate Y variates\n')
    }
    ##
    Xv <- colnames(X)
    if (!all(Xv %in% auxmetadata$name)) {
        stop('unknown X variates\n')
    }
    if (length(unique(Xv)) != length(Xv)) {
        stop('duplicate X variates\n')
    }
    ##
    if (length(intersect(Yv, Xv)) > 0) {
        stop('overlap in Y and X variates\n')
    }

#### Subsample and get ncomponents and nsamples
    ## source('mcsubset.R')
    if (!missing(subsamples) &&
        (is.numeric(subsamples) || (is.character(subsamples)
            && length(subsamples) == 1))) {
        if (is.character(subsamples)) {
            subsamples <- round(seq(1, ncol(mcoutput$W),
                length.out = as.numeric(subsamples)
            ))
        }
        mcoutput <- mcsubset(mcoutput, subsamples)
    }

    nsamples <- ncol(mcoutput$W)
    ncomponents <- nrow(mcoutput$W)

### Guide to indices:
    ## .i. = order in X/Y corresponding to appearance in vnames
    ## .t. = vnames present in X/Y, kept in their vnames-order

#### Type R
    vnames <- auxmetadata[auxmetadata$mcmctype == 'R', 'name']
    XiR <- match(vnames, Xv)
    XtR <- which(!is.na(XiR))
    XiR <- XiR[XtR]
    XnR <- length(XiR)
    ##
    YiR <- match(vnames, Yv)
    YtR <- which(!is.na(YiR))
    YiR <- YiR[YtR]
    YnR <- length(YiR)
    if (YnR > 0 || XnR > 0) {
        mcoutput$Rvar <- sqrt(mcoutput$Rvar)
    }

#### Type C
    vnames <- auxmetadata[auxmetadata$mcmctype == 'C', 'name']
    XiC <- match(vnames, Xv)
    XtC <- which(!is.na(XiC))
    XiC <- XiC[XtC]
    XnC <- length(XiC)
    ##
    YiC <- match(vnames, Yv)
    YtC <- which(!is.na(YiC))
    YiC <- YiC[YtC]
    YnC <- length(YiC)
    if (YnC > 0 || XnC > 0) {
        mcoutput$Cvar <- sqrt(mcoutput$Cvar)
        Cbounds <- cbind(
            auxmetadata[match(vnames, auxmetadata$name), 'tdomainmin'],
            ## sign is important here:
            ## for upper tail, take opposite mean and value
            -auxmetadata[match(vnames, auxmetadata$name), 'tdomainmax']
        )
    }

#### Type D
    vnames <- auxmetadata[auxmetadata$mcmctype == 'D', 'name']
    XiD <- match(vnames, Xv)
    XtD <- which(!is.na(XiD))
    XiD <- XiD[XtD]
    XnD <- length(XiD)
    ##
    YiD <- match(vnames, Yv)
    YtD <- which(!is.na(YiD))
    YiD <- YiD[YtD]
    YnD <- length(YiD)
    if (YnD > 0 || XnD > 0) {
        mcoutput$Dvar <- sqrt(mcoutput$Dvar)
        Dsteps <- auxmetadata[match(vnames, auxmetadata$name), 'halfstep'] /
            auxmetadata[match(vnames, auxmetadata$name), 'tscale']
        Dlefts <- auxmetadata[match(vnames, auxmetadata$name), 'tleftbound']
        Drights <- auxmetadata[match(vnames, auxmetadata$name), 'trightbound']
    }

#### Type O
    vnames <- auxmetadata[auxmetadata$mcmctype == 'O', 'name']
    XiO <- match(vnames, Xv)
    XtO <- which(!is.na(XiO))
    XiO <- XiO[XtO]
    XnO <- length(XiO)
    ##
    YiO <- match(vnames, Yv)
    YtO <- which(!is.na(YiO))
    YiO <- YiO[YtO]
    YnO <- length(YiO)

#### Type N
    vnames <- auxmetadata[auxmetadata$mcmctype == 'N', 'name']
    XiN <- match(vnames, Xv)
    XtN <- which(!is.na(XiN))
    XiN <- XiN[XtN]
    XnN <- length(XiN)
    ##
    YiN <- match(vnames, Yv)
    YtN <- which(!is.na(YiN))
    YiN <- YiN[YtN]
    YnN <- length(YiN)

#### Type B
    vnames <- auxmetadata[auxmetadata$mcmctype == 'B', 'name']
    XiB <- match(vnames, Xv)
    XtB <- which(!is.na(XiB))
    XiB <- XiB[XtB]
    XnB <- length(XiB)
    ##
    YiB <- match(vnames, Yv)
    YtB <- which(!is.na(YiB))
    YiB <- YiB[YtB]
    YnB <- length(YiB)


#### transformation of inputs
    Y2 <- as.matrix(vtransform(Y, auxmetadata = auxmetadata,
        Rout = 'normalized',
        Cout = 'boundisinf',
        Dout = 'normalized',
        Oout = 'numeric',
        Nout = 'numeric',
        Bout = 'numeric'))

    if (!is.null(X)) {
        X2 <- as.matrix(vtransform(X, auxmetadata = auxmetadata,
            Rout = 'normalized',
            Cout = 'boundisinf',
            Dout = 'normalized',
            Oout = 'numeric',
            Nout = 'numeric',
            Bout = 'numeric'))
        if (nrow(X2) < nrow(Y2)) {
            message('*Note: X has fewer data than Y. Recycling*')
            X2 <- t(matrix(rep(t(X2), ceiling(nrow(Y2) / nrow(X2))),
                nrow = ncol(X2),
                dimnames = list(colnames(X2), NULL))
            )[seq_len(nrow(Y2)), , drop = FALSE]
        }
        if (nrow(X2) > nrow(Y2)) {
            message('*Note: X has more data than Y. Recycling*')
            Y2 <- t(matrix(rep(t(Y2), ceiling(nrow(X2) / nrow(Y2))),
                nrow = ncol(Y2),
                dimnames = list(colnames(Y2), NULL))
            )[seq_len(nrow(X2)), , drop = FALSE]
            Y <- t(matrix(rep(t(Y), ceiling(nrow(X2) / nrow(Y))),
                nrow = ncol(Y),
                dimnames = list(colnames(Y), NULL))
            )[seq_len(nrow(X2)), , drop = FALSE]
        }
    } else {
        X2 <- sapply(seq_len(nrow(Y2)), function(x) NA)
    }
    ## ndata <- nrow(Y2)

    foreach(y = t(Y2), x = t(X2),
        .combine = combine,
        .inorder = TRUE) %dochains% {

#### the loop is over the columns of y and x
#### each instance is a 1-column vector
            if (all(is.na(x))) {
                probX <- log(mcoutput$W)
            } else {
                ## rows: components, cols: samples
                probX <- log(mcoutput$W) +
                    (if (XnR > 0) { # continuous
                        colSums(
                            dnorm(
                                x = x[XiR, ],
                                mean = mcoutput$Rmean[XtR, , , drop = FALSE],
                                sd = mcoutput$Rvar[XtR, , , drop = FALSE],
                                log = TRUE
                            ), na.rm = TRUE)
                    } else {
                        0
                    }) +
                    (if (XnC > 0) { # censored
                        isfin <- is.finite(x[XiC, ])
                        indf <- which(isfin)
                        indi <- which(!isfin)
                        (if (length(indf) > 0) {
                            colSums(
                                dnorm(
                                    x = x[XiC[indf], ],
                                    mean = mcoutput$Cmean[XtC[indf], , , drop = FALSE],
                                    sd = mcoutput$Cvar[XtC[indf], , , drop = FALSE],
                                    log = TRUE
                                ), na.rm = TRUE)
                        } else {
                            0
                        }) +
                            (if (length(indi) > 0) {
                                v2 <- XtC[indi]
                                v1 <- -sign(x[XiC[indi], ])
                                ## for upper tail, take opposite mean and value
                                colSums(
                                    pnorm(
                                        q = Cbounds[cbind(v2, 1.5 - 0.5 * v1)],
                                        mean = v1 * mcoutput$Cmean[v2, , , drop = FALSE],
                                        sd = mcoutput$Cvar[v2, , , drop = FALSE],
                                        log.p = TRUE
                                    ), na.rm = TRUE)
                            } else {
                                0
                            })
                    } else {
                        0
                    }) +
                    (if (XnD > 0) { # discretized
                        vrights <- pmax(x[XiD, ] + Dsteps[XtD], Dlefts[XtD])
                        vrights[vrights >= Drights[XtD]] <- +Inf
                        vlefts <- pmin(x[XiD, ] - Dsteps[XtD], Drights[XtD])
                        vlefts[vlefts <= Dlefts[XtD]] <- -Inf
                        colSums(log(
                            pnorm(
                                q = vrights,
                                mean = mcoutput$Dmean[XtD, , , drop = FALSE],
                                sd = mcoutput$Dvar[XtD, , , drop = FALSE]
                            ) -
                            pnorm(
                                q = vlefts,
                                mean = mcoutput$Dmean[XtD, , , drop = FALSE],
                                sd = mcoutput$Dvar[XtD, , , drop = FALSE]
                            )
                        ), na.rm = TRUE)
                    } else {
                        0
                    }) +
                    (if (XnO > 0) { # nominal
                        colSums(log(
                            aperm(
                                vapply(seq_len(XnO), function(v) {
                                    mcoutput$Oprob[XtO[v], , x[XiO[v], ], ]
                                }, mcoutput$W),
                                c(3, 1, 2))
                            ), na.rm = TRUE)
                    } else {
                        0
                    }) +
                    (if (XnN > 0) { # nominal
                        colSums(log(
                            aperm(
                                vapply(seq_len(XnN), function(v) {
                                    mcoutput$Nprob[XtN[v], , x[XiN[v], ], ]
                                }, mcoutput$W),
                                c(3, 1, 2))
                            ), na.rm = TRUE)
                    } else {
                        0
                    }) +
                    (if (XnB > 0) { # binary
                        colSums(log(
                            x[XiB, ] * mcoutput$Bprob[XtB, , , drop = FALSE] +
                            (1L - x[XiB, ]) *
                            (1 - mcoutput$Bprob[XtB, , , drop = FALSE])
                        ), na.rm = TRUE)
                    } else {
                        0
                    })
            } # end probX
            probX <- apply(probX, 2, function(xx) {
                xx - max(xx[is.finite(xx)])
            })
            ##
            ##
            if (all(is.na(y))) {
                probY <- array(NA, dim = dim(mcoutput$W))
            } else {
                probY <- (if (YnR > 0) { # continuous
                    colSums(
                        dnorm(
                            x = y[YiR, ],
                            mean = mcoutput$Rmean[YtR, , , drop = FALSE],
                            sd = mcoutput$Rvar[YtR, , , drop = FALSE],
                            log = TRUE
                        ), na.rm = TRUE)
                } else {
                    0
                }) +
                    (if (YnC > 0) { # censored
                        isfin <- is.finite(y[YiC, ])
                        indf <- which(isfin)
                        indi <- which(!isfin)
                        (if (length(indf) > 0) {
                            colSums(
                                dnorm(
                                    x = y[YiC[indf], ],
                                    mean = mcoutput$Cmean[YtC[indf], , , drop = FALSE],
                                    sd = mcoutput$Cvar[YtC[indf], , , drop = FALSE],
                                    log = TRUE
                                ), na.rm = TRUE)
                        } else {
                            0
                        }) +
                            (if (length(indi) > 0) {
                                v2 <- YtC[indi]
                                v1 <- -sign(y[YiC[indi], ])
                                ## for upper tail, take opposite mean and value
                                colSums(
                                    pnorm(
                                        q = Cbounds[cbind(v2, 1.5 - 0.5 * v1)],
                                        mean = v1 * mcoutput$Cmean[v2, , , drop = FALSE],
                                        sd = mcoutput$Cvar[v2, , , drop = FALSE],
                                        log.p = TRUE
                                    ), na.rm = TRUE)
                            } else {
                                0
                            })
                    } else {
                        0
                    }) +
                    (if (YnD > 0) { # discretized
                        vrights <- y[YiD, ] + Dsteps[YtD]
                        vrights[vrights >= Drights[YtD]] <- +Inf
                        vlefts <- y[YiD, ] - Dsteps[YtD]
                        vlefts[vlefts <= Dlefts[YtD]] <- -Inf
                        colSums(log(
                            pnorm(
                                q = vrights,
                                mean = mcoutput$Dmean[YtD, , , drop = FALSE],
                                sd = mcoutput$Dvar[YtD, , , drop = FALSE]
                            ) -
                            pnorm(
                                q = vlefts,
                                mean = mcoutput$Dmean[YtD, , , drop = FALSE],
                                sd = mcoutput$Dvar[YtD, , , drop = FALSE]
                            )
                        ), na.rm = TRUE)
                    } else {
                        0
                    }) +
                    (if (YnO > 0) { # nominal
                        colSums(log(
                            aperm(
                                vapply(seq_len(YnO), function(v) {
                                    mcoutput$Oprob[YtO[v], , y[YiO[v], ], ]
                                }, mcoutput$W),
                                c(3, 1, 2))
                            ), na.rm = TRUE)
                    } else {
                        0
                    }) +
                    (if (YnN > 0) { # nominal
                        colSums(log(
                            aperm(
                                vapply(seq_len(YnN), function(v) {
                                    mcoutput$Nprob[YtN[v], , y[YiN[v], ], ]
                                }, mcoutput$W),
                                c(3, 1, 2))
                            ), na.rm = TRUE)
                    } else {
                        0
                    }) +
                    (if (YnB > 0) { # binary
                        colSums(log(
                        (y[YiB, ] * mcoutput$Bprob[YtB, , , drop = FALSE]) +
                        ((1 - y[YiB, ]) *
                         (1 - mcoutput$Bprob[YtB, , , drop = FALSE]))
                        ), na.rm = TRUE)
                    } else {
                        0
                    })
            }
#### Output: rows=components, columns=samples
            ##
            ## if(all(is.na(x))){
            ##     out <- rowSums(exp(probX+probY))
            ## }else{
            ## str(probX)
            ## print(any(is.na(probX)))
            ## print(apply(probX, 1, max, na.rm=TRUE))
            ## str(probX)
            ## print(any(is.na(probX)))
            ## str(probY)
            ## print(any(is.na(probY)))
            ## }
            ## rbind(log(W),probX,probY)
            ## using if + NULL condition is faster than using fn=identity
            if(is.null(fn)) {
                colSums(exp(probX + probY)) / colSums(exp(probX))
            } else {
                fn(colSums(exp(probX + probY)) / colSums(exp(probX)))
            }
        } *
    (if (jacobian) {
        exp(-rowSums(
            log(vtransform(Y,
                auxmetadata = auxmetadata,
                invjacobian = TRUE)),
            na.rm = TRUE
        ))
    } else {
        1L
    })
}

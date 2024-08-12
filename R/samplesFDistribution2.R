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
#' @return The frequencies F(Y|X) corresponding to the Monte Carlo samples
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
                if(!silent) {
                    cat('\nClosing connections to cores.\n')
                }
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

    if (ncores < 2) {
        `%dochains%` <- `%do%`
    } else {
        `%dochains%` <- `%dopar%`
    }

    ## Extract Monte Carlo output & auxmetadata
    ## If mcoutput is a string, check if it's a folder name or file name
    if (is.character(mcoutput)) {
        ## Check if 'mcoutput' is a folder containing Fdistribution.rds
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
                stop('The argument "mcoutput" must be a folder containing Fdistribution.rds, or the path to an rds-file containing the output from "inferpopulation".')
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
        Clefts <- auxmetadata[match(vnames, auxmetadata$name), 'tleftbound']
        Crights <- auxmetadata[match(vnames, auxmetadata$name), 'trightbound']
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
                lprobX <- log(mcoutput$W)
            } else {
                ## rows: components, cols: samples
                lprobX <- log(mcoutput$W) +
                    util_lprob(
                        x = x,
                        mcoutput = mcoutput,
                        nR = XnR, iR = XiR, tR = XtR,
                        nC = XnC, iC = XiC, tC = XtC,
                        Clefts = Clefts, Crights = Crights,
                        nD = XnD, iD = XiD, tD = XtD,
                        Dsteps = Dsteps, Dlefts = Dlefts, Drights = Drights,
                        nO = XnO, iO = XiO, tO = XtO,
                        nN = XnN, iN = XiN, tN = XtN,
                        nB = XnB, iB = XiB, tB = XtB
                    )
            } # end lprobX
            ##
            ##
            if (all(is.na(y))) {
                lprobY <- array(NA, dim = dim(mcoutput$W))
            } else {
                lprobY <- util_lprob(
                        x = y,
                        mcoutput = mcoutput,
                        nR = YnR, iR = YiR, tR = YtR,
                        nC = YnC, iC = YiC, tC = YtC,
                        Clefts = Clefts, Crights = Crights,
                        nD = YnD, iD = YiD, tD = YtD,
                        Dsteps = Dsteps, Dlefts = Dlefts, Drights = Drights,
                        nO = YnO, iO = YiO, tO = YtO,
                        nN = YnN, iN = YiN, tN = YtN,
                        nB = YnB, iB = YiB, tB = YtB
                    )
            }
#### Output: rows=components, columns=samples

            lprobX <- apply(lprobX, 2, function(xx) {
                xx - max(xx[is.finite(xx)])
            })

            if(is.null(fn)) {
                colSums(exp(lprobX + lprobY)) / colSums(exp(lprobX))
            } else {
                fn(colSums(exp(lprobX + lprobY)) / colSums(exp(lprobX)))
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
        1
    })
}

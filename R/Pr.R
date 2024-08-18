#' Calculate posterior probabilities
#'
#' This function calculates the probability P(Y | X, data), where Y and X are two (non overlapping) sets of joint variates. The function also gives quantiles about the possible variability of the probability P(Y | X, newdata, data) that we could have if more learning data were provided, as well as a number of samples of the possible values of such probabilities. If several joint values are given for Y or X, the function will create a 2D grid of results for all possible compbinations of the given Y and X values.
#'
#' @param Y matrix or data.table: set of values of variates of which we want
#'   the joint probability of. One variate per column, one set of values per row.
#' @param X matrix or data.table or `NULL`: set of values of variates on which we want to condition the joint probability of `Y`. If `NULL` (default), no conditioning is made (except for conditioning on the learning dataset and prior assumptions). One variate per column, one set of values per row.
#' @param learned Either a string with the name of a directory or full
#'   path for an 'learned.rds' object, or such an object itself
#' @param quantiles numeric vector, between 0 and 1, or `NULL`: desired quantiles of the variability of the probability for `Y`. Default `c(0.05, 0.95)` or the 5% and 95% quantiles.
#' @param nsamples integer or `NULL`: desired number of samples of the variability of the probability for `Y`. Default `100`.
#' @param parallel logical or integer: whether to use pre-existing parallel
#'   workers, or how many to create and use. Default `TRUE`.
#' @param silent logical: give warnings or updates in the computation?
#'   Default `FALSE`.
#' @param usememory logical: save partial results to disc, to avoid crashes?
#'   Default `TRUE`.
#'
#' @return A list of: (1) a matrix with the probabilities P(Y|X,data,assumptions), for all combinations of values of `Y` (rows) and `X` (columns); (2) an array with the variability quantiles (3rd dimension of the array) for such probabilities; (3) an array with the variability samples (3rd dimension of the array) for such probabilities.
#'
#' @import parallel foreach doParallel
#'
#' @export
Pr <- function(
    Y,
    X = NULL,
    learned,
    quantiles = c(5, 95)/100,
    nsamples = 100L,
    parallel = TRUE,
    silent = FALSE,
    usememory = TRUE
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

    ## Extract Monte Carlo output & auxmetadata
    ## If learned is a string, check if it's a folder name or file name
    if (is.character(learned)) {
        ## Check if 'learned' is a folder containing learned.rds
        if (file_test('-d', learned) &&
                file.exists(file.path(learned, 'learned.rds'))) {
            learned <- readRDS(file.path(learned, 'learned.rds'))
        } else {
            ## Assume 'learned' the full path of learned.rds
            ## possibly without the file extension '.rds'
            learned <- paste0(sub('.rds$', '', learned), '.rds')
            if (file.exists(learned)) {
                learned <- readRDS(learned)
            } else {
                stop('The argument "learned" must be a folder containing learned.rds, or the path to an rds-file containing the output from "learn()".')
            }
        }
    }
    ## Add check to see that learned is correct type of object?
    auxmetadata <- learned$auxmetadata
    learned$auxmetadata <- NULL
    learned$auxinfo <- NULL
    ncomponents <- nrow(learned$W)
    nmcsamples <- ncol(learned$W)

    ## Consistency checks
    if (length(dim(Y)) != 2) {
        stop('Y must have two dimensions')
    }
    if (!is.null(X) && length(dim(X)) != 2) {
        stop('X must be NULL or have two dimensions')
    }
    ##
    if (!is.null(X) && ncol(X) == 0) {
        stop('X must be NULL or have two dimensions')
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

    nY <- nrow(Y)
    nX <- max(nrow(X), 1L)

    ## determine if parallel computation is possible and needed
    if (ncores < 2) {
        `%dox%` <- `%do%`
        `%doy%` <- `%do%`
    } else {
        if(nX >= 2 * ncores) {
            `%dox%` <- `%dopar%`
        } else {
            ## no point using more parallel cores than X values
            `%dox%` <- `%do%`
        }
        ##
        if(nY * nX >= 2 * ncores) {
            `%doy%` <- `%dopar%`
        } else {
            ## no point using more parallel cores than Y values
            `%doy%` <- `%do%`
        }
    }



    ## #### Subsample and get ncomponents and nsamples
    ##     ## source('mcsubset.R')
    ##     if (!missing(subsamples) &&
    ##             (is.numeric(subsamples) || (is.character(subsamples)
    ##                 && length(subsamples) == 1))) {
    ##         if (is.character(subsamples)) {
    ##             subsamples <- round(seq(1, ncol(learned$W),
    ##                 length.out = as.numeric(subsamples)
    ##             ))
    ##         }
    ##         learned <- mcsubset(learned, subsamples)
    ##     }


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
        learned$Rvar <- sqrt(learned$Rvar)
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
        learned$Cvar <- sqrt(learned$Cvar)
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
        learned$Dvar <- sqrt(learned$Dvar)
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



#### First calculate and save arrays for X values:
    if (is.null(X)) {
        lprobX <- log(learned$W)
        usememory <- FALSE
    } else {
        X2 <- as.matrix(vtransform(X, auxmetadata = auxmetadata,
            Rout = 'normalized',
            Cout = 'boundisinf',
            Dout = 'normalized',
            Oout = 'numeric',
            Nout = 'numeric',
            Bout = 'numeric'))

        temporarydir <- tempdir() # where to save X objects
        ##
        todelete <- foreach(jj = seq_len(nX), x = t(X2),
            .combine = `c`,
            .inorder = TRUE) %dox% {
                lprobX <- c(log(learned$W)) +
                    util_lprob(
                        x = x,
                        learned = learned,
                        nR = XnR, iR = XiR, tR = XtR,
                        nC = XnC, iC = XiC, tC = XtC,
                        Clefts = Clefts, Crights = Crights,
                        nD = XnD, iD = XiD, tD = XtD,
                        Dsteps = Dsteps, Dlefts = Dlefts, Drights = Drights,
                        nO = XnO, iO = XiO, tO = XtO,
                        nN = XnN, iN = XiN, tN = XtN,
                        nB = XnB, iB = XiB, tB = XtB
                    ) # rows=components, columns=samples

                lprobX <- apply(lprobX, 2, function(xx) {
                    xx - max(xx[is.finite(xx)])
                })

                saveRDS(lprobX,
                    file.path(temporarydir,
                        paste0('__X', jj, '__.rds'))
                )
                NULL
            }
    }

#### Now calculate for each Y value, combining with each X value
    ## transformation of inputs
    Y2 <- as.matrix(vtransform(Y, auxmetadata = auxmetadata,
        Rout = 'normalized',
        Cout = 'boundisinf',
        Dout = 'normalized',
        Oout = 'numeric',
        Nout = 'numeric',
        Bout = 'numeric'))

    ## jacobians <- exp(-rowSums(
    ##     log(vtransform(Y,
    ##         auxmetadata = auxmetadata,
    ##         invjacobian = TRUE)),
    ##     na.rm = TRUE
    ## ))

    sampleseq <- round(seq(1, nmcsamples, length.out = nsamples))

    keys <- c('values',
        if(!is.null(quantiles)){'quantiles'},
        if(!is.null(nsamples)) {'samples'}
    )
    ##
    combfnr <- function(...){setNames(do.call(mapply,
        c(FUN = `rbind`, lapply(list(...), `[`, keys))),
        keys)}
    ## combfnc <- function(...){setNames(do.call(mapply, c(FUN=cbind, lapply(list(...), `[`, keys))), keys)}

    out <- foreach(y = t(Y2),
        .combine = `combfnr`, .inorder = TRUE) %:%
        foreach(jj = seq_len(nX),
            .combine = `combfnr`, .inorder = TRUE) %doy% {
                if (all(is.na(y))) {
                    lprobY <- array(NA, dim = c(ncomponents, nmcsamples))
                } else {
                    lprobY <- util_lprob(
                        x = y,
                        learned = learned,
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

                if(usememory) {
                    lprobX <- readRDS(file.path(temporarydir,
                        paste0('__X', jj, '__.rds')
                    ))
                }

                FF <- colSums(exp(lprobX + lprobY)) / colSums(exp(lprobX))
                FF <- FF[!is.na(FF)]

                list(
                    values = rbind(mean(FF, na.rm = TRUE)),
                    ##
                    quantiles = rbind(if(!is.null(quantiles)) {
                        quantile(FF, probs = quantiles, na.rm = TRUE, type = 6)
                    }),
                    ##
                    samples = rbind(if(!is.null(nsamples)) {
                        FF[sampleseq]
                    })
                    ##
                    ## error = sd(FF, na.rm = TRUE)/sqrt(nmcsamples)
                )
            }
    ## in the output-list elements the Y & X values are the rows
    if(!is.null(nsamples)){
        colnames(out$samples) <- sampleseq
    }

    jacobians <- exp(-rowSums(
        log(as.matrix(vtransform(Y,
            auxmetadata = auxmetadata,
            invjacobian = TRUE))),
        na.rm = TRUE
    ))

    ## cat('\nnext\n')
    ## transform each element into a Y-X grid
    lapply(out, function(xx){
        temp <- colnames(xx)
        if(!is.null(temp)){
            dim(xx) <- c(nY, nX, ncol(xx))
            dimnames(xx) <- list( Y = NULL, X = NULL, temp)
        } else {
            dim(xx) <- c(nY, nX)
            dimnames(xx) <- list( Y = NULL, X = NULL)
        }
        jacobians * xx
    })
}

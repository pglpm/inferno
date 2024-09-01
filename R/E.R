#' Calculate posterior probabilities
#'
#' This function calculates the probability P(Y | X, data), where Y and X are two (non overlapping) sets of joint variates. The function also gives quantiles about the possible variability of the probability P(Y | X, newdata, data) that we could have if more learning data were provided, as well as a number of samples of the possible values of such probabilities. If several joint values are given for Y or X, the function will create a 2D grid of results for all possible compbinations of the given Y and X values.
#'
#' @param Y1names String vector: joint variates
#'   the joint probability of. One variate per column, one set of values per row.
#' @param X matrix or data.table or `NULL`: set of values of variates on which we want to condition the joint probability of `Y`. If `NULL` (default), no conditioning is made (except for conditioning on the learning dataset and prior assumptions). One variate per column, one set of values per row.
#' @param learnt Either a string with the name of a directory or full
#'   path for an 'learnt.rds' object, or such an object itself
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
E <- function(
    Ynames,
    X = NULL,
    learnt,
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
    workers <- setupParallel(parallel, silent)
    ncores <- workers$ncores
    
    if (!is.logical(workers$cluster)) {
        on.exit(closeCoresOnExit(workers$cluster, silent))
    }

    ## Extract Monte Carlo output & auxmetadata
    ## If learnt is a string, check if it's a folder name or file name
    if (is.character(learnt)) {
        ## Check if 'learnt' is a folder containing learnt.rds
        if (file_test('-d', learnt) &&
                file.exists(file.path(learnt, 'learnt.rds'))) {
            learnt <- readRDS(file.path(learnt, 'learnt.rds'))
        } else {
            ## Assume 'learnt' the full path of learnt.rds
            ## possibly without the file extension '.rds'
            learnt <- paste0(sub('.rds$', '', learnt), '.rds')
            if (file.exists(learnt)) {
                learnt <- readRDS(learnt)
            } else {
                stop('The argument "learnt" must be a folder containing learnt.rds, or the path to an rds-file containing the output from "learn()".')
            }
        }
    }
    ## Add check to see that learnt is correct type of object?
    auxmetadata <- learnt$auxmetadata
    learnt$auxmetadata <- NULL
    learnt$auxinfo <- NULL
    ncomponents <- nrow(learnt$W)
    nmcsamples <- ncol(learnt$W)

    ## Consistency checks
    if(!is.character(Ynames) || any(is.na(Ynames))){
        stop('Ynames must be a vector of variate names')
    }
    if (!is.null(X) && length(dim(X)) != 2) {
        stop('X must be NULL or have two dimensions')
    }
    ##
    if (!is.null(X) && ncol(X) == 0) {
        stop('X must be NULL or have two dimensions')
    }

    ## More consistency checks
    if(!all(Ynames %in% auxmetadata$name)) {
        stop('unknown Y variates\n')
    }
    if(length(unique(Ynames)) != length(Ynames)) {
        stop('duplicate Y variates\n')
    }
    ##
    Xnames <- colnames(X)
    if (!all(Xnames %in% auxmetadata$name)) {
        stop('unknown X variates\n')
    }
    if (length(unique(Xnames)) != length(Xnames)) {
        stop('duplicate X variates\n')
    }
    ##
    if (length(intersect(Ynames, Xnames)) > 0) {
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
    ##             subsamples <- round(seq(1, ncol(learnt$W),
    ##                 length.out = as.numeric(subsamples)
    ##             ))
    ##         }
    ##         learnt <- mcsubset(learnt, subsamples)
    ##     }


### Guide to indices:
    ## .i. = order in X/Y corresponding to appearance in vnames
    ## .t. = vnames present in X/Y, kept in their vnames-order

#### Type R
    vnames <- auxmetadata[auxmetadata$mcmctype == 'R', 'name']
    XiR <- match(vnames, Xnames)
    XtR <- which(!is.na(XiR))
    XiR <- XiR[XtR]
    XnR <- length(XiR)
    ##
    YiR <- match(vnames, Ynames)
    YtR <- which(!is.na(YiR))
    YiR <- YiR[YtR]
    YnR <- length(YiR)
    if (YnR > 0 || XnR > 0) {
        learnt$Rvar <- sqrt(learnt$Rvar)
    }

#### Type C
    vnames <- auxmetadata[auxmetadata$mcmctype == 'C', 'name']
    XiC <- match(vnames, Xnames)
    XtC <- which(!is.na(XiC))
    XiC <- XiC[XtC]
    XnC <- length(XiC)
    ##
    YiC <- match(vnames, Ynames)
    YtC <- which(!is.na(YiC))
    YiC <- YiC[YtC]
    YnC <- length(YiC)
    if (YnC > 0 || XnC > 0) {
        learnt$Cvar <- sqrt(learnt$Cvar)
        Clefts <- auxmetadata[match(vnames, auxmetadata$name), 'tdomainminplushs']
        Crights <- auxmetadata[match(vnames, auxmetadata$name), 'tdomainmaxminushs']
    }

#### Type D
    vnames <- auxmetadata[auxmetadata$mcmctype == 'D', 'name']
    XiD <- match(vnames, Xnames)
    XtD <- which(!is.na(XiD))
    XiD <- XiD[XtD]
    XnD <- length(XiD)
    ##
    YiD <- match(vnames, Ynames)
    YtD <- which(!is.na(YiD))
    YiD <- YiD[YtD]
    YnD <- length(YiD)
    if (YnD > 0 || XnD > 0) {
        learnt$Dvar <- sqrt(learnt$Dvar)
        Dsteps <- auxmetadata[match(vnames, auxmetadata$name), 'halfstep'] /
            auxmetadata[match(vnames, auxmetadata$name), 'tscale']
        Dlefts <- auxmetadata[match(vnames, auxmetadata$name), 'tdomainminplushs']
        Drights <- auxmetadata[match(vnames, auxmetadata$name), 'tdomainmaxminushs']
    }

#### Type O
    vnames <- auxmetadata[auxmetadata$mcmctype == 'O', 'name']
    XiO <- match(vnames, Xnames)
    XtO <- which(!is.na(XiO))
    XiO <- XiO[XtO]
    XnO <- length(XiO)
    ##
    YiO <- match(vnames, Ynames)
    YtO <- which(!is.na(YiO))
    YiO <- YiO[YtO]
    YnO <- length(YiO)

#### Type N
    vnames <- auxmetadata[auxmetadata$mcmctype == 'N', 'name']
    XiN <- match(vnames, Xnames)
    XtN <- which(!is.na(XiN))
    XiN <- XiN[XtN]
    XnN <- length(XiN)
    ##
    YiN <- match(vnames, Ynames)
    YtN <- which(!is.na(YiN))
    YiN <- YiN[YtN]
    YnN <- length(YiN)

#### Type B
    vnames <- auxmetadata[auxmetadata$mcmctype == 'B', 'name']
    XiB <- match(vnames, Xnames)
    XtB <- which(!is.na(XiB))
    XiB <- XiB[XtB]
    XnB <- length(XiB)
    ##
    YiB <- match(vnames, Ynames)
    YtB <- which(!is.na(YiB))
    YiB <- YiB[YtB]
    YnB <- length(YiB)



#### First calculate and save arrays for X values:
    if (is.null(X)) {
        lprobX <- log(learnt$W)
        usememory <- FALSE
    } else {
        X2 <- as.matrix(vtransform(X, auxmetadata = auxmetadata,
            Rout = 'normalized',
            Cout = 'boundisinf',
            Dout = 'normalized',
            Oout = 'numeric',
            Nout = 'numeric',
            Bout = 'numeric',
            logjacobianOr = NULL))

        temporarydir <- tempdir() # where to save X objects
        ##
        todelete <- foreach(jj = seq_len(nX), x = t(X2),
            .combine = `c`,
            .inorder = TRUE) %dox% {
                lprobX <- c(log(learnt$W)) +
                    util_lprob(
                        x = x,
                        learnt = learnt,
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
        Bout = 'numeric',
        logjacobianOr = NULL))

    ## jacobians <- exp(-rowSums(
    ##     log(vtransform(Y,
    ##         auxmetadata = auxmetadata,
    ##         invjacobian = TRUE)),
    ##     na.rm = TRUE
    ## ))

    if(!is.null(nsamples)){
        sampleseq <- round(seq(1, nmcsamples, length.out = nsamples))
    }

    keys <- c('values',
        if(!is.null(quantiles)){'quantiles'},
        if(!is.null(nsamples)) {'samples'}
    )
    ##
    combfnr <- function(...){setNames(do.call(mapply,
        c(FUN = `rbind`, lapply(list(...), `[`, keys))),
        keys)}
    ## combfnc <- function(...){setNames(do.call(mapply, c(FUN=cbind, lapply(list(...), `[`, keys))), keys)}

    out <- foreach(jj = seq_len(nX),
        .combine = `combfnr`, .inorder = TRUE) %:%
        foreach(y = t(Y2),
            .combine = `combfnr`, .inorder = TRUE) %doy% {
                if (all(is.na(y))) {
                    lprobY <- array(NA, dim = c(ncomponents, nmcsamples))
                } else {
                    lprobY <- util_lprob(
                        x = y,
                        learnt = learnt,
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
                    values = mean(FF, na.rm = TRUE),
                    ##
                    quantiles = (if(!is.null(quantiles)) {
                        quantile(FF, probs = quantiles, type = 6,
                            na.rm = TRUE, names = FALSE)
                    }),
                    ##
                    samples = (if(!is.null(nsamples)) {
                        FF[sampleseq]
                    })
                    ##
                    ## error = sd(FF, na.rm = TRUE)/sqrt(nmcsamples)
                )
            } # End foreach over Y and X

    jacobians <- exp(rowSums(
        as.matrix(vtransform(Y,
            auxmetadata = auxmetadata,
            logjacobianOr = TRUE)),
        na.rm = TRUE
    ))

    ## multiply by jacobian factors
    out$values <- out$values * jacobians

    ## transform to grid
    ## in the output-list elements the Y & X values are the rows
    dim(out$values) <- c(nY, nX)
    dimnames(out$values) <- list(Y = NULL, X = NULL)

    if(!is.null(quantiles)){
        out$quantiles <- out$quantiles * jacobians
        temp <- names(quantile(1, probs = quantiles, names = TRUE))
        ## transform to grid
        dim(out$quantiles) <- c(nY, nX, length(quantiles))
        dimnames(out$quantiles) <- list(Y = NULL, X = NULL, temp)
    }

    if(!is.null(nsamples)){
    ## transform to grid
        out$samples <- out$samples * jacobians
        dim(out$samples) <- c(nY, nX, nsamples)
        dimnames(out$samples) <- list(Y = NULL, X = NULL, sampleseq)
    }

    out
}

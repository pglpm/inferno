#' (Obsolete) Calculate cumulative posterior probabilities
#'
#' This function calculates the probability `P(Y <= y | X, data)`, where Y and X are two (non overlapping) sets of joint variates (if the `lower.tail` argument is `FALSE`, then `P(Y > y | X, data)` is calculated). The function also gives quantiles about the possible variability of the probability `P(Y <= y | X, newdata, data)` that we could have if more learning data were provided, as well as a number of samples of the possible values of such probabilities. If several joint values are given for Y or X, the function will create a 2D grid of results for all possible compbinations of the given Y and X values.
#'
#' @param Y Matrix or data.table: set of values of variates of which we want
#'   the joint probability of. One variate per column, one set of values per row.
#' @param X Matrix or data.table or `NULL`: set of values of variates on which we want to condition the joint probability of `Y`. If `NULL` (default), no conditioning is made (except for conditioning on the learning dataset and prior assumptions). One variate per column, one set of values per row.
#' @param learnt Either a character with the name of a directory or full path for a 'learnt.rds' object, produced by the [learn()] function, or such an object itself.
#' @param nsamples Integer or `NULL` or `"all"`: desired number of samples of the variability of the probability for `Y`. If `NULL`, no samples are reported. If `"all"` (or `Inf`), all samples obtained by the [learn()] function are used. Default `"all"`.
#' @param quantiles Numeric vector, between 0 and 1, or `NULL`: desired quantiles of the variability of the probability for `Y`. Default `c(0.055, 0.25, 0.75, 0.945)`, that is, the 5.5%, 25%, 75%, 94.5% quantiles (these are typical quantile values in the Bayesian literature: they give 50% and 89% credibility intervals, which correspond to 1 shannons and 0.5 shannons of uncertainty). If `NULL`, no quantiles are calculated.
#' @param parallel Logical or `NULL` or positive integer: `TRUE`: use roughly half of available cores; `FALSE`: use serial computation; `NULL`: don't do anything (use pre-registered condition); integer: use this many cores. Default `NULL`
#' @param eq Logical: include `Y = y` in the cumulative probability? Default `TRUE`.
#' @param lower.tail Logical: calculate `P(Y < y)`? (`TRUE`, default) Or `P(Y > y)`? (`FALSE`).
#' @param silent Logical: give warnings or updates in the computation?
#'   Default `FALSE`.
#' @param usememory Logical: save partial results to disc, to avoid crashes?
#'   Default `TRUE`.
#' @param keepYX Logical, default `TRUE`: keep a copy of the `Y` and `X` arguments in the output? This is used for the plot method.
#'
#' @return A list of class `probability`, consisting of the elements `values`,  `quantiles` (possibly `NULL`), `samples` (possibly `NULL`), `Y`, `X`. Element `values`: a matrix with the probabilities P(Y|X,data,assumptions), for all combinations of values of `Y` (rows) and `X` (columns). Element `quantiles`: an array with the variability quantiles (3rd dimension of the array) for such probabilities. Element `samples`: an array with the variability samples (3rd dimension of the array) for such probabilities. Elements `Y`, `X`: copies of the `Y` and `X` arguments.
#'
#' @import parallel foreach doParallel
#'
#' @export
oldtailPr <- function(
    Y,
    X = NULL,
    learnt,
    nsamples = 'all',
    quantiles = c(0.055, 0.25, 0.75, 0.945),
    parallel = NULL,
    eq = TRUE,
    lower.tail = TRUE,
    silent = FALSE,
    usememory = TRUE,
    keepYX = TRUE
) {
    if (!silent) {
        cat('\n')
    }

#### Requested parallel processing
    ## NB: doesn't make sense to have more cores than chains
    closeexit <- FALSE
    if (isTRUE(parallel)) {
        ## user wants us to register a parallel backend
        ## and to choose number of cores
        ncores <- max(1,
            floor(parallel::detectCores() / 2))
        cl <- parallel::makeCluster(ncores)
        doParallel::registerDoParallel(cl)
        closeexit <- TRUE
        if(!silent){cat('Registered', foreach::getDoParName(),
            'with', foreach::getDoParWorkers(), 'workers\n')}
    } else if (isFALSE(parallel)) {
        ## user wants us not to use parallel cores
        ncores <- 1
        foreach::registerDoSEQ()
    } else if (is.null(parallel)) {
        ## user wants us not to do anything
        ncores <- foreach::getDoParWorkers()
    } else if (is.finite(parallel) && parallel >= 1) {
        ## user wants us to register 'parallal' # of cores
        ncores <- parallel
        cl <- parallel::makeCluster(ncores)
        doParallel::registerDoParallel(cl)
        closeexit <- TRUE
        if(!silent){cat('Registered', foreach::getDoParName(),
            'with', foreach::getDoParWorkers(), 'workers\n')}
    } else {
        stop("Unknown value of argument 'parallel'")
    }

    ## Close parallel connections if any were opened
    if(closeexit) {
        closecoresonexit <- function(){
            if(!silent){cat('\nClosing connections to cores.\n')}
            foreach::registerDoSEQ()
            parallel::stopCluster(cl)
            ## parallel::setDefaultCluster(NULL)
            env <- foreach:::.foreachGlobals
            rm(list=ls(name=env), pos=env)
        }
        on.exit(closecoresonexit())
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

    if(is.numeric(nsamples)){
        if(is.na(nsamples) || nsamples < 1) {
            nsamples <- NULL
        } else if(!is.finite(nsamples)) {
            nsamples <- nmcsamples
        }
    } else if (is.character(nsamples) && nsamples == 'all'){
        nsamples <- nmcsamples
    }

    Y <- as.data.frame(Y)
    if(!is.null(X)){X <- as.data.frame(X)}

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

    doquantiles <- !is.null(quantiles)
    dosamples <- !is.null(nsamples)

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
        learnt$Rvar <- sqrt(learnt$Rvar)
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
        learnt$Cvar <- sqrt(learnt$Cvar)
        Clefts <- auxmetadata[match(vnames, auxmetadata$name), 'tdomainmin']
        Crights <- auxmetadata[match(vnames, auxmetadata$name), 'tdomainmax']
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
        learnt$Dvar <- sqrt(learnt$Dvar)
        Dsteps <- auxmetadata[match(vnames, auxmetadata$name), 'halfstep'] /
            auxmetadata[match(vnames, auxmetadata$name), 'tscale']
        Dlefts <- auxmetadata[match(vnames, auxmetadata$name), 'tdomainmin']
        Drights <- auxmetadata[match(vnames, auxmetadata$name), 'tdomainmax']
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
    if(any(Yv %in% vnames)) {
        stop('It does not make sense to ask for the cumulative probability of a nominal variate.')
    }

#### Type B
    vnames <- auxmetadata[auxmetadata$mcmctype == 'B', 'name']
    XiB <- match(vnames, Xv)
    XtB <- which(!is.na(XiB))
    XiB <- XiB[XtB]
    XnB <- length(XiB)
    ##
    if(any(Yv %in% vnames)) {
        stop('It does not make sense to ask for the cumulative probability of a two-valued variate.')
    }


#### First calculate and save arrays for X values:
    if (is.null(X)) {
        lprobX <- log(learnt$W)
        usememory <- FALSE
    } else {
        X2 <- as.matrix(vtransform(X, auxmetadata = auxmetadata,
            Rout = 'normalized',
            Cout = 'boundisinf',
            Dout = 'boundnormalized',
            Oout = 'index',
            Nout = 'index',
            Bout = 'numeric',
            logjacobianOr = NULL))

        ## create unique dir where to save X objects
        temporarydir <- tempdir()
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

                ## ## seems to lead to garbage for extreme values
                ## lprobX <- apply(lprobX, 2, function(xx) {
                ##     xx - max(xx[is.finite(xx)])
                ## })

                saveRDS(lprobX,
                    file.path(temporarydir,
                        paste0('__X', jj, '__.rds'))
                )
                NULL
            }
    }

#### Now calculate for each Y value, combining with each X value
    ## transformation of inputs
    if(eq){Cout <- 'boundisinf'}else{Cout <- 'boundnormalized'}
    Y2 <- as.matrix(vtransform(Y, auxmetadata = auxmetadata,
        Rout = 'normalized',
        Cout = Cout,
        Dout = 'boundnormalized',
        Oout = 'index',
        Nout = 'index',
        Bout = 'numeric',
        logjacobianOr = NULL))

    if((eq && lower.tail) || (!eq && !lower.tail)){
        if(YnD > 0){
            Y2[ , YiD] <- Y2[ , YiD] + Dsteps[YtD]
        }
    } else {
        if(YnO > 0){
            Y2[ , YiO] <- Y2[ , YiO] - 1
        }
        if(YnD > 0){
            Y2[ , YiD] <- Y2[ , YiD] - Dsteps[YtD]
        }
    }
    ## jacobians <- exp(-rowSums(
    ##     log(vtransform(Y,
    ##         auxmetadata = auxmetadata,
    ##         invjacobian = TRUE)),
    ##     na.rm = TRUE
    ## ))

    keys <- c('values', 'samples', 'quantiles')
    ##
    combfnr <- function(...){setNames(do.call(mapply,
        c(FUN = `rbind`, lapply(list(...), `[`, keys, drop = FALSE))),
        keys)}
    ## combfnc <- function(...){setNames(do.call(mapply, c(FUN=cbind, lapply(list(...), `[`, keys))), keys)}

    out <- foreach(jj = seq_len(nX),
        .combine = `combfnr`, .inorder = TRUE) %:%
        foreach(y = t(Y2),
            .combine = `combfnr`, .inorder = TRUE) %doy% {
                if (all(is.na(y))) {
                    lprobY <- array(NA, dim = c(ncomponents, nmcsamples))
                } else {
                    lprobY <- util_ltailprob(
                        x = y,
                        learnt = learnt,
                        nR = YnR, iR = YiR, tR = YtR,
                        nC = YnC, iC = YiC, tC = YtC,
                        Clefts = Clefts, Crights = Crights,
                        nD = YnD, iD = YiD, tD = YtD,
                        Dsteps = Dsteps, Dlefts = Dlefts, Drights = Drights,
                        nO = YnO, iO = YiO, tO = YtO,
                        ## nN = YnN, iN = YiN, tN = YtN,
                        ## nB = YnB, iB = YiB, tB = YtB,
                        lower.tail = lower.tail
                    )
                }

                if(usememory) {
                    lprobX <- readRDS(file.path(temporarydir,
                        paste0('__X', jj, '__.rds')
                    ))
                }

                FF <- colSums(exp(lprobX + lprobY), na.rm = TRUE) /
                    colSums(exp(lprobX), na.rm = TRUE)

                list(
                    values = mean(FF, na.rm = TRUE),
                    ##
                    samples = (if(dosamples) {
                        FF <- FF[!is.na(FF)]
                        FF[round(seq(1, length(FF), length.out = nsamples))]
                    }),
                    ##
                    quantiles = (if(doquantiles) {
                        quantile(FF, probs = quantiles, type = 6,
                            na.rm = TRUE, names = FALSE)
                    })
                    ##
                    ## error = sd(FF, na.rm = TRUE)/sqrt(nmcsamples)
                )
            } # End foreach over Y and X

    ## transform to grid
    ## in the output-list elements the Y & X values are the rows

    dim(out$values) <- c(nY, nX)

    ## if(ncol(Y) == 1){Ynames <- Y[, 1]} else {Ynames <- NULL}
    Ynames <- apply(Y, 1, paste0, collapse=',')
    ## if(isTRUE(ncol(X) == 1)){Xnames <- X[, 1]} else {Xnames <- NULL}
    if(!is.null(X)){
        Xnames <- apply(X, 1, paste0, collapse=',')
    } else {
        Xnames <- NULL
    }
    dimnames(out$values) <- list(Y = Ynames, X = Xnames)

    if(dosamples){
    ## transform to grid
        dim(out$samples) <- c(nY, nX, nsamples)
        dimnames(out$samples) <- list(Y = Ynames, X = Xnames,
            round(seq(1, nmcsamples, length.out = nsamples)))
    }

    if(doquantiles){
        out$quantiles <- out$quantiles
        temp <- names(quantile(1, probs = quantiles, names = TRUE))
        ## transform to grid
        dim(out$quantiles) <- c(nY, nX, length(quantiles))
        dimnames(out$quantiles) <- list(Y = Ynames, X = Xnames, temp)
    }

    if(isTRUE(keepYX)){
    ## save Y and X values in the output; useful for plotting methods
        out$Y <- Y
        out$X <- X
    }

    out$lowertail <- lower.tail

    class(out) <- 'probability'
    out
}

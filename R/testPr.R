#' Calculate posterior probabilities
#'
#' This function calculates the posterior probability `Pr(Y | X, data)`, where `Y` and `X` are two (non overlapping) sets of joint variate values. If `X` is omitted or `NULL`, then the posterior probability `Pr(Y | data)` is calculated. The function also gives quantiles about the possible variability of the probability `Pr(Y | X, newdata, data)` that we could have if more learning data were provided, as well as a number of samples of the possible values of such probabilities. If several joint values are given for `Y` or `X`, the function will create a 2D grid of results for all possible compbinations of the given `Y` and `X` values. This function also allows for base-rate or other prior-probability corrections: If a prior (for instance, a base rate) for `Y` is given, the function will calculate the `Pr(Y | X, data, prior)` from `Pr(X | Y, data)` and the prior by means of Bayes's theorem. Each variate in each argument `Y`, `X` can be specified either as a point-value `Y = y` or as a left-open interval `Y ≤ y` or as a right-open interval `Y ≥ y`, through the argument `tails`.
#'
#' @param Y Matrix or data.table: set of values of variates of which we want
#'   the joint probability of. One variate per column, one set of values per row.
#' @param X Matrix or data.table or `NULL` (default): set of values of variates on which we want to condition the joint probability of `Y`. If `NULL`, no conditioning is made (except for conditioning on the learning dataset and prior assumptions). One variate per column, one set of values per row.
#' @param learnt Either a character with the name of a directory or full path for a 'learnt.rds' object, produced by the [learn()] function, or such an object itself.
#' @param tails Named vector or list, or `NULL` (default). The names must match some or all of the variates in arguments `Y` and `X`. For variates in this list, the probability arguments are understood in an semi-open interval sense: `Y ≤ y` or `Y ≥ y`, an so on. This is true for variates on the left and on the right of the conditional sign `\|`. A left-open interval `Y ≤ y` is indicated by the values `'<='` or `'left'` or `-1`; a right-open interval `Y ≥ y` is indicated by the values `'>='` or `'right'` or `+1`. Values `NULL`, `'=='`, `0` indicate that a point value `Y = y` (not an interval) should be calculated. **NB**: the semi-open intervals *always* include the given value; this is important for ordinal or rounded variates. For instance, if `Y` is an integer variate, then to calculate  `P(Y < 3)` you should require `P(Y <= 2)`; for this reason we also have that `P(Y <= 2)` and  `P(Y >= 2)` generally add up to *more* than 1.
#' @param priorY Numeric vector with the same length as the rows of `Y`, or `TRUE`, or `NULL` (default): prior probabilities or base rates for the `Y` values. If `TRUE`, the prior probabilities are assumed to be all equal.
#' @param nsamples Integer or `NULL` or `'all'` (default): desired number of samples of the variability of the probability for `Y`. If `NULL`, no samples are reported. If `'all'` (or `Inf`), all samples obtained by the [learn()] function are used.
#' @param quantiles Numeric vector, between 0 and 1, or `NULL`: desired quantiles of the variability of the probability for `Y`. Default `c(0.055, 0.25, 0.75, 0.945)`, that is, the 5.5%, 25%, 75%, 94.5% quantiles (these are typical quantile values in the Bayesian literature: they give 50% and 89% credibility intervals, which correspond to 1 shannons and 0.5 shannons of uncertainty). If `NULL`, no quantiles are calculated.
#' @param parallel Logical or `NULL` (default) or positive integer: `TRUE`: use roughly half of available cores; `FALSE`: use serial computation; `NULL`: don't do anything (use pre-registered condition); integer: use this many cores.
#' @param silent Logical, default `FALSE`: give warnings or updates in the computation?
#' @param usememory Logical, default `TRUE`: save partial results to disc, to avoid crashes?
#' @param keepYX Logical, default `TRUE`: keep a copy of the `Y` and `X` arguments in the output? This is used for the plot method.
#'
#' @return A list of class `probability`, consisting of the elements `values`,  `quantiles` (possibly `NULL`), `samples` (possibly `NULL`), `values.MCaccuracy`, `quantiles.MCaccuracy` (possibly `NULL`), `Y`, `X`. Element `values`: a matrix with the probabilities P(Y|X,data,assumptions), for all combinations of values of `Y` (rows) and `X` (columns). Element `quantiles`: an array with the variability quantiles (3rd dimension of the array) for such probabilities. Element `samples`: an array with the variability samples (3rd dimension of the array) for such probabilities. Elements `values.MCaccuracy` and `quantiles.MCaccuracy`: arrays with the numerical accuracies (roughly speaking a standard deviation) of the Monte Carlo calculations for the `values` and `quantiles` elements. Elements `Y`, `X`: copies of the `Y` and `X` arguments.
#'
#' @import parallel foreach doParallel
#'
#' @export
testPr <- function(
    Y,
    X = NULL,
    learnt,
    tails = NULL,
    priorY = NULL,
    nsamples = 'all',
    quantiles = c(0.055, 0.25, 0.75, 0.945),
    parallel = NULL,
    silent = FALSE,
    usememory = TRUE,
    keepYX = TRUE
) {
    Qerror <- pnorm(c(-1, 1))

    ## if (!silent) {
    ##     cat('\n')
    ## }

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
    Yv <- colnames(Y)

    if(all(is.na(X))){X <- NULL}
    if(!is.null(X)){X <- as.data.frame(X)}
    Xv <- colnames(X)

    if(!is.null(tails)){
        tails <- as.list(tails)
        if(is.null(names(tails))) {
            stop('Missing variate names in "tails"')
        }
    }
    tailsv <- names(tails)
    tailscentre <- list('==', 0, '0', NULL)
    tailsleft <- list('<=', -1, '-1', 'left')
    tailsright <- list('>=', 1, '+1', 'right')
    tailsvalues <- c(tailscentre, tailsleft, tailsright)

    ## Consistency checks

    if (!all(Yv %in% auxmetadata$name)) {
        stop('unknown Y variate ',
            paste0(Yv[!(Yv %in% auxmetadata$name)], collapse = ' '),
            '\n')
    }
    if (length(unique(Yv)) != length(Yv)) {
        stop('duplicate Y variates\n')
    }

    if (!all(Xv %in% auxmetadata$name)) {
        stop('unknown X variate ',
            paste0(Xv[!(Xv %in% auxmetadata$name)], collapse = ' '),
            '\n')
    }
    if (length(unique(Xv)) != length(Xv)) {
        stop('duplicate X variates\n')
    }

    if (any(Yv %in% Xv)) {
        stop('overlap in Y and X variates\n')
    }

    if (!all(tailsv %in% c(Yv, Xv))) {
        warning('variate ',
            paste0(tailsv[!(tailsv %in% c(Yv, Xv))], collapse = ' '),
            ' not among Y and X; ignored\n')
    }
    if (length(unique(tailsv)) != length(tailsv)) {
        stop('duplicate "tails" variates\n')
    }
    if(!all(tails %in% tailsvalues)) {
        stop('"tails" values must be ',
            paste0(tailsvalues, collapse = ' '), '\n')
    }

    ## transform 'tails' to -1, +1
    ## +1: '<=',    -1: '>='
    ## this is opposite of the argument convention because
    ## interval probabilities are calculated with `lower.tail = TRUE`
    ## eg:
    ## pnorm(x, mean, sd, lower.tail = FALSE) ==
    ##     pnorm(-x, -mean, sd, lower.tail = TRUE)
    tails[tails %in% tailscentre] <- NULL
    cleft <- tails %in% tailsleft
    cright <- tails %in% tailsright
    tails[cleft] <- +1
    tails[cright] <- -1
    tails <- unlist(tails)

    ## Check if a prior for Y is given, in that case Y and X will be swapped
    if (isFALSE(priorY) || is.null(priorY)) {
        priorY <- NULL
        doquantiles <- !is.null(quantiles)
    } else {
        if(is.null(X)){ stop('X must be non-null if priorY is given') }

        if(!isTRUE(priorY) && length(priorY) != nrow(Y)) {
            stop('priorY must have as many elements a the rows of Y')
        }

        doquantiles <- FALSE
        if(!is.null(quantiles)) {
            ## message('For the moment, ',
            ##     'computation of quantiles is affected by a larger error',
            ##     'if "priorY" is specified.')
            nsamples0 <- nsamples
            nsamples <- max(nsamples0, 1200L)
        }

        ## if priorY is TRUE, the user wants a uniform prior distribution
        if(isTRUE(priorY)){
            priorY <- 1 + numeric(nrow(Y))
        } else {
            ## Check for invalid probabilities
            if (!is.numeric(priorY) || any(priorY < 0)) {
                stop('priorY contains invalid probabilities')
            }
        }

        ## Swap X and Y, to use Bayes's theorem
        . <- Y
        Y <- X
        X <- .
        rm(.)
        . <- Yv
        Yv <- Xv
        Xv <- .
        rm(.)
    }

    nY <- nrow(Y)
    nX <- max(nrow(X), 1L)

    ## determine if parallel computation is possible and needed
    if (ncores < 2) {
        `%doyx%` <- `%do%`
    } else {
        `%doyx%` <- `%dopar%`
    }
    dosamples <- !is.null(nsamples)


#### tmp dir where to save X and Y objects
    temporarydir <- tempdir()


#### Calculate and save arrays for X values:
    if (is.null(X)) {
        lprobX <- log(learnt$W)
        usememory <- FALSE
    } else {
        ## Construction of the arguments for util_lprobs, X argument
        lpargs <- util_lprobsargsyx(
            x = X,
            auxmetadata = auxmetadata,
            learnt = learnt,
            tails = tails
        )

        invisible(parLapply(cl = cl,
            X = lpargs$xVs,
            fun = util_lprobsave,
            params = lpargs$params,
            logW = c(log(learnt$W)),
            temporarydir = temporarydir,
            lab = '__X'
        ))
    }


#### Calculate and save arrays for Y values:

    ## Construction of the arguments for util_lprobs, Y argument
    ## jacobians <- exp(-rowSums(
    ##     log(vtransform(Y,
    ##         auxmetadata = auxmetadata,
    ##         invjacobian = TRUE)),
    ##     na.rm = TRUE
    ## ))

    lpargs <- util_lprobsargs(
        x = Y,
        auxmetadata = auxmetadata,
        learnt = learnt,
        tails = tails
    )

    invisible(parLapply(cl = cl,
        X = lpargs$xVs,
        fun = util_lprobsave,
        params = lpargs$params,
        logW = 0,
        temporarydir = temporarydir,
        lab = '__Y'
    ))

    keys <- c('values', 'quantiles', 'samples', 'values.MCaccuracy', 'quantiles.MCaccuracy')
    ##
    combfnr <- function(...){setNames(do.call(mapply,
        c(FUN = `rbind`, lapply(X = list(...), FUN = `[`, keys, drop = FALSE))),
        keys)}
    ## combfnc <- function(...){setNames(do.call(mapply, c(FUN=cbind, lapply(list(...), `[`, keys))), keys)}

    out <- foreach(
        jx = seq_len(nX),
        .combine = `combfnr`,
        .inorder = TRUE
    ) %:% foreach(
        jy = seq_len(nY),
        .combine = `combfnr`,
        .inorder = TRUE
    ) %doyx% {

        if(usememory) {
            lprobX <- readRDS(file.path(temporarydir,
                paste0('__X', jx, '__.rds')
            ))
            lprobY <- readRDS(file.path(temporarydir,
                paste0('__Y', jy, '__.rds')
            ))
        }

        FF <- colSums(x = exp(lprobX + lprobY), na.rm = TRUE) /
            colSums(x = exp(lprobX), na.rm = TRUE)

        list(
            values = mean(x = FF, na.rm = TRUE),
            ##
            quantiles = if(doquantiles) {
                quantile(x = FF, probs = quantiles, type = 6,
                    na.rm = TRUE, names = FALSE)
            },
            ##
            samples = if(dosamples) {
                FF <- FF[!is.na(FF)]
                FF[round(seq(1, length(FF), length.out = nsamples))]
            },
            ##
            values.MCaccuracy = funMCSELD(x = FF),
            ##
            quantiles.MCaccuracy = if(doquantiles) {
                temp <- funMCEQ(x = FF, prob = quantiles, Qpair = Qerror)
                (temp[2, ] - temp[1, ]) / 2
            }
            ##
            ## error = sd(FF, na.rm = TRUE)/sqrt(nmcsamples)
        )
    } # End foreach over Y and X

    ## clean temp files
    if(usememory) {
        unlink(x = c(
            sapply(seq_len(nX), function(jx){
                file.path(temporarydir, paste0('__X', jx, '__.rds'))
            }),
            sapply(seq_len(nY), function(jy){
                file.path(temporarydir, paste0('__Y', jy, '__.rds'))
            })
        ))
    }


    if(is.null(priorY)){
        y <- Y
        y[, colnames(Y) %in% tailsv] <- NA
        jacobians <- exp(rowSums(
            as.matrix(vtransform(y,
                auxmetadata = auxmetadata,
                logjacobianOr = TRUE)),
            na.rm = TRUE
        ))
        rm(y)
    }

    ## transform to grid
    ## in the output-list elements the Y & X values are the rows
    dim(out$values.MCaccuracy) <- dim(out$values) <- c(nY, nX)

    if(is.null(priorY)){
        ## multiply by jacobian factors
        out$values <- out$values * jacobians
        out$values.MCaccuracy <- signif(x = out$values.MCaccuracy * jacobians, digits = 2)

        ## if(ncol(Y) == 1){Ynames <- Y[, 1]} else {Ynames <- NULL}
        Ynames <- apply(X = Y, MARGIN = 1, FUN = paste0, collapse=',',
            simplify = TRUE)

        if(!is.null(X)){
            Xnames <- apply(X = X, MARGIN = 1, FUN = paste0, collapse=',',
                simplify = TRUE)
        } else {
            Xnames <- NULL
        }
        dimnames(out$values) <- list(Y = Ynames, X = Xnames)
        dimnames(out$values.MCaccuracy) <- dimnames(out$values)
    } else {
        ## Bayes's theorem
        out$values <- t(priorY * t(out$values))
        out$values.MCaccuracy <- NULL
        normf <- rowSums(out$values, na.rm = TRUE)
        out$values <- t(out$values/normf)

        ## if(ncol(X) == 1){Ynames <- X[, 1]} else {Ynames <- NULL}
        Ynames <- apply(X = X, MARGIN = 1, FUN = paste0, collapse=',',
            simplify = TRUE)

        if(!is.null(Y)){
            Xnames <- apply(X = Y, MARGIN = 1, FUN = paste0, collapse=',',
                simplify = TRUE)
        } else {
            Xnames <- NULL
        }
        dimnames(out$values) <- list(Y = Ynames, X = Xnames)
    }

    if(dosamples){
        ## transform to grid
        dim(out$samples) <- c(nY, nX, nsamples)

        if(is.null(priorY)){
            ## multiply by jacobian factors
            out$samples <- out$samples * jacobians
        } else {
            ## Bayes's theorem
            out$samples <- priorY * aperm(a = out$samples, perm = c(2, 1, 3),
                resize = FALSE)
            normf <- c(t(colSums(out$samples, na.rm = TRUE)))
            out$samples <- aperm(a = aperm(a = out$samples, perm = NULL,
                resize = FALSE) / normf, perm = NULL,
                resize = FALSE)
        }

        dimnames(out$samples) <- list(Y = Ynames, X = Xnames,
            round(seq(1, nmcsamples, length.out = nsamples)))
    }

    if(doquantiles){
        if(is.null(priorY)){
            ## transform to grid
            dim(out$quantiles) <- c(nY, nX, length(quantiles))
            dim(out$quantiles.MCaccuracy) <- c(nY, nX, length(quantiles))
            ## multiply by jacobian factors
            out$quantiles <- out$quantiles * jacobians
            out$quantiles.MCaccuracy <- signif(x = out$quantiles.MCaccuracy * jacobians,
                digits = 2)
        } else {
            ## calculate quantiles from samples
            out$quantiles <- aperm(
                a = apply(X = out$samples, MARGIN = c(1, 2), FUN = quantile,
                    probs = quantiles, type = 6,
                    na.rm = TRUE, names = FALSE,
                    simplify = TRUE),
                perm = c(2, 3, 1), resize = FALSE)

            ## adjust number of samples as originally requested
            if(is.null(nsamples0)) {
                out$samples <- NULL
            } else if(nsamples0 < nsamples) {
                out$samples <-out$samples[ , ,
                    round(seq(1, nsamples, length.out = nsamples0))]
            }
        }

        temp <- names(quantile(1, probs = quantiles, names = TRUE))
        dimnames(out$quantiles) <- list(Y = Ynames, X = Xnames, temp)
        dimnames(out$quantiles.MCaccuracy) <- dimnames(out$quantiles)
    }

    if(isTRUE(keepYX)){
        ## save Y and X values in the output; useful for plotting methods
        if(is.null(priorY)){
            out$Y <- Y
            out$X <- X
        } else {
            out$Y <- X
            out$X <- Y
            }
    }

    class(out) <- 'probability'
    out
}

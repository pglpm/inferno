#' Calculate quantiles
#'
#' This function calculates the posterior probability `Pr(Y | X, data)`, where `Y` and `X` are two (non overlapping) sets of joint variate values. If `X` is omitted or `NULL`, then the posterior probability `Pr(Y | data)` is calculated. The function also gives quantiles about the possible variability of the probability `Pr(Y | X, newdata, data)` that we could have if more learning data were provided, as well as a number of samples of the possible values of such probabilities. If several joint values are given for `Y` or `X`, the function will create a 2D grid of results for all possible compbinations of the given `Y` and `X` values. This function also allows for base-rate or other prior-probability corrections: If a prior (for instance, a base rate) for `Y` is given, the function will calculate the `Pr(Y | X, data, prior)` from `Pr(X | Y, data)` and the prior by means of Bayes's theorem. Each variate in each argument `Y`, `X` can be specified either as a point-value `Y = y` or as a left-open interval `Y ≤ y` or as a right-open interval `Y ≥ y`, through the argument `tails`.
#'
#' @param pY Named numerical vector or one-element list: set of probabilities of a variate for which we want to find the quantiles. The variate is given as the name of the vector or list. Default: `c(0.055, 0.25, 0.75, 0.945)`, that is, the 5.5%, 25%, 75%, 94.5% quantiles.
#' @param X Matrix or data.table or `NULL` (default): set of values of variates on which we want to condition. If `NULL`, no conditioning is made (except for conditioning on the learning dataset and prior assumptions). One variate per column, one set of values per row.
#' @param learnt Either a character with the name of a directory or full path for a 'learnt.rds' object, produced by the [learn()] function, or such an object itself.
#' @param tails Named vector or list, or `NULL` (default). The names must match some or all of the variates in argument `X`. For variates in this list, the probability arguments are understood in an semi-open interval sense: `X ≤ x` or `X ≥ x`, an so on. A left-open interval `X ≤ x` is indicated by the values `'<='` or `'left'` or `-1`; a right-open interval `X ≥ x` is indicated by the values `'>='` or `'right'` or `+1`. Values `NULL`, `'=='`, `0` indicate that a point value `Y = y` (not an interval) should be calculated. **NB**: the semi-open intervals *always* include the given value; this is important for ordinal or rounded variates. For instance, if `X` is an integer variate, then to condition on  `X < 3` we should require `X <= 2`.
#' @param priorY Numeric vector with the same length as the rows of `Y`, or `TRUE`, or `NULL` (default): prior probabilities or base rates for the `Y` values. If `TRUE`, the prior probabilities are assumed to be all equal. For the moment only the value `NULL` is accepted.
#' @param nsamples Integer or `NULL` or `'all'` (default): desired number of samples of the variability of the quantile for `Y`. If `NULL`, no samples are reported. If `'all'` (or `Inf`), all samples obtained by the [learn()] function are used.
#' @param quantiles Numeric vector, between 0 and 1, or `NULL`: desired quantiles of the variability of the quantile for `Y`. Default `c(0.055, 0.25, 0.75, 0.945)`, that is, the 5.5%, 25%, 75%, 94.5% quantiles (these are typical quantile values in the Bayesian literature: they give 50% and 89% credibility intervals, which correspond to 1 shannons and 0.5 shannons of uncertainty). If `NULL`, no quantiles are calculated.
#' @param parallel Logical or positive integer or cluster object. `TRUE`: use roughly half of available cores; `FALSE`: use serial computation; integer: use this many cores. It can also be a cluster object previously created with [parallel::makeCluster()]; in this case the parallel computation will use this object.
#' @param silent Logical, default `FALSE`: give warnings or updates in the computation?
#' @param keepYX Logical, default `TRUE`: keep a copy of the `Y` and `X` arguments in the output? This is used for the plot method.
#'
#' @return A list of the elements `values`,  `quantiles` (possibly `NULL`), `samples` (possibly `NULL`), `values.MCaccuracy`, `quantiles.MCaccuracy` (possibly `NULL`), `Y`, `X`. Element `values`: a matrix with the requested `Y`-quantiles conditional on the requested `X`-values, for all combinations of `pY` (rows) and `X` (columns). Element `quantiles`: an array with the variability quantiles (3rd dimension of the array). Element `samples`: an array with the variability samples (3rd dimension of the array). Elements `values.MCaccuracy` and `quantiles.MCaccuracy`: arrays with the numerical accuracies (roughly speaking a standard deviation) of the Monte Carlo calculations for the `values` and `quantiles` elements. Elements `pY`, `X`: copies of the `pY` and `X` arguments.
#'
#' @import parallel
#' 
##  #' @export
qPr <- function(
    pY,
    X = NULL,
    learnt,
    tails = NULL,
    priorY = NULL,
    nsamples = 'all',
    quantiles = c(0.055, 0.5, 0.945),
    parallel = NULL,
    silent = FALSE,
    keepYX = TRUE
) {
    ## #' @param usememory Logical, default `TRUE`: save partial results to disc, to avoid excessive RAM use. (For the moment only possible value is `TRUE`.)
    usememory <- TRUE

    Qerror <- pnorm(c(-1, 1))

    ## if (!silent) {
    ##     cat('\n')
    ## }

#### Requested parallel processing
    ## NB: doesn't make sense to have more cores than chains
    closeexit <- FALSE
    if ('cluster' %in% class(parallel)){
        ## user provides a cluster object
        cl <- parallel
    } else if (isTRUE(parallel)) {
        ## user wants us to register a parallel backend
        ## and to choose number of cores
        ncores <- max(1,
            floor(parallel::detectCores() / 2))
        cl <- parallel::makeCluster(ncores)
        ## doParallel::registerDoParallel(cl)
        closeexit <- TRUE
        message('Registered ', capture.output(print(cl)), '.')
    } else if (isFALSE(parallel)) {
        ## user wants us not to use parallel cores
        ncores <- 1
        cl <- parallel::makeCluster(ncores)
    } else if (is.numeric(parallel) &&
                   is.finite(parallel) && parallel >= 1) {
        ## user wants us to register 'parallel' # of cores
        ncores <- parallel
        cl <- parallel::makeCluster(ncores)
        closeexit <- TRUE
        message('Registered ', capture.output(print(cl)), '.')
    } else {
        stop("Unknown value of argument 'parallel'.")
    }

    ## Close parallel connections if any were opened
    if(closeexit) {
        closecoresonexit <- function(){
            message('Closing connections to cores.')
            parallel::stopCluster(cl)
            ## parallel::setDefaultCluster(NULL)
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


    pY <- as.list(pY)
    if(length(pY) > 1){stop('Specify only one variate in "pY".')}
    Yv <- names(pY)

    if(all(is.na(X))){X <- NULL}
    if(!is.null(X)){X <- as.data.frame(X)}
    Xv <- colnames(X)

    if(!is.null(tails)){
        tails <- as.list(tails)
        if(is.null(names(tails))) {
            stop('Missing variate names in "tails"')
        }
    }

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
    if (auxmetadata[auxmetadata$name == Yv, 'mcmctype'] %in% c('B', 'N')){
        stop('quantiles are undefined for binary and nominal variates.')
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

    tailsv <- names(tails)
    if (!all(tailsv %in% Xv)) {
        warning('"tails" variate ',
            paste0(tailsv[!(tailsv %in% Xv)], collapse = ' '),
            ' not among X; ignored\n')
        tails <- tails[(tailsv %in% Xv)]
        tailsv <- names(tails)
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

    nY <- length(pY[[1]])
    Y <- pY[[1]]
    nX <- max(nrow(X), 1L)
    auxY <- auxmetadata[auxmetadata$name == Yv, ]

    dosamples <- !is.null(nsamples)


#### tmp dir where to save X and Y objects
    temporarydir <- tempdir()


#### Calculate and save arrays for X values:
    if (is.null(X)) {
        lprobX <- log(learnt$W)
        saveRDS(lprobX,
            file.path(temporarydir,
                paste0('__X', 1, '__.rds'))
        )
    } else {
        ## Construction of the arguments for util_lprobs, X argument
        lpargs <- util_lprobsargsyx(
            x = X,
            auxmetadata = auxmetadata,
            learnt = learnt,
            tails = tails
        )

        invisible(parallel::parLapply(cl = cl,
            X = lpargs$xVs,
            fun = util_lprobsbase,
            params = lpargs$params,
            logW = c(log(learnt$W)),
            temporarydir = temporarydir,
            lab = '__X'
        ))
    }

#### Determine the type of Y variate, set params accordingly
    if(auxY$mcmctype == 'O'){
        params1 = learnt$Oprob[auxY$indexpos + seq_len(auxY$Nvalues), ,]
        params2 = NULL
        util_qYX <- util_qYXdiscr
    } else if(auxY$mcmctype == 'R'){
        params1 = learnt$Rmean[auxY$id, ,]
        params2 = sqrt(learnt$Rvar[auxY$id, ,])
        util_qYX <- util_qYXcont
    } else if(auxY$mcmctype == 'D'){
        params1 = learnt$Dmean[auxY$id, ,]
        params2 = sqrt(learnt$Dvar[auxY$id, ,])
        util_qYX <- util_qYXcont
    } else if(auxY$mcmctype == 'C'){
        params1 = learnt$Cmean[auxY$id, ,]
        params2 = sqrt(learnt$Cvar[auxY$id, ,])
        util_qYX <- util_qYXcont
    } else {
        stop('type of Y not found')
    }

#### Calculation with all pY and X combinations
    ## keys <- c('values', 'quantiles', 'samples', 'values.MCaccuracy', 'quantiles.MCaccuracy')
    keys <- c('values', 'quantiles', 'samples')
    ##
    combfnr <- function(...){setNames(do.call(mapply,
        c(FUN = `rbind`, lapply(X = ..., FUN = `[`, keys, drop = FALSE))),
        keys)}
    ## combfnc <- function(...){setNames(do.call(mapply, c(FUN=cbind, lapply(list(...), `[`, keys))), keys)}

    out <- combfnr(parallel::parApply(cl = cl,
            X = expand.grid(pY = Y, jx = seq_len(nX)),
            MARGIN = 1,
            FUN = util_qYX,
            auxmetadata = auxY,
            params1 = params1, params2 = params2,
            temporarydir = temporarydir, usememory = usememory,
            doquantiles = doquantiles, quantiles = quantiles,
            dosamples = dosamples, nsamples = nsamples,
            Qerror = Qerror,
            eps = .Machine$double.eps * 10
    ))

    ## clean temp files
    if(usememory) {
        unlink(x = sapply(seq_len(nX), function(jx){
            file.path(temporarydir, paste0('__X', jx, '__.rds'))
            }))
    }


    ## if(is.null(priorY)){
    ##     y <- Y
    ##     y[, colnames(Y) %in% tailsv] <- NA
    ##     jacobians <- exp(rowSums(
    ##         as.matrix(vtransform(y,
    ##             auxmetadata = auxmetadata,
    ##             logjacobianOr = TRUE)),
    ##         na.rm = TRUE
    ##     ))
    ##     rm(y)
    ## }

    ## transform to grid
    ## in the output-list elements the Y & X values are the rows
    ## dim(out$values.MCaccuracy) <- dim(out$values) <- c(nY, nX)
    if(is.null(priorY)){

        ## if(ncol(Y) == 1){Ynames <- Y[, 1]} else {Ynames <- NULL}
        Ynames <- pY

        if(!is.null(X)){
            Xnames <- apply(X = X, MARGIN = 1, FUN = paste0, collapse=',',
                simplify = TRUE)
        } else {
            Xnames <- NULL
        }
        dimnames(out$values) <- c(Ynames, list(X = Xnames))
        ## dimnames(out$values.MCaccuracy) <- dimnames(out$values)
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
        dimnames(out$values) <- c(Ynames, list(X = Xnames))
    }

    if(dosamples){
        ## transform to grid
        dim(out$samples) <- c(nY, nX, nsamples)

        if(is.null(priorY)){
            ## multiply by jacobian factors
            ## out$samples <- out$samples * jacobians
        } else {
            ## Bayes's theorem
            out$samples <- priorY * aperm(a = out$samples, perm = c(2, 1, 3),
                resize = FALSE)
            normf <- c(t(colSums(out$samples, na.rm = TRUE)))
            out$samples <- aperm(a = aperm(a = out$samples, perm = NULL,
                resize = FALSE) / normf, perm = NULL,
                resize = FALSE)
        }

        dimnames(out$samples) <- c(Ynames, list(X = Xnames,
            round(seq(1, nmcsamples, length.out = nsamples))))
    }

    if(doquantiles){
        if(is.null(priorY)){
            ## transform to grid
            dim(out$quantiles) <- c(nY, nX, length(quantiles))
            ## dim(out$quantiles.MCaccuracy) <- c(nY, nX, length(quantiles))
            ## multiply by jacobian factors
            ## out$quantiles <- out$quantiles * jacobians
            ## out$quantiles.MCaccuracy <- signif(x = out$quantiles.MCaccuracy * jacobians,
            ##     digits = 2)
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
        dimnames(out$quantiles) <- c(Ynames, list(X = Xnames, temp))
        ## dimnames(out$quantiles.MCaccuracy) <- dimnames(out$quantiles)
    }

    if(isTRUE(keepYX)){
        ## save Y and X values in the output; useful for plotting methods
        if(is.null(priorY)){
            out$pY <- pY
            out$X <- X
        } else {
            out$Y <- X
            out$X <- Y
            }
    }

    class(out) <- 'probability'
    out
}





#' Calculate quantiles for Y by bisection
#'
#' @keywords internal
util_qYXcont <- function(
    iyx,
    params1, params2,
    auxmetadata,
    temporarydir, usememory = TRUE,
    doquantiles, quantiles,
    dosamples, nsamples,
    Qerror,
    eps = .Machine$double.eps * 3
) {
    pY <- iyx['pY']

    if(usememory) {
        lprobX <- readRDS(file.path(temporarydir,
            paste0('__X', iyx['jx'], '__.rds')
        ))
    }
    sumlpX <- colSums(exp(lprobX), na.rm = TRUE)

    maxsamples <- ncol(params1)

    Yvals <- .Machine$double.xmax * c(-0.125, 0.125)

#### Calculate quantile for the posterior probability distribution

    values <- (Yvals[1] + Yvals[2]) / 2

    FF <- mean(colSums(exp(
        lprobX + pnorm(q = values,
            mean = params1, sd = params2,
            lower.tail = TRUE, log.p = TRUE)
    ), na.rm = TRUE) / sumlpX) - pY

    while(abs(FF) > eps && Yvals[2] - Yvals[1] > eps){
        Yvals[(FF > 0) + 1L] <- values
        values <- (Yvals[1] + Yvals[2]) / 2
        FF <- mean(colSums(exp(
            lprobX + pnorm(q = values,
                mean = params1, sd = params2,
                lower.tail = TRUE, log.p = TRUE)
        ), na.rm = TRUE) / sumlpX) - pY
    }
    values <- unname(unlist(vtransform(values,
        auxmetadata = auxmetadata,
        Rout = 'original',
        Cout = 'original',
        Dout = 'original',
        Oout = 'original',
        Nout = 'original',
        Bout = 'original',
        variates <- auxmetadata$name,
        logjacobianOr = NULL)))


#### Calculate quantile for the frequency samples
    if(doquantiles){
        nmaxsamples <- ncol(params1)
        selsamples <- TRUE
    } else if(dosamples) {
        nmaxsamples <- nsamples
        selsamples <- round(seq(1, ncol(params1), length.out = nsamples))
    }

    if(doquantiles || dosamples) {
        params1 <- t(params1[, selsamples])
        params2 <- t(params2[, selsamples])
        lprobX <- t(lprobX[, selsamples])

        Yvals <- rep(.Machine$double.xmax * c(-0.125, 0.125), each = nmaxsamples)
        dim(Yvals) <- c(nmaxsamples, 2)
        sampleseq <- 1:nmaxsamples

        samples <- (Yvals[, 1] + Yvals[, 2]) / 2
        FF <- rowSums(exp(
            lprobX + pnorm(q = samples,
                mean = params1, sd = params2,
                lower.tail = TRUE, log.p = TRUE)
        ), na.rm = TRUE) / sumlpX - pY

        tocheck <- abs(FF) > eps
        while(any(tocheck)) {
            choose <- c(sampleseq[tocheck], (FF[tocheck] > 0) + 1L)
            dim(choose) <- c(sum(tocheck), 2)
            Yvals[choose] <- samples[tocheck]
            samples[tocheck] <- (Yvals[tocheck, 1] + Yvals[tocheck, 2]) / 2
            FF <- rowSums(exp(
                lprobX + pnorm(q = samples,
                    mean = params1, sd = params2,
                    lower.tail = TRUE, log.p = TRUE)
            ), na.rm = TRUE) / sumlpX - pY
        tocheck <- abs(FF) > eps & Yvals[, 2] - Yvals[, 1] > eps
        }
        samples <- unname(unlist(vtransform(samples,
            auxmetadata = auxmetadata,
            Rout = 'original',
            Cout = 'original',
            Dout = 'original',
            Oout = 'original',
            Nout = 'original',
            Bout = 'original',
            variates <- auxmetadata$name,
            logjacobianOr = NULL)))
    }

    list(
        values = values,
        ##
        quantiles = if(doquantiles) {
            quantile(x = samples, probs = quantiles, type = 6,
                na.rm = TRUE, names = FALSE)
        },
        ##
        samples = if(dosamples) {
            samples[round(seq(1, length(samples), length.out = nsamples))]
        }
        ## values.MCaccuracy
        ## quantiles.MCaccuracy
)
}

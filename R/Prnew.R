#' Calculate posterior probabilities
#'
#' This function calculates the probability P(Y | X, data), where Y and X are two (non overlapping) sets of joint variates. The function also gives quantiles about the possible variability of the probability P(Y | X, newdata, data) that we could have if more learning data were provided, as well as a number of samples of the possible values of such probabilities. If several joint values are given for Y or X, the function will create a 2D grid of results for all possible compbinations of the given Y and X values. This function also allows for base-rate or other prior-probability corrections: If a prior (for instance, a base rate) for Y is given, the function will calculate the P(Y | X, data, prior) from P(X | Y, data) and the prior by means of Bayes's theorem.
#'
#' @param Y matrix or data.table: set of values of variates of which we want
#'   the joint probability of. One variate per column, one set of values per row.
#' @param X matrix or data.table or `NULL`: set of values of variates on which we want to condition the joint probability of `Y`. If `NULL` (default), no conditioning is made (except for conditioning on the learning dataset and prior assumptions). One variate per column, one set of values per row.
#' @param learnt either a string with the name of a directory or full path for a 'learnt.rds' object, produced by the \code{\link{learn}} function, or such an object itself.
#' @param cumul named vector or list, or `NULL` (default). The names must match some or all of the variates in arguments `Y` and `X`. The values can be one of `'<='` or `'>='` or `'=='` or `NULL`.
#' @param priorY numeric vector with the same length as the rows of `Y`, or `TRUE`, or `NULL` (default): prior probabilities or base rates for the `Y` values. If `TRUE`, the prior probabilities are assumed to be all equal.
#' @param nsamples integer or `NULL` or `"all"`: desired number of samples of the variability of the probability for `Y`. If `NULL`, no samples are reported. If `"all"` (or `Inf`), all samples obtained by the \code{\link{learn}} function are used. Default `"all"`.
#' @param quantiles numeric vector, between 0 and 1, or `NULL`: desired quantiles of the variability of the probability for `Y`. Default `c(0.055, 0.25, 0.75, 0.945)`, that is, the 5.5%, 25%, 75%, 94.5% quantiles (these are typical quantile values in the Bayesian literature: they give 50% and 89% credibility intervals, which correspond to 1 shannons and 0.5 shannons of uncertainty). If `NULL`, no quantiles are calculated.
#' @param parallel Logical or `NULL` or positive integer: `TRUE`: use roughly half of available cores; `FALSE`: use serial computation; `NULL`: don't do anything (use pre-registered condition); integer: use this many cores. Default `NULL`
#' @param silent logical: give warnings or updates in the computation?
#'   Default `FALSE`.
#' @param usememory logical: save partial results to disc, to avoid crashes?
#'   Default `TRUE`.
#' @param keepYX logical, default `TRUE`: keep a copy of the `Y` and `X` arguments in the output? This is used for the plot method.
#'
#' @return A list of class `probability`, consisting of the elements `values`,  `quantiles` (possibly `NULL`), `samples` (possibly `NULL`), `Y`, `X`. Element `values`: a matrix with the probabilities P(Y|X,data,assumptions), for all combinations of values of `Y` (rows) and `X` (columns). Element `quantiles`: an array with the variability quantiles (3rd dimension of the array) for such probabilities. Element `samples`: an array with the variability samples (3rd dimension of the array) for such probabilities. Elements `Y`, `X`: copies of the `Y` and `X` arguments.
#'
#' @import parallel foreach doParallel
#'
#' @export
Pr2 <- function(
    Y,
    X = NULL,
    learnt,
    cumul = NULL,
    priorY = NULL,
    nsamples = 'all',
    quantiles = c(0.055, 0.25, 0.75, 0.945),
    parallel = NULL,
    silent = FALSE,
    usememory = TRUE,
    keepYX = TRUE
) {
    cumulcentre <- list('==', 0, NULL)
    cumulleft <- list('<=', -1, 'left')
    cumulright <- list('>=', +1, 'right')

    cumulvalues <- c(cumulcentre, cumulleft, cumulright)

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
    Yv <- colnames(Y)

    if(all(is.na(X))){X <- NULL}
    if(!is.null(X)){X <- as.data.frame(X)}
    Xv <- colnames(X)

    if(!is.null(cumul)){cumul <- as.list(cumul)}
    cumulv <- names(cumul)


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

    if (!all(cumulv %in% c(Yv, Xv))) {
        warning('variate ',
            paste0(cumulv[!(cumulv %in% c(Yv, Xv))], collapse = ' '),
            ' not among Y and X; ignored\n')
    }
    if (length(unique(cumulv)) != length(cumulv)) {
        stop('duplicate "cumul" variates\n')
    }
    if(!all(cumul %in% cumulvalues)) {
        stop('"cumul" values must be ',
            paste0(cumulvalues, collapse = ' '), '\n')
    }

    ## transform 'cumul' to -1, +1
    ## +1: '<',    -1: '>'
    ## this is opposite of the argument convention because
    ## interval probabilities are calculated with `lower.tail = TRUE`
    ## eg:
    ## pnorm(x, mean, sd, lower.tail = FALSE) ==
    ##     pnorm(-x, -mean, sd, lower.tail = TRUE)
    cumul[cumul %in% cumulcentre] <- NULL
    cleft <- cumul %in% cumulleft
    cright <- cumul %in% cumulright
    cumul[cleft] <- +1
    cumul[cright] <- -1
    cumul <- unlist(cumul)

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

    dosamples <- !is.null(nsamples)



#### Construction of the arguments for util_lprobs, X argument

    lprobsargs <- util_lprobsargs(
        x = X,
        cumul = cumul,
        auxmetadata = auxmetadata,
        learnt = learnt
    )

#### First calculate and save arrays for X values:
    if (is.null(X)) {
        lprobX <- log(learnt$W)
        usememory <- FALSE
    } else {
        ## create unique dir where to save X objects
        temporarydir <- tempdir()

        invisible(foreach(jj = seq_len(nX),
            xV0 = lprobsargs$xV0,
            xV1 = lprobsargs$xV1,
            xV2 = lprobsargs$xV2,
            xVN = lprobsargs$xVn,
            xVB = lprobsargs$xVb,
            .combine = `c`,
            .inorder = TRUE) %dox% {
                lprobX <- c(log(learnt$W)) +
                   util_lprobs(
                       nV0 = lprobsargs$nV0,
                       V0mean = lprobsargs$V0mean,
                       V0sd = lprobsargs$V0sd,
                       xV0 = xV0,
                       nV1 = lprobsargs$nV1,
                       V1mean = lprobsargs$V1mean,
                       V1sd = lprobsargs$V1sd,
                       xV1 = xV1,
                       nV2 = lprobsargs$nV2,
                       V2mean = lprobsargs$V2mean,
                       V2sd = lprobsargs$V2sd,
                       V2steps = lprobsargs$V2steps,
                       xV2 = xV2,
                       nVN = lprobsargs$nVN,
                       VNprobs = lprobsargs$VNprobs,
                       xVN = xVN,
                       nVB = lprobsargs$nVB,
                       VBprobs = lprobsargs$VBprobs,
                       xVB = c(xVB)
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
            })
    }

#### Now calculate for each Y value, combining with each X value
    ## transformation of inputs
    Y2 <- as.matrix(vtransform(Y, auxmetadata = auxmetadata,
        Rout = 'normalized',
        Cout = 'boundisinf',
        Dout = 'normalized',
        Oout = 'index',
        Nout = 'index',
        Bout = 'numeric',
        logjacobianOr = NULL))

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

    if(is.null(priorY)){
        jacobians <- exp(rowSums(
            as.matrix(vtransform(Y,
                auxmetadata = auxmetadata,
                logjacobianOr = TRUE)),
            na.rm = TRUE
        ))
    }

    ## transform to grid
    ## in the output-list elements the Y & X values are the rows
    dim(out$values) <- c(nY, nX)

    if(is.null(priorY)){
        ## multiply by jacobian factors
        out$values <- out$values * jacobians

        ## if(ncol(Y) == 1){Ynames <- Y[, 1]} else {Ynames <- NULL}
        Ynames <- apply(Y, 1, paste0, collapse=',')

        if(!is.null(X)){
        Xnames <- apply(X, 1, paste0, collapse=',')
        } else {
            Xnames <- NULL
        }
    } else {
        ## Bayes's theorem
        out$values <- t(priorY * t(out$values))
        normf <- rowSums(out$values, na.rm = TRUE)
        out$values <- t(out$values/normf)

        ## if(ncol(X) == 1){Ynames <- X[, 1]} else {Ynames <- NULL}
        Ynames <- apply(X, 1, paste0, collapse=',')

        if(!is.null(Y)){
            Xnames <- apply(Y, 1, paste0, collapse=',')
        } else {
            Xnames <- NULL
        }
    }
    dimnames(out$values) <- list(Y = Ynames, X = Xnames)

    if(dosamples){
        ## transform to grid
        dim(out$samples) <- c(nY, nX, nsamples)

        if(is.null(priorY)){
            ## multiply by jacobian factors
            out$samples <- out$samples * jacobians
        } else {
            ## Bayes's theorem
            out$samples <- priorY * aperm(out$samples, c(2,1,3))
            normf <- c(t(colSums(out$samples, na.rm = TRUE)))
            out$samples <- aperm(aperm(out$samples)/normf)
        }

        dimnames(out$samples) <- list(Y = Ynames, X = Xnames,
            round(seq(1, nmcsamples, length.out = nsamples)))
    }

    if(doquantiles){
        if(is.null(priorY)){
            ## transform to grid
            dim(out$quantiles) <- c(nY, nX, length(quantiles))
            ## multiply by jacobian factors
            out$quantiles <- out$quantiles * jacobians
        } else {
            ## calculate quantiles from samples
            out$quantiles <- aperm(
                apply(out$samples, c(1, 2), quantile,
                    probs = quantiles, type = 6,
                    na.rm = TRUE, names = FALSE),
                c(2,3,1) )

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

    out$lowertail = NA

    class(out) <- 'probability'
    out
}

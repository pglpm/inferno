#' Generate datapoints
#'
#' This function generate datapoints according to the posterior probability `Pr(Y | X, data)` calculated with [learn()], for the variates specified in the argument `Y`, and conditional on the variate values specified in the argument `X`. If `X` is omitted or `NULL`, then the posterior probability `Pr(Y | data)` is used. Each variate in the argument `X` can be specified either as a point-value `X = x` or as a left-open interval `X ≤ x` or as a right-open interval `X ≥ x`, through the argument `tails`.
#'
#' @param n Positive integer: number of samples to draw.
#' @param Ynames Character vector: names of variates to draw jointly
#' @param X List or data.table or `NULL`: set of values of variates on which we want to condition the joint probability for `Y`. If `NULL` (default), no conditioning is made. Any rows beyond the first are discarded
#' @param learnt Either a character with the name of a directory or full path for a 'learnt.rds' object, produced by the [learn()] function, or such an object itself.
#' @param tails Named vector or list, or `NULL` (default). The names must match some or all of the variates in arguments `X`. For variates in this list, the probability conditional is understood in an semi-open interval sense: `X ≤ x` or `X ≥ x`, an so on. See analogous argument in [Pr()].
#' @param mcsamples Vector of integers, or `'all'`, or `NULL` (default): which Monte Carlo samples calculated by the [learn()] function should be used to draw the variate values. The default is to choose a random subset if `n` is smaller than their number, otherwise to recycle them as necessary.
#'
#' @return A data frame of joint draws of the variates `Ynames` from the posterior distribution, conditional on `X`. The row names of the data frame report the Monte Carlo sample (from [learn()]) used for that draw, and the total number of draws from that sample so far.
#'
#' @importFrom extraDistr rcat
#' @importFrom extraDistr rbern
#'
#' @export
rPr <- function(
    n,
    Ynames,
    X = NULL,
    learnt,
    tails = NULL,
    mcsamples = NULL
) {

    ## Extract Monte Carlo output & aux-metadata
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
                stop("The argument 'learnt' must be a folder containing learnt.rds, or the path to an rds-file containing the output from 'learn()'.")
            }
        }
    }
    ## Add check to see that learnt is correct type of object?
    auxmetadata <- learnt$auxmetadata
    learnt$auxmetadata <- NULL
    learnt$auxinfo <- NULL
    ncomponents <- nrow(learnt$W)
    nmcsamples <- ncol(learnt$W)

    if(is.null(mcsamples) ||
           (is.character(mcsamples) && mcsamples == 'all') ||
           isTRUE(mcsamples)) {
        mcsamples <- seq_len(nmcsamples)
    } else if (any(!is.finite(mcsamples)) || any(mcsamples < 1)) {
        stop("'mcsamples' should be a list of positive integers or NULL")
    } else {
        mcsamples <- round(mcsamples[mcsamples <= nmcsamples])
    }

    nmcs <- length(mcsamples)
    if(n <= nmcs) {
        sseq <- mcsamples[sort.int(sample.int(nmcs, n))]
    } else {
        sseq <- c(rep(x = mcsamples, times = n %/% nmcs),
            mcsamples[sort.int(sample.int(nmcs, n %% nmcs))])
    }

    if(all(is.na(X))){X <- NULL}
    if(!is.null(X)){
        X <- as.data.frame(X)
        if (nrow(X) > 1) {
            message('Only the first row of X is considered')
            X <- X[1, , drop = FALSE]
        }
    }
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
    if(!is.character(Ynames) || any(is.na(Ynames))){
        stop('Ynames must be a vector of variate names')
    }
    if(!all(Ynames %in% auxmetadata$name)) {
        stop('unknown Y variates\n')
    }
    if(length(unique(Ynames)) != length(Ynames)) {
        stop('duplicate Y variates\n')
    }

    if (!all(Xv %in% auxmetadata$name)) {
        stop('unknown X variates\n')
    }
    if (length(unique(Xv)) != length(Xv)) {
        stop('duplicate X variates\n')
    }
    if(any(Ynames %in% Xv)) {
        stop('overlap in Y and X variates\n')
    }

    if (!all(tailsv %in% Xv)) {
        warning('"tails" variate ',
            paste0(tailsv[!(tailsv %in% Xv)], collapse = ' '),
            ' not among X; ignored\n')
        tails <- tails[tailsv %in% Xv]
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


#### Adjust component weights W for conditioning on X
    if(is.null(X)){
        lW <- log(learnt$W)
    } else {
        lpargs <- util_lprobsargsyx(
            x = X,
            auxmetadata = auxmetadata,
            learnt = learnt,
            tails = tails
        )

        lW <- util_lprobsbase(
            xVs = lpargs$xVs[[1]],
            params = lpargs$params,
            logW =  log(learnt$W)
        ) # rows=components, columns=samples

    } # end definition of lW if non-null X


    ## Utility function to avoid finite-precision errors
    ## now externally defined
    ## denorm <- function(lprob) {
    ##     apply(X = lprob, MARGIN = 2, FUN = function(xx) {
    ##         xx - max(xx[is.finite(xx)])
    ##     }, simplify = TRUE)
    ## }

#### Draw samples of Ynames

    lWnorm <- util_denorm(lW[, sseq, drop = FALSE])
    Ws <- extraDistr::rcat(n = n, prob = t(
        apply(X = lWnorm, MARGIN = 2, FUN = function(xx){
            xx <- exp(xx)
            xx/sum(xx, na.rm = TRUE)
        }, simplify = TRUE)
    ))

    Yout <- NULL
    vYout <- NULL

    ## R
    toselect <- which((auxmetadata$name %in% Ynames) &
                          (auxmetadata$mcmctype == 'R'))
    nvrt <- length(toselect)
    if(nvrt > 0){
        aux <- auxmetadata[toselect, ]
        vYout <- c(vYout, aux$name)
        totake <- cbind(rep.int(x = aux$id, times = rep(n, nvrt)), Ws, sseq)
        Yout <- c(Yout,
            rnorm(n = n * nvrt,
                mean = learnt$Rmean[totake],
                sd = learnt$Rsd[totake] )
        )
    }

    ## C
    toselect <- which((auxmetadata$name %in% Ynames) &
                          (auxmetadata$mcmctype == 'C'))
    nvrt <- length(toselect)
    if(nvrt > 0){
        aux <- auxmetadata[toselect, ]
        vYout <- c(vYout, aux$name)
        totake <- cbind(rep.int(x = aux$id, times = rep(n, nvrt)), Ws, sseq)
        Yout <- c(Yout,
            rnorm(n = n * nvrt,
                mean = learnt$Cmean[totake],
                sd = learnt$Csd[totake] )
        )
    }

    ## D
    toselect <- which((auxmetadata$name %in% Ynames) &
                          (auxmetadata$mcmctype == 'D'))
    nvrt <- length(toselect)
    if(nvrt > 0){
        aux <- auxmetadata[toselect, ]
        vYout <- c(vYout, aux$name)
        totake <- cbind(rep.int(x = aux$id, times = rep(n, nvrt)), Ws, sseq)
        Yout <- c(Yout,
            rnorm(n = n * nvrt,
                mean = learnt$Dmean[totake],
                sd = learnt$Dsd[totake] )
        )
    }

    ## O
    toselect <- which((auxmetadata$name %in% Ynames) &
                          (auxmetadata$mcmctype == 'O'))
    nvrt <- length(toselect)
    if(nvrt > 0){
        vYout <- c(vYout, auxmetadata$name[toselect])
        for(i in toselect) {
            aux <- auxmetadata[i, ]
            totake <- cbind(Ws, sseq)
            Yout <- c(Yout,
                extraDistr::rcat(n = n,
                    prob = apply(
                        X = learnt$Oprob[aux$indexpos + seq_len(aux$Nvalues), ,],
                        MARGIN = 1, FUN = `[`, totake,
                        simplify = TRUE) )
            )
        }
    }

    ## N
    toselect <- which((auxmetadata$name %in% Ynames) &
                          (auxmetadata$mcmctype == 'N'))
    nvrt <- length(toselect)
    if(nvrt > 0){
        vYout <- c(vYout, auxmetadata$name[toselect])
        for(i in toselect) {
            aux <- auxmetadata[i, ]
            totake <- cbind(Ws, sseq)
            Yout <- c(Yout,
                extraDistr::rcat(n = n,
                    prob = apply(
                        X = learnt$Nprob[aux$indexpos + seq_len(aux$Nvalues), ,],
                        MARGIN = 1, FUN = `[`, totake,
                        simplify = TRUE) )
            )
        }
    }

    ## B
    toselect <- which((auxmetadata$name %in% Ynames) &
                          (auxmetadata$mcmctype == 'B'))
    nvrt <- length(toselect)
    if(nvrt > 0){
        aux <- auxmetadata[toselect, ]
        vYout <- c(vYout, aux$name)
        totake <- cbind(rep.int(x = aux$id, times = rep(n, nvrt)), Ws, sseq)
        Yout <- c(Yout,
            extraDistr::rbern(n = n * nvrt,
                prob = learnt$Bprob[totake])
        )
    }

    dim(Yout) <- c(n, length(Ynames))
    Yout <- Yout[, match(Ynames, vYout), drop = FALSE]
    colnames(Yout) <- Ynames

    Yout <- vtransform(Yout,
        auxmetadata = auxmetadata,
        Rout = 'original',
        Cout = 'original',
        Dout = 'original',
        Oout = 'original',
        Nout = 'original',
        Bout = 'original',
        logjacobianOr = NULL)

    ## row-name scheme: 'mcsample.draw'
    rownames(Yout) <- paste0(sseq, '_', ((seq_len(n) - 1L) %/% nmcs) + 1L)

    Yout
}

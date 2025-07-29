#' Calculate mutual information between groups of joint variates
#'
#' This function calculates various entropic information measures of two variates (each variate may consist of joint variates): the mutual information, the conditional entropies, and the entropies.
#'
#' @param Y1names String vector: first group of joint variates
#' @param Y2names String vector or NULL: second group of joint variates
#' @param X matrix or data.frame or NULL: values of some variates conditional on which we want the probabilities.
#' @param learnt Either a string with the name of a directory or full path
#'   for an 'learnt.rds' object, or such an object itself.
#' @param tails named vector or list, or `NULL` (default). The names must match some or all of the variates in arguments `X`. For variates in this list, the probability conditional is understood in an semi-open interval sense: `X ≤ x` or `X ≥ x`, an so on. See analogous argument in \code{\link{Pr()}}.
#' @param n integer or `NULL` (default): number of samples from which to approximately calculate the mutual information. Default as many as Monte Carlo samples in `learnt`.
#' @param unit Either one of 'Sh' for *shannon* (default), 'Hart' for *hartley*, 'nat' for *natural unit*, or a positive real indicating the base of the logarithms to be used.
#' @param parallel Logical or `NULL` or positive integer: `TRUE`: use roughly half of available cores; `FALSE`: use serial computation; `NULL`: don't do anything (use pre-registered condition); integer: use this many cores. Default `NULL`
#' @param silent logical: give warnings or updates in the computation?
#'
#' @return A list consisting of the elements `MI`, `CondEn12`, `CondEn21`, `En1`, `En2`, `MImax`, `unit`, `Y1names`, `Y1names`. All elements except `unit`, `Y1names`, `Y2names` are a vector of `value` and `error`. Element `MI` is the mutual information between (joint) variates `Y1names` and (joint) variates `Y2names`. Element`CondEn12` is the conditional entropy of the first variate given the second, and vice versa for `CondEn21`. Elements `En1` and `En1` are the (differential) entropies of the first and second variates. Element `MImax` is the maximum possible value of the mutual information. Elements `unit`, `Y1names`, `Y2names` are identical to the same inputs.
#'
#' @import parallel foreach doParallel
#'
#' @export
mutualinfo <- function(
    Y1names,
    Y2names,
    X = NULL,
    learnt,
    tails = NULL,
    n = NULL,
    unit = 'Sh',
    parallel = NULL,
    silent = FALSE
){

#### Mutual information and conditional entropy between Y2 and Y1
#### conditional on X, data, prior information
#### are calculated by Monte Carlo integration:
#### 0. adjusted component weights are calculated for conditioning on X:
####     all probabilities below are conditional on X & data & prior
#### 1. joint samples of Y1_i, Y2_i are drawn
#### 2. probabilities p(Y1|Y2) are calculated for each sample
####    the conditional entropy is Monte-Carlo approximated by
####    H(Y1|Y2) = - sum_{i} log p(Y1_i | Y2_i)
#### 3. probabilities p(Y1) are calculated for each sample
####    the entropy is Monte-Carlo approximated by
####    H(Y1) = - sum_{i} log p(Y1_i)
#### 4. the mutual info is Monte-Carlo approximated by
####    I(Y1|Y2) = sum_{i} [log p(Y1_i | Y2_i) - log p(Y1_i)]
####           = -H(Y1|Y2) + H(Y1)
####
#### For these computations it is not necessary to transform the Y1,Y2 variates
#### from the internal Monte Carlo representation to the original one

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

    ## determine if parallel computation is possible and needed
    if (ncores < 2) {
        `%dochains%` <- `%do%`
    } else {
        `%dochains%` <- `%dopar%`
    }

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
    nmcs <- ncol(learnt$W)

    if(is.null(n) || n == 0) {n <- 1}
    if(n > 0){
        n <- n * nmcs
    } else if(n < 0) {
        n <- -n
    }

    if(n <= nmcs) {
        sseq <- sort(sample.int(nmcs, n))
    } else {
        sseq <- c(rep(x = seq_len(nmcs), times = n %/% nmcs),
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
    if (unit == 'Sh') {
        base <- 2
    } else if (unit == 'Hart') {
        base <- 10
    } else if (unit == 'nat') {
        base <- exp(1)
    } else if (is.numeric(unit) && unit > 0) {
        base <- unit
    } else {
        stop("unit must be 'Sh', 'Hart', 'nat', or a positive real")
    }

    if(!is.character(Y1names) || any(is.na(Y1names))){
        stop('Y1names must be a vector of variate names')
    }
    if(!is.null(Y2names) && (!is.character(Y2names) || any(is.na(Y2names)))){
        stop('Y2names must be NULL or a vector of variate names')
    }

    ## More consistency checks
    if(!all(Y1names %in% auxmetadata$name)) {
        stop('unknown Y1 variates\n')
    }
    if(length(unique(Y1names)) != length(Y1names)) {
        stop('duplicate Y1 variates\n')
    }
    ##
    if(!is.null(Y2names) && !all(Y2names %in% auxmetadata$name)){
        stop('unknown Y2 variates\n')
    }
    if(length(unique(Y2names)) != length(Y2names)) {
        stop('duplicate Y2 variates\n')
    }


    if (!all(Xv %in% auxmetadata$name)) {
        stop('unknown X variates\n')
    }
    if (length(unique(Xv)) != length(Xv)) {
        stop('duplicate X variates\n')
    }
    ##
    if(any(Y1names %in% Y2names)) {
        stop('overlap in Y1 and Y2 variates\n')
    }
    if(any(Y1names %in% Xv)) {
        stop('overlap in Y1 and X variates\n')
    }
    if(any(Y2names %in% Xv)) {
        stop('overlap in Y2 and X variates\n')
    }

    if (!all(tailsv %in% Xv)) {
        warning('variate ',
            paste0(tailsv[!(tailsv %in% Xv)], collapse = ' '),
            ' not among X; ignored\n')
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


#### Step 0. Adjust component weights W for conditioning on X
    if(is.null(X)){
        lW <- log(learnt$W)
    } else {
        lpargs <- util_lprobsargs(
            x = X,
            auxmetadata = auxmetadata,
            learnt = learnt,
            tails = tails
        )

        lW <- log(learnt$W) +
            util_lprobs(
                nV0 = lpargs$nV0,
                V0mean = lpargs$V0mean,
                V0sd = lpargs$V0sd,
                xV0 = lpargs$xV0,
                nV1 = lpargs$nV1,
                V1mean = lpargs$V1mean,
                V1sd = lpargs$V1sd,
                xV1 = lpargs$xV1,
                nV2 = lpargs$nV2,
                V2mean = lpargs$V2mean,
                V2sd = lpargs$V2sd,
                V2steps = lpargs$V2steps,
                xV2 = lpargs$xV2,
                nVN = lpargs$nVN,
                VNprobs = lpargs$VNprobs,
                xVN = lpargs$xVN,
                nVB = lpargs$nVB,
                VBprobs = lpargs$VBprobs,
                xVB = c(lpargs$xVB)
            ) # rows=components, columns=samples
    } # end definition of lW if non-null X


    ## Utility function to avoid finite-precision errors
    denorm <- function(lprob) {
        apply(X = lprob, MARGIN = 2, FUN = function(xx) {
            xx - max(xx[is.finite(xx)])
        }, simplify = TRUE)
    }


#### Combine Y1,Y2 into single Y for speed
    Ynames <- c(Y1names, Y2names)

#### STEP 1. Draw samples of Ynames (that is, Y1names,Y2names)

    lWnorm <- denorm(lW)
    Ws <- extraDistr::rcat(n = n, prob = t(
        apply(X = lWnorm[, sseq, drop = FALSE], MARGIN = 2, FUN = function(xx){
            xx <- exp(xx)
            xx / sum(xx, na.rm = TRUE)
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
                sd = sqrt(learnt$Rvar[totake]) )
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
                sd = sqrt(learnt$Cvar[totake]) )
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
                sd = sqrt(learnt$Dvar[totake]) )
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

    Y1transf <- Yout[, Y1names, drop = FALSE]
    Y2transf <- Yout[, Y2names, drop = FALSE]
    rm(Yout)
    gc()


#### STEP 2. Calculate sum_i log2_p(Y1|Y2) for all samples
    lpargs1 <- util_lprobsargs(
        x = Y1transf,
        auxmetadata = auxmetadata,
        learnt = learnt,
        tails = NULL
    )

    lpargs2 <- util_lprobsargs(
        x = Y2transf,
        auxmetadata = auxmetadata,
        learnt = learnt,
        tails = NULL
    )

    out <- foreach(
            x1V0 = lpargs1$xV0,
            x1V1 = lpargs1$xV1,
            x1V2 = lpargs1$xV2,
            x1VN = lpargs1$xVN,
            x1VB = lpargs1$xVB,
            ##
            x2V0 = lpargs2$xV0,
            x2V1 = lpargs2$xV1,
            x2V2 = lpargs2$xV2,
            x2VN = lpargs2$xVN,
            x2VB = lpargs2$xVB,
            .combine = rbind,
            .inorder = TRUE
        ) %dochains% {
### lprobY2
            lprobY2 <- util_lprobs(
                    nV0 = lpargs2$nV0,
                    V0mean = lpargs2$V0mean,
                    V0sd = lpargs2$V0sd,
                    xV0 = x2V0,
                    nV1 = lpargs2$nV1,
                    V1mean = lpargs2$V1mean,
                    V1sd = lpargs2$V1sd,
                    xV1 = x2V1,
                    nV2 = lpargs2$nV2,
                    V2mean = lpargs2$V2mean,
                    V2sd = lpargs2$V2sd,
                    V2steps = lpargs2$V2steps,
                    xV2 = x2V2,
                    nVN = lpargs2$nVN,
                    VNprobs = lpargs2$VNprobs,
                    xVN = x2VN,
                    nVB = lpargs2$nVB,
                    VBprobs = lpargs2$VBprobs,
                    xVB = c(x2VB)
            ) # rows=components, columns=samples

### lprobY1
            lprobY1 <- util_lprobs(
                    nV0 = lpargs1$nV0,
                    V0mean = lpargs1$V0mean,
                    V0sd = lpargs1$V0sd,
                    xV0 = x1V0,
                    nV1 = lpargs1$nV1,
                    V1mean = lpargs1$V1mean,
                    V1sd = lpargs1$V1sd,
                    xV1 = x1V1,
                    nV2 = lpargs1$nV2,
                    V2mean = lpargs1$V2mean,
                    V2sd = lpargs1$V2sd,
                    V2steps = lpargs1$V2steps,
                    xV2 = x1V2,
                    nVN = lpargs1$nVN,
                    VNprobs = lpargs1$VNprobs,
                    xVN = x1VN,
                    nVB = lpargs1$nVB,
                    VBprobs = lpargs1$VBprobs,
                    xVB = c(x1VB)
                ) # rows=components, columns=samples

            celWnorm <- colSums(exp(lWnorm))

### Construct probabilities from lprobY1, lprobY2
            lpY1and2 <- log(mean(
                colSums(exp(lprobY1 + lprobY2 + lWnorm)) / celWnorm,
                na.rm = TRUE), base = base)

            lpY1 <- log(mean(
                colSums(exp(lprobY1 + lWnorm)) / celWnorm,
                na.rm = TRUE), base = base)

            lpY2 <- log(mean(
                colSums(exp(lprobY2 + lWnorm)) / celWnorm,
                na.rm = TRUE), base = base)

            lprobnorm <- denorm(lprobY2 + lW)
            lpY1given2 <- log(mean(
                colSums(exp(lprobY1 + lprobnorm)) / colSums(exp(lprobnorm)),
                na.rm = TRUE), base = base)

            lprobnorm <- denorm(lprobY1 + lW)
            lpY2given1 <- log(mean(
                colSums(exp(lprobY2 + lprobnorm)) / colSums(exp(lprobnorm)),
                na.rm = TRUE), base = base)

            c(
                MI = lpY1and2 - lpY1 - lpY2,
                ## mi1 = lpY1given2 - lpY1,
                ## mi2 = lpY2given1 - lpY2,
                CondEn12 = -lpY1given2,
                CondEn21 = -lpY2given1,
                En1 = -lpY1,
                En2 = -lpY2
            )
        } # End foreach loop

        ## Jacobian factors
    logjacobians1 <- rowSums(
        as.matrix(vtransform(Y1transf,
            auxmetadata = auxmetadata,
            logjacobianOr = FALSE)),
        na.rm = TRUE)

    logjacobians2 <- rowSums(
        as.matrix(vtransform(Y2transf,
            auxmetadata = auxmetadata,
            logjacobianOr = FALSE)),
        na.rm = TRUE)

    out[, -1] <- out[, -1] - c(logjacobians1, logjacobians2)/log(base)

    out <- unlist(apply(X = rbind(
        value = colMeans(out, na.rm = TRUE),
        error = signif(x = apply(
            X = out, MARGIN = 2, FUN = sd, na.rm = TRUE, simplify = TRUE
        )/sqrt(n), digits = 2)
    ), MARGIN = 2, FUN = list, simplify = TRUE), recursive = FALSE)

    ## ## generally there's no MI maximum for continous variates
    ## mmax <- paste0('En', which.min(c(out$En1['value'], out$En2['value'])) )

    c(out,
        list(#MImax = out[[mmax]],
            unit = unit, Y1names = Y1names, Y2names = Y2names
        )
    )
}

#' Calculate mutual information between groups of joint variates
#'
#' @param Y1names String vector: first group of joint variates
#' @param Y2names String vector or NULL: second group of joint variates
#' @param X matrix or data.frame or NULL: values of some variates conditional on
#'   which we want the probabilities
#' @param learnt Either a string with the name of a directory or full path
#'   for an 'learnt.rds' object, or such an object itself
#' @param nsamples numeric: number of samples from which to approximately
#'   calculate the mutual information. Default 3600
#' @param unit Either one of 'Sh' (default), 'Hart', 'nat', or a positive real
#' @param parallel, logical or numeric: whether to use pre-existing parallel
#'   workers, or how many to create and use
#' @param silent logical: give warnings or updates in the computation
#'
#' @return A list with the mutual information, its error, and its unit
#'
#' @export
mutualinfo <- function(
    Y1names,
    Y2names,
    X = NULL,
    learnt,
    nsamples = 3600,
    unit = 'Sh',
    parallel = TRUE,
    silent = TRUE
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

    if (!silent) {
        cat('\n')
    }

    ## Utility function to avoid finite-precision errors
    denorm <- function(lprob) {
        apply(lprob, 2, function(xx) {
            xx - max(xx[is.finite(xx)])
        })
    }
    ## denorm <- function(lprob) {
    ##     apply(lprob, 2, function(xx) {
    ##         xx <- exp(xx - max(xx[is.finite(xx)]))
    ##         xx/sum(xx)
    ##     })
    ## }

#### Determine the status of parallel processing
    workers <- setupParallel(parallel, silent)
    ncores <- workers$ncores
    
    if (!is.logical(workers$cluster)) {
        on.exit(closeCoresOnExit(workers$cluster, silent))
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

    nMCsamples <- ncol(learnt$W)
    ncomponents <- nrow(learnt$W)

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
    if (!is.null(X) && length(dim(X)) != 2) {
        stop('X must be NULL or have two dimensions')
    }
    if (!is.null(X) && dim(X)[1] > 1) {
        message('Only the first row of X is considered')
        X <- X[1, , drop = FALSE]
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
    ##
    Xnames <- colnames(X)
    if (!all(Xnames %in% auxmetadata$name)) {
        stop('unknown X variates\n')
    }
    if (length(unique(Xnames)) != length(Xnames)) {
        stop('duplicate X variates\n')
    }
    ##
    if(length(intersect(Y1names, Y2names)) > 0) {
        stop('overlap in Y1 and Y2 variates\n')
    }
    if(length(intersect(Y1names, Xnames)) > 0) {
        stop('overlap in Y1 and X variates\n')
    }
    if(length(intersect(Y2names, Xnames)) > 0) {
        stop('overlap in Y2 and X variates\n')
    }


#### Calculate how many samples per MC sample
    ## ***todo: account for the case where nsamples < nMCsamples
    if(nsamples < 0) {
        nsamples <- -nsamples * nMCsamples
    } else if (nsamples == 0 || !is.finite(nsamples)) {
        stop("'nsamples' cannot be zero")
    }
    n <- ceiling(nsamples/nMCsamples)*nMCsamples
    sseq <- seq_len(nMCsamples)
    ## source('vtransform.R')

#### Combine Y1,Y2 into single Y for speed
    Ynames <- c(Y1names, Y2names)

### Guide to indices:
    ## .i. = order in X/Y corresponding to appearance in vnames
    ## .t. = vnames present in X/Y, kept in their vnames-order

#### Type R
    vnames <- auxmetadata[auxmetadata$mcmctype == 'R', 'name']
    Y2iR <- match(vnames, Y2names)
    Y2tR <- which(!is.na(Y2iR))
    Y2iR <- Y2iR[Y2tR]
    Y2nR <- length(Y2iR)
    ##
    Y1iR <- match(vnames, Y1names)
    Y1tR <- which(!is.na(Y1iR))
    Y1iR <- Y1iR[Y1tR]
    Y1nR <- length(Y1iR)
    ##
    XiR <- match(vnames, Xnames)
    XtR <- which(!is.na(XiR))
    XiR <- XiR[XtR]
    XnR <- length(XiR)
    if (Y1nR > 0 || Y2nR > 0 || XnR > 0) {
        learnt$Rvar <- sqrt(learnt$Rvar)
    }
    ##
    YiR <- match(vnames, Ynames)
    YtR <- which(!is.na(YiR))
    YiR <- YiR[YtR]
    YnR <- length(YiR)

#### Type C
    vnames <- auxmetadata[auxmetadata$mcmctype == 'C', 'name']
    Y2iC <- match(vnames, Y2names)
    Y2tC <- which(!is.na(Y2iC))
    Y2iC <- Y2iC[Y2tC]
    Y2nC <- length(Y2iC)
    ##
    Y1iC <- match(vnames, Y1names)
    Y1tC <- which(!is.na(Y1iC))
    Y1iC <- Y1iC[Y1tC]
    Y1nC <- length(Y1iC)
    ##
    XiC <- match(vnames, Xnames)
    XtC <- which(!is.na(XiC))
    XiC <- XiC[XtC]
    XnC <- length(XiC)
    if (Y1nC > 0 || Y2nC > 0 || XnC > 0) {
        learnt$Cvar <- sqrt(learnt$Cvar)
        Clefts <- auxmetadata[match(vnames, auxmetadata$name), 'tdomainminplushs']
        Crights <- auxmetadata[match(vnames, auxmetadata$name), 'tdomainmaxminushs']
    }
    ##
    YiC <- match(vnames, Ynames)
    YtC <- which(!is.na(YiC))
    YiC <- YiC[YtC]
    YnC <- length(YiC)

#### Type D
    vnames <- auxmetadata[auxmetadata$mcmctype == 'D', 'name']
    Y2iD <- match(vnames, Y2names)
    Y2tD <- which(!is.na(Y2iD))
    Y2iD <- Y2iD[Y2tD]
    Y2nD <- length(Y2iD)
    ##
    Y1iD <- match(vnames, Y1names)
    Y1tD <- which(!is.na(Y1iD))
    Y1iD <- Y1iD[Y1tD]
    Y1nD <- length(Y1iD)
    ##
    XiD <- match(vnames, Xnames)
    XtD <- which(!is.na(XiD))
    XiD <- XiD[XtD]
    XnD <- length(XiD)
    if (Y1nD > 0 || Y2nD > 0 || XnD > 0) {
        learnt$Dvar <- sqrt(learnt$Dvar)
        Dsteps <- auxmetadata[match(vnames, auxmetadata$name), 'halfstep'] /
            auxmetadata[match(vnames, auxmetadata$name), 'tscale']
        Dlefts <- auxmetadata[match(vnames, auxmetadata$name), 'tdomainminplushs']
        Drights <- auxmetadata[match(vnames, auxmetadata$name), 'tdomainmaxminushs']
    }
    ##
    YiD <- match(vnames, Ynames)
    YtD <- which(!is.na(YiD))
    YiD <- YiD[YtD]
    YnD <- length(YiD)

#### Type O
    vnames <- auxmetadata[auxmetadata$mcmctype == 'O', 'name']
    Y2iO <- match(vnames, Y2names)
    Y2tO <- which(!is.na(Y2iO))
    Y2iO <- Y2iO[Y2tO]
    Y2nO <- length(Y2iO)
    ##
    Y1iO <- match(vnames, Y1names)
    Y1tO <- which(!is.na(Y1iO))
    Y1iO <- Y1iO[Y1tO]
    Y1nO <- length(Y1iO)
    ##
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
    Y2iN <- match(vnames, Y2names)
    Y2tN <- which(!is.na(Y2iN))
    Y2iN <- Y2iN[Y2tN]
    Y2nN <- length(Y2iN)
    ##
    Y1iN <- match(vnames, Y1names)
    Y1tN <- which(!is.na(Y1iN))
    Y1iN <- Y1iN[Y1tN]
    Y1nN <- length(Y1iN)
    ##
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
    Y2iB <- match(vnames, Y2names)
    Y2tB <- which(!is.na(Y2iB))
    Y2iB <- Y2iB[Y2tB]
    Y2nB <- length(Y2iB)
    ##
    Y1iB <- match(vnames, Y1names)
    Y1tB <- which(!is.na(Y1iB))
    Y1iB <- Y1iB[Y1tB]
    Y1nB <- length(Y1iB)
    ##
    XiB <- match(vnames, Xnames)
    XtB <- which(!is.na(XiB))
    XiB <- XiB[XtB]
    XnB <- length(XiB)
    ##
    YiB <- match(vnames, Ynames)
    YtB <- which(!is.na(YiB))
    YiB <- YiB[YtB]
    YnB <- length(YiB)


#### STEP 0. Adjust component weights W for conditioning on X
    if(is.null(X)){
        lW <- log(learnt$W)
    } else {
        x <- t(as.matrix(vtransform(X, auxmetadata = auxmetadata,
            Rout = 'normalized',
            Cout = 'boundisinf',
            Dout = 'normalized',
            Oout = 'numeric',
            Nout = 'numeric',
            Bout = 'numeric',
            logjacobianOr = NULL)))
        ##
        lW <- log(learnt$W) +
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
            )
    } # end definition of lW if non-null X


#### STEP 1. Draw samples of Y (that is, Y1,Y2)

    ## Y is drawn as follows, for each MC sample:
    ## 1. draw a component, according to its probability
    ## 2. draw from the appropriate kernel distributions
    ## using the parameters of that component

    lWnorm <- denorm(lW)
    Ws <- extraDistr::rcat(n=n, prob=t(
        apply(lWnorm, 2, function(xx){
            xx <- exp(xx)
            xx/sum(xx, na.rm = TRUE)})
    ))
    Yout <- c(
    (if(YnR > 0){# continuous
        totake <- cbind(rep(YtR,each=n), Ws, sseq)
        rnorm(n=n*YnR,
            mean=learnt$Rmean[totake],
            sd=learnt$Rvar[totake]
        )
    }else{NULL}),
    (if(YnC > 0){# censored
        totake <- cbind(rep(YtC,each=n), Ws, sseq)
        rnorm(n=n*YnC,
            mean=learnt$Cmean[totake],
            sd=learnt$Cvar[totake]
        )
    }else{NULL}),
    (if(YnD > 0){# continuous discretized
        totake <- cbind(rep(YtD,each=n), Ws, sseq)
        rnorm(n=n*YnD,
            mean=learnt$Dmean[totake],
            sd=learnt$Dvar[totake]
        )
    }else{NULL}),
    (if(YnO > 0){# nominal
        totake <- cbind(rep(YtO,each=n), Ws, sseq)
        extraDistr::rcat(n=n*YnO,
            prob=apply(learnt$Oprob,3,`[`,totake))
    }else{NULL}),
    (if(YnN > 0){# nominal
        totake <- cbind(rep(YtN,each=n), Ws, sseq)
        extraDistr::rcat(n=n*YnN,
            prob=apply(learnt$Nprob,3,`[`,totake))
    }else{NULL}),
    (if(YnB > 0){# binary
        totake <- cbind(rep(YtB,each=n), Ws, sseq)
        extraDistr::rbern(n=n*YnB,
            prob=learnt$Bprob[totake])
    }else{NULL})
    )

    ## rows: n samples, cols: Y variates
    dim(Yout) <- c(n, length(Ynames))
    ## Match to original order of Ynames
    Yout <- Yout[, match(Ynames, Ynames[c(YiR, YiC, YiD, YiO, YiN, YiB)]),
        drop = FALSE]
    colnames(Yout) <- Ynames
    Yout <- as.matrix(vtransform(Yout,
        auxmetadata = auxmetadata,
        Rout = 'mi',
        Cout = 'mi',
        Dout = 'mi',
        Oout = 'mi',
        Nout = 'mi',
        Bout = 'mi',
        logjacobianOr = NULL))

    Y1transf <- Yout[, Y1names, drop = FALSE]
    Y2transf <- Yout[, Y2names, drop = FALSE]
    rm(Yout)
    gc()


#### STEP 2. Calculate sum_i log2_p(Y1|Y2) for all samples

    ## from samplesFDistribution.R with some modifications
    out <- foreach(y1 = t(Y1transf), y2 = t(Y2transf),
        .combine = `rbind`,
        .inorder = TRUE) %dochains% {
#### the loop is over the columns of y1 and y2
#### each instance is a 1-column vector

### lprobY2
            ## rows: components, cols: samples
            if (all(is.na(y2))) {
                lprobY2 <- array(0, dim = dim(lW))
            } else {
                lprobY2 <- util_lprob(
                    x = y2,
                    learnt = learnt,
                    nR = Y2nR, iR = Y2iR, tR = Y2tR,
                    nC = Y2nC, iC = Y2iC, tC = Y2tC,
                    Clefts = Clefts, Crights = Crights,
                    nD = Y2nD, iD = Y2iD, tD = Y2tD,
                    Dsteps = Dsteps, Dlefts = Dlefts, Drights = Drights,
                    nO = Y2nO, iO = Y2iO, tO = Y2tO,
                    nN = Y2nN, iN = Y2iN, tN = Y2tN,
                    nB = Y2nB, iB = Y2iB, tB = Y2tB
                )

            }# end lprobY2

### lprobY1
            lprobY1 <- util_lprob(
                x = y1,
                learnt = learnt,
                nR = Y1nR, iR = Y1iR, tR = Y1tR,
                nC = Y1nC, iC = Y1iC, tC = Y1tC,
                Clefts = Clefts, Crights = Crights,
                nD = Y1nD, iD = Y1iD, tD = Y1tD,
                Dsteps = Dsteps, Dlefts = Dlefts, Drights = Drights,
                nO = Y1nO, iO = Y1iO, tO = Y1tO,
                nN = Y1nN, iN = Y1iN, tN = Y1tN,
                nB = Y1nB, iB = Y1iB, tB = Y1tB
            )
            ## Output: rows=components, columns=samples

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
        } # End loop through generated samples

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

    c(
        unlist(apply(rbind(
            value = colMeans(out, na.rm = TRUE),
            error = apply(out, 2, sd, na.rm = TRUE)/sqrt(n)
        ), 2, list), recursive = FALSE),
        list(Y1names = Y1names, Y2names = Y2names, unit = unit)
    )
}

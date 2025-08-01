#' Calculate posterior probabilities
#'
#' This function calculates the probability P(Y | X, data), where Y and X are two (non overlapping) sets of joint variates. The function also gives quantiles about the possible variability of the probability P(Y | X, newdata, data) that we could have if more learning data were provided, as well as a number of samples of the possible values of such probabilities. If several joint values are given for Y or X, the function will create a 2D grid of results for all possible compbinations of the given Y and X values. This function also allows for base-rate or other prior-probability corrections: If a prior (for instance, a base rate) for Y is given, the function will calculate the P(Y | X, data, prior) from P(X | Y, data) and the prior by means of Bayes's theorem.
#'
#' @param Y matrix or data.table: set of values of variates of which we want
#'   the joint probability of. One variate per column, one set of values per row.
#' @param X matrix or data.table or `NULL`: set of values of variates on which we want to condition the joint probability of `Y`. If `NULL` (default), no conditioning is made (except for conditioning on the learning dataset and prior assumptions). One variate per column, one set of values per row.
#' @param learnt either a string with the name of a directory or full path for a 'learnt.rds' object, produced by the \code{\link{learn}} function, or such an object itself.
#' @param less to be documented
#' @param greater to be documented
#' @param priorY numeric vector with the same length as the rows of `Y`, or `TRUE`, or `NULL` (default): prior probabilities or base rates for the `Y` values. If `TRUE`, the prior probabilities are assumed to be all equal.
#' @param nsamples integer or `NULL` or `"all"`: desired number of samples of the variability of the probability for `Y`. If `NULL`, no samples are reported. If `"all"` (or `Inf`), all samples obtained by the \code{\link{learn}} function are used. Default `100`.
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
Pr2 <- function(
    Y,
    X = NULL,
    learnt,
    less = NULL,
    greater = NULL,
    priorY = NULL,
    nsamples = 100L,
    quantiles = c(0.055, 0.25, 0.75, 0.945),
    parallel = NULL,
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
    if(all(is.na(X))){X <- NULL}
    if(!is.null(X)){X <- as.data.frame(X)}

    ## Consistency checks
    if (length(dim(Y)) != 2) {
        stop('Y must have two dimensions')
    }
    if (!is.null(X) && length(dim(X)) != 2) {
        stop('X must be NULL or NA or have two dimensions')
    }
    ##
    if (!is.null(X) && ncol(X) == 0) {
        stop('X must be NULL or NA or have two dimensions')
    }

    ## More consistency checks
    Yv0 <- colnames(Y)
    if (!all(Yv0 %in% auxmetadata$name)) {
        stop('unknown Y variates\n')
    }
    if (length(unique(Yv0)) != length(Yv0)) {
        stop('duplicate Y variates\n')
    }
    ##
    Xv0 <- colnames(X)
    if (!all(Xv0 %in% auxmetadata$name)) {
        stop('unknown X variates\n')
    }
    if (length(unique(Xv0)) != length(Xv0)) {
        stop('duplicate X variates\n')
    }
    ##
    if (length(intersect(Yv0, Xv0)) > 0) {
        stop('overlap in Y and X variates\n')
    }


    ## Consistency checks for 'less' and 'greater' args
    if (all(is.na(less))){less <- NULL}
    if (all(is.na(greater))){greater <- NULL}

    if ((isTRUE(less) && !is.null(greater)) ||
           (isTRUE(greater) && !is.null(less))){
        stop("if either 'less' or 'greater' is true, the other must be NULL\n")
    }

    ## meaningless to require < or > for nominal and binary variates
    cumulvnames <- auxmetadata[
        auxmetadata$mcmctype %in% c('R', 'C', 'D', 'O'),
        'name']
    if (isTRUE(less)){
        less <- cumulvnames
    }
    if (isTRUE(greater)){
        greater <- cumulvnames
    }

    if (!is.null(less) && !all(less %in% cumulvnames)) {
        stop("unknown or impossible variates in 'less'\n")
    }
    if (!is.null(greater) && !all(greater %in% cumulvnames)) {
        stop("unknown or impossible variates in 'greater'\n")
    }
    if (length(intersect(less, greater)) > 0){
        stop("overlap in 'less' and 'greater' variates\n")
    }

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
        . <- Yv0
        Yv0 <- Xv0
        Xv0 <- .
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

    ## division into 'less' and 'greater' groups
    XvL <- Xv0[Xv0 %in% less]
    XvU <- Xv0[Xv0 %in% greater]
    Xv <- setdiff(Xv0, c(XvL, XvU))
    YvL <- Yv0[Yv0 %in% less]
    YvU <- Yv0[Yv0 %in% greater]
    Yv <- setdiff(Yv0, c(YvL, YvU))


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
    XiLR <- match(vnames, XvL)
    XtLR <- which(!is.na(XiLR))
    XiLR <- XiLR[XtLR]
    XnLR <- length(XiLR)
    ##
    XiUR <- match(vnames, XvU)
    XtUR <- which(!is.na(XiUR))
    XiUR <- XiUR[XtUR]
    XnUR <- length(XiUR)
    ##
    YiR <- match(vnames, Yv)
    YtR <- which(!is.na(YiR))
    YiR <- YiR[YtR]
    YnR <- length(YiR)
    ##
    YiLR <- match(vnames, YvL)
    YtLR <- which(!is.na(YiLR))
    YiLR <- YiLR[YtLR]
    YnLR <- length(YiLR)
    ##
    YiUR <- match(vnames, YvU)
    YtUR <- which(!is.na(YiUR))
    YiUR <- YiUR[YtUR]
    YnUR <- length(YiUR)
    if (YnR > 0 || YnLR > 0 || YnUR > 0 ||
            XnR > 0 || XnLR > 0 || XnUR > 0) {
        learnt$Rvar <- sqrt(learnt$Rvar)
    }

#### Type C
    vnames <- auxmetadata[auxmetadata$mcmctype == 'C', 'name']
    XiC <- match(vnames, Xv)
    XtC <- which(!is.na(XiC))
    XiC <- XiC[XtC]
    XnC <- length(XiC)
    ##
    XiLC <- match(vnames, XvL)
    XtLC <- which(!is.na(XiLC))
    XiLC <- XiLC[XtLC]
    XnLC <- length(XiLC)
    ##
    XiUC <- match(vnames, XvU)
    XtUC <- which(!is.na(XiUC))
    XiUC <- XiUC[XtUC]
    XnUC <- length(XiUC)
    ##
    YiC <- match(vnames, Yv)
    YtC <- which(!is.na(YiC))
    YiC <- YiC[YtC]
    YnC <- length(YiC)
    ##
    YiLC <- match(vnames, YvL)
    YtLC <- which(!is.na(YiLC))
    YiLC <- YiLC[YtLC]
    YnLC <- length(YiLC)
    ##
    YiUC <- match(vnames, YvU)
    YtUC <- which(!is.na(YiUC))
    YiUC <- YiUC[YtUC]
    YnUC <- length(YiUC)
    if (YnC > 0 || YnLC > 0 || YnUC > 0 ||
            XnC > 0 || XnLC > 0 || XnUC > 0) {
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
    XiLD <- match(vnames, XvL)
    XtLD <- which(!is.na(XiLD))
    XiLD <- XiLD[XtLD]
    XnLD <- length(XiLD)
    ##
    XiUD <- match(vnames, XvU)
    XtUD <- which(!is.na(XiUD))
    XiUD <- XiUD[XtUD]
    XnUD <- length(XiUD)
    ##
    YiD <- match(vnames, Yv)
    YtD <- which(!is.na(YiD))
    YiD <- YiD[YtD]
    YnD <- length(YiD)
    ##
    YiLD <- match(vnames, YvL)
    YtLD <- which(!is.na(YiLD))
    YiLD <- YiLD[YtLD]
    YnLD <- length(YiLD)
    ##
    YiUD <- match(vnames, YvU)
    YtUD <- which(!is.na(YiUD))
    YiUD <- YiUD[YtUD]
    YnUD <- length(YiUD)
    if (YnD > 0 || YnLD > 0 || YnUD > 0 ||
            XnD > 0 || XnLD > 0 || XnUD > 0) {
        learnt$Dvar <- sqrt(learnt$Dvar)
        Dsteps <- auxmetadata[match(vnames, auxmetadata$name), 'halfstep'] /
            auxmetadata[match(vnames, auxmetadata$name), 'tscale']
        Dlefts <- auxmetadata[match(vnames, auxmetadata$name), 'tdomainminplushs']
        Drights <- auxmetadata[match(vnames, auxmetadata$name), 'tdomainmaxminushs']
    }

#### Type O
    vnames <- auxmetadata[auxmetadata$mcmctype == 'O', 'name']
    XiO <- match(vnames, Xv)
    XtO <- which(!is.na(XiO))
    XiO <- XiO[XtO]
    XnO <- length(XiO)
    ##
    XiLO <- match(vnames, XvL)
    XtLO <- which(!is.na(XiLO))
    XiLO <- XiLO[XtLO]
    XnLO <- length(XiLO)
    ##
    XiUO <- match(vnames, XvU)
    XtUO <- which(!is.na(XiUO))
    XiUO <- XiUO[XtUO]
    XnUO <- length(XiUO)
    ##
    YiO <- match(vnames, Yv)
    YtO <- which(!is.na(YiO))
    YiO <- YiO[YtO]
    YnO <- length(YiO)
    ##
    YiLO <- match(vnames, YvL)
    YtLO <- which(!is.na(YiLO))
    YiLO <- YiLO[YtLO]
    YnLO <- length(YiLO)
    ##
    YiUO <- match(vnames, YvU)
    YtUO <- which(!is.na(YiUO))
    YiUO <- YiUO[YtUO]
    YnUO <- length(YiUO)

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

        ## adjust for cumulative-probability requests
        if (XnLC > 0) {
            tochoose <- X2[ , XiLC] == -Inf
            X2[tochoose, XiLC] <- Clefts[XtLC]
        }
        if(XnLD > 0){
            X2[ , XiLD] <- X2[ , XiLD] + Dsteps[XtLD]
            tochoose <- X2[ , XiLC] > Drights[XtLD]
            X2[tochoose, XiLD] <- +Inf
        }
        ##
        if (XnUC > 0) {
            tochoose <- X2[ , XiUC] == Inf
            X2[tochoose, XiUC] <- Crights[XtUC]
        }
        if(XnUD > 0){
            X2[ , XiUD] <- X2[ , XiUD] - Dsteps[XtUD]
            tochoose <- X2[ , XiUC] < Dlefts[XtUD]
            X2[tochoose, XiUD] <- -Inf
        }
        if(XnUO > 0){
            X2[ , XiUO] <- X2[ , XiUO] - 1
        }

        temporarydir <- tempdir() # where to save X objects
        ##
        todelete <- foreach(jj = seq_len(nX), x = t(X2),
            .combine = `c`,
            .inorder = TRUE) %dox% {
                lprobX <- c(log(learnt$W)) +
                    util_lprobs(
                        x = x,
                        learnt = learnt,
                        nR = XnR, iR = XiR, tR = XtR,
                        nLR = XnLR, iLR = XiLR, tLR = XtLR,
                        nUR = XnUR, iUR = XiUR, tUR = XtUR,
                        ##
                        nC = XnC, iC = XiC, tC = XtC,
                        nLC = XnLC, iLC = XiLC, tLC = XtLC,
                        nUC = XnUC, iUC = XiUC, tUC = XtUC,
                        Clefts = Clefts, Crights = Crights,
                        ##
                        nD = XnD, iD = XiD, tD = XtD,
                        nLD = XnLD, iLD = XiLD, tLD = XtLD,
                        nUD = XnUD, iUD = XiUD, tUD = XtUD,
                        Dsteps = Dsteps, Dlefts = Dlefts, Drights = Drights,
                        ##
                        nO = XnO, iO = XiO, tO = XtO,
                        nLO = XnLO, iLO = XiLO, tLO = XtLO,
                        nUO = XnUO, iUO = XiUO, tUO = XtUO,
                        ##
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

    ## adjust for cumulative-probability requests
        if (YnLC > 0) {
            tochoose <- Y2[ , YiLC] == -Inf
            Y2[tochoose, YiLC] <- Clefts[YtLC]
        }
        if(YnLD > 0){
            Y2[ , YiLD] <- Y2[ , YiLD] + Dsteps[YtLD]
            tochoose <- Y2[ , YiLC] > Drights[YtLD]
            Y2[tochoose, YiLD] <- +Inf
        }
        ##
        if (YnUC > 0) {
            tochoose <- Y2[ , YiUC] == Inf
            Y2[tochoose, YiUC] <- Crights[YtUC]
        }
        if(YnUD > 0){
            Y2[ , YiUD] <- Y2[ , YiUD] - Dsteps[YtUD]
            tochoose <- Y2[ , YiUC] < Dlefts[YtUD]
            Y2[tochoose, YiUD] <- -Inf
        }
        if(YnUO > 0){
            Y2[ , YiUO] <- Y2[ , YiUO] - 1
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
                    lprobY <- util_lprobs(
                        x = y,
                        learnt = learnt,
                        nR = YnR, iR = YiR, tR = YtR,
                        nLR = YnLR, iLR = YiLR, tLR = YtLR,
                        nUR = YnUR, iUR = YiUR, tUR = YtUR,
                        ##
                        nC = YnC, iC = YiC, tC = YtC,
                        nLC = YnLC, iLC = YiLC, tLC = YtLC,
                        nUC = YnUC, iUC = YiUC, tUC = YtUC,
                        Clefts = Clefts, Crights = Crights,
                        ##
                        nD = YnD, iD = YiD, tD = YtD,
                        nLD = YnLD, iLD = YiLD, tLD = YtLD,
                        nUD = YnUD, iUD = YiUD, tUD = YtUD,
                        Dsteps = Dsteps, Dlefts = Dlefts, Drights = Drights,
                        ##
                        nO = YnO, iO = YiO, tO = YtO,
                        nLO = YnLO, iLO = YiLO, tLO = YtLO,
                        nUO = YnUO, iUO = YiUO, tUO = YtUO,
                        ##
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
    } else {
        ## Bayes's theorem
        out$values <- t(priorY * t(out$values))
        normf <- rowSums(out$values, na.rm = TRUE)
        out$values <- t(out$values/normf)
    }
    dimnames(out$values) <- list(Y = NULL, X = NULL)

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

        dimnames(out$samples) <- list(Y = NULL, X = NULL,
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
        dimnames(out$quantiles) <- list(Y = NULL, X = NULL, temp)
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

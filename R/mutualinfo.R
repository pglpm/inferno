#' Calculate mutual information between groups of joint variates
#'
#' @param Y1names String vector: first group of joint variates
#' @param Y2names String vector: second group of joint variates
#' @param X matrix or data.frame or NULL: values of some variates conditional on
#'   which we want the probabilities
#' @param mcoutput Either a string with the name of a directory or full path
#'   for a 'FDistribution.rds' object, or such an object itself
#' @param nsamples numeric: number of samples from which to approximately
#'   calculate the mutual information. Default 3600
#' @param unit Either one of 'Sh' (default), 'Hart', 'nat', or a positive real
#' @param parallel, logical or numeric: whether to use pre-existing parallel
#'   workers, or how many to create and use
#' @param silent logical: give warnings or updates in the computation
#'
#' @return A list with the mutual information, its error, and its unit
#'
#' @import parallel foreach doParallel
#' @importFrom extraDistr rcat rbern
#'
#' @export
mutualinfo <- function(
    Y1names,
    Y2names,
    X = NULL,
    mcoutput,
    nsamples = 3600,
    unit = 'Sh',
    parallel = TRUE,
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

    if (!silent) {
        cat('\n')
    }

    ## Utility function to avoid finite-precision errors
    denorm <- function(lprob) {
        apply(lprob, 2, function(xx) {
            xx <- exp(xx - max(xx[is.finite(xx)]))
            xx/sum(xx)
        })
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

    if (ncores < 2) {
        `%dochains%` <- `%do%`
    } else {
        `%dochains%` <- `%dopar%`
    }

    ## Extract Monte Carlo output & aux-metadata
    ## If mcoutput is a string, check if it's a folder name or file name
    if (is.character(mcoutput)) {
        ## Check if 'mcoutput' is a folder containing Fdistribution.rds
        if (file_test('-d', mcoutput) &&
            file.exists(file.path(mcoutput, 'Fdistribution.rds'))) {
            mcoutput <- readRDS(file.path(mcoutput, 'Fdistribution.rds'))
        } else {
            ## Assume 'mcoutput' the full path of Fdistributions.rds
            ## possibly without the file extension '.rds'
            mcoutput <- paste0(sub('.rds$', '', mcoutput), '.rds')
            if (file.exists(mcoutput)) {
                mcoutput <- readRDS(mcoutput)
            } else {
                stop("The argument 'mcoutput' must be a folder containing Fdistribution.rds, or the path to an rds-file containing the output from 'inferpopulation'.")
            }
        }
    }
    ## Add check to see that mcoutput is correct type of object?
    auxmetadata <- mcoutput$auxmetadata
    mcoutput$auxmetadata <- NULL
    mcoutput$auxinfo <- NULL

    nMCsamples <- ncol(mcoutput$W)
    ncomponents <- nrow(mcoutput$W)

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
    if(!is.character(Y2names) || any(is.na(Y2names))){
        stop('Y2names must be a vector of variate names')
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
    if(!all(Y2names %in% auxmetadata$name)){
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
        mcoutput$Rvar <- sqrt(mcoutput$Rvar)
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
        mcoutput$Cvar <- sqrt(mcoutput$Cvar)
        Clefts <- auxmetadata[match(vnames, auxmetadata$name), 'tleftbound']
        Crights <- auxmetadata[match(vnames, auxmetadata$name), 'trightbound']
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
        mcoutput$Dvar <- sqrt(mcoutput$Dvar)
        Dsteps <- auxmetadata[match(vnames, auxmetadata$name), 'halfstep'] /
            auxmetadata[match(vnames, auxmetadata$name), 'tscale']
        Dlefts <- auxmetadata[match(vnames, auxmetadata$name), 'tleftbound']
        Drights <- auxmetadata[match(vnames, auxmetadata$name), 'trightbound']
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
        lW <- log(mcoutput$W)
    } else {
        x <- t(as.matrix(vtransform(X, auxmetadata = auxmetadata,
            Rout = 'normalized',
            Cout = 'boundisinf',
            Dout = 'normalized',
            Oout = 'numeric',
            Nout = 'numeric',
            Bout = 'numeric')))
        ##
        lW <- log(mcoutput$W) +
            (if (XnR > 0) { # continuous
                colSums(
                    dnorm(
                        x = x[XiR, ],
                        mean = mcoutput$Rmean[XtR, , , drop = FALSE],
                        sd = mcoutput$Rvar[XtR, , , drop = FALSE],
                        log = TRUE
                    ),na.rm = TRUE)
            } else {
                0
            }) +
            (if (XnC > 0) { # censored
                isfin <- is.finite(x[XiC, ])
                indf <- which(isfin)
                indi <- which(!isfin)
                (if (length(indf) > 0) {
                    colSums(
                        dnorm(
                            x = x[XiC[indf], ],
                            mean = mcoutput$Cmean[XtC[indf], , , drop = FALSE],
                            sd = mcoutput$Cvar[XtC[indf], , , drop = FALSE],
                            log = TRUE
                        ), na.rm = TRUE)
                } else {
                    0
                }) +
                    (if (length(indi) > 0) {
                        vt <- XtC[indi]
                        vx <- pmin(
                            pmax(Clefts[vt], x[XiC[indi], ]),
                            Crights[vt])
                        ## for upper tail, take opposite mean and value
                        colSums(
                            pnorm(
                                q = -abs(vx),
                                mean = -sign(vx) *
                                    mcoutput$Cmean[vt, , , drop = FALSE],
                                sd = mcoutput$Cvar[vt, , , drop = FALSE],
                                log.p = TRUE
                            ), na.rm = TRUE)
                    } else {
                        0
                    })
            } else {
                0
            }) +
            (if (XnD > 0) { # continuous discretized
                vrights <- pmax(x[XiD, ] + Dsteps[XtD], Dlefts[XtD])
                vrights[vrights >= Drights[XtD]] <- +Inf
                vlefts <- pmin(x[XiD, ] - Dsteps[XtD], Drights[XtD])
                vlefts[vlefts <= Dlefts[XtD]] <- -Inf
                colSums(log(
                    pnorm(
                        q = vrights,
                        mean = mcoutput$Dmean[XtD, , , drop = FALSE],
                        sd = mcoutput$Dvar[XtD, , , drop = FALSE]
                    ) -
                    pnorm(
                        q = vlefts,
                        mean = mcoutput$Dmean[XtD, , , drop = FALSE],
                        sd = mcoutput$Dvar[XtD, , , drop = FALSE]
                    )
                ), na.rm = TRUE)
            } else {
                0
            }) +
            (if (XnO > 0) { # nominal
                colSums(log(
                    aperm(
                        vapply(seq_len(XnO), function(v) {
                            mcoutput$Oprob[XtO[v], , x[XiO[v], ], ]
                        }, mcoutput$W),
                        c(3, 1, 2))
                ), na.rm = TRUE)
            } else {
                0
            }) +
            (if (XnN > 0) { # nominal
                colSums(log(
                    aperm(
                        vapply(seq_len(XnN), function(v) {
                            mcoutput$Nprob[XtN[v], , x[XiN[v], ], ]
                        }, mcoutput$W),
                        c(3, 1, 2))
                ), na.rm = TRUE)
            } else {
                0
            }) +
            (if (XnB > 0) { # binary
                colSums(log(
                    x[XiB, ] * mcoutput$Bprob[XtB, , , drop = FALSE] +
                    (1L - x[XiB, ]) *
                    (1 - mcoutput$Bprob[XtB, , , drop = FALSE])
                ), na.rm = TRUE)
            } else {
                0
            })
        ##
    } # end definition of lW if non-null X


#### STEP 1. Draw samples of Y (that is, Y1,Y2)

    ## Y is drawn as follows, for each MC sample:
    ## 1. draw a component, according to its probability
    ## 2. draw from the appropriate kernel distributions
    ## using the parameters of that component

    lWnorm <- denorm(lW)

    Ws <- extraDistr::rcat(n=n, prob=t(lWnorm))
    Yout <- c(
    (if(YnR > 0){# continuous
        totake <- cbind(rep(YtR,each=n), Ws, sseq)
        rnorm(n=n*YnR,
            mean=mcoutput$Rmean[totake],
            sd=mcoutput$Rvar[totake]
        )
    }else{NULL}),
    (if(YnC > 0){# censored
        totake <- cbind(rep(YtC,each=n), Ws, sseq)
        rnorm(n=n*YnC,
            mean=mcoutput$Cmean[totake],
            sd=mcoutput$Cvar[totake]
        )
    }else{NULL}),
    (if(YnD > 0){# continuous discretized
        totake <- cbind(rep(YtD,each=n), Ws, sseq)
        rnorm(n=n*YnD,
            mean=mcoutput$Dmean[totake],
            sd=mcoutput$Dvar[totake]
        )
    }else{NULL}),
    (if(YnO > 0){# nominal
        totake <- cbind(rep(YtO,each=n), Ws, sseq)
        extraDistr::rcat(n=n*YnO,
            prob=apply(mcoutput$Oprob,3,`[`,totake))
    }else{NULL}),
    (if(YnN > 0){# nominal
        totake <- cbind(rep(YtN,each=n), Ws, sseq)
        extraDistr::rcat(n=n*YnN,
            prob=apply(mcoutput$Nprob,3,`[`,totake))
    }else{NULL}),
    (if(YnB > 0){# binary
        totake <- cbind(rep(YtB,each=n), Ws, sseq)
        extraDistr::rbern(n=n*YnB,
            prob=mcoutput$Bprob[totake])
    }else{NULL})
    )

    ## rows: n samples, cols: Y variates
    dim(Yout) <- c(n, length(Ynames))
    ## Match to original order of Ynames
    Yout <- Yout[, match(Ynames, Ynames[c(YiR, YiC, YiD, YiO, YiN, YiB)])]
    colnames(Yout) <- Ynames
    Yout <- as.matrix(vtransform(Yout,
        auxmetadata = auxmetadata,
        Rout = 'mi',
        Cout = 'mi',
        Dout = 'mi',
        Oout = 'mi',
        Nout = 'mi',
        Bout = 'mi'))

    Y1transf <- Yout[, Y1names]
    Y2transf <- Yout[, Y2names]
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
            lprobY2 <- (
                (if (Y2nR > 0) { # continuous
                    colSums(
                        dnorm(
                            x = y2[Y2iR, ],
                            mean = mcoutput$Rmean[Y2tR, , , drop = FALSE],
                            sd = mcoutput$Rvar[Y2tR, , , drop = FALSE],
                            log = TRUE
                        ), na.rm = TRUE)
                } else {
                    0
                }) +
                (if (Y2nC > 0) { # censored
                    isfin <- is.finite(y2[Y2iC, ])
                    indf <- which(isfin)
                    indi <- which(!isfin)
                    (if (length(indf) > 0) {
                        colSums(
                            dnorm(
                                x = y2[Y2iC[indf], ],
                                mean = mcoutput$Cmean[Y2tC[indf], , , drop = FALSE],
                                sd = mcoutput$Cvar[Y2tC[indf], , , drop = FALSE],
                                log = TRUE
                            ), na.rm = TRUE)
                    } else {
                        0
                    }) +
                        (if (length(indi) > 0) {
                            vt <- Y2tC[indi]
                            vx <- pmin(
                                pmax(Clefts[vt], y2[Y2iC[indi], ]),
                                Crights[vt])
                            ## for upper tail, take opposite mean and value
                            colSums(
                                pnorm(
                                    q = -abs(vx),
                                    mean = -sign(vx) *
                                        mcoutput$Cmean[vt, , , drop = FALSE],
                                    sd = mcoutput$Cvar[vt, , , drop = FALSE],
                                    log.p = TRUE
                                ), na.rm = TRUE)
                        } else {
                            0
                        })
                } else {
                    0
                }) +
                (if (Y2nD > 0) { # discretized
                    vrights <- pmax(y2[Y2iD, ] + Dsteps[Y2tD], Dlefts[Y2tD])
                    vrights[vrights >= Drights[Y2tD]] <- +Inf
                    vlefts <- pmin(y2[Y2iD, ] - Dsteps[Y2tD], Drights[Y2tD])
                    vlefts[vlefts <= Dlefts[Y2tD]] <- -Inf
                    colSums(log(
                        pnorm(
                            q = vrights,
                            mean = mcoutput$Dmean[Y2tD, , , drop = FALSE],
                            sd = mcoutput$Dvar[Y2tD, , , drop = FALSE]
                        ) -
                        pnorm(
                            q = vlefts,
                            mean = mcoutput$Dmean[Y2tD, , , drop = FALSE],
                            sd = mcoutput$Dvar[Y2tD, , , drop = FALSE]
                        )
                    ), na.rm = TRUE)
                } else {
                    0
                }) +
                (if (Y2nO > 0) { # nominal
                    colSums(log(
                        aperm(
                            vapply(seq_len(Y2nO), function(v) {
                                mcoutput$Oprob[Y2tO[v], , y2[Y2iO[v], ], ]
                            }, mcoutput$W),
                            c(3, 1, 2))
                    ), na.rm = TRUE)
                } else {
                    0
                }) +
                (if (Y2nN > 0) { # nominal
                    colSums(log(
                        aperm(
                            vapply(seq_len(Y2nN), function(v) {
                                mcoutput$Nprob[Y2tN[v], , y2[Y2iN[v], ], ]
                            }, mcoutput$W),
                            c(3, 1, 2))
                    ), na.rm = TRUE)
                } else {
                    0
                }) +
                (if (Y2nB > 0) { # binary
                    colSums(log(
                        y2[Y2iB, ] * mcoutput$Bprob[Y2tB, , , drop = FALSE] +
                        (1L - y2[Y2iB, ]) *
                        (1 - mcoutput$Bprob[Y2tB, , , drop = FALSE])
                    ), na.rm = TRUE)
                } else {
                    0
                })
            ) # end lprobY2

### lprobY1
            ##
            if (all(is.na(y1))) {
                lprobY1 <- array(NA, dim = dim(W))
            } else {
                lprobY1 <- (
                    (if (Y1nR > 0) { # continuous
                    colSums(
                        dnorm(
                            x = y1[Y1iR, ],
                            mean = mcoutput$Rmean[Y1tR, , , drop = FALSE],
                            sd = mcoutput$Rvar[Y1tR, , , drop = FALSE],
                            log = TRUE
                        ), na.rm = TRUE)
                } else {
                    0
                }) +
                    (if (Y1nC > 0) { # censored
                        isfin <- is.finite(y1[Y1iC, ])
                        indf <- which(isfin)
                        indi <- which(!isfin)
                        (if (length(indf) > 0) {
                            colSums(
                                dnorm(
                                    x = y1[Y1iC[indf], ],
                                    mean = mcoutput$Cmean[Y1tC[indf], , , drop = FALSE],
                                    sd = mcoutput$Cvar[Y1tC[indf], , , drop = FALSE],
                                    log = TRUE
                                ), na.rm = TRUE)
                        } else {
                            0
                        }) +
                            (if (length(indi) > 0) {
                                vt <- Y1tC[indi]
                                vx <- pmin(
                                    pmax(Clefts[vt], y1[Y1iC[indi], ]),
                                    Crights[vt])
                                ## for upper tail, take opposite mean and value
                                colSums(
                                    pnorm(
                                        q = -abs(vx),
                                        mean = -sign(vx) *
                                            mcoutput$Cmean[vt, , , drop = FALSE],
                                        sd = mcoutput$Cvar[vt, , , drop = FALSE],
                                        log.p = TRUE
                                    ), na.rm = TRUE)
                            } else {
                                0
                            })
                    } else {
                        0
                    }) +
                    (if (Y1nD > 0) { # discretized
                        vrights <- y1[Y1iD, ] + Dsteps[Y1tD]
                        vrights[vrights >= Drights[Y1tD]] <- +Inf
                        vlefts <- y1[Y1iD, ] - Dsteps[Y1tD]
                        vlefts[vlefts <= Dlefts[Y1tD]] <- -Inf
                        colSums(log(
                            pnorm(
                                q = vrights,
                                mean = mcoutput$Dmean[Y1tD, , , drop = FALSE],
                                sd = mcoutput$Dvar[Y1tD, , , drop = FALSE]
                            ) -
                            pnorm(
                                q = vlefts,
                                mean = mcoutput$Dmean[Y1tD, , , drop = FALSE],
                                sd = mcoutput$Dvar[Y1tD, , , drop = FALSE]
                            )
                        ), na.rm = TRUE)
                    } else {
                        0
                    }) +
                    (if (Y1nO > 0) { # nominal
                        colSums(log(
                            aperm(
                                vapply(seq_len(Y1nO), function(v) {
                                    mcoutput$Oprob[Y1tO[v], , y1[Y1iO[v], ], ]
                                }, mcoutput$W),
                                c(3, 1, 2))
                        ), na.rm = TRUE)
                    } else {
                        0
                    }) +
                    (if (Y1nN > 0) { # nominal
                        colSums(log(
                            aperm(
                                vapply(seq_len(Y1nN), function(v) {
                                    mcoutput$Nprob[Y1tN[v], , y1[Y1iN[v], ], ]
                                }, mcoutput$W),
                                c(3, 1, 2))
                        ), na.rm = TRUE)
                    } else {
                        0
                    }) +
                    (if (Y1nB > 0) { # binary
                        colSums(log(
                        (y1[Y1iB, ] * mcoutput$Bprob[Y1tB, , , drop = FALSE]) +
                        ((1 - y1[Y1iB, ]) *
                         (1 - mcoutput$Bprob[Y1tB, , , drop = FALSE]))
                        ), na.rm = TRUE)
                    } else {
                        0
                    })
                )
            } # end lprobYY
            ## Output: rows=components, columns=samples


### Construct probabilities from lprobY1, lprobY2


            lpY1and2 <- log(mean(
                colSums(exp(lprobY1 + lprobY2 + lWnorm)) /
                colSums(exp(lWnorm))
            ), base = base)

            lpY1 <- log(mean(
                colSums(exp(lprobY1 + lWnorm)) /
                colSums(exp(lWnorm))
            ), base = base)

            lpY2 <- log(mean(
                colSums(exp(lprobY2 + lWnorm)) /
                colSums(exp(lWnorm))
            ), base = base)

            lprobnorm <- denorm(lprobY2 + lW)
            lpY1given2 <- log(mean(
                colSums(exp(lprobY1 + lprobnorm)) /
                colSums(exp(lprobnorm))
            ), base = base)

            lprobnorm <- denorm(lprobY1 + lW)
            lpY2given1 <- log(mean(
                colSums(exp(lprobY2 + lprobnorm)) /
                colSums(exp(lprobnorm))
            ), base = base)

            c(
                mi0 = lpY1and2 - lpY1 - lpY2,
                mi1 = lpY1given2 - lpY1,
                mi2 = lpY2given1 - lpY2,
                ce12 = lpY1given2,
                ce21 = lpY2given1,
                cond12a = lpY1given2,
                cond12b = lpY1and2 - lpY2
                )
        }

    str(out)

    list(
        MI = colMeans(out, na.rm = TRUE),
        error = apply(out, 2, sd, na.rm = TRUE),
        conds = exp(out[1,c('cond12a', 'cond12b')]),
        unit = unit
    )
}

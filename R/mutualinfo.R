#' Calculate mutual information between groups of joint variates
#'
#' @param Yvrt String vector: first group of joint variates
#' @param Xvrt String vector: second group of joint variates
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
mutualinfo <- function(Yvrt, Xvrt, mcoutput,
                       nsamples = 3600,
                       unit = 'Sh',
                       parallel = TRUE,
                       silent = FALSE){

#### Mutual information and conditional entropy between X and Y
#### are calculated by Monte Carlo integration:
#### 1. joint samples of Y_i, X_i are drawn
#### 2. probabilities p(Y|X) are calculated for each sample
####    the conditional entropy is Monte-Carlo approximated by
####    H(Y|X) = - sum_{i} log p(Y_i | X_i)
#### 3. probabilities p(Y|X) are calculated for each sample
####    the entropy is Monte-Carlo approximated by
####    H(Y) = - sum_{i} log p(Y_i)
#### 4. the mutual info is Monte-Carlo approximated by
####    I(Y|X) = sum_{i} [log p(Y_i | X_i) - log p(Y_i)]
####           = -H(Y|X) + H(Y)
####
#### For these computations it is not necessary to transform the Y,X variates
#### from the internal Monte Carlo representation to the original one

  if (!silent) {
    cat('\n')
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
                                        # Check if 'mcoutput' is a folder containing Fdistribution.rds
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

  nMCsamples <- ncol(mcoutput$W)
  nclusters <- nrow(mcoutput$W)

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

  
  if(!is.character(Yvrt) || any(is.na(Yvrt))){
    stop('Yvrt must be a vector of variate names')
  }
  if(!is.character(Xvrt) || any(is.na(Xvrt))){
    stop('Xvrt must be a vector of variate names')
  }

  ## More consistency checks
  if(!all(Yvrt %in% auxmetadata$name)){stop('unknown Y variates\n')}
  if(length(unique(Yvrt)) != length(Yvrt)){stop('duplicate Y variates\n')}
  ##
  if(!all(Xvrt %in% auxmetadata$name)){stop('unknown X variates\n')}
  if(length(unique(Xvrt)) != length(Xvrt)){stop('duplicate X variates\n')}
  ##
  if(length(intersect(Yvrt, Xvrt)) > 0){stop('overlap in Y and X variates\n')}


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

#### Combine Y,X into single Z for speed
  Zvrt <- c(Yvrt, Xvrt)

#### Type R
  vnames <- auxmetadata[auxmetadata$mcmctype == 'R', 'name']
  XiR <- match(vnames, Xvrt)
  XtR <- which(!is.na(XiR))
  XiR <- XiR[XtR]
  XnR <- length(XiR)
  ##
  YiR <- match(vnames, Yvrt)
  YtR <- which(!is.na(YiR))
  YiR <- YiR[YtR]
  YnR <- length(YiR)
  if (YnR > 0 || XnR > 0) {
    mcoutput$Rvar <- sqrt(mcoutput$Rvar)
  }
  ##
  ZiR <- match(vnames, Zvrt)
  ZtR <- which(!is.na(ZiR))
  ZiR <- ZiR[ZtR]
  ZnR <- length(ZiR)

#### Type C
  vnames <- auxmetadata[auxmetadata$mcmctype == 'C', 'name']
  XiC <- match(vnames, Xvrt)
  XtC <- which(!is.na(XiC))
  XiC <- XiC[XtC]
  XnC <- length(XiC)
  ##
  YiC <- match(vnames, Yvrt)
  YtC <- which(!is.na(YiC))
  YiC <- YiC[YtC]
  YnC <- length(YiC)
  if (YnC > 0 || XnC > 0) {
    mcoutput$Cvar <- sqrt(mcoutput$Cvar)
    Cbounds <- cbind(
      auxmetadata[match(vnames, auxmetadata$name), 'tcensormin'],
      ## sign is important here:
      ## for upper tail, take opposite mean and value
      -auxmetadata[match(vnames, auxmetadata$name), 'tcensormax']
      )
    ## Cbounds <- cbind(
    ##   c(vtransform(
    ##     x = matrix(NA,
    ##                nrow = 1, ncol = length(vnames),
    ##                dimnames = list(NULL, vnames)
    ##                ),
    ##     auxmetadata = auxmetadata,
    ##     Cout = 'sleft'
    ##   )),
    ##   ## sign is important here:
    ##   ## for upper tail, take opposite mean and value
    ##   -c(vtransform(
    ##      x = matrix(NA,
    ##                 nrow = 1, ncol = length(vnames),
    ##                 dimnames = list(NULL, vnames)
    ##                 ),
    ##      auxmetadata = auxmetadata,
    ##      Cout = 'sright'
    ##    ))
    ## )
  }
  ##
  ZiC <- match(vnames, Zvrt)
  ZtC <- which(!is.na(ZiC))
  ZiC <- ZiC[ZtC]
  ZnC <- length(ZiC)

#### Type D
  vnames <- auxmetadata[auxmetadata$mcmctype == 'D', 'name']
  XiD <- match(vnames, Xvrt)
  XtD <- which(!is.na(XiD))
  XiD <- XiD[XtD]
  XnD <- length(XiD)
  ##
  YiD <- match(vnames, Yvrt)
  YtD <- which(!is.na(YiD))
  YiD <- YiD[YtD]
  YnD <- length(YiD)
  if (YnD > 0 || XnD > 0) {
    mcoutput$Dvar <- sqrt(mcoutput$Dvar)
    Dbounds <- cbind(
      auxmetadata[match(vnames, auxmetadata$name), 'tcensormin'],
      ## sign is important here:
      ## for upper tail, take opposite mean and value
      -auxmetadata[match(vnames, auxmetadata$name), 'tcensormax']
      )
    ## Dbounds <- cbind(
    ##   c(vtransform(
    ##     x = matrix(NA,
    ##                nrow = 1, ncol = length(vnames),
    ##                dimnames = list(NULL, vnames)
    ##                ),
    ##     auxmetadata = auxmetadata,
    ##     Dout = 'sleft'
    ##   )),
    ##   ## sign is important here:
    ##   ## for upper tail, take opposite mean and value
    ##   -c(vtransform(
    ##      x = matrix(NA,
    ##                 nrow = 1, ncol = length(vnames),
    ##                 dimnames = list(NULL, vnames)
    ##                 ),
    ##      auxmetadata = auxmetadata,
    ##      Dout = 'sright'
    ##    ))
    ## )
  }
  ##
  ZiD <- match(vnames, Zvrt)
  ZtD <- which(!is.na(ZiD))
  ZiD <- ZiD[ZtD]
  ZnD <- length(ZiD)

#### Type L
  vnames <- auxmetadata[auxmetadata$mcmctype == 'L', 'name']
  XiL <- match(vnames, Xvrt)
  XtL <- which(!is.na(XiL))
  XiL <- XiL[XtL]
  XnL <- length(XiL)
  ##
  YiL <- match(vnames, Yvrt)
  YtL <- which(!is.na(YiL))
  YiL <- YiL[YtL]
  YnL <- length(YiL)
  if (YnL > 0 || XnL > 0) {
    mcoutput$Lvar <- sqrt(mcoutput$Lvar)
    ##
    Lmaxn <- max(auxmetadata[auxmetadata$name %in% vnames, 'Nvalues'])
    Lleft <- t(sapply(vnames, function(avar) {
      nn <- auxmetadata[auxmetadata$name == avar, 'Nvalues']
      seqs <- seq(auxmetadata[auxmetadata$name == avar, 'domainmin'],
                  auxmetadata[auxmetadata$name == avar, 'domainmax'],
                  length.out = nn
                  )
      c(
        vtransform(seqs, auxmetadata = auxmetadata,
                   Lout = 'left',
                   variates = avar),
        rep(NA, Lmaxn - nn)
      )
    }))
    Lright <- t(sapply(vnames, function(avar) {
      nn <- auxmetadata[auxmetadata$name == avar, 'Nvalues']
      seqs <- seq(auxmetadata[auxmetadata$name == avar, 'domainmin'],
                  auxmetadata[auxmetadata$name == avar, 'domainmax'],
                  length.out = nn
                  )
      c(
        vtransform(seqs, auxmetadata = auxmetadata,
                   Lout = 'right',
                   variates = avar),
        rep(NA, Lmaxn - nn)
      )
    }))
  }
  ##
  ZiL <- match(vnames, Zvrt)
  ZtL <- which(!is.na(ZiL))
  ZiL <- ZiL[ZtL]
  ZnL <- length(ZiL)


#### Type O
  vnames <- auxmetadata[auxmetadata$mcmctype == 'O', 'name']
  XiO <- match(vnames, Xvrt)
  XtO <- which(!is.na(XiO))
  XiO <- XiO[XtO]
  XnO <- length(XiO)
  ##
  YiO <- match(vnames, Yvrt)
  YtO <- which(!is.na(YiO))
  YiO <- YiO[YtO]
  YnO <- length(YiO)
  ##
  ZiO <- match(vnames, Zvrt)
  ZtO <- which(!is.na(ZiO))
  ZiO <- ZiO[ZtO]
  ZnO <- length(ZiO)

#### Type N
  vnames <- auxmetadata[auxmetadata$mcmctype == 'N', 'name']
  XiN <- match(vnames, Xvrt)
  XtN <- which(!is.na(XiN))
  XiN <- XiN[XtN]
  XnN <- length(XiN)
  ##
  YiN <- match(vnames, Yvrt)
  YtN <- which(!is.na(YiN))
  YiN <- YiN[YtN]
  YnN <- length(YiN)
  ##
  ZiN <- match(vnames, Zvrt)
  ZtN <- which(!is.na(ZiN))
  ZiN <- ZiN[ZtN]
  ZnN <- length(ZiN)

#### Type B
  vnames <- auxmetadata[auxmetadata$mcmctype == 'B', 'name']
  XiB <- match(vnames, Xvrt)
  XtB <- which(!is.na(XiB))
  XiB <- XiB[XtB]
  XnB <- length(XiB)
  ##
  YiB <- match(vnames, Yvrt)
  YtB <- which(!is.na(YiB))
  YiB <- YiB[YtB]
  YnB <- length(YiB)
  ##
  ZiB <- match(vnames, Zvrt)
  ZtB <- which(!is.na(ZiB))
  ZiB <- ZiB[ZtB]
  ZnB <- length(ZiB)

#### STEP 1. Draw samples of Z (that is, Y,X)

  ## Z is drawn as follows, for each MC sample:
  ## 1. draw a cluster, according to its probability
  ## 2. draw from the appropriate kernel distributions
  ## using the parameters of that cluster
  Ws <- extraDistr::rcat(n=n, prob=t(mcoutput$W))
  Zout <- c(
  (if(ZnR > 0){# continuous
     totake <- cbind(rep(ZtR,each=n), Ws, sseq)
     rnorm(n=n*ZnR,
           mean=mcoutput$Rmean[totake],
           sd=mcoutput$Rvar[totake]
           )
   }else{NULL}),
  (if(ZnC > 0){# censored
     totake <- cbind(rep(ZtC,each=n), Ws, sseq)
     rnorm(n=n*ZnC,
           mean=mcoutput$Rmean[totake],
           sd=mcoutput$Rvar[totake]
           )
   }else{NULL}),
  (if(ZnD > 0){# continuous discretized
     totake <- cbind(rep(ZtD,each=n), Ws, sseq)
     rnorm(n=n*ZnD,
           mean=mcoutput$Rmean[totake],
           sd=mcoutput$Rvar[totake]
           )
   }else{NULL}),
  (if(ZnL > 0){# ordinal
     totake <- cbind(rep(ZtL,each=n), Ws, sseq)
     rnorm(n=n*ZnL,
           mean=mcoutput$Rmean[totake],
           sd=mcoutput$Rvar[totake]
           )
   }else{NULL}),
  (if(ZnO > 0){# nominal
     totake <- cbind(rep(ZtO,each=n), Ws, sseq)
     extraDistr::rcat(n=n*ZnO,
                      prob=apply(mcoutput$Oprob,3,`[`,totake))
   }else{NULL}),
  (if(ZnN > 0){# nominal
     totake <- cbind(rep(ZtN,each=n), Ws, sseq)
     extraDistr::rcat(n=n*ZnN,
                      prob=apply(mcoutput$Nprob,3,`[`,totake))
   }else{NULL}),
  (if(ZnB > 0){# binary
     totake <- cbind(rep(ZtB,each=n), Ws, sseq)
     extraDistr::rbern(n=n*ZnB,
                       prob=mcoutput$Bprob[totake])
   }else{NULL})
  )

  ## rows: n samples, cols: Z variates
  dim(Zout) <- c(n, length(Zvrt))
  ## Match to original order of Zvrt
  Zout <- Zout[, match(Zvrt, Zvrt[c(ZiR, ZiC, ZiD, ZiL, ZiO, ZiN, ZiB)])]
  colnames(Zout) <- Zvrt
  Zout <- vtransform(Zout,
                     auxmetadata = auxmetadata,
                     Rout = 'mi',
                     Cout = 'mi',
                     Dout = 'mi',
                     Lout = 'mi',
                     Oout = 'mi',
                     Nout = 'mi',
                     Bout = 'mi')

  Y2 <- Zout[, Yvrt]
  X2 <- Zout[, Xvrt]
  rm(Zout)
  gc()

#### STEP 2. Calculate sum_i log2_p(Y|X) for all samples

  ## from samplesFDistribution.R with some modifications
  lpYX <- log(foreach(y = t(Y2), x = t(X2),
                           .combine = c,
                           .inorder = TRUE) %dochains% {
#### the loop is over the columns of y and x
#### each instance is a 1-column vector
                             ## rows: clusters, cols: samples
                             probX <- log(mcoutput$W) +
                               (if (XnR > 0) { # continuous
                                  colSums(
                                    dnorm(
                                      x = x[XiR, ],
                                      mean = mcoutput$Rmean[XtR, , , drop = FALSE],
                                      sd = mcoutput$Rvar[XtR, , , drop = FALSE],
                                      log = TRUE
                                    ),
                                    na.rm = TRUE
                                  )
                                } else {
                                  0
                                }) +
                               (if (XnC > 0) { # censored
                                  indf <- which(is.finite(x[XiC, ]))
                                  indi <- which(!is.finite(x[XiC, ]))
                                  (if (length(indf) > 0) {
                                     colSums(
                                       dnorm(
                                         x = x[XiC[indf], ],
                                         mean = mcoutput$Cmean[XtC[indf], , , drop = FALSE],
                                         sd = mcoutput$Cvar[XtC[indf], , , drop = FALSE],
                                         log = TRUE
                                       ),
                                       na.rm = TRUE
                                     )
                                   } else {
                                     0
                                   }) +
                                    (if (length(indi) > 0) {
                                       v2 <- XtC[indi]
                                       v1 <- -sign(x[XiC[indi], ])
                                       ## for upper tail, take opposite mean and value
                                       colSums(
                                         pnorm(
                                           q = Cbounds[cbind(v2, 1.5 - 0.5 * v1)],
                                           mean = v1 * mcoutput$Cmean[v2, , , drop = FALSE],
                                           sd = mcoutput$Cvar[v2, , , drop = FALSE],
                                           log.p = TRUE
                                         ),
                                         na.rm = TRUE
                                       )
                                     } else {
                                       0
                                     })
                                } else {
                                  0
                                }) +
                               (if (XnD > 0) { # continuous discretized
                                  indf <- which(is.finite(x[XiD, ]))
                                  indi <- which(!is.finite(x[XiD, ]))
                                  (if (length(indf) > 0) {
                                     colSums(
                                       dnorm(
                                         x = x[XiD[indf], ],
                                         mean = mcoutput$Dmean[XtD[indf], , , drop = FALSE],
                                         sd = mcoutput$Dvar[XtD[indf], , , drop = FALSE],
                                         log = TRUE
                                       ),
                                       na.rm = TRUE
                                     )
                                   } else {
                                     0
                                   }) +
                                    (if (length(indi) > 0) {
                                       v2 <- XtD[indi]
                                       v1 <- -sign(x[XiD[indi], ])
                                       ## for upper tail, take opposite mean and value
                                       colSums(
                                         pnorm(
                                           q = Dbounds[v2, 1.5 - 0.5 * v1],
                                           mean = v1 * mcoutput$Dmean[v2, , , drop = FALSE],
                                           sd = mcoutput$Dvar[v2, , , drop = FALSE],
                                           log.p = TRUE
                                         ),
                                         na.rm = TRUE
                                       )
                                     } else {
                                       0
                                     })
                                } else {
                                  0
                                }) +
                               (if (XnL > 0) { # ordinal
                                  v2 <- cbind(XtL, x[XiL, ])
                                  colSums(
                                    log(
                                      pnorm(
                                        q = Lright[v2],
                                        mean = mcoutput$Lmean[XtL, , , drop = FALSE],
                                        sd = mcoutput$Lvar[XtL, , , drop = FALSE]
                                      ) -
                                      pnorm(
                                        q = Lleft[v2],
                                        mean = mcoutput$Lmean[XtL, , , drop = FALSE],
                                        sd = mcoutput$Lvar[XtL, , , drop = FALSE]
                                      )
                                    ),
                                    na.rm = TRUE
                                  )
                                } else {
                                  0
                                }) +
                               (if (XnO > 0) { # nominal
                                  colSums(
                                    log(aperm(
                                      vapply(seq_len(XnO), function(v) {
                                        mcoutput$Oprob[XtO[v], , x[XiO[v], ], ]
                                      }, mcoutput$W),
                                      c(3, 1, 2)
                                    )),
                                    na.rm = TRUE
                                  )
                                } else {
                                  0
                                }) +
                               (if (XnN > 0) { # nominal
                                  colSums(
                                    log(aperm(
                                      vapply(seq_len(XnN), function(v) {
                                        mcoutput$Nprob[XtN[v], , x[XiN[v], ], ]
                                      }, mcoutput$W),
                                      c(3, 1, 2)
                                    )),
                                    na.rm = TRUE
                                  )
                                } else {
                                  0
                                }) +
                               (if (XnB > 0) { # binary
                                  colSums(
                                    log(x[XiB, ] * mcoutput$Bprob[XtB, , , drop = FALSE] +
                                        (1L - x[XiB, ]) * (1 - mcoutput$Bprob[XtB, , , drop = FALSE])),
                                    na.rm = TRUE
                                  )
                                } else {
                                  0
                                })
                             ## end probX
                             probX <- apply(probX, 2, function(xx) {
                               xx - max(xx[is.finite(xx)])
                             })
                             ##
                             ##
                             if (all(is.na(y))) {
                               probY <- array(NA, dim = dim(mcoutput$W))
                             } else {
                               probY <- (if (YnR > 0) { # continuous
                                           colSums(
                                             dnorm(
                                               x = y[YiR, ],
                                               mean = mcoutput$Rmean[YtR, , , drop = FALSE],
                                               sd = mcoutput$Rvar[YtR, , , drop = FALSE],
                                               log = TRUE
                                             ),
                                             na.rm = TRUE
                                           )
                                         } else {
                                           0
                                         }) +
                                 (if (YnC > 0) { # censored
                                    indf <- which(is.finite(y[YiC, ]))
                                    indi <- which(!is.finite(y[YiC, ]))
                                    (if (length(indf) > 0) {
                                       colSums(
                                         dnorm(
                                           x = y[YiC[indf], ],
                                           mean = mcoutput$Cmean[YtC[indf], , , drop = FALSE],
                                           sd = mcoutput$Cvar[YtC[indf], , , drop = FALSE],
                                           log = TRUE
                                         ),
                                         na.rm = TRUE
                                       )
                                     } else {
                                       0
                                     }) +
                                      (if (length(indi) > 0) {
                                         v2 <- YtC[indi]
                                         v1 <- -sign(y[YiC[indi], ])
                                         ## for upper tail, take opposite mean and value
                                         colSums(
                                           pnorm(
                                             q = Cbounds[cbind(v2, 1.5 - 0.5 * v1)],
                                             mean = v1 * mcoutput$Cmean[v2, , , drop = FALSE],
                                             sd = mcoutput$Cvar[v2, , , drop = FALSE],
                                             log.p = TRUE
                                           ),
                                           na.rm = TRUE
                                         )
                                       } else {
                                         0
                                       })
                                  } else {
                                    0
                                  }) +
                                 (if (YnD > 0) { # continuous discretized
                                    indf <- which(is.finite(y[YiD, ]))
                                    indi <- which(!is.finite(y[YiD, ]))
                                    (if (length(indf) > 0) {
                                       colSums(
                                         dnorm(
                                           x = y[YiD[indf], ],
                                           mean = mcoutput$Dmean[YtD[indf], , , drop = FALSE],
                                           sd = mcoutput$Dvar[YtD[indf], , , drop = FALSE],
                                           log = TRUE
                                         ),
                                         na.rm = TRUE
                                       )
                                     } else {
                                       0
                                     }) +
                                      (if (length(indi) > 0) {
                                         v2 <- YtD[indi]
                                         v1 <- -sign(y[YiD[indi], ])
                                         ## for upper tail, take opposite mean and value
                                         colSums(
                                           pnorm(
                                             q = Dbounds[v2, 1.5 - 0.5 * v1],
                                             mean = v1 * mcoutput$Dmean[v2, , , drop = FALSE],
                                             sd = mcoutput$Dvar[v2, , , drop = FALSE],
                                             log.p = TRUE
                                           ),
                                           na.rm = TRUE
                                         )
                                       } else {
                                         0
                                       })
                                  } else {
                                    0
                                  }) +
                                 (if (YnL > 0) { # ordinal
                                    v2 <- cbind(YtL, y[YiL, ])
                                    colSums(
                                      log(
                                        pnorm(
                                          q = Lright[v2],
                                          mean = mcoutput$Lmean[YtL, , , drop = FALSE],
                                          sd = mcoutput$Lvar[YtL, , , drop = FALSE]
                                        ) -
                                        pnorm(
                                          q = Lleft[v2],
                                          mean = mcoutput$Lmean[YtL, , , drop = FALSE],
                                          sd = mcoutput$Lvar[YtL, , , drop = FALSE]
                                        )
                                      ),
                                      na.rm = TRUE
                                    )
                                  } else {
                                    0
                                  }) +
                                 (if (YnO > 0) { # nominal
                                    colSums(
                                      log(aperm(
                                        vapply(seq_len(YnO), function(v) {
                                          mcoutput$Oprob[YtO[v], , y[YiO[v], ], ]
                                        }, mcoutput$W),
                                        c(3, 1, 2)
                                      )),
                                      na.rm = TRUE
                                    )
                                  } else {
                                    0
                                  }) +
                                 (if (YnN > 0) { # nominal
                                    colSums(
                                      log(aperm(
                                        vapply(seq_len(YnN), function(v) {
                                          mcoutput$Nprob[YtN[v], , y[YiN[v], ], ]
                                        }, mcoutput$W),
                                        c(3, 1, 2)
                                      )),
                                      na.rm = TRUE
                                    )
                                  } else {
                                    0
                                  }) +
                                 (if (YnB > 0) { # binary
                                    colSums(
                                      log((y[YiB, ] * mcoutput$Bprob[YtB, , , drop = FALSE]) +
                                          ((1 - y[YiB, ]) * (1 - mcoutput$Bprob[YtB, , , drop = FALSE]))),
                                      na.rm = TRUE
                                    )
                                  } else {
                                    0
                                  })
                             }
                             ## Output: rows=clusters, columns=samples
                            ## temp <- colSums(exp(probX + probY))/colSums(exp(probX))
                            ## c(mean(temp), 1+sd(temp)/length(temp)/mean(temp))
                            mean(colSums(exp(probX + probY))/colSums(exp(probX)))
                           }, base = base)

#### STEP 3. Calculate sum_i log2_p(Y) for all samples

  ## from samplesFDistribution.R with some modifications
  lpY <- log(foreach(y = t(Y2),
                          .combine = c,
                          .inorder = TRUE) %dochains% {
#### the loop is over the columns of y and x
#### each instance is a 1-column vector
                            probX <- apply(log(mcoutput$W), 2, function(xx) {
                              xx - max(xx[is.finite(xx)])
                            })
                            ##
                            ##
                            if (all(is.na(y))) {
                              probY <- array(NA, dim = dim(mcoutput$W))
                            } else {
                              probY <- (if (YnR > 0) { # continuous
                                          colSums(
                                            dnorm(
                                              x = y[YiR, ],
                                              mean = mcoutput$Rmean[YtR, , , drop = FALSE],
                                              sd = mcoutput$Rvar[YtR, , , drop = FALSE],
                                              log = TRUE
                                            ),
                                            na.rm = TRUE
                                          )
                                        } else {
                                          0
                                        }) +
                                (if (YnC > 0) { # censored
                                   indf <- which(is.finite(y[YiC, ]))
                                   indi <- which(!is.finite(y[YiC, ]))
                                   (if (length(indf) > 0) {
                                      colSums(
                                        dnorm(
                                          x = y[YiC[indf], ],
                                          mean = mcoutput$Cmean[YtC[indf], , , drop = FALSE],
                                          sd = mcoutput$Cvar[YtC[indf], , , drop = FALSE],
                                          log = TRUE
                                        ),
                                        na.rm = TRUE
                                      )
                                    } else {
                                      0
                                    }) +
                                     (if (length(indi) > 0) {
                                        v2 <- YtC[indi]
                                        v1 <- -sign(y[YiC[indi], ])
                                        ## for upper tail, take opposite mean and value
                                        colSums(
                                          pnorm(
                                            q = Cbounds[cbind(v2, 1.5 - 0.5 * v1)],
                                            mean = v1 * mcoutput$Cmean[v2, , , drop = FALSE],
                                            sd = mcoutput$Cvar[v2, , , drop = FALSE],
                                            log.p = TRUE
                                          ),
                                          na.rm = TRUE
                                        )
                                      } else {
                                        0
                                      })
                                 } else {
                                   0
                                 }) +
                                (if (YnD > 0) { # continuous discretized
                                   indf <- which(is.finite(y[YiD, ]))
                                   indi <- which(!is.finite(y[YiD, ]))
                                   (if (length(indf) > 0) {
                                      colSums(
                                        dnorm(
                                          x = y[YiD[indf], ],
                                          mean = mcoutput$Dmean[YtD[indf], , , drop = FALSE],
                                          sd = mcoutput$Dvar[YtD[indf], , , drop = FALSE],
                                          log = TRUE
                                        ),
                                        na.rm = TRUE
                                      )
                                    } else {
                                      0
                                    }) +
                                     (if (length(indi) > 0) {
                                        v2 <- YtD[indi]
                                        v1 <- -sign(y[YiD[indi], ])
                                        ## for upper tail, take opposite mean and value
                                        colSums(
                                          pnorm(
                                            q = Dbounds[v2, 1.5 - 0.5 * v1],
                                            mean = v1 * mcoutput$Dmean[v2, , , drop = FALSE],
                                            sd = mcoutput$Dvar[v2, , , drop = FALSE],
                                            log.p = TRUE
                                          ),
                                          na.rm = TRUE
                                        )
                                      } else {
                                        0
                                      })
                                 } else {
                                   0
                                 }) +
                                (if (YnL > 0) { # ordinal
                                   v2 <- cbind(YtL, y[YiL, ])
                                   colSums(
                                     log(
                                       pnorm(
                                         q = Lright[v2],
                                         mean = mcoutput$Lmean[YtL, , , drop = FALSE],
                                         sd = mcoutput$Lvar[YtL, , , drop = FALSE]
                                       ) -
                                       pnorm(
                                         q = Lleft[v2],
                                         mean = mcoutput$Lmean[YtL, , , drop = FALSE],
                                         sd = mcoutput$Lvar[YtL, , , drop = FALSE]
                                       )
                                     ),
                                     na.rm = TRUE
                                   )
                                 } else {
                                   0
                                 }) +
                                (if (YnO > 0) { # nominal
                                   colSums(
                                     log(aperm(
                                       vapply(seq_len(YnO), function(v) {
                                         mcoutput$Oprob[YtO[v], , y[YiO[v], ], ]
                                       }, mcoutput$W),
                                       c(3, 1, 2)
                                     )),
                                     na.rm = TRUE
                                   )
                                 } else {
                                   0
                                 }) +
                                (if (YnN > 0) { # nominal
                                   colSums(
                                     log(aperm(
                                       vapply(seq_len(YnN), function(v) {
                                         mcoutput$Nprob[YtN[v], , y[YiN[v], ], ]
                                       }, mcoutput$W),
                                       c(3, 1, 2)
                                     )),
                                     na.rm = TRUE
                                   )
                                 } else {
                                   0
                                 }) +
                                (if (YnB > 0) { # binary
                                   colSums(
                                     log((y[YiB, ] * mcoutput$Bprob[YtB, , , drop = FALSE]) +
                                         ((1 - y[YiB, ]) * (1 - mcoutput$Bprob[YtB, , , drop = FALSE]))),
                                     na.rm = TRUE
                                   )
                                 } else {
                                   0
                                 })
                            }
                            ## Output: rows=clusters, columns=samples
                            ## temp <- colSums(exp(probX + probY))/colSums(exp(probX))
                            ## c(mean(temp), 1+sd(temp)/length(temp)/mean(temp))
                            mean(colSums(exp(probX + probY))/colSums(exp(probX)))
                          }, base = base)

  if (!silent && exists('cl')) {
    cat('\nClosing connections to cores.\n')
    parallel::stopCluster(cl)
  }

  ## error <- sd(lpYX[,1] - lpY[,1])/sqrt(nrow(lpYX))
  ## error2 <- mean(lpYX[,2] + lpY[,2])

  ## ## ***todo: transform variates to original space & calculate Jacobian***
  ## ljac <- -rowSums(
  ##           log(vtransform(Y, auxmetadata = auxmetadata, invjacobian = TRUE)),
  ##           na.rm = TRUE
  ##          )


  ## ## ***todo: also output H(Y|X) and H(Y); these needs the Jacobian***
  ## list(MI = lpYX - lpY,
  ##      condH = -lpYX,
  ##      H = -lpY)

  list(MI = mean(lpYX - lpY),
       error = sd(lpYX - lpY)/sqrt(n),
       unit = unit)
}

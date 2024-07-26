#' Calculate mutual information between groups of joint variates
#'
#' @param Y1vrt String vector: first group of joint variates
#' @param Y2vrt String vector: second group of joint variates
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
mutualinfo <- function(Y1vrt, Y2vrt, X = NULL, mcoutput,
                       nsamples = 3600,
                       unit = 'Sh',
                       parallel = TRUE,
                       silent = FALSE){

#### Mutual information and conditional entropy between Y2 and Y1
#### conditional on X, data, prior information
#### are calculated by Monte Carlo integration:
#### 0. adjusted cluster weights are calculated for conditioning on X
#### 1. joint samples of Y1_i, Y2_i are drawn
#### 2. probabilities p(Y1|Y2) are calculated for each sample
####    the conditional entropy is Monte-Carlo approximated by
####    H(Y1|Y2) = - sum_{i} log p(Y1_i | Y2_i)
#### 3. probabilities p(Y1|Y2) are calculated for each sample
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


  if(!is.character(Y1vrt) || any(is.na(Y1vrt))){
    stop('Y1vrt must be a vector of variate names')
  }
  if(!is.character(Y2vrt) || any(is.na(Y2vrt))){
    stop('Y2vrt must be a vector of variate names')
  }
  if (!is.null(X) && length(dim(X)) != 2) {
    stop('X must be NULL or have two dimensions')
  }
  if (!is.null(X) && dim(X)[1] > 1) {
    message('Only the first row of X is considered')
    X <- X[1, , drop = FALSE]
  }

  ## More consistency checks
  if(!all(Y1vrt %in% auxmetadata$name)) {
    stop('unknown Y1 variates\n')
  }
  if(length(unique(Y1vrt)) != length(Y1vrt)) {
    stop('duplicate Y1 variates\n')
  }
  ##
  if(!all(Y2vrt %in% auxmetadata$name)){
    stop('unknown Y2 variates\n')
  }
  if(length(unique(Y2vrt)) != length(Y2vrt)) {
    stop('duplicate Y2 variates\n')
  }
  ##
  Xvrt <- colnames(X)
  if (!all(Xvrt %in% auxmetadata$name)) {
    stop('unknown X variates\n')
  }
  if (length(unique(Xvrt)) != length(Xvrt)) {
    stop('duplicate X variates\n')
  }
  ##
  if(length(intersect(Y1vrt, Y2vrt)) > 0) {
    stop('overlap in Y1 and Y2 variates\n')
  }
  if(length(intersect(Y1vrt, Xvrt)) > 0) {
    stop('overlap in Y1 and X variates\n')
  }
  if(length(intersect(Y2vrt, Xvrt)) > 0) {
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

#### Combine Y1,Y2 into single Z for speed
  Zvrt <- c(Y1vrt, Y2vrt)

#### Type R
  vnames <- auxmetadata[auxmetadata$mcmctype == 'R', 'name']
  Y2iR <- match(vnames, Y2vrt)
  Y2tR <- which(!is.na(Y2iR))
  Y2iR <- Y2iR[Y2tR]
  Y2nR <- length(Y2iR)
  ##
  Y1iR <- match(vnames, Y1vrt)
  Y1tR <- which(!is.na(Y1iR))
  Y1iR <- Y1iR[Y1tR]
  Y1nR <- length(Y1iR)
  ##
  XiR <- match(vnames, Xvrt)
  XtR <- which(!is.na(XiR))
  XiR <- XiR[XtR]
  XnR <- length(XiR)
  if (Y1nR > 0 || Y2nR > 0 || XnR > 0) {
    mcoutput$Rvar <- sqrt(mcoutput$Rvar)
  }
  ##
  ZiR <- match(vnames, Zvrt)
  ZtR <- which(!is.na(ZiR))
  ZiR <- ZiR[ZtR]
  ZnR <- length(ZiR)

#### Type C
  vnames <- auxmetadata[auxmetadata$mcmctype == 'C', 'name']
  Y2iC <- match(vnames, Y2vrt)
  Y2tC <- which(!is.na(Y2iC))
  Y2iC <- Y2iC[Y2tC]
  Y2nC <- length(Y2iC)
  ##
  Y1iC <- match(vnames, Y1vrt)
  Y1tC <- which(!is.na(Y1iC))
  Y1iC <- Y1iC[Y1tC]
  Y1nC <- length(Y1iC)
  ##
  XiC <- match(vnames, Xvrt)
  XtC <- which(!is.na(XiC))
  XiC <- XiC[XtC]
  XnC <- length(XiC)
  if (Y1nC > 0 || Y2nC > 0 || XnC > 0) {
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
  Y2iD <- match(vnames, Y2vrt)
  Y2tD <- which(!is.na(Y2iD))
  Y2iD <- Y2iD[Y2tD]
  Y2nD <- length(Y2iD)
  ##
  Y1iD <- match(vnames, Y1vrt)
  Y1tD <- which(!is.na(Y1iD))
  Y1iD <- Y1iD[Y1tD]
  Y1nD <- length(Y1iD)
  ##
  XiD <- match(vnames, Xvrt)
  XtD <- which(!is.na(XiD))
  XiD <- XiD[XtD]
  XnD <- length(XiD)
  if (Y1nD > 0 || Y2nD > 0 || XnD > 0) {
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
  Y2iL <- match(vnames, Y2vrt)
  Y2tL <- which(!is.na(Y2iL))
  Y2iL <- Y2iL[Y2tL]
  Y2nL <- length(Y2iL)
  ##
  Y1iL <- match(vnames, Y1vrt)
  Y1tL <- which(!is.na(Y1iL))
  Y1iL <- Y1iL[Y1tL]
  Y1nL <- length(Y1iL)
  ##
  XiL <- match(vnames, Xvrt)
  XtL <- which(!is.na(XiL))
  XiL <- XiL[XtL]
  XnL <- length(XiL)
  if (Y1nL > 0 || Y2nL > 0 || XnL > 0) {
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
  Y2iO <- match(vnames, Y2vrt)
  Y2tO <- which(!is.na(Y2iO))
  Y2iO <- Y2iO[Y2tO]
  Y2nO <- length(Y2iO)
  ##
  Y1iO <- match(vnames, Y1vrt)
  Y1tO <- which(!is.na(Y1iO))
  Y1iO <- Y1iO[Y1tO]
  Y1nO <- length(Y1iO)
  ##
  XiO <- match(vnames, Xvrt)
  XtO <- which(!is.na(XiO))
  XiO <- XiO[XtO]
  XnO <- length(XiO)
  ##
  ZiO <- match(vnames, Zvrt)
  ZtO <- which(!is.na(ZiO))
  ZiO <- ZiO[ZtO]
  ZnO <- length(ZiO)

#### Type N
  vnames <- auxmetadata[auxmetadata$mcmctype == 'N', 'name']
  Y2iN <- match(vnames, Y2vrt)
  Y2tN <- which(!is.na(Y2iN))
  Y2iN <- Y2iN[Y2tN]
  Y2nN <- length(Y2iN)
  ##
  Y1iN <- match(vnames, Y1vrt)
  Y1tN <- which(!is.na(Y1iN))
  Y1iN <- Y1iN[Y1tN]
  Y1nN <- length(Y1iN)
  ##
  XiN <- match(vnames, Xvrt)
  XtN <- which(!is.na(XiN))
  XiN <- XiN[XtN]
  XnN <- length(XiN)
  ##
  ZiN <- match(vnames, Zvrt)
  ZtN <- which(!is.na(ZiN))
  ZiN <- ZiN[ZtN]
  ZnN <- length(ZiN)

#### Type B
  vnames <- auxmetadata[auxmetadata$mcmctype == 'B', 'name']
  Y2iB <- match(vnames, Y2vrt)
  Y2tB <- which(!is.na(Y2iB))
  Y2iB <- Y2iB[Y2tB]
  Y2nB <- length(Y2iB)
  ##
  Y1iB <- match(vnames, Y1vrt)
  Y1tB <- which(!is.na(Y1iB))
  Y1iB <- Y1iB[Y1tB]
  Y1nB <- length(Y1iB)
  ##
  XiB <- match(vnames, Xvrt)
  XtB <- which(!is.na(XiB))
  XiB <- XiB[XtB]
  XnB <- length(XiB)
  ##
  ZiB <- match(vnames, Zvrt)
  ZtB <- which(!is.na(ZiB))
  ZiB <- ZiB[ZtB]
  ZnB <- length(ZiB)

#### STEP 0. Adjust cluster weights W for conditioning on X
  if(is.null(X)){
    W <- mcoutput$W
  } else {
    X <- t(vtransform(X, auxmetadata = auxmetadata,
                     Rout = 'normalized',
                     Cout = 'boundisinf',
                     Dout = 'boundisinf',
                     Lout = 'normalized',
                     Oout = 'numeric',
                     Nout = 'numeric',
                     Bout = 'numeric'))
    ##
    W <- log(mcoutput$W) +
      (if (XnR > 0) { # continuous
         colSums(
           dnorm(
             x = X[XiR, ],
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
         indf <- which(is.finite(X[XiC, ]))
         indi <- which(!is.finite(X[XiC, ]))
         (if (length(indf) > 0) {
            colSums(
              dnorm(
                x = X[XiC[indf], ],
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
              v1 <- -sign(X[XiC[indi], ])
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
         indf <- which(is.finite(X[XiD, ]))
         indi <- which(!is.finite(X[XiD, ]))
         (if (length(indf) > 0) {
            colSums(
              dnorm(
                x = X[XiD[indf], ],
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
              v1 <- -sign(X[XiD[indi], ])
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
         v2 <- cbind(XtL, X[XiL, ])
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
               mcoutput$Oprob[XtO[v], , X[XiO[v], ], ]
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
               mcoutput$Nprob[XtN[v], , X[XiN[v], ], ]
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
           log(X[XiB, ] * mcoutput$Bprob[XtB, , , drop = FALSE] +
               (1L - X[XiB, ]) * (1 - mcoutput$Bprob[XtB, , , drop = FALSE])),
           na.rm = TRUE
         )
       } else {
         0
       })
    ##
    W <- apply(W, 2, function(xx) {
      xx <- xx - max(xx[is.finite(xx)])
      exp(xx)/sum(exp(xx))
    })
  } # end definition of W if non-null X

#### STEP 1. Draw samples of Z (that is, Y1,Y2)

  ## Z is drawn as follows, for each MC sample:
  ## 1. draw a cluster, according to its probability
  ## 2. draw from the appropriate kernel distributions
  ## using the parameters of that cluster

  Ws <- extraDistr::rcat(n=n, prob=t(W))
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

  Y1transf <- Zout[, Y1vrt]
  Y2transf <- Zout[, Y2vrt]
  rm(Zout)
  gc()

#### STEP 2. Calculate sum_i log2_p(Y1|Y2) for all samples

  ## from samplesFDistribution.R with some modifications
  lpYX <- log(foreach(y = t(Y1transf), x = t(Y2transf),
                      .combine = c,
                      .inorder = TRUE) %dochains% {
#### the loop is over the columns of y and x
#### each instance is a 1-column vector
                        ## rows: clusters, cols: samples
                        probX <- log(W) +
                          (if (Y2nR > 0) { # continuous
                             colSums(
                               dnorm(
                                 x = x[Y2iR, ],
                                 mean = mcoutput$Rmean[Y2tR, , , drop = FALSE],
                                 sd = mcoutput$Rvar[Y2tR, , , drop = FALSE],
                                 log = TRUE
                               ),
                               na.rm = TRUE
                             )
                           } else {
                             0
                           }) +
                          (if (Y2nC > 0) { # censored
                             indf <- which(is.finite(x[Y2iC, ]))
                             indi <- which(!is.finite(x[Y2iC, ]))
                             (if (length(indf) > 0) {
                                colSums(
                                  dnorm(
                                    x = x[Y2iC[indf], ],
                                    mean = mcoutput$Cmean[Y2tC[indf], , , drop = FALSE],
                                    sd = mcoutput$Cvar[Y2tC[indf], , , drop = FALSE],
                                    log = TRUE
                                  ),
                                  na.rm = TRUE
                                )
                              } else {
                                0
                              }) +
                               (if (length(indi) > 0) {
                                  v2 <- Y2tC[indi]
                                  v1 <- -sign(x[Y2iC[indi], ])
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
                          (if (Y2nD > 0) { # continuous discretized
                             indf <- which(is.finite(x[Y2iD, ]))
                             indi <- which(!is.finite(x[Y2iD, ]))
                             (if (length(indf) > 0) {
                                colSums(
                                  dnorm(
                                    x = x[Y2iD[indf], ],
                                    mean = mcoutput$Dmean[Y2tD[indf], , , drop = FALSE],
                                    sd = mcoutput$Dvar[Y2tD[indf], , , drop = FALSE],
                                    log = TRUE
                                  ),
                                  na.rm = TRUE
                                )
                              } else {
                                0
                              }) +
                               (if (length(indi) > 0) {
                                  v2 <- Y2tD[indi]
                                  v1 <- -sign(x[Y2iD[indi], ])
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
                          (if (Y2nL > 0) { # ordinal
                             v2 <- cbind(Y2tL, x[Y2iL, ])
                             colSums(
                               log(
                                 pnorm(
                                   q = Lright[v2],
                                   mean = mcoutput$Lmean[Y2tL, , , drop = FALSE],
                                   sd = mcoutput$Lvar[Y2tL, , , drop = FALSE]
                                 ) -
                                 pnorm(
                                   q = Lleft[v2],
                                   mean = mcoutput$Lmean[Y2tL, , , drop = FALSE],
                                   sd = mcoutput$Lvar[Y2tL, , , drop = FALSE]
                                 )
                               ),
                               na.rm = TRUE
                             )
                           } else {
                             0
                           }) +
                          (if (Y2nO > 0) { # nominal
                             colSums(
                               log(aperm(
                                 vapply(seq_len(Y2nO), function(v) {
                                   mcoutput$Oprob[Y2tO[v], , x[Y2iO[v], ], ]
                                 }, W),
                                 c(3, 1, 2)
                               )),
                               na.rm = TRUE
                             )
                           } else {
                             0
                           }) +
                          (if (Y2nN > 0) { # nominal
                             colSums(
                               log(aperm(
                                 vapply(seq_len(Y2nN), function(v) {
                                   mcoutput$Nprob[Y2tN[v], , x[Y2iN[v], ], ]
                                 }, W),
                                 c(3, 1, 2)
                               )),
                               na.rm = TRUE
                             )
                           } else {
                             0
                           }) +
                          (if (Y2nB > 0) { # binary
                             colSums(
                               log(x[Y2iB, ] * mcoutput$Bprob[Y2tB, , , drop = FALSE] +
                                   (1L - x[Y2iB, ]) * (1 - mcoutput$Bprob[Y2tB, , , drop = FALSE])),
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
                          probY <- array(NA, dim = dim(W))
                        } else {
                          probY <- (if (Y1nR > 0) { # continuous
                                      colSums(
                                        dnorm(
                                          x = y[Y1iR, ],
                                          mean = mcoutput$Rmean[Y1tR, , , drop = FALSE],
                                          sd = mcoutput$Rvar[Y1tR, , , drop = FALSE],
                                          log = TRUE
                                        ),
                                        na.rm = TRUE
                                      )
                                    } else {
                                      0
                                    }) +
                            (if (Y1nC > 0) { # censored
                               indf <- which(is.finite(y[Y1iC, ]))
                               indi <- which(!is.finite(y[Y1iC, ]))
                               (if (length(indf) > 0) {
                                  colSums(
                                    dnorm(
                                      x = y[Y1iC[indf], ],
                                      mean = mcoutput$Cmean[Y1tC[indf], , , drop = FALSE],
                                      sd = mcoutput$Cvar[Y1tC[indf], , , drop = FALSE],
                                      log = TRUE
                                    ),
                                    na.rm = TRUE
                                  )
                                } else {
                                  0
                                }) +
                                 (if (length(indi) > 0) {
                                    v2 <- Y1tC[indi]
                                    v1 <- -sign(y[Y1iC[indi], ])
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
                            (if (Y1nD > 0) { # continuous discretized
                               indf <- which(is.finite(y[Y1iD, ]))
                               indi <- which(!is.finite(y[Y1iD, ]))
                               (if (length(indf) > 0) {
                                  colSums(
                                    dnorm(
                                      x = y[Y1iD[indf], ],
                                      mean = mcoutput$Dmean[Y1tD[indf], , , drop = FALSE],
                                      sd = mcoutput$Dvar[Y1tD[indf], , , drop = FALSE],
                                      log = TRUE
                                    ),
                                    na.rm = TRUE
                                  )
                                } else {
                                  0
                                }) +
                                 (if (length(indi) > 0) {
                                    v2 <- Y1tD[indi]
                                    v1 <- -sign(y[Y1iD[indi], ])
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
                            (if (Y1nL > 0) { # ordinal
                               v2 <- cbind(Y1tL, y[Y1iL, ])
                               colSums(
                                 log(
                                   pnorm(
                                     q = Lright[v2],
                                     mean = mcoutput$Lmean[Y1tL, , , drop = FALSE],
                                     sd = mcoutput$Lvar[Y1tL, , , drop = FALSE]
                                   ) -
                                   pnorm(
                                     q = Lleft[v2],
                                     mean = mcoutput$Lmean[Y1tL, , , drop = FALSE],
                                     sd = mcoutput$Lvar[Y1tL, , , drop = FALSE]
                                   )
                                 ),
                                 na.rm = TRUE
                               )
                             } else {
                               0
                             }) +
                            (if (Y1nO > 0) { # nominal
                               colSums(
                                 log(aperm(
                                   vapply(seq_len(Y1nO), function(v) {
                                     mcoutput$Oprob[Y1tO[v], , y[Y1iO[v], ], ]
                                   }, W),
                                   c(3, 1, 2)
                                 )),
                                 na.rm = TRUE
                               )
                             } else {
                               0
                             }) +
                            (if (Y1nN > 0) { # nominal
                               colSums(
                                 log(aperm(
                                   vapply(seq_len(Y1nN), function(v) {
                                     mcoutput$Nprob[Y1tN[v], , y[Y1iN[v], ], ]
                                   }, W),
                                   c(3, 1, 2)
                                 )),
                                 na.rm = TRUE
                               )
                             } else {
                               0
                             }) +
                            (if (Y1nB > 0) { # binary
                               colSums(
                                 log((y[Y1iB, ] * mcoutput$Bprob[Y1tB, , , drop = FALSE]) +
                                     ((1 - y[Y1iB, ]) * (1 - mcoutput$Bprob[Y1tB, , , drop = FALSE]))),
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

#### STEP 3. Calculate sum_i log2_p(Y1) for all samples

  ## from samplesFDistribution.R with some modifications
  lpY <- log(foreach(y = t(Y1transf),
                     .combine = c,
                     .inorder = TRUE) %dochains% {
#### the loop is over the columns of y and x
#### each instance is a 1-column vector
                       probX <- apply(log(W), 2, function(xx) {
                         xx - max(xx[is.finite(xx)])
                       })
                       ##
                       ##
                       if (all(is.na(y))) {
                         probY <- array(NA, dim = dim(W))
                       } else {
                         probY <- (if (Y1nR > 0) { # continuous
                                     colSums(
                                       dnorm(
                                         x = y[Y1iR, ],
                                         mean = mcoutput$Rmean[Y1tR, , , drop = FALSE],
                                         sd = mcoutput$Rvar[Y1tR, , , drop = FALSE],
                                         log = TRUE
                                       ),
                                       na.rm = TRUE
                                     )
                                   } else {
                                     0
                                   }) +
                           (if (Y1nC > 0) { # censored
                              indf <- which(is.finite(y[Y1iC, ]))
                              indi <- which(!is.finite(y[Y1iC, ]))
                              (if (length(indf) > 0) {
                                 colSums(
                                   dnorm(
                                     x = y[Y1iC[indf], ],
                                     mean = mcoutput$Cmean[Y1tC[indf], , , drop = FALSE],
                                     sd = mcoutput$Cvar[Y1tC[indf], , , drop = FALSE],
                                     log = TRUE
                                   ),
                                   na.rm = TRUE
                                 )
                               } else {
                                 0
                               }) +
                                (if (length(indi) > 0) {
                                   v2 <- Y1tC[indi]
                                   v1 <- -sign(y[Y1iC[indi], ])
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
                           (if (Y1nD > 0) { # continuous discretized
                              indf <- which(is.finite(y[Y1iD, ]))
                              indi <- which(!is.finite(y[Y1iD, ]))
                              (if (length(indf) > 0) {
                                 colSums(
                                   dnorm(
                                     x = y[Y1iD[indf], ],
                                     mean = mcoutput$Dmean[Y1tD[indf], , , drop = FALSE],
                                     sd = mcoutput$Dvar[Y1tD[indf], , , drop = FALSE],
                                     log = TRUE
                                   ),
                                   na.rm = TRUE
                                 )
                               } else {
                                 0
                               }) +
                                (if (length(indi) > 0) {
                                   v2 <- Y1tD[indi]
                                   v1 <- -sign(y[Y1iD[indi], ])
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
                           (if (Y1nL > 0) { # ordinal
                              v2 <- cbind(Y1tL, y[Y1iL, ])
                              colSums(
                                log(
                                  pnorm(
                                    q = Lright[v2],
                                    mean = mcoutput$Lmean[Y1tL, , , drop = FALSE],
                                    sd = mcoutput$Lvar[Y1tL, , , drop = FALSE]
                                  ) -
                                  pnorm(
                                    q = Lleft[v2],
                                    mean = mcoutput$Lmean[Y1tL, , , drop = FALSE],
                                    sd = mcoutput$Lvar[Y1tL, , , drop = FALSE]
                                  )
                                ),
                                na.rm = TRUE
                              )
                            } else {
                              0
                            }) +
                           (if (Y1nO > 0) { # nominal
                              colSums(
                                log(aperm(
                                  vapply(seq_len(Y1nO), function(v) {
                                    mcoutput$Oprob[Y1tO[v], , y[Y1iO[v], ], ]
                                  }, W),
                                  c(3, 1, 2)
                                )),
                                na.rm = TRUE
                              )
                            } else {
                              0
                            }) +
                           (if (Y1nN > 0) { # nominal
                              colSums(
                                log(aperm(
                                  vapply(seq_len(Y1nN), function(v) {
                                    mcoutput$Nprob[Y1tN[v], , y[Y1iN[v], ], ]
                                  }, W),
                                  c(3, 1, 2)
                                )),
                                na.rm = TRUE
                              )
                            } else {
                              0
                            }) +
                           (if (Y1nB > 0) { # binary
                              colSums(
                                log((y[Y1iB, ] * mcoutput$Bprob[Y1tB, , , drop = FALSE]) +
                                    ((1 - y[Y1iB, ]) * (1 - mcoutput$Bprob[Y1tB, , , drop = FALSE]))),
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
  ##           log(vtransform(Y1, auxmetadata = auxmetadata, invjacobian = TRUE)),
  ##           na.rm = TRUE
  ##          )


  ## ## ***todo: also output H(Y1|Y2) and H(Y1); these needs the Jacobian***
  ## list(MI = lpYX - lpY,
  ##      condH = -lpYX,
  ##      H = -lpY)

  list(MI = mean(lpYX - lpY),
       error = sd(lpYX - lpY)/sqrt(n),
       unit = unit)
}

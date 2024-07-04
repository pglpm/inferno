#### In progress
mutualinfo <- function(Yvrt, Xvrt, mcoutput, n=NULL, useOquantiles=TRUE, parallel=TRUE, silent=FALSE){

#### Mutual information and conditional entropy between X and Y
#### are calculated by Monte Carlo integration:
#### 1. joint samples of Y_i, X_i are drawn
#### 2. probabilities p(Y|X) and p(Y) are calculated for each sample
#### 3. the conditional entropy is Monte-Carlo approximated by
####    H(Y|X) = sum_{i} log p(Y_i | X_i)
#### 4. the mutual info is Monte-Carlo approximated by
####    I(Y|X) = sum_{i} [log p(Y_i | X_i) - log p(Y_i)]
####           = H(Y|X) - sum_{i} log p(Y_i)
#### (we also obtain the entropy of Y for free.)
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
      cl <- makeCluster(parallel)
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
    `%dochains%` <- `%dorng%`
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
  n <- ceiling(nsamples/nMCsamples)
  sseq <- seq_len(n)

  ## source('vtransform.R')

#### Combine Y,X into single Z for speed
  Zvrt <- c(Yvrt, Xvrt)

  ## Type R
  vnames <- auxmetadata[mcmctype == 'R', name]
  XiR <- match(vnames, Xvrt)
  XtR <- which(!is.na(XiR))
  XiR <- XiR[XtR]
  XnR <- length(XiR)
  ##
  YiR <- match(vnames, Yvrt)
  YtR <- which(!is.na(YiR))
  YiR <- YiR[YtR]
  YnR <- length(YiR)
  if(YnR > 0 || XnR > 0){
    mcoutput$Rvar <- sqrt(mcoutput$Rvar)
  }
  ##
  ZiR <- match(vnames, Zvrt)
  ZtR <- which(!is.na(ZiR))
  ZiR <- ZiR[ZtR]
  ZnR <- length(ZiR)

  ## Type C
  vnames <- auxmetadata[mcmctype == 'C', name]
  XiC <- match(vnames, Xvrt)
  XtC <- which(!is.na(XiC))
  XiC <- XiC[XtC]
  XnC <- length(XiC)
  ##
  YiC <- match(vnames, Yvrt)
  YtC <- which(!is.na(YiC))
  YiC <- YiC[YtC]
  YnC <- length(YiC)
  if(YnC > 0 || XnC > 0){
    mcoutput$Cvar <- sqrt(mcoutput$Cvar)
    Cbounds <- cbind(
      c(vtransform(x=matrix(NA,nrow=1,ncol=length(vnames),
                            dimnames=list(NULL,vnames)),
                   auxmetadata=auxmetadata, Cout='sleft',
                   useOquantiles=useOquantiles)),
      ## sign is important here:
      ## for upper tail, take opposite mean and value
      -c(vtransform(x=matrix(NA,nrow=1,ncol=length(vnames),
                             dimnames=list(NULL,vnames)),
                    auxmetadata=auxmetadata, Cout='sright',
                    useOquantiles=useOquantiles))
    )
  }
  ##
  ZiC <- match(vnames, Zvrt)
  ZtC <- which(!is.na(ZiC))
  ZiC <- ZiC[ZtC]
  ZnC <- length(ZiC)

  ## Type D
  vnames <- auxmetadata[mcmctype == 'D', name]
  XiD <- match(vnames, Xvrt)
  XtD <- which(!is.na(XiD))
  XiD <- XiD[XtD]
  XnD <- length(XiD)
  ##
  YiD <- match(vnames, Yvrt)
  YtD <- which(!is.na(YiD))
  YiD <- YiD[YtD]
  YnD <- length(YiD)
  if(YnD > 0 || XnD > 0){
    mcoutput$Dvar <- sqrt(mcoutput$Dvar)
    Dbounds <- cbind(
      c(vtransform(x=matrix(NA,nrow=1,ncol=length(vnames),
                            dimnames=list(NULL,vnames)),
                   auxmetadata=auxmetadata, Dout='sleft',
                   useOquantiles=useOquantiles)),
      ## sign is important here:
      ## for upper tail, take opposite mean and value
      -c(vtransform(x=matrix(NA,nrow=1,ncol=length(vnames),
                             dimnames=list(NULL,vnames)),
                    auxmetadata=auxmetadata, Dout='sright',
                    useOquantiles=useOquantiles))
    )
  }
  ##
  ZiD <- match(vnames, Zvrt)
  ZtD <- which(!is.na(ZiD))
  ZiD <- ZiD[ZtD]
  ZnD <- length(ZiD)

  ## Type O
  vnames <- auxmetadata[mcmctype == 'O', name]
  XiO <- match(vnames, Xvrt)
  XtO <- which(!is.na(XiO))
  XiO <- XiO[XtO]
  XnO <- length(XiO)
  ##
  YiO <- match(vnames, Yvrt)
  YtO <- which(!is.na(YiO))
  YiO <- YiO[YtO]
  YnO <- length(YiO)
  if(YnO > 0 || XnO > 0){
    mcoutput$Ovar <- sqrt(mcoutput$Ovar)
    ##
    Omaxn <- max(auxmetadata[name %in% vnames, Nvalues])
    Oleft <- t(sapply(vnames, function(avar){
      nn <- auxmetadata[name == avar, Nvalues]
      seqs <- seq( auxmetadata[name == avar, domainmin],
                  auxmetadata[name == avar, domainmax],
                  length.out=nn )
      c(vtransform(seqs, auxmetadata, Oout='left', variates=avar,
                   useOquantiles=useOquantiles),
        rep(NA,Omaxn-nn))
    }))
    Oright <- t(sapply(vnames, function(avar){
      nn <- auxmetadata[name == avar, Nvalues]
      seqs <- seq( auxmetadata[name == avar, domainmin],
                  auxmetadata[name == avar, domainmax],
                  length.out=nn )
      c(vtransform(seqs, auxmetadata, Oout='right', variates=avar,
                   useOquantiles=useOquantiles),
        rep(NA,Omaxn-nn))
    }))
  }
  ##
  ZiO <- match(vnames, Zvrt)
  ZtO <- which(!is.na(ZiO))
  ZiO <- ZiO[ZtO]
  ZnO <- length(ZiO)


  ## Type N
  vnames <- auxmetadata[mcmctype == 'N', name]
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

  ## Type B
  vnames <- auxmetadata[mcmctype == 'B', name]
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
    (if(ZnO > 0){# ordinal
       totake <- cbind(rep(ZtO,each=n), Ws, sseq)
       rnorm(n=n*ZnO,
             mean=mcoutput$Rmean[totake],
             sd=mcoutput$Rvar[totake]
             )
     }else{NULL}),
    (if(ZnN > 0){# nominal
       totake <- cbind(rep(ZtN,each=n), Ws, sseq)
       extraDistr::rcat(n=n*ZnO,
                        prob=apply(mcoutput$Nprob,3,`[`,totake))
     }else{NULL}),
    (if(ZnB > 0){# binary
       totake <- cbind(rep(ZtB,each=n), Ws, sseq)
       extraDistr::rbern(n=n*ZnO,
                         prob=mcoutput$Bprob[totake])
     }else{NULL})
    )

  ## rows: n samples, cols: Z variates
  dim(Zout) <- c(n, length(Zvrt))
  ## Match to original order of Zvrt
  Zout <- vtransform(
    Zout[, match(Zvrt, Zvrt[c(ZiR, ZiC, ZiD, ZiO, ZiN, ZiB)])],
    auxmetadata = auxmetadata,
    Rout = 'id',
    Cout = 'idboundinf',
    Dout = 'idboundinf',
    Oout = '',
    Nout = '',
    Bout = '',
    useOquantiles = useOquantiles)
    
    
  colnames(Zout) <- Zvrt

  

  Yout <- Zout[, Yvrt]
  Xout <- Zout[, Xvrt]
  rm(Zout)
  gc()

#### STEP 2a. Calculate p(X) for all samples




}


set.seed(1)
system.time(for(i in 1:1000){
                out <- rnorm(n=3*4*10000,mean=test1,sd=test1)
                dim(out) <- c(3,4*10000)
                })
##
set.seed(1)
system.time(for(i in 1:1000){
                out2 <- matrix(rnorm(n=3*4*10000,mean=test1,sd=test1), nrow=3, ncol=4*10000)
                })
identical(out,out2)

set.seed(1)
system.time(for(i in 1:1000){
                out2 <- matrix(rnorm(n=3*4*10000,mean=test1,sd=test1), nrow=3, ncol=4*10000)
                })
##
set.seed(1)
system.time(for(i in 1:1000){
                out <- rnorm(n=3*4*10000,mean=test1,sd=test1)
                dim(out) <- c(3,4*10000)
                })
identical(out,out2)

#### Unfinished
generateVariates <- function(n=1, Yvrt, X=NULL, mcoutput, combine='rbind', useOquantiles=TRUE, parallel=TRUE, silent=FALSE){

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

  nsamples <- ncol(mcoutput$W)
  nclusters <- nrow(mcoutput$W)

  ## Consistency checks
  if(!is.character(Yvrt) || any(is.na(Yvrt))){
    stop('Yvrt must be a vector of variate names')
  }
  if(!is.null(X) && length(dim(X)) != 2){
    stop('X must be NULL or have two dimensions')
  }

  if(!all(Yvrt %in% auxmetadata$name)){stop('unknown Y variates\n')}
  if(length(unique(Yvrt)) != length(Yvrt)){stop('duplicate Y variates\n')}
  ##
  Xv <- colnames(X)
  if(!all(Xv %in% auxmetadata$name)){stop('unknown X variates\n')}
  if(length(unique(Xv)) != length(Xv)){stop('duplicate X variates\n')}
  ##
  if(length(intersect(Yvrt, Xv)) > 0){stop('overlap in Y and X variates\n')}


#### Subsample and get nclusters and nsamples
  ## source('mcsubset.R')
  if((is.logical(n) && n) || (is.character(n) && n=='all')){
    n <- nsamples
  }
  subsamples <- c(rep(1:nsamples, n %/% nsamples),
                  sort(sample(1:nsamples, n %% nsamples, replace=F)))
  mcoutput <- mcsubset(mcoutput, subsamples)
  sseq <- seq_len(n)

  ## source('vtransform.R')

#### Type R
  vnames <- auxmetadata[mcmctype == 'R', name]
  XiR <- match(vnames, Xv)
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

#### Type C
  vnames <- auxmetadata[mcmctype == 'C', name]
  XiC <- match(vnames, Xv)
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

#### Type D
  vnames <- auxmetadata[mcmctype == 'D', name]
  XiD <- match(vnames, Xv)
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

#### Type O
  vnames <- auxmetadata[mcmctype == 'O', name]
  XiO <- match(vnames, Xv)
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

#### Type N
  vnames <- auxmetadata[mcmctype == 'N', name]
  XiN <- match(vnames, Xv)
  XtN <- which(!is.na(XiN))
  XiN <- XiN[XtN]
  XnN <- length(XiN)
  ##
  YiN <- match(vnames, Yvrt)
  YtN <- which(!is.na(YiN))
  YiN <- YiN[YtN]
  YnN <- length(YiN)
  ## if(YnN > 0 || XnN > 0){
  ##     }

#### Type B
  vnames <- auxmetadata[mcmctype == 'B', name]
  XiB <- match(vnames, Xv)
  XtB <- which(!is.na(XiB))
  XiB <- XiB[XtB]
  XnB <- length(XiB)
  ##
  YiB <- match(vnames, Yvrt)
  YtB <- which(!is.na(YiB))
  YiB <- YiB[YtB]
  YnB <- length(YiB)

#### Match original order of Yvrt with the order from foreach loop
  Yorder <- match(Yvrt, Yvrt[c(YiR, YiC, YiD, YiO, YiN, YiB)])
  Yn <- length(Yvrt)


#### transformation of inputs
  if(!is.null(X)){
    X2 <- vtransform(X, auxmetadata,
                     Cout='boundisinf',
                     Dout='boundisinf',
                     Oout='',
                     Nout='numeric',
                     Bout='numeric',
                     useOquantiles=useOquantiles)
  }else{X2 <- t(NA)}

  out <- foreach(x=t(X2), .combine=rbind, .inorder=TRUE)%dochains%{

    if(all(is.na(x))){
      probX <- mcoutput$W
    }else{
      ## rows: clusters, cols: samples
      probX <- log(mcoutput$W) +
        (if(XnR > 0){# continuous
           colSums(
             dnorm(x=x[XiR,],
                   mean=mcoutput$Rmean[XtR,,,drop=FALSE],
                   sd=mcoutput$Rvar[XtR,,,drop=FALSE],
                   log=TRUE),
             na.rm=TRUE)
         }else{0}) +
        (if(XnC > 0){# censored
           indf <- which(is.finite(x[XiC,]))
           indi <- which(!is.finite(x[XiC,]))
           (if(length(indf) > 0){
              colSums(
                dnorm(x=x[XiC[indf],],
                      mean=mcoutput$Cmean[XtC[indf],,,drop=FALSE],
                      sd=mcoutput$Cvar[XtC[indf],,,drop=FALSE],
                      log=TRUE),
                na.rm=TRUE)
            }else{0}) +
             (if(length(indi) > 0){
                v2 <- XtC[indi]
                v1 <- -sign(x[XiC[indi],])
                ## for upper tail, take opposite mean and value
                colSums(
                  pnorm(q=Cbounds[cbind(v2, 1.5-0.5*v1)],
                        mean=v1*mcoutput$Cmean[v2,,,drop=FALSE],
                        sd=mcoutput$Cvar[v2,,,drop=FALSE],
                        log.p=TRUE),
                  na.rm=TRUE)
              }else{0})
         }else{0}) +
        (if(XnD > 0){# continuous discretized
           indf <- which(is.finite(x[XiD,]))
           indi <- which(!is.finite(x[XiD,]))
           (if(length(indf) > 0){
              colSums(
                dnorm(x=x[XiD[indf],],
                      mean=mcoutput$Dmean[XtD[indf],,,drop=FALSE],
                      sd=mcoutput$Dvar[XtD[indf],,,drop=FALSE],
                      log=TRUE),
                na.rm=TRUE)
            }else{0}) +
             (if(length(indi) > 0){
                v2 <- XtD[indi]
                v1 <- -sign(x[XiD[indi],])
                ## for upper tail, take opposite mean and value
                colSums(
                  pnorm(q=Dbounds[v2, 1.5-0.5*v1],
                        mean=v1*mcoutput$Dmean[v2,,,drop=FALSE],
                        sd=mcoutput$Dvar[v2,,,drop=FALSE],
                        log.p=TRUE),
                  na.rm=TRUE)
              }else{0})
         }else{0}) +
        (if(XnO > 0){# ordinal
           v2 <- cbind(XtO,x[XiO,])
           colSums(
             log(
               pnorm(q=Oright[v2],
                     mean=mcoutput$Omean[XtO,,,drop=FALSE],
                     sd=mcoutput$Ovar[XtO,,,drop=FALSE]) -
               pnorm(q=Oleft[v2],
                     mean=mcoutput$Omean[XtO,,,drop=FALSE],
                     sd=mcoutput$Ovar[XtO,,,drop=FALSE])
             ),
             na.rm=TRUE)
         }else{0}) +
        (if(XnN > 0){# nominal
           colSums(
             log( aperm(
               vapply(seq_len(XnN), function(v){
                 mcoutput$Nprob[XtN[v],,x[XiN[v],],]
               }, mcoutput$W),
               c(3,1,2)) ),
             na.rm=TRUE)
           ## colSums(
           ##      log(array(
           ##          t(sapply(seq_len(XnN), function(v){
           ##              mcoutput$Nprob[XtN[v],,x[XiN[v],],]
           ##          })),
           ##          dim=c(XnN, nclusters, nsamples))),
           ##      na.rm=TRUE)
           ## temp <- apply(mcoutput$Nprob, c(2,4), function(xx){
           ##     xx[cbind(XtN, x[XiN,])]
           ## })
           ## dim(temp) <- c(XnN, nclusters, nsamples)
           ## ##
           ## colSums(log(temp), na.rm=TRUE)
         }else{0}) +
        (if(XnB > 0){# binary
           colSums(
             log( x[XiB,]*mcoutput$Bprob[XtB,,,drop=FALSE] +
                  (1L-x[XiB,])*(1-mcoutput$Bprob[XtB,,,drop=FALSE]) ),
             na.rm=TRUE)
         }else{0})
    } # end probX
    probX <- t(exp(apply(probX, 2, function(xx){xx - max(xx[is.finite(xx)])})))
    ## probX: each row is a MC sample
    ## the row elements are the probabilities of the clusters
    ## for that MC sample

    ## Y is drawn as follows, for each MC sample:
    ## 1. draw a cluster, according to its probability
    ## 2. draw from the appropriate kernel distributions
    ## using the parameters of that cluster
    Ws <- extraDistr::rcat(n=n, prob=probX)
    Yout <- c( # rows: n samples, cols: Y variates
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
             mean=mcoutput$Rmean[totake],
             sd=mcoutput$Rvar[totake]
             )
     }else{NULL}),
    (if(YnD > 0){# continuous discretized
       totake <- cbind(rep(YtD,each=n), Ws, sseq)
       rnorm(n=n*YnD,
             mean=mcoutput$Rmean[totake],
             sd=mcoutput$Rvar[totake]
             )
     }else{NULL}),
    (if(YnO > 0){# ordinal
       totake <- cbind(rep(YtO,each=n), Ws, sseq)
       rnorm(n=n*YnO,
             mean=mcoutput$Rmean[totake],
             sd=mcoutput$Rvar[totake]
             )
     }else{NULL}),
    (if(YnN > 0){# nominal
       totake <- cbind(rep(YtN,each=n), Ws, sseq)
       extraDistr::rcat(n=n*YnO,
                        prob=apply(mcoutput$Nprob,3,`[`,totake))
     }else{NULL}),
    (if(YnB > 0){# binary
       totake <- cbind(rep(YtB,each=n), Ws, sseq)
       extraDistr::rbern(n=n*YnO,
                         prob=mcoutput$Bprob[totake])
     }else{NULL})
    )
    dim(Yout) <- c(n, Yn)
  } # end foreach

  ## Remove doRNG information
  attr(out,'rng') <- NULL
  attr(out,'doRNG_version') <- NULL

  ## Reorder Y outputs (columns) according to original Yvrt
  Yout <- Yout[,Yorder]
  



}


## set.seed(1)
## system.time(for(i in 1:1000){
##                 out <- rnorm(n=3*4*10000,mean=test1,sd=test1)
##                 dim(out) <- c(3,4*10000)
##                 })
## ##
## set.seed(1)
## system.time(for(i in 1:1000){
##                 out2 <- matrix(rnorm(n=3*4*10000,mean=test1,sd=test1), nrow=3, ncol=4*10000)
##                 })
## identical(out,out2)
## 
## set.seed(1)
## system.time(for(i in 1:1000){
##                 out2 <- matrix(rnorm(n=3*4*10000,mean=test1,sd=test1), nrow=3, ncol=4*10000)
##                 })
## ##
## set.seed(1)
## system.time(for(i in 1:1000){
##                 out <- rnorm(n=3*4*10000,mean=test1,sd=test1)
##                 dim(out) <- c(3,4*10000)
##                 })
## identical(out,out2)

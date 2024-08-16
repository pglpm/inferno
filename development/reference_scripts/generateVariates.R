#### Unfinished
generateVariates <- function(n=1, Yvrt, X=NULL, agent, combine='rbind', useLquantiles=TRUE, parallel=TRUE, silent=FALSE){

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
    `%dochains%` <- `%dorng%`
  }

  ## Extract Monte Carlo output & aux-metadata
  ## If agent is a string, check if it's a folder name or file name
  if (is.character(agent)) {
                                        # Check if 'agent' is a folder containing agent.rds
    if (file_test('-d', agent) &&
        file.exists(file.path(agent, 'agent.rds'))) {
      agent <- readRDS(file.path(agent, 'agent.rds'))
    } else {
      ## Assume 'agent' the full path of agent.rds
      ## possibly without the file extension '.rds'
      agent <- paste0(sub('.rds$', '', agent), '.rds')
      if (file.exists(agent)) {
        agent <- readRDS(agent)
      } else {
        stop("The argument 'agent' must be a folder containing agent.rds, or the path to an rds-file containing the output from learn().")
      }
    }
  }
  ## Add check to see that agent is correct type of object?
  auxmetadata <- agent$auxmetadata
  agent$auxmetadata <- NULL

  nsamples <- ncol(agent$W)
  ncomponents <- nrow(agent$W)

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


#### Subsample and get ncomponents and nsamples
  ## source('mcsubset.R')
  if((is.logical(n) && n) || (is.character(n) && n=='all')){
    n <- nsamples
  }
  subsamples <- c(rep(1:nsamples, n %/% nsamples),
                  sort(sample(1:nsamples, n %% nsamples, replace=F)))
  agent <- mcsubset(agent, subsamples)
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
    agent$Rvar <- sqrt(agent$Rvar)
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
    agent$Cvar <- sqrt(agent$Cvar)
    Cbounds <- cbind(
      c(vtransform(x=matrix(NA,nrow=1,ncol=length(vnames),
                            dimnames=list(NULL,vnames)),
                   auxmetadata=auxmetadata, Cout='sleft')),
      ## sign is important here:
      ## for upper tail, take opposite mean and value
      -c(vtransform(x=matrix(NA,nrow=1,ncol=length(vnames),
                             dimnames=list(NULL,vnames)),
                    auxmetadata=auxmetadata, Cout='sright'))
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
    agent$Dvar <- sqrt(agent$Dvar)
    Dbounds <- cbind(
      c(vtransform(x=matrix(NA,nrow=1,ncol=length(vnames),
                            dimnames=list(NULL,vnames)),
                   auxmetadata=auxmetadata, Dout='sleft')),
      ## sign is important here:
      ## for upper tail, take opposite mean and value
      -c(vtransform(x=matrix(NA,nrow=1,ncol=length(vnames),
                             dimnames=list(NULL,vnames)),
                    auxmetadata=auxmetadata, Dout='sright'))
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
    agent$Ovar <- sqrt(agent$Ovar)
    ##
    Omaxn <- max(auxmetadata[name %in% vnames, Nvalues])
    Oleft <- t(sapply(vnames, function(avar){
      nn <- auxmetadata[name == avar, Nvalues]
      seqs <- seq( auxmetadata[name == avar, domainmin],
                  auxmetadata[name == avar, domainmax],
                  length.out=nn )
      c(vtransform(seqs, auxmetadata = auxmetadata, Oout='left', variates=avar,
                   useLquantiles=useLquantiles),
        rep(NA,Omaxn-nn))
    }))
    Oright <- t(sapply(vnames, function(avar){
      nn <- auxmetadata[name == avar, Nvalues]
      seqs <- seq( auxmetadata[name == avar, domainmin],
                  auxmetadata[name == avar, domainmax],
                  length.out=nn )
      c(vtransform(seqs, auxmetadata = auxmetadata, Oout='right', variates=avar,
                   useLquantiles=useLquantiles),
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
    X2 <- vtransform(X, auxmetadata = auxmetadata,
                     Cout='boundisinf',
                     Dout='boundisinf',
                     Oout='',
                     Nout='numeric',
                     Bout='numeric',
                     useLquantiles=useLquantiles)
  }else{X2 <- t(NA)}

  out <- foreach(x=t(X2), .combine=rbind, .inorder=TRUE)%dochains%{

    if(all(is.na(x))){
      probX <- agent$W
    }else{
      ## rows: components, cols: samples
      probX <- log(agent$W) +
        (if(XnR > 0){# continuous
           colSums(
             dnorm(x=x[XiR,],
                   mean=agent$Rmean[XtR,,,drop=FALSE],
                   sd=agent$Rvar[XtR,,,drop=FALSE],
                   log=TRUE),
             na.rm=TRUE)
         }else{0}) +
        (if(XnC > 0){# censored
           indf <- which(is.finite(x[XiC,]))
           indi <- which(!is.finite(x[XiC,]))
           (if(length(indf) > 0){
              colSums(
                dnorm(x=x[XiC[indf],],
                      mean=agent$Cmean[XtC[indf],,,drop=FALSE],
                      sd=agent$Cvar[XtC[indf],,,drop=FALSE],
                      log=TRUE),
                na.rm=TRUE)
            }else{0}) +
             (if(length(indi) > 0){
                v2 <- XtC[indi]
                v1 <- -sign(x[XiC[indi],])
                ## for upper tail, take opposite mean and value
                colSums(
                  pnorm(q=Cbounds[cbind(v2, 1.5-0.5*v1)],
                        mean=v1*agent$Cmean[v2,,,drop=FALSE],
                        sd=agent$Cvar[v2,,,drop=FALSE],
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
                      mean=agent$Dmean[XtD[indf],,,drop=FALSE],
                      sd=agent$Dvar[XtD[indf],,,drop=FALSE],
                      log=TRUE),
                na.rm=TRUE)
            }else{0}) +
             (if(length(indi) > 0){
                v2 <- XtD[indi]
                v1 <- -sign(x[XiD[indi],])
                ## for upper tail, take opposite mean and value
                colSums(
                  pnorm(q=Dbounds[v2, 1.5-0.5*v1],
                        mean=v1*agent$Dmean[v2,,,drop=FALSE],
                        sd=agent$Dvar[v2,,,drop=FALSE],
                        log.p=TRUE),
                  na.rm=TRUE)
              }else{0})
         }else{0}) +
        (if(XnO > 0){# ordinal
           v2 <- cbind(XtO,x[XiO,])
           colSums(
             log(
               pnorm(q=Oright[v2],
                     mean=agent$Omean[XtO,,,drop=FALSE],
                     sd=agent$Ovar[XtO,,,drop=FALSE]) -
               pnorm(q=Oleft[v2],
                     mean=agent$Omean[XtO,,,drop=FALSE],
                     sd=agent$Ovar[XtO,,,drop=FALSE])
             ),
             na.rm=TRUE)
         }else{0}) +
        (if(XnN > 0){# nominal
           colSums(
             log( aperm(
               vapply(seq_len(XnN), function(v){
                 agent$Nprob[XtN[v],,x[XiN[v],],]
               }, agent$W),
               c(3,1,2)) ),
             na.rm=TRUE)
           ## colSums(
           ##      log(array(
           ##          t(sapply(seq_len(XnN), function(v){
           ##              agent$Nprob[XtN[v],,x[XiN[v],],]
           ##          })),
           ##          dim=c(XnN, ncomponents, nsamples))),
           ##      na.rm=TRUE)
           ## temp <- apply(agent$Nprob, c(2,4), function(xx){
           ##     xx[cbind(XtN, x[XiN,])]
           ## })
           ## dim(temp) <- c(XnN, ncomponents, nsamples)
           ## ##
           ## colSums(log(temp), na.rm=TRUE)
         }else{0}) +
        (if(XnB > 0){# binary
           colSums(
             log( x[XiB,]*agent$Bprob[XtB,,,drop=FALSE] +
                  (1L-x[XiB,])*(1-agent$Bprob[XtB,,,drop=FALSE]) ),
             na.rm=TRUE)
         }else{0})
    } # end probX
    probX <- t(exp(apply(probX, 2, function(xx){xx - max(xx[is.finite(xx)])})))
    ## probX: each row is a MC sample
    ## the row elements are the probabilities of the components
    ## for that MC sample

    ## Y is drawn as follows, for each MC sample:
    ## 1. draw a component, according to its probability
    ## 2. draw from the appropriate kernel distributions
    ## using the parameters of that component
    Ws <- extraDistr::rcat(n=n, prob=probX)
    Yout <- c( # rows: n samples, cols: Y variates
    (if(YnR > 0){# continuous
       totake <- cbind(rep(YtR,each=n), Ws, sseq)
       rnorm(n=n*YnR,
             mean=agent$Rmean[totake],
             sd=agent$Rvar[totake]
             )
     }else{NULL}),
    (if(YnC > 0){# censored
       totake <- cbind(rep(YtC,each=n), Ws, sseq)
       rnorm(n=n*YnC,
             mean=agent$Rmean[totake],
             sd=agent$Rvar[totake]
             )
     }else{NULL}),
    (if(YnD > 0){# continuous discretized
       totake <- cbind(rep(YtD,each=n), Ws, sseq)
       rnorm(n=n*YnD,
             mean=agent$Rmean[totake],
             sd=agent$Rvar[totake]
             )
     }else{NULL}),
    (if(YnO > 0){# ordinal
       totake <- cbind(rep(YtO,each=n), Ws, sseq)
       rnorm(n=n*YnO,
             mean=agent$Rmean[totake],
             sd=agent$Rvar[totake]
             )
     }else{NULL}),
    (if(YnN > 0){# nominal
       totake <- cbind(rep(YtN,each=n), Ws, sseq)
       extraDistr::rcat(n=n*YnO,
                        prob=apply(agent$Nprob,3,`[`,totake))
     }else{NULL}),
    (if(YnB > 0){# binary
       totake <- cbind(rep(YtB,each=n), Ws, sseq)
       extraDistr::rbern(n=n*YnO,
                         prob=agent$Bprob[totake])
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

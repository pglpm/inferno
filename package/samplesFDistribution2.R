samplesFDistribution2 <- function(Y, X, mcsamples, auxmetadata, subsamples, jacobian=TRUE, fn=identity, parallel=TRUE, useOquantiles=TRUE){
    ## Consistency checks
    if(length(dim(Y)) != 2){stop('Y must have two dimensions')}
    if(missing(X)){X <- NULL}
    if(!is.null(X) && length(dim(X)) != 2){stop('X must be NULL or have two dimensions')}
    ##
    if(!is.null(X) && ncol(X) == 0){X <- NULL}
    ## auxmetadata
    if(is.character(auxmetadata) && file.exists(auxmetadata)){
        auxmetadata <- readRDS(auxmetadata)
    }

    ## More consistency checks
    Yv <- colnames(Y)
    if(!all(Yv %in% auxmetadata$name)){stop('unknown Y variates\n')}
    if(length(unique(Yv)) != length(Yv)){stop('duplicate Y variates\n')}
    ##
    Xv <- colnames(X)
    if(!all(Xv %in% auxmetadata$name)){stop('unknown X variates\n')}
    if(length(unique(Xv)) != length(Xv)){stop('duplicate X variates\n')}
    ##
    if(length(intersect(Yv, Xv)) > 0){stop('overlap in Y and X variates\n')}

    ## mcsamples and subsamples
    if(is.character(mcsamples) && file.exists(mcsamples)){
        mcsamples <- readRDS(mcsamples)
    }
    ##
    if(missing(subsamples) || is.null(subsamples) || (is.logical(subsamples) && !subsamples)){
        subsamples <- 1:ncol(mcsamples$W)
    }else if(is.character(subsamples)){
        subsamples <- round(seq(1,nrow(subsamples),length.out=as.numeric(subsamples)))
    }
    ##
    mcsamples <- lapply(mcsamples,function(xx){
        do.call('[',c(list(xx),rep(TRUE,length(dim(xx))-1), list(subsamples), list(drop=FALSE)) )
    })

    W <- mcsamples$W
    mcsamples$W <- NULL
    nsamples <- ncol(W)

    nclusters <- nrow(W)

#### Type R
    vnames <- auxmetadata[mcmctype == 'R', name]
    XiR <- match(vnames, Xv)
    XtR <- which(!is.na(XiR))
    XiR <- XiR[XtR]
    XnR <- length(XiR)
    ##
    YiR <- match(vnames, Yv)
    YtR <- which(!is.na(YiR))
    YiR <- YiR[YtR]
    YnR <- length(YiR)
    if(YnR > 0 || XnR > 0){
        Rmean <- mcsamples$Rmean
        mcsamples$Rmean <- NULL
        Rvar <- sqrt(mcsamples$Rvar)
        mcsamples$Rvar <- NULL
    }

#### Type C
    vnames <- auxmetadata[mcmctype == 'C', name]
    XiC <- match(vnames, Xv)
    XtC <- which(!is.na(XiC))
    XiC <- XiC[XtC]
    XnC <- length(XiC)
    ##
    YiC <- match(vnames, Yv)
    YtC <- which(!is.na(YiC))
    YiC <- YiC[YtC]
    YnC <- length(YiC)
    if(YnC > 0 || XnC > 0){
        Cmean <- mcsamples$Cmean
        mcsamples$Cmean <- NULL
        Cvar <- sqrt(mcsamples$Cvar)
        mcsamples$Cvar <- NULL
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
    YiD <- match(vnames, Yv)
    YtD <- which(!is.na(YiD))
    YiD <- YiD[YtD]
    YnD <- length(YiD)
    if(YnD > 0 || XnD > 0){
        Dmean <- mcsamples$Dmean
        mcsamples$Dmean <- NULL
        Dvar <- sqrt(mcsamples$Dvar)
        mcsamples$Dvar <- NULL
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
    YiO <- match(vnames, Yv)
    YtO <- which(!is.na(YiO))
    YiO <- YiO[YtO]
    YnO <- length(YiO)
    if(YnO > 0 || XnO > 0){
        Omean <- mcsamples$Omean
        mcsamples$Omean <- NULL
        Ovar <- sqrt(mcsamples$Ovar)
        mcsamples$Ovar <- NULL
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
    YiN <- match(vnames, Yv)
    YtN <- which(!is.na(YiN))
    YiN <- YiN[YtN]
    YnN <- length(YiN)
    if(YnN > 0 || XnN > 0){
        Nprob <- mcsamples$Nprob
        mcsamples$Nprob <- NULL
    }

#### Type B
    vnames <- auxmetadata[mcmctype == 'B', name]
    XiB <- match(vnames, Xv)
    XtB <- which(!is.na(XiB))
    XiB <- XiB[XtB]
    XnB <- length(XiB)
    ##
    YiB <- match(vnames, Yv)
    YtB <- which(!is.na(YiB))
    YiB <- YiB[YtB]
    YnB <- length(YiB)
    if(YnB > 0 || XnB > 0){
        Bprob <- mcsamples$Bprob
        mcsamples$Bprob <- NULL
    }

#### transformation of inputs
    Y2 <- vtransform(Y, auxmetadata, Cout='boundisinf', Dout='boundisinf', Oout='', Nout='numeric', Bout='numeric', useOquantiles=useOquantiles)
    if(!is.null(X)){
        X2 <- vtransform(X, auxmetadata, Cout='boundisinf', Dout='boundisinf', Oout='', Nout='numeric', Bout='numeric', useOquantiles=useOquantiles)
        if(nrow(X2) < nrow(Y2)){
            warning('*Note: X has fewer data than Y. Recycling*')
            X2 <- t(matrix(rep(t(X2), ceiling(nrow(Y2)/nrow(X2))), nrow=ncol(X2), dimnames=list(colnames(X2),NULL)))[1:nrow(Y2),,drop=FALSE]
        }
        if(nrow(X2) > nrow(Y2)){
            warning('*Note: X has more data than Y. Recycling*')
            Y2 <- t(matrix(rep(t(Y2), ceiling(nrow(X2)/nrow(Y2))), nrow=ncol(Y2), dimnames=list(colnames(Y2),NULL)))[1:nrow(X2),,drop=FALSE]
            Y <- t(matrix(rep(t(Y), ceiling(nrow(X2)/nrow(Y))), nrow=ncol(Y), dimnames=list(colnames(Y),NULL)))[1:nrow(X2),,drop=FALSE]
        }
    }else{X2 <- lapply(seq_len(nrow(Y2)),function(x)NA)}
    ## ndata <- nrow(Y2)
    ##
    ##

    if(parallel){`%thisdo%` <- `%dopar%`}else{`%thisdo%` <- `%do%`}
    foreach(y=t(Y2), x=t(X2), .combine=rbind, .inorder=T)%thisdo%{
#### the loop is over the columns of y and x
#### each instance is a 1-column vector
        ##
        if(all(is.na(x))){
            probX <- log(W) 
        }else{
            ## rows: clusters, cols: samples
            probX <- log(W) + 
                (if(XnR > 0){# continuous
                     colSums(
                         dnorm(x=x[XiR,],
                               mean=Rmean[XtR,,,drop=FALSE],
                               sd=Rvar[XtR,,,drop=FALSE],
                               log=T),
                         na.rm=T)
                 }else{0}) +
                (if(XnC > 0){# censored
                     indf <- which(is.finite(x[XiC,]))
                     indi <- which(!is.finite(x[XiC,]))
                     (if(length(indf) > 0){
                          colSums(
                              dnorm(x=x[XiC[indf],],
                                    mean=Cmean[XtC[indf],,,drop=FALSE],
                                    sd=Cvar[XtC[indf],,,drop=FALSE],
                                    log=T),
                              na.rm=T)
                      }else{0}) +
                         (if(length(indi) > 0){
                              v2 <- XtC[indi]
                              v1 <- -sign(x[XiC[indi],])
                              ## for upper tail, take opposite mean and value
                              colSums(
                                  pnorm(q=Cbounds[cbind(v2, 1.5-0.5*v1)],
                                        mean=v1*Cmean[v2,,,drop=FALSE],
                                        sd=Cvar[v2,,,drop=FALSE],
                                        log.p=T),
                                  na.rm=T)
                          }else{0})
                 }else{0}) +
                (if(XnD > 0){# continuous discretized
                     indf <- which(is.finite(x[XiD,]))
                     indi <- which(!is.finite(x[XiD,]))
                     (if(length(indf) > 0){
                          colSums(
                              dnorm(x=x[XiD[indf],],
                                    mean=Dmean[XtD[indf],,,drop=FALSE],
                                    sd=Dvar[XtD[indf],,,drop=FALSE],
                                    log=T),
                              na.rm=T)
                      }else{0}) +
                         (if(length(indi) > 0){
                              v2 <- XtD[indi]
                              v1 <- -sign(x[XiD[indi],])
                              ## for upper tail, take opposite mean and value
                              colSums(
                                  pnorm(q=Dbounds[v2, 1.5-0.5*v1],
                                        mean=v1*Dmean[v2,,,drop=FALSE],
                                        sd=Dvar[v2,,,drop=FALSE],
                                        log.p=T),
                                  na.rm=T)
                          }else{0})
                 }else{0}) +
                (if(XnO > 0){# ordinal
                     v2 <- cbind(XtO,x[XiO,])
                     colSums(
                         log(
                             pnorm(q=Oright[v2],
                                   mean=Omean[XtO,,,drop=FALSE],
                                   sd=Ovar[XtO,,,drop=FALSE]) -
                             pnorm(q=Oleft[v2],
                                   mean=Omean[XtO,,,drop=FALSE],
                                   sd=Ovar[XtO,,,drop=FALSE])
                         ),
                         na.rm=T)
                 }else{0}) +
                (if(XnN > 0){# nominal
                     temp <- apply(Nprob, c(2,4), function(xx){
                         xx[cbind(XtN, x[XiN,])]
                     })
                     dim(temp) <- c(XnN, nclusters, nsamples)
                     ##
                     colSums(log(temp), na.rm=T)
                 }else{0}) +
                (if(XnB > 0){# binary
                     colSums(
                         log( x[XiB,]*Bprob[XtB,,,drop=FALSE] +
                              (1L-x[XiB,])*(1-Bprob[XtB,,,drop=FALSE]) ),
                         na.rm=T)
                 }else{0})
        } # end probX
        probX <- apply(probX, 2, function(xx){xx - max(xx[is.finite(xx)])})
        ##
        ##
        if(all(is.na(y))){
            probY <- array(NA, dim=dim(W))
        }else{
            probY <- (if(YnR > 0){# continuous
                          colSums(
                              dnorm(x=y[YiR,],
                                    mean=Rmean[YtR,,,drop=FALSE],
                                    sd=Rvar[YtR,,,drop=FALSE],
                                    log=T),
                              na.rm=T)
                      }else{0}) +
                (if(YnC > 0){# censored
                     indf <- which(is.finite(y[YiC,]))
                     indi <- which(!is.finite(y[YiC,]))
                     (if(length(indf) > 0){
                          colSums(
                              dnorm(x=y[YiC[indf],],
                                    mean=Cmean[YtC[indf],,,drop=FALSE],
                                    sd=Cvar[YtC[indf],,,drop=FALSE],
                                    log=T),
                              na.rm=T)
                      }else{0}) +
                         (if(length(indi) > 0){
                              v2 <- YtC[indi]
                              v1 <- -sign(y[YiC[indi],])
                              ## for upper tail, take opposite mean and value
                              colSums(
                                  pnorm(q=Cbounds[cbind(v2, 1.5-0.5*v1)],
                                        mean=v1*Cmean[v2,,,drop=FALSE],
                                        sd=Cvar[v2,,,drop=FALSE],
                                        log.p=T),
                                  na.rm=T)
                          }else{0})
                 }else{0}) +
                (if(YnD > 0){# continuous discretized
                     indf <- which(is.finite(y[YiD,]))
                     indi <- which(!is.finite(y[YiD,]))
                     (if(length(indf) > 0){
                          colSums(
                              dnorm(x=y[YiD[indf],],
                                    mean=Dmean[YtD[indf],,,drop=FALSE],
                                    sd=Dvar[YtD[indf],,,drop=FALSE],
                                    log=T),
                              na.rm=T)
                      }else{0}) +
                         (if(length(indi) > 0){
                              v2 <- YtD[indi]
                              v1 <- -sign(y[YiD[indi],])
                              ## for upper tail, take opposite mean and value
                              colSums(
                                  pnorm(q=Dbounds[v2, 1.5-0.5*v1],
                                        mean=v1*Dmean[v2,,,drop=FALSE],
                                        sd=Dvar[v2,,,drop=FALSE],
                                        log.p=T),
                                  na.rm=T)
                          }else{0})
                 }else{0}) +
                (if(YnO > 0){# ordinal
                     v2 <- cbind(YtO,y[YiO,])
                     colSums(
                         log(
                             pnorm(q=Oright[v2],
                                   mean=Omean[YtO,,,drop=FALSE],
                                   sd=Ovar[YtO,,,drop=FALSE]) -
                             pnorm(q=Oleft[v2],
                                   mean=Omean[YtO,,,drop=FALSE],
                                   sd=Ovar[YtO,,,drop=FALSE])
                         ),
                         na.rm=T)
                 }else{0}) +
                (if(YnN > 0){# nominal
                     temp <- apply(Nprob, c(2,4), function(xx){
                         xx[cbind(YtN, y[YiN,])]
                     })
                     dim(temp) <- c(YnN, nclusters, nsamples)
                     ##
                     colSums(log(temp), na.rm=T)
                 }else{0}) +
                (if(YnB > 0){# binary
                     colSums(
                         log( y[YiB,]*Bprob[YtB,,,drop=FALSE] +
                              (1-y[YiB,])*(1-Bprob[YtB,,,drop=FALSE]) ),
                         na.rm=T)
                 }else{0})
        }
#### Output: rows=clusters, columns=samples
        ##
        ## if(all(is.na(x))){
        ##     out <- rowSums(exp(probX+probY)) 
        ## }else{
        ## str(probX)
        ## print(any(is.na(probX)))
        ## print(apply(probX, 1, max, na.rm=T))
        ## str(probX)
        ## print(any(is.na(probX)))
        ## str(probY)
        ## print(any(is.na(probY)))
        ## }
        fn( colSums(exp(probX+probY))/colSums(exp(probX)) )
        ## rbind(log(W),probX,probY)
    } *
    (if(jacobian){
         exp(-rowSums(
                  log(vtransform(Y, auxmetadata=auxmetadata, invjacobian=TRUE, useOquantiles=useOquantiles)),
                  na.rm=T))}else{1L})
}

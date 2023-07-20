samplesFDistribution3 <- function(Y, X, mcsamples, auxmetadata, subsamples, jacobian=TRUE, fn=identity, parallel=TRUE, useOquantiles=TRUE){
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
    nsamples <- ncol(mcsamples$W)

    nclusters <- nrow(mcsamples$W)

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
        mcsamples$Rvar <- sqrt(mcsamples$Rvar)
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
        mcsamples$Cvar <- sqrt(mcsamples$Cvar)
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
        YseqC <- 1:YnC
        XseqC <- 1:XnC
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
        mcsamples$Dvar <- sqrt(mcsamples$Dvar)
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
        YseqD <- 1:YnD
        XseqD <- 1:XnD
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
        mcsamples$Ovar <- sqrt(mcsamples$Ovar)
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
    YseqN <- 1:YnN
    XseqN <- 1:XnN

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
            probX <- log(mcsamples$W) 
        }else{
            ## rows: clusters, cols: samples
            probX <- log(mcsamples$W) + 
                (if(XnR > 0){# continuous
                     colSums(
                         array(dnorm(x=x[XiR,],
                                     mean=mcsamples$Rmean[XtR,,],
                                     sd=mcsamples$Rvar[XtR,,],log=T),
                               dim=c(XnR, nclusters, nsamples)),
                         na.rm=T)
                 }else{0}) +
                (if(XnC > 0){# censored
                     colSums(
                         array(
                             t(sapply(XseqC, function(v){
                                 v1 <- x[XiC[v],]
                                 v2 <- XtC[v]
                                 if(is.finite(v1)){
                                     (dnorm(x=v1,
                                            mean=mcsamples$Cmean[v2,,],
                                            sd=mcsamples$Cvar[v2,,],log=T))
                                 }else{
                                     (pnorm(q=Cbounds[v2, 2L-(v1 < 0)],
                                            mean=mcsamples$Cmean[v2,,],
                                            sd=mcsamples$Cvar[v2,,],
                                            lower.tail=(v1 < 0),
                                            log.p=T))
                                 }
                             })),
                             dim=c(XnC, nclusters, nsamples)),
                         na.rm=T)
                 }else{0}) +
                (if(XnD > 0){# continuous discretized
                     colSums(
                         array(
                             t(sapply(XseqD, function(v){
                                 v1 <- x[XiD[v],]
                                 v2 <- XtD[v]
                                 if(is.finite(v1)){
                                     (dnorm(x=v1,
                                            mean=mcsamples$Dmean[v2,,],
                                            sd=mcsamples$Dvar[v2,,],log=T))
                                 }else{
                                     (pnorm(q=Dbounds[v2, 2L-(v1 < 0)],
                                            mean=mcsamples$Dmean[v2,,],
                                            sd=mcsamples$Dvar[v2,,],
                                            lower.tail=(v1 < 0),
                                            log.p=T))
                                 }
                             })),
                             dim=c(XnD, nclusters, nsamples)),
                         na.rm=T)
                 }else{0}) +
                ## (if(XnD > 0){# continuous
                ##      colSums(
                ##          array(dnorm(x=x[XiD,],
                ##                      mean=mcsamples$Dmean[XtD,,],
                ##                      sd=mcsamples$Dvar[XtD,,],log=T),
                ##                dim=c(XnD, nclusters, nsamples)),
                ##          na.rm=T)
                ##  }else{0}) +
                (if(XnO > 0){
                     v2 <- cbind(XtO,x[XiO,])
                     colSums(
                         log(array(
                             pnorm(q=Oright[v2],
                                   mean=mcsamples$Omean[XtO,,],
                                   sd=mcsamples$Ovar[XtO,,]) -
                             pnorm(q=Oleft[v2],
                                   mean=mcsamples$Omean[XtO,,],
                                   sd=mcsamples$Ovar[XtO,,]),
                             dim=c(XnO, nclusters, nsamples))),
                         na.rm=T)
                 }else{0}) +
                (if(XnN > 0){
                     colSums(
                         log(array(
                             t(sapply(XseqN, function(v){
                                 mcsamples$Nprob[XtN[v],,x[XiN[v],],]
                             })),
                             dim=c(XnN, nclusters, nsamples))),
                         na.rm=T)
                 }else{0}) +
                (if(XnB > 0){
                     colSums(
                         log(array(x[XiB,]*mcsamples$Bprob[XtB,,] +
                                   (1-x[XiB,])*(1-mcsamples$Bprob[XtB,,]),
                                   dim=c(XnB, nclusters, nsamples))),
                         na.rm=T)
                 }else{0})
            }
        probX <- apply(probX, 2, function(xx){xx - max(xx[is.finite(xx)])})
        ##
        ##
        if(all(is.na(y))){
            probY <- array(NA, dim=dim(W))
        }else{
            probY <- (if(YnR > 0){# continuous
                     colSums(
                         array(dnorm(x=y[YiR,],
                                     mean=mcsamples$Rmean[YtR,,],
                                     sd=mcsamples$Rvar[YtR,,],log=T),
                               dim=c(YnR, nclusters, nsamples)),
                         na.rm=T)
                 }else{0}) +
                (if(YnC > 0){# censored
                     colSums(
                         array(
                             t(sapply(YseqC, function(v){
                                 v1 <- y[YiC[v],]
                                 v2 <- YtC[v]
                                 if(is.finite(v1)){
                                     (dnorm(x=v1,
                                            mean=mcsamples$Cmean[v2,,],
                                            sd=mcsamples$Cvar[v2,,],log=T))
                                 }else{
                                     (pnorm(q=Cbounds[v2, 2L-(v1 < 0)],
                                            mean=mcsamples$Cmean[v2,,],
                                            sd=mcsamples$Cvar[v2,,],
                                            lower.tail=(v1 < 0),
                                            log.p=T))
                                 }
                             })),
                             dim=c(YnC, nclusters, nsamples)),
                         na.rm=T)
                 }else{0}) +
                (if(YnD > 0){# continuous discretized
                     colSums(
                         array(
                             t(sapply(YseqD, function(v){
                                 v1 <- y[YiD[v],]
                                 v2 <- YtD[v]
                                 if(is.finite(v1)){
                                     (dnorm(x=v1,
                                            mean=mcsamples$Dmean[v2,,],
                                            sd=mcsamples$Dvar[v2,,],log=T))
                                 }else{
                                     (pnorm(q=Dbounds[v2, 2L-(v1 < 0)],
                                            mean=mcsamples$Dmean[v2,,],
                                            sd=mcsamples$Dvar[v2,,],
                                            lower.tail=(v1 < 0),
                                            log.p=T))
                                 }
                             })),
                             dim=c(YnD, nclusters, nsamples)),
                         na.rm=T)
                 }else{0}) +
                ## (if(YnD > 0){# continuous
                ##      colSums(
                ##          array(dnorm(x=y[YiD,],
                ##                      mean=mcsamples$Dmean[YtD,,],
                ##                      sd=mcsamples$Dvar[YtD,,],log=T),
                ##                dim=c(YnD, nclusters, nsamples)),
                ##          na.rm=T)
                ##  }else{0}) +
                (if(YnO > 0){
                     v2 <- cbind(YtO,y[YiO,])
                     colSums(
                         log(array(
                             pnorm(q=Oright[v2],
                                   mean=mcsamples$Omean[YtO,,],
                                   sd=mcsamples$Ovar[YtO,,]) -
                             pnorm(q=Oleft[v2],
                                   mean=mcsamples$Omean[YtO,,],
                                   sd=mcsamples$Ovar[YtO,,]),
                             dim=c(YnO, nclusters, nsamples))),
                         na.rm=T)
                 }else{0}) +
                (if(YnN > 0){
                     colSums(
                         log(array(
                             t(sapply(YseqN, function(v){
                                 mcsamples$Nprob[YtN[v],,y[YiN[v],],]
                             })),
                             dim=c(YnN, nclusters, nsamples))),
                         na.rm=T)
                 }else{0}) +
                (if(YnB > 0){
                     colSums(
                         log(array(y[YiB,]*mcsamples$Bprob[YtB,,] +
                                   (1-y[YiB,])*(1-mcsamples$Bprob[YtB,,]),
                                   dim=c(YnB, nclusters, nsamples))),
                         na.rm=T)
                 }else{0})
            ## end probY
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

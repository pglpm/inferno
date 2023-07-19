samplesFDistribution <- function(Y, X, mcsamples, auxmetadata, subsamples, jacobian=TRUE, fn=identity, parallel=TRUE, useOquantiles=TRUE){
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

    ##
    allv <- union(Yv, Xv)
    vn <- vnames <- vindices <- list()
    ##    auxmetadata <- auxmetadata[name %in% allv]

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
            c(vtransform(x=matrix(NA,nrow=1,ncol=vn$C,
                                  dimnames=list(NULL,vnames$C)),
                         auxmetadata=auxmetadata,Cout='sleft')),
            ## sign is important here:
            ## for upper tail, take opposite mean and value
            -c(vtransform(x=matrix(NA,nrow=1,ncol=vn$C,
                                   dimnames=list(NULL,vnames$C)),
                          auxmetadata=auxmetadata,Cout='sright'))
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
        mcsamples$Dvar <- sqrt(mcsamples$Dvar)
        Dbounds <- cbind(
            c(vtransform(x=matrix(NA,nrow=1,ncol=vn$D,
                                  dimnames=list(NULL,vnames$D)),
                         auxmetadata=auxmetadata,Dout='sleft')),
            ## sign is important here:
            ## for upper tail, take opposite mean and value
            -c(vtransform(x=matrix(NA,nrow=1,ncol=vn$D,
                                   dimnames=list(NULL,vnames$D)),
                          auxmetadata=auxmetadata,Dout='sright'))
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
        mcsamples$Ovar <- sqrt(mcsamples$Ovar)
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


    
    for(atype in c('R','C','D','O','N','B')){
        ## 
        ## To save memory, we'll extract from mcsamples
        ## only sample parameters belonging to vnames
        vnames[[atype]] <- allv[allv %in% auxmetadata[mcmctype == atype, name]]
        vn[[atype]] <- length(vnames[[atype]])
        vindices[[atype]] <- sapply(vnames[[atype]], function(xx)auxmetadata[name == xx, id])
        ordering <- order(vindices[[atype]])
        vnames[[atype]] <- vnames[[atype]][ordering]
        vindices[[atype]] <- vindices[[atype]][ordering]
    }

    allparams <- colnames(mcsamples)

    ## W
    W <- mcsamples[,grep('^W', allparams),drop=F]
    nclusters <- nrow(mcsamples$W)

    if(vn$R > 0){# continuous
        inds <- paste0(vindices$R,collapse='|')
        Rmean <- array(t(mcsamples[,grep(paste0('^Rmean\\[(',inds,')'), allparams),drop=F]),
                       dim=c(vn$R,nclusters,nsamples), dimnames=NULL)
        Rvarsd <- sqrt(array(t(mcsamples[,grep(paste0('^Rvar\\[(',inds,')'), allparams),drop=F]),
                      dim=c(vn$R,nclusters,nsamples), dimnames=NULL))
        ##
        totake <- intersect(vnames$R, Yv)
YnR <- length(totake)
        YiR <- unname(sapply(totake, function(xx){which(Yv==xx)}))
        YtR <- unname(sapply(totake, function(xx){which(vnames$R==xx)}))
        ##
        totake <- intersect(vnames$R, Xv)
XnR <- length(totake)
        XiR <- unname(sapply(totake, function(xx){which(Xv==xx)}))
        XtR <- unname(sapply(totake, function(xx){which(vnames$R==xx)}))
    }else{
        YnR <- XnR <- 0
    }
    if(vn$C > 0){# censored
        inds <- paste0(vindices$C,collapse='|')
        Cmean <- array(t(mcsamples[,grep(paste0('^Cmean\\[(',inds,')'), allparams),drop=F]),
                       dim=c(vn$C,nclusters,nsamples), dimnames=NULL)
        Cvarsd <- sqrt(array(t(mcsamples[,grep(paste0('^Cvar\\[(',inds,')'), allparams),drop=F]),
                      dim=c(vn$C,nclusters,nsamples), dimnames=NULL))
        Cbounds <- cbind(
            c(vtransform(x=matrix(NA,nrow=1,ncol=vn$C,dimnames=list(NULL,vnames$C)),
                         auxmetadata=auxmetadata,Cout='sleft')),
            ## sign is important here:
            ## for upper tail, take opposite mean and value
            -c(vtransform(x=matrix(NA,nrow=1,ncol=vn$C,dimnames=list(NULL,vnames$C)),
                         auxmetadata=auxmetadata,Cout='sright'))
        )
        ##
        totake <- intersect(vnames$C, Yv)
        YnC <- length(totake)
        YiC <- unname(sapply(totake, function(xx){which(Yv==xx)}))
        YtC <- unname(sapply(totake, function(xx){which(vnames$C==xx)}))
        ##
        totake <- intersect(vnames$C, Xv)
        XnC <- length(totake)
        XiC <- unname(sapply(totake, function(xx){which(Xv==xx)}))
        XtC <- unname(sapply(totake, function(xx){which(vnames$C==xx)}))
        ##
        YseqC <- 1:YnC
        XseqC <- 1:XnC
    }else{
        YnC <- XnC <- 0
    }
    if(vn$D > 0){## discretized
        inds <- paste0(vindices$D,collapse='|')
        Dmean <- array(t(mcsamples[,grep(paste0('^Dmean\\[(',inds,')'), allparams),drop=F]),
                       dim=c(vn$D,nclusters,nsamples), dimnames=NULL)
        Dvarsd <- sqrt(array(t(mcsamples[,grep(paste0('^Dvar\\[(',inds,')'), allparams),drop=F]),
                      dim=c(vn$D,nclusters,nsamples), dimnames=NULL))
        Dbounds <- cbind(
            c(vtransform(x=matrix(NA,nrow=1,ncol=vn$D,dimnames=list(NULL,vnames$D)),
                         auxmetadata=auxmetadata,Dout='sleft')),
            ## sign is important here:
            ## for upper tail, take opposite mean and value
            -c(vtransform(x=matrix(NA,nrow=1,ncol=vn$D,dimnames=list(NULL,vnames$D)),
                         auxmetadata=auxmetadata,Dout='sright'))
        )
        ##
        totake <- intersect(vnames$D, Yv)
        YnD <- length(totake)
        YiD <- unname(sapply(totake, function(xx){which(Yv==xx)}))
        YtD <- unname(sapply(totake, function(xx){which(vnames$D==xx)}))
        ##
        totake <- intersect(vnames$D, Xv)
        XnD <- length(totake)
        XiD <- unname(sapply(totake, function(xx){which(Xv==xx)}))
        XtD <- unname(sapply(totake, function(xx){which(vnames$D==xx)}))
        ##
        YseqD <- 1:YnD
        XseqD <- 1:XnD
    }else{
        YnD <- XnD <- 0
    }
    if(vn$O > 0){# ordinal
        Qfunction <- readRDS('Qfunction8192.rds')
        inds <- paste0(vindices$O,collapse='|')
        Omean <- array(t(mcsamples[,grep(paste0('^Omean\\[(',inds,')'), allparams),drop=F]),
                       dim=c(vn$O,nclusters,nsamples), dimnames=NULL)
        Ovarsd <- sqrt(array(t(mcsamples[,grep(paste0('^Ovar\\[(',inds,')'), allparams),drop=F]),
                      dim=c(vn$O,nclusters,nsamples), dimnames=NULL))
        ##
        Omaxn <- max(auxmetadata[name %in% vnames$O, Nvalues])
        Oseq <- 1:vn$O
        ##
        Oleft <- t(sapply(vnames$O, function(avar){
            nn <- auxmetadata[name == avar, Nvalues]
            ## c(Qfunction((0:(nn-1))/nn), rep(NA,Omaxn-nn))
            c(vtransform(seq(auxmetadata[name == avar, domainmin],
                             auxmetadata[name == avar, domainmax],
                             length.out=auxmetadata[name == avar, Nvalues]
                             ),
                         auxmetadata, Oout='left', variates=avar,
                             useOquantiles=useOquantiles),
              rep(NA,Omaxn-nn))
        }))
        Oright <- t(sapply(vnames$O, function(avar){
            nn <- auxmetadata[name == avar, Nvalues]
            c(vtransform(seq(auxmetadata[name == avar, domainmin],
                             auxmetadata[name == avar, domainmax],
                             length.out=auxmetadata[name == avar, Nvalues]
                             ),
                         auxmetadata, Oout='right', variates=avar,
                             useOquantiles=useOquantiles),
                rep(NA,Omaxn-nn))
        }))
        ##
        totake <- intersect(vnames$O, Yv)
YnO <- length(totake)
        YiO <- unname(sapply(totake, function(xx){which(Yv==xx)}))
        YtO <- unname(sapply(totake, function(xx){which(vnames$O==xx)}))
        ##
        totake <- intersect(vnames$O, Xv)
XnO <- length(totake)
        XiO <- unname(sapply(totake, function(xx){which(Xv==xx)}))
        XtO <- unname(sapply(totake, function(xx){which(vnames$O==xx)}))
    }else{
        YnO <- XnO <- 0
    }
    if(vn$N > 0){# nominal
        Nmaxn <- max(auxmetadata[name %in% vnames$N, Nvalues])
        inds <- paste0(vindices$N,collapse='|')
        indn <- paste0(1:Nmaxn,collapse='|')
        Nprob <- array(t(mcsamples[,grep(paste0('^Nprob\\[(',inds,'), .*, (',indn,')\\]'), allparams),drop=F]),
                       dim=c(vn$N,nclusters,Nmaxn,nsamples), dimnames=NULL)
        ##
        totake <- intersect(vnames$N, Yv)
YnN <- length(totake)
        YiN <- unname(sapply(totake, function(xx){which(Yv==xx)}))
        YtN <- unname(sapply(totake, function(xx){which(vnames$N==xx)}))
        ##
        totake <- intersect(vnames$N, Xv)
XnN <- length(totake)
        XiN <- unname(sapply(totake, function(xx){which(Xv==xx)}))
        XtN <- unname(sapply(totake, function(xx){which(vnames$N==xx)}))
        ##
        YseqN <- 1:YnN
        XseqN <- 1:XnN
    }else{
        YnN <- XnN <- 0
    }
    if(vn$B > 0){## binary
        inds <- paste0(vindices$B,collapse='|')
        Bprob <- array(t(mcsamples[,grep(paste0('^Bprob\\[(',inds,')'), allparams),drop=F]),
                       dim=c(vn$B,nclusters,nsamples), dimnames=NULL)
        ##
        totake <- intersect(vnames$B, Yv)
YnB <- length(totake)
        YiB <- unname(sapply(totake, function(xx){which(Yv==xx)}))
        YtB <- unname(sapply(totake, function(xx){which(vnames$B==xx)}))
        ##
        totake <- intersect(vnames$B, Xv)
XnB <- length(totake)
        XiB <- unname(sapply(totake, function(xx){which(Xv==xx)}))
        XtB <- unname(sapply(totake, function(xx){which(vnames$B==xx)}))
    }else{
        YnB <- XnB <- 0
    }
    rm(mcsamples)

    ##
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
        ## ## for debugging
        ## for(iii in 1:nrow(Y2)){
        ## print(iii)
        ## iii <- iii+1
        ##     y <- t(Y2)[,1,drop=F]
        ##     x <- t(X2)[,1,drop=F]
        ##
        ##            
        ##
        if(all(is.na(x))){
            probX <- log(mcsamples$W) 
        }else{
            ## rows: clusters, cols: samples
            probX <- log(mcsamples$W) + 
            (if(XnR > 0){# continuous
                 colSums(
                     dnorm(x=x[XiR,],
                           mean=mcsamples$Rmean[XtR,,drop=F],
                           sd=mcsamples$Rvar[XtR,,drop=F],
                           log=T),
                     na.rm=T)
             }else{0}) +
            (if(XnC > 0){# censored
                 indf <- which(is.finite(x[XiC,]))
                 indi <- which(!is.finite(x[XiC,]))
                 (if(length(indf) > 0){
                      colSums(
                          dnorm(x=x[XiC[indf],],
                                mean=mcsamples$Cmean[XtC[indf],,drop=F],
                                sd=mcsamples$Cvar[XtC[indf],,drop=F],
                                log=T),
                          na.rm=T)
                  }else{0}) +
                     (if(length(indi) > 0){
                          v2 <- XtC[indi]
                          v1 <- -sign(x[XiC[indi],])
                          ## for upper tail, take opposite mean and value
                          colSums(
                              pnorm(q=Cbounds[v2, 1.5-0.5*v1],
                                    mean=v1*mcsamples$Cmean[v2,,drop=F],
                                    sd=mcsamples$Cvar[v2,,drop=F],
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
                                mean=mcsamples$Dmean[XtD[indf],,drop=F],
                                sd=mcsamples$Dvar[XtD[indf],,drop=F],
                                log=T),
                          na.rm=T)
                  }else{0}) +
                     (if(length(indi) > 0){
                          v2 <- XtD[indi]
                          v1 <- -sign(x[XiD[indi],])
                          ## for upper tail, take opposite mean and value
                          colSums(
                              pnorm(q=Dbounds[v2, 1.5-0.5*v1],
                                    mean=v1*mcsamples$Dmean[v2,,drop=F],
                                    sd=mcsamples$Dvar[v2,,drop=F],
                                    log.p=T),
                              na.rm=T)
                      }else{0})
             }else{0}) +
                (if(XnO > 0){# ordinal
                     v2 <- cbind(XtO,x[XiO,])
                     colSums(
                         log(
                             pnorm(q=Oright[v2],
                                   mean=mcsamples$Omean[XtO,,drop=F],
                                   sd=mcsamples$Ovarsd[XtO,,drop=F]) -
                             pnorm(q=Oleft[v2],
                                   mean=mcsamples$Omean[XtO,,drop=F],
                                   sd=mcsamples$Ovar[XtO,,drop=F])
                         ),
                         na.rm=T)
                 }else{0}) +
                (if(XnN > 0){# nominal
                     temp <- apply(mcsamples$Nprob, c(2,4), function(xx){
                         xx[cbind(XtN, x[XiN,])]
                     })
                     dim(temp) <- c(XnN, nclusters, nsamples)
                     ##
                     colSums(log(temp), na.rm=T)
                 }else{0}) +
                (if(XnB > 0){# binary
                     colSums(
                         log( x[XiB,]*mcsamples$Bprob[XtB,,drop=F] +
                              (1-x[XiB,])*(1-mcsamples$Bprob[XtB,,drop=F]) ),
                         na.rm=T)
                 }else{0})
        } # end probX
        ##
        if(all(is.na(y))){
            probY <- array(NA, dim=dim(W))
        }else{
            (if(YnR > 0){# continuous
                 colSums(
                     dnorm(x=y[YiR,],
                           mean=mcsamples$Rmean[YtR,,drop=F],
                           sd=mcsamples$Rvar[YtR,,drop=F],
                           log=T),
                     na.rm=T)
             }else{0}) +
            (if(YnC > 0){# censored
                 indf <- which(is.finite(y[YiC,]))
                 indi <- which(!is.finite(y[YiC,]))
                 (if(length(indf) > 0){
                      colSums(
                          dnorm(x=y[YiC[indf],],
                                mean=mcsamples$Cmean[YtC[indf],,drop=F],
                                sd=mcsamples$Cvar[YtC[indf],,drop=F],
                                log=T),
                          na.rm=T)
                  }else{0}) +
                     (if(length(indi) > 0){
                          v2 <- YtC[indi]
                          v1 <- -sign(y[YiC[indi],])
                          ## for upper tail, take opposite mean and value
                          colSums(
                              pnorm(q=Cbounds[v2, 1.5-0.5*v1],
                                    mean=v1*mcsamples$Cmean[v2,,drop=F],
                                    sd=mcsamples$Cvar[v2,,drop=F],
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
                                mean=mcsamples$Dmean[YtD[indf],,drop=F],
                                sd=mcsamples$Dvar[YtD[indf],,drop=F],
                                log=T),
                          na.rm=T)
                  }else{0}) +
                     (if(length(indi) > 0){
                          v2 <- YtD[indi]
                          v1 <- -sign(y[YiD[indi],])
                          ## for upper tail, take opposite mean and value
                          colSums(
                              pnorm(q=Dbounds[v2, 1.5-0.5*v1],
                                    mean=v1*mcsamples$Dmean[v2,,drop=F],
                                    sd=mcsamples$Dvar[v2,,drop=F],
                                    log.p=T),
                              na.rm=T)
                      }else{0})
             }else{0}) +
                (if(YnO > 0){# ordinal
                     v2 <- cbind(YtO,y[YiO,])
                     colSums(
                         log(
                             pnorm(q=Oright[v2],
                                   mean=mcsamples$Omean[YtO,,drop=F],
                                   sd=mcsamples$Ovarsd[YtO,,drop=F]) -
                             pnorm(q=Oleft[v2],
                                   mean=mcsamples$Omean[YtO,,drop=F],
                                   sd=mcsamples$Ovar[YtO,,drop=F])
                         ),
                         na.rm=T)
                 }else{0}) +
                (if(YnN > 0){# nominal
                     temp <- apply(mcsamples$Nprob, c(2,4), function(xx){
                         xx[cbind(YtN, y[YiN,])]
                     })
                     dim(temp) <- c(YnN, nclusters, nsamples)
                     ##
                     colSums(log(temp), na.rm=T)
                 }else{0}) +
                (if(YnB > 0){# binary
                     colSums(
                         log( y[YiB,]*mcsamples$Bprob[YtB,,drop=F] +
                              (1-y[YiB,])*(1-mcsamples$Bprob[YtB,,drop=F]) ),
                         na.rm=T)
                 }else{0})








            probY <- t( # rows: MCsamples, cols: clusters
                (if(YnR > 0){# continuous
                     colSums(
                         dnorm(x=y[YiR,],
                                     mean=mcsamples$Rmean[YtR,,drop=F],
                                     sd=mcsamples$Rvar[YtR,,drop=F]),log=T,
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
                                            mean=Cmean[v2,,],
                                            sd=Cvarsd[v2,,],log=T))
                                 }else{
                                     (pnorm(q=Cbounds[v2, 2L-(v1 < 0)],
                                            mean=Cmean[v2,,],
                                            sd=Cvarsd[v2,,],
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
                                            mean=Dmean[v2,,],
                                            sd=Dvarsd[v2,,],log=T))
                                 }else{
                                     (pnorm(q=Dbounds[v2, 2L-(v1 < 0)],
                                            mean=Dmean[v2,,],
                                            sd=Dvarsd[v2,,],
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
                ##                      mean=Dmean[YtD,,],
                ##                      sd=Dvarsd[YtD,,],log=T),
                ##                dim=c(YnD, nclusters, nsamples)),
                ##          na.rm=T)
                ##  }else{0}) +
                (if(YnO > 0){
                     v2 <- cbind(YtO,y[YiO,])
                     colSums(
                         log(array(
                             pnorm(q=Oright[v2],
                                   mean=Omean[YtO,,],
                                   sd=Ovarsd[YtO,,]) -
                             pnorm(q=Oleft[v2],
                                   mean=Omean[YtO,,],
                                   sd=Ovarsd[YtO,,]),
                             dim=c(YnO, nclusters, nsamples))),
                         na.rm=T)
                 }else{0}) +
                (if(YnN > 0){
                     colSums(
                         log(array(
                             t(sapply(YseqN, function(v){
                                 Nprob[YtN[v],,y[YiN[v],],]
                             })),
                             dim=c(YnN, nclusters, nsamples))),
                         na.rm=T)
                 }else{0}) +
                (if(YnB > 0){
                     colSums(
                         log(array(y[YiB,]*Bprob[YtB,,] +
                                   (1-y[YiB,])*(1-Bprob[YtB,,]),
                                   dim=c(YnB, nclusters, nsamples))),
                         na.rm=T)
                 }else{0})
                ) # end probY
        }
#### Output: rows=samples, columns=clusters
        ##
        ## if(all(is.na(x))){
        ##     out <- rowSums(exp(probX+probY)) 
        ## }else{
        ## str(probX)
        ## print(any(is.na(probX)))
        ## print(apply(probX, 1, max, na.rm=T))
        probX <- probX - apply(probX, 1, function(xx){max(xx[is.finite(xx)])})
        ## str(probX)
        ## print(any(is.na(probX)))
        ## str(probY)
        ## print(any(is.na(probY)))
        ## }
        fn( rowSums(exp(probX+probY))/rowSums(exp(probX)) )
        ## rbind(log(W),probX,probY)
    } *
        (if(jacobian){
             exp(-rowSums(
                      log(vtransform(Y, auxmetadata=auxmetadata, invjacobian=TRUE, useOquantiles=useOquantiles)),
                      na.rm=T))}else{1L})
}

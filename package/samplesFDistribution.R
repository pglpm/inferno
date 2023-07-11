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
    Xv <- colnames(X)
    if(!all(Xv %in% auxmetadata$name)){stop('unknown X variates\n')}
    if(length(intersect(Yv, Xv)) > 0){stop('overlap in Y and X variates\n')}

    ## mcsamples and subsamples
    if(is.character(mcsamples) && file.exists(mcsamples)){
        mcsamples <- readRDS(mcsamples)
    }
    if(missing(subsamples) || is.null(subsamples) || (is.logical(subsamples) && !subsamples)){
        subsamples <- 1:nrow(mcsamples)
    }else if(is.character(subsamples)){
        subsamples <- round(seq(1,nrow(subsamples),length.out=as.numeric(subsamples)))
    }

    mcsamples <- mcsamples[subsamples,,drop=F]
    nsamples <- nrow(mcsamples)

    ##
    allv <- union(Yv, Xv)
    vn <- vnames <- vindices <- list()
    ##    auxmetadata <- auxmetadata[name %in% allv]

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
    nclusters <- ncol(W)

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
            c(vtransform(x=matrix(NA,nrow=1,ncol=vn$C,dimnames=list(NULL,vnames$C)),
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
            c(vtransform(x=matrix(NA,nrow=1,ncol=vn$D,dimnames=list(NULL,vnames$D)),
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
            probX <- log(W) 
        }else{
            probX <- 0*log(W) +
                t( # rows: MCsamples, cols: clusters
                (if(XnR > 0){# continuous
                     colSums(
                         array(dnorm(x=x[XiR,],
                                     mean=Rmean[XtR,,],
                                     sd=Rvarsd[XtR,,],log=T),
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
                             dim=c(XnD, nclusters, nsamples)),
                         na.rm=T)
                 }else{0}) +
                ## (if(XnD > 0){# continuous
                ##      colSums(
                ##          array(dnorm(x=x[XiD,],
                ##                      mean=Dmean[XtD,,],
                ##                      sd=Dvarsd[XtD,,],log=T),
                ##                dim=c(XnD, nclusters, nsamples)),
                ##          na.rm=T)
                ##  }else{0}) +
                (if(XnO > 0){
                     v2 <- cbind(XtO,x[XiO,])
                     colSums(
                         log(array(
                             pnorm(q=Oright[v2],
                                   mean=Omean[XtO,,],
                                   sd=Ovarsd[XtO,,]) -
                             pnorm(q=Oleft[v2],
                                   mean=Omean[XtO,,],
                                   sd=Ovarsd[XtO,,]),
                             dim=c(XnO, nclusters, nsamples))),
                         na.rm=T)
                 }else{0}) +
                (if(XnN > 0){
                     colSums(
                         log(array(
                             t(sapply(XseqN, function(v){
                                 Nprob[XtN[v],,x[XiN[v],],]
                             })),
                             dim=c(XnN, nclusters, nsamples))),
                         na.rm=T)
                 }else{0}) +
                (if(XnB > 0){
                     colSums(
                         log(array(x[XiB,]*Bprob[XtB,,] +
                                   (1-x[XiB,])*(1-Bprob[XtB,,]),
                                   dim=c(XnB, nclusters, nsamples))),
                         na.rm=T)
                 }else{0})
                ) +
                log(W)
        } # end probX
        ##
        if(all(is.na(y))){
            probY <- array(NA, dim=dim(W))
        }else{
            probY <- t( # rows: MCsamples, cols: clusters
                (if(YnR > 0){# continuous
                     colSums(
                         array(dnorm(x=y[YiR,],
                                     mean=Rmean[YtR,,],
                                     sd=Rvarsd[YtR,,],log=T),
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

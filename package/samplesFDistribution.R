samplesFDistribution <- function(Y, X=NULL, mcsamples, varinfoaux, subsamples=NULL, jacobian=TRUE, fn=identity){
    ## Consistency checks
    if(length(dim(Y)) != 2){stop('Y must have two dimensions')}
    if(!is.null(X) && length(dim(X)) != 2){stop('X must be NULL or have two dimensions')}
    ##
    if(!is.null(X) && ncol(X) == 0){X <- NULL}
    ## varinfoaux
    if(is.character(varinfoaux) && file.exists(varinfoaux)){
        varinfoaux <- readRDS(varinfoaux)
    }

    ## More consistency checks
    Yv <- colnames(Y)
    if(!all(Yv %in% varinfoaux$name)){stop('unknown Y variates\n')}
    Xv <- colnames(X)
    if(!all(Xv %in% varinfoaux$name)){stop('unknown X variates\n')}
    if(length(intersect(Yv, Xv)) > 0){stop('overlap in Y and X variates\n')}

    ## mcsamples and subsamples
    if(is.character(mcsamples) && file.exists(mcsamples)){
        mcsamples <- readRDS(mcsamples)
    }
    if(is.null(subsamples) || (is.logical(subsamples) && !subsamples)){
        subsamples <- 1:nrow(mcsamples)
    }else if(is.character(subsamples)){
        subsamples <- round(seq(1,nrow(subsamples),length.out=as.numeric(subsamples)))
    }

    mcsamples <- mcsamples[subsamples,,drop=F]
    nsamples <- nrow(mcsamples)

    ##
    allv <- union(Yv, Xv)
    vn <- vnames <- vindices <- list()
    ##    varinfoaux <- varinfoaux[name %in% allv]
    for(atype in c('R','C','D','O','N','B')){
        ## 
        ## To save memory, we'll extract from mcsamples
        ## only sample parameters belonging to vnames
        vnames[[atype]] <- allv[allv %in% varinfoaux[mcmctype == atype, name]]
        vn[[atype]] <- length(vnames[[atype]])
        vindices[[atype]] <- sapply(vnames[[atype]], function(xx)varinfoaux[name == xx, id])
        ordering <- order(vindices[[atype]])
        vnames[[atype]] <- vnames[[atype]][ordering]
        vindices[[atype]] <- vindices[[atype]][ordering]
    }

    allparams <- colnames(mcsamples)

    ## W
    W <- mcsamples[,grep('^W', allparams)]
    nclusters <- ncol(W)

    if(vn$R > 0){# continuous
        inds <- paste0(vindices$R,collapse='|')
        Rmean <- array(t(mcsamples[,grep(paste0('^Rmean\\[(',inds,')'), allparams),drop=F]),
                       dim=c(vn$R,nclusters,nsamples), dimnames=list(vnames$R,NULL))
        Rvar <- array(t(mcsamples[,grep(paste0('^Rvar\\[(',inds,')'), allparams),drop=F]),
                      dim=c(vn$R,nclusters,nsamples), dimnames=list(vnames$R,NULL))
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
                       dim=c(vn$C,nclusters,nsamples), dimnames=list(vnames$C,NULL))
        Cvar <- array(t(mcsamples[,grep(paste0('^Cvar\\[(',inds,')'), allparams),drop=F]),
                      dim=c(vn$C,nclusters,nsamples), dimnames=list(vnames$C,NULL))
        Cbounds <- cbind(
            c(vtransform(x=matrix(NA,nrow=1,ncol=vn$C,dimnames=list(NULL,vnames$C)),
                         varinfoaux=varinfoaux,variates=vnames$C,Cout='sleft')),
            c(vtransform(x=matrix(NA,nrow=1,ncol=vn$C,dimnames=list(NULL,vnames$C)),
                         varinfoaux=varinfoaux,variates=vnames$C,Cout='sright'))
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
                       dim=c(vn$D,nclusters,nsamples), dimnames=list(vnames$D,NULL))
        Dvar <- array(t(mcsamples[,grep(paste0('^Dvar\\[(',inds,')'), allparams),drop=F]),
                      dim=c(vn$D,nclusters,nsamples), dimnames=list(vnames$D,NULL))
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
    }else{
        YnD <- XnD <- 0
    }
    if(vn$O > 0){# ordinal
        Qfunction <- readRDS('Qfunction512.rds')
        inds <- paste0(vindices$O,collapse='|')
        Omean <- array(t(mcsamples[,grep(paste0('^Omean\\[(',inds,')'), allparams),drop=F]),
                       dim=c(vn$O,nclusters,nsamples), dimnames=list(vnames$O,NULL))
        Ovar <- array(t(mcsamples[,grep(paste0('^Ovar\\[(',inds,')'), allparams),drop=F]),
                      dim=c(vn$O,nclusters,nsamples), dimnames=list(vnames$O,NULL))
        ##
        Omaxn <- max(varinfoaux[name %in% vnames$O, Nvalues])
        Oseq <- 1:vn$O
        ##
        Oleft <- t(sapply(vnames$O, function(avar){
            nn <- varinfoaux[name == avar, Nvalues]
            c(Qfunction((0:(nn-1))/nn), rep(NA,Omaxn-nn))
        }))
        Oright <- t(sapply(vnames$O, function(avar){
            nn <- varinfoaux[name == avar, Nvalues]
            c(Qfunction((1:nn)/nn), rep(NA,Omaxn-nn))
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
        Nmaxn <- max(varinfoaux[name %in% vnames$N, Nvalues])
        inds <- paste0(vindices$N,collapse='|')
        indn <- paste0(1:Nmaxn,collapse='|')
        Nprob <- array(t(mcsamples[,grep(paste0('^Nprob\\[(',inds,'), .*, (',indn,')\\]'), allparams),drop=F]),
                       dim=c(vn$N,nclusters,Nmaxn,nsamples), dimnames=list(vnames$N,NULL,NULL,NULL))
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
                       dim=c(vn$B,nclusters,nsamples), dimnames=list(vnames$B,NULL))
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
    Y2 <- vtransform(Y, varinfoaux, Cout='index', Dout='', Oout='', Nout='numeric', Bout='numeric')
    if(!is.null(X)){
        X2 <- vtransform(X, varinfoaux, Cout='index', Dout='', Oout='', Nout='numeric', Bout='numeric')
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
    foreach(y=t(Y2), x=t(X2), .combine=rbind, .inorder=T)%do%{
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
            probX <- 0*log(W) 
        }else{
            probX <- 0*log(W) +
                t( # rows: MCsamples, cols: clusters
                (if(XnR > 0){# continuous
                     colSums(
                         array(dnorm(x=x[XiR,],
                                     mean=Rmean[XtR,,],
                                     sd=sqrt(Rvar[XtR,,]),log=T),
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
                                            sd=sqrt(Cvar[v2,,]),log=T))
                                 }else{
                                     (pnorm(q=Cbounds[v2, 2L-(v1 < 0)],
                                            mean=Cmean[v2,,],
                                            sd=sqrt(Cvar[v2,,]),
                                            lower.tail=(v1 < 0),
                                            log.p=T))
                                 }
                             })),
                             dim=c(XnC, nclusters, nsamples)),
                         na.rm=T)
                 }else{0}) +
                (if(XnD > 0){# continuous
                     colSums(
                         array(dnorm(x=x[XiD,],
                                     mean=Dmean[XtD,,],
                                     sd=sqrt(Dvar[XtD,,]),log=T),
                               dim=c(XnD, nclusters, nsamples)),
                         na.rm=T)
                 }else{0}) +
                (if(XnO > 0){
                     v2 <- cbind(XtO,x[XiO,])
                     colSums(
                         log(array(
                             pnorm(q=Oright[v2],
                                   mean=Omean[XtO,,],
                                   sd=sqrt(Ovar[XtO,,])) -
                             pnorm(q=Oleft[v2],
                                   mean=Omean[XtO,,],
                                   sd=sqrt(Ovar[XtO,,])),
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
                )
        }# end probX
        ##
        if(all(is.na(y))){
            probY <- NA
        }else{
            probY <- t( # rows: MCsamples, cols: clusters
                (if(YnR > 0){# continuous
                     colSums(
                         array(dnorm(x=y[YiR,],
                                     mean=Rmean[YtR,,],
                                     sd=sqrt(Rvar[YtR,,]),log=T),
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
                                            sd=sqrt(Cvar[v2,,]),log=T))
                                 }else{
                                     (pnorm(q=Cbounds[v2, 2L-(v1 < 0)],
                                            mean=Cmean[v2,,],
                                            sd=sqrt(Cvar[v2,,]),
                                            lower.tail=(v1 < 0),
                                            log.p=T))
                                 }
                             })),
                             dim=c(YnC, nclusters, nsamples)),
                         na.rm=T)
                 }else{0}) +
                (if(YnD > 0){# continuous
                     colSums(
                         array(dnorm(x=y[YiD,],
                                     mean=Dmean[YtD,,],
                                     sd=sqrt(Dvar[YtD,,]),log=T),
                               dim=c(YnD, nclusters, nsamples)),
                         na.rm=T)
                 }else{0}) +
                (if(YnO > 0){
                     v2 <- cbind(YtO,y[YiO,])
                     colSums(
                         log(array(
                             pnorm(q=Oright[v2],
                                   mean=Omean[YtO,,],
                                   sd=sqrt(Ovar[YtO,,])) -
                             pnorm(q=Oleft[v2],
                                   mean=Omean[YtO,,],
                                   sd=sqrt(Ovar[YtO,,])),
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
        ##probX <- probX - apply(probX, 1, function(xx){max(xx[is.finite(xx)])})
        ## str(probX)
        ## print(any(is.na(probX)))
        ## str(probY)
        ## print(any(is.na(probY)))
        ## }
        ##        fn( rowSums(exp(probX+probY))/rowSums(exp(probX)) )
        rbind(log(W),probX,probY)
    } *
        (if(jacobian){
             exp(-rowSums(
                      log(vtransform(Y, varinfoaux=varinfoaux, invjacobian=TRUE)),
                      na.rm=T))}else{1L})
}

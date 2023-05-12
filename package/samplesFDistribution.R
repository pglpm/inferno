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

    allparams <- colnames(mcsamples)
    ##
    vn <- vnames <- Yt <- Xt <- Yi <- Xi <- Yn <- Xn <- Yseq <- Xseq <- vseq <- list()
    varinfoaux <- varinfoaux[name %in% c(Yv,Xv)]
    for(atype in c('R','C','D','O','N','B')){
        vn[[atype]] <- length(varinfoaux[mcmctype == atype, name])
        vnames[[atype]] <- varinfoaux[mcmctype == atype, name]
        ##
        common <- Yv[Yv %in% vnames$R]
        Yt[[atype]] <- which(Yv %in% common)
        Yi[[atype]] <- sapply(common, function(xx)varinfoaux[name == xx, id])
        Yn[[atype]] <- length(Yi[[atype]])
        Yseq[[atype]] <- 1:Yn[[atype]]
        ##
        Xt[[atype]] <- sapply(intersect(vnames[[atype]], Xv),
                              function(xx)which(Xv==xx))
        Xi[[atype]] <- sapply(intersect(vnames[[atype]], Xv),
                              function(xx)varinfoaux[name==xx,id])
        Xn[[atype]] <- length(Xi[[atype]])
        Xseq[[atype]] <- 1:Xn[[atype]]

        print(atype)
        print(Xn[[atype]])
        print(Xv[Xt[[atype]]])
        print(varinfoaux[mcmctype==atype,name][Xi[[atype]]])
        print(Yn[[atype]])
        print(Yv[Yt[[atype]]])
        print(varinfoaux[mcmctype==atype,name][Yi[[atype]]])
    }

    
    ##    
    if(vn$N > 0){
        Nmaxn <- max(varinfoaux[mcmctype == 'N', Nvalues])
    }
    ##    
    if(vn$O > 0){
        Omaxn <- max(varinfoaux[mcmctype == 'O', Nvalues])
        Oseq <- 1:vn$O
    }

    ## W
    W <- mcsamples[,grep('^W', allparams)]
    nclusters <- ncol(W)

    if(vn$R > 0){# continuous
        Rmean <- array(t(mcsamples[,grep('^Rmean', allparams),drop=F]),
                       dim=c(vn$R,nclusters,nsamples), dimnames=list(vnames$R,NULL))
        Rvar <- array(t(mcsamples[,grep('^Rvar', allparams),drop=F]),
                      dim=c(vn$R,nclusters,nsamples), dimnames=list(vnames$R,NULL))
    }
    if(vn$C > 0){# censored
        Cmean <- array(t(mcsamples[,grep('^Cmean', allparams),drop=F]),
                       dim=c(vn$C,nclusters,nsamples), dimnames=list(vnames$C,NULL))
        Cvar <- array(t(mcsamples[,grep('^Cvar', allparams),drop=F]),
                      dim=c(vn$C,nclusters,nsamples), dimnames=list(vnames$C,NULL))
        Cbounds <- cbind(
            c(vtransform(x=rbind(rep(NA,vn$C)),varinfoaux=varinfoaux,variates=vnames$C,Cout='sleft')),
            c(vtransform(x=rbind(rep(NA,vn$C)),varinfoaux=varinfoaux,variates=vnames$C,Cout='sright'))
        )
    }
    if(vn$D > 0){## discretized
        Dmean <- array(t(mcsamples[,grep('^Dmean', allparams),drop=F]),
                       dim=c(vn$D,nclusters,nsamples), dimnames=list(vnames$D,NULL))
        Dvar <- array(t(mcsamples[,grep('^Dvar', allparams),drop=F]),
                      dim=c(vn$D,nclusters,nsamples), dimnames=list(vnames$D,NULL))
    }
    if(vn$O > 0){# ordinal
        Qfunction <- readRDS('Qfunction512.rds')
        Omean <- array(t(mcsamples[,grep('^Omean', allparams),drop=F]),
                       dim=c(vn$O,nclusters,nsamples), dimnames=list(vnames$O,NULL))
        Ovar <- array(t(mcsamples[,grep('^Ovar', allparams),drop=F]),
                      dim=c(vn$O,nclusters,nsamples), dimnames=list(vnames$O,NULL))
        Oleft <- t(sapply(vnames$O, function(avar){
            nn <- varinfoaux[name == avar, Nvalues]
            c(Qfunction((0:(nn-1))/nn), rep(NA,Omaxn-nn))
        }))
        Oright <- t(sapply(vnames$O, function(avar){
            nn <- varinfoaux[name == avar, Nvalues]
            c(Qfunction((1:nn)/nn), rep(NA,Omaxn-nn))
        }))
    }
    if(vn$N > 0){# nominal
        Nprob <- array(t(mcsamples[,grep('^Nprob', allparams),drop=F]),
                       dim=c(vn$N,nclusters,Nmaxn,nsamples), dimnames=list(vnames$N,NULL))
    }
    if(vn$B > 0){## binary
        Bprob <- array(t(mcsamples[,grep('^Bprob', allparams),drop=F]),
                       dim=c(vn$B,nclusters,nsamples), dimnames=list(vnames$B,NULL))
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
    foreach(y=t(Y2), x=t(X2), .combine=rbind, .inorder=T)%dopar%{
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
            probX <- log(W) +
                t( # rows: MCsamples, cols: clusters
                (if(Xn$R > 0){# continuous
                     colSums(
                         array(dnorm(x=x[Xt$R,],
                                     mean=Rmean[Xi$R,,],
                                     sd=sqrt(Rvar[Xi$R,,]),log=T),
                               dim=c(Xn$R, nclusters, nsamples)),
                         na.rm=T)
                 }else{0}) +
                (if(Xn$C > 0){# censored
                     colSums(
                         array(
                             t(sapply(Xseq$C, function(v){
                                 v1 <- x[Xt$C[v],]
                                 v2 <- Xi$C[v]
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
                             dim=c(Xn$C, nclusters, nsamples)),
                         na.rm=T)
                 }else{0}) +
                (if(Xn$D > 0){# continuous
                     colSums(
                         array(dnorm(x=x[Xt$D,],
                                     mean=Dmean[Xi$D,,],
                                     sd=sqrt(Dvar[Xi$D,,]),log=T),
                               dim=c(Xn$D, nclusters, nsamples)),
                         na.rm=T)
                 }else{0}) +
                (if(Xn$O > 0){
                     v2 <- cbind(Oseq,x[Xt$O,])
                     colSums(
                         log(array(
                             pnorm(q=Oright[v2],
                                   mean=Omean[Xi$O,,],
                                   sd=sqrt(Ovar[Xi$O,,])) -
                             pnorm(q=Oleft[v2],
                                   mean=Omean[Xi$O,,],
                                   sd=sqrt(Ovar[Xi$O,,])),
                             dim=c(Xn$O, nclusters, nsamples))),
                         na.rm=T)
                 }else{0}) +
                (if(Xn$N > 0){
                     colSums(
                         log(array(
                             t(sapply(Xseq$N, function(v){
                                 Nprob[Xi$N[v],,x[Xt$N[v],],]
                             })),
                             dim=c(Xn$N, nclusters, nsamples))),
                         na.rm=T)
                 }else{0}) +
                (if(Xn$B > 0){
                     colSums(
                         log(array(x[Xt$B,]*Bprob[Xi$B,,] +
                                   (1-x[Xt$B,])*(1-Bprob[Xi$B,,]),
                                   dim=c(Xn$B, nclusters, nsamples))),
                         na.rm=T)
                 }else{0})
                )
        }# end probX
        ##
        if(all(is.na(y))){
            probY <- NA
        }else{
            probY <- t( # rows: MCsamples, cols: clusters
                (if(Yn$R > 0){# continuous
                     colSums(
                         array(dnorm(x=y[Yt$R,],
                                     mean=Rmean[Yi$R,,],
                                     sd=sqrt(Rvar[Yi$R,,]),log=T),
                               dim=c(Yn$R, nclusters, nsamples)),
                         na.rm=T)
                 }else{0}) +
                (if(Yn$C > 0){# censored
                     colSums(
                         array(
                             t(sapply(Yseq$C, function(v){
                                 v1 <- y[Yt$C[v],]
                                 v2 <- Yi$C[v]
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
                             dim=c(Yn$C, nclusters, nsamples)),
                         na.rm=T)
                 }else{0}) +
                (if(Yn$D > 0){# continuous
                     colSums(
                         array(dnorm(x=y[Yt$D,],
                                     mean=Dmean[Yi$D,,],
                                     sd=sqrt(Dvar[Yi$D,,]),log=T),
                               dim=c(Yn$D, nclusters, nsamples)),
                         na.rm=T)
                 }else{0}) +
                (if(Yn$O > 0){
                     v2 <- cbind(Oseq,y[Yt$O,])
                     colSums(
                         log(array(
                             pnorm(q=Oright[v2],
                                   mean=Omean[Yi$O,,],
                                   sd=sqrt(Ovar[Yi$O,,])) -
                             pnorm(q=Oleft[v2],
                                   mean=Omean[Yi$O,,],
                                   sd=sqrt(Ovar[Yi$O,,])),
                             dim=c(Yn$O, nclusters, nsamples))),
                         na.rm=T)
                 }else{0}) +
                (if(Yn$N > 0){
                     colSums(
                         log(array(
                             t(sapply(Yseq$N, function(v){
                                 Nprob[Yi$N[v],,y[Yt$N[v],],]
                             })),
                             dim=c(Yn$N, nclusters, nsamples))),
                         na.rm=T)
                 }else{0}) +
                (if(Yn$B > 0){
                     colSums(
                         log(array(y[Yt$B,]*Bprob[Yi$B,,] +
                                   (1-y[Yt$B,])*(1-Bprob[Yi$B,,]),
                                   dim=c(Yn$B, nclusters, nsamples))),
                         na.rm=T)
                 }else{0})
                ) # end probY
        }
        ##
        ## if(all(is.na(x))){
        ##     out <- rowSums(exp(probX+probY)) 
        ## }else{
        probX <- probX - apply(probX, 1, max, na.rm=T)
        ## }
        fn( rowSums(exp(probX+probY))/rowSums(exp(probX)) )
    } *
        (if(jacobian){
             exp(-rowSums(
                      log(vtransform(Y, varinfoaux=varinfoaux, invjacobian=TRUE)),
                      na.rm=T))}else{1L})
}

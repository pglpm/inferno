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
    varinfoaux <- varinfoaux[name %in% c(Yv,Xv)]

    

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
    vn <- vnames <- Yt <- Xt <- list()
    vnames <- list()
    for(atype in c('R','C','D','O','N','B')){
        vn[[atype]] <- length(varinfoaux[mcmctype == atype, name])
        vseq[[atype]] <- 1:vn[[atype]]
        vnames[[atype]] <- varinfoaux[mcmctype == atype, name]
        ## these keep the column indices in Y,X
        ## of the variates 'atype', in the order of varinfoaux
        Yt[[atype]] <- sapply(vnames[[atype]],function(xx)which(Yv==xx))
        Xt[[atype]] <- sapply(vnames[[atype]],function(xx)which(Xv==xx))
        Yn[[atype]] <- length(Yt[[atype]])
        Xn[[atype]] <- length(Xt[[atype]])
        ## these keep the column indices in the list of variates 'atype'
        ## present in Y,X
        Yi[[atype]] <- which(vnames[[atype]] %in% Yv)
        Xi[[atype]] <- which(vnames[[atype]] %in% Xv)
    }
    ##    
    if(vn$N > 0){
        Nmaxn <- max(varinfoaux[mcmctype == 'N', Nvalues])
    }
    ##    
    if(vn$O > 0){
        Omaxn <- max(varinfoaux[mcmctype == 'O', Nvalues])
        Oseq <- 1:Omaxn
    }

    ## W
    W <- mcsamples[,grep('^W', allparams)]
    nclusters <- length(W)

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
    if(vn$B > 0){## binary
        Bprob <- array(t(mcsamples[,grep('^Bprob', allparams),drop=F]),
                       dim=c(vn$B,nclusters,nsamples), dimnames=list(vnames$B,NULL))
    }





    test <- array(grep('^Rvar', colnames(mcsamples)),
                  dim=c(vn$R,nclusters), dimnames=list(vnames$R,NULL))



    mcsamples <- t(mcsamples[subsamples,,drop=FALSE])
    ##



    Yvn <- length(Yv)
    Xv <- colnames(X)
    Xvn <- length(Xv)
    if(length(intersect(Yv, Xv)) > 0){cat('WARNING: overlap in Y and X variates\n')}
    if(!all(Xv %in% Vv)){cat('Warning: unknown X variates\n')}
    ##
    variate <- lapply(variatetypes, function(x)names(varinfo[['type']])[varinfo[['type']]==x])
    len <- lapply(variate,length)
    names(variate) <- names(len) <- variatetypes
    ##
    Wi <- grep('W',rownames(mcsamples))
    nclusters <- length(Wi)
    W <- mcsamples[Wi,]
    ## seqclusters <- seq_len(nclusters)
    ##
    Yv <- lapply(variatetypes, function(zzz){out <- Yv[varinfo[['type']][Yv]==zzz]
        names(out) <- NULL
        out})
    Yn <- lapply(Yv,length)
    names(Yv) <- names(Yn) <- variatetypes
    ##
    Xv <- lapply(variatetypes, function(zzz){out <- Xv[varinfo[['type']][Xv]==zzz]
        names(out) <- NULL
        out})
    Xn <- lapply(Xv,length)
    names(Xv) <- names(Xn) <- variatetypes
    ##
    Vv <- lapply(variatetypes, function(zzz){out <- Vv[varinfo[['type']][Vv]==zzz]
        names(out) <- NULL
        out})
    Vn <- lapply(Vv,length)
    names(Vv) <- names(Vn) <- variatetypes
    ##
    ##
    if(len$C > 0){## categorical to be redone
        Cprob <- sapply(paste0('Cprob\\['),
                        grep,rownames(mcsamples))
        ncategories <- length(Cprob)/len$C/nclusters
        seqcategories <- seq_len(ncategories)
        Cprob <- aperm(array(Cprob,
                             dim=c(nclusters,ncategories,length(totake)),
                             dimnames=list(NULL,NULL,Vv$C)),
                       c(2,1,3))
    }
    if(len$B > 0){## binary
        Bprob <- array(grep('Bprob\\[', rownames(mcsamples)),
                       dim=c(len$B,nclusters), dimnames=list(variate$B,NULL))
    }
    if(len$I > 0){## integer ordinal
        Imean <- array(grep('Imean\\[', rownames(mcsamples)),
                       dim=c(len$I,nclusters), dimnames=list(variate$I,NULL))
        Ivar <- array(grep('Ivar\\[', rownames(mcsamples)),
                      dim=c(len$I,nclusters), dimnames=list(variate$I,NULL))
        Ilefts <- lapply(variate$I,function(v){ qnorm((0:(varinfo[['n']][v]-1L))/varinfo[['n']][v]) })
        Irights <- lapply(variate$I,function(v){ qnorm((1:(varinfo[['n']][v]))/varinfo[['n']][v]) })
        names(Ilefts) <- names(Irights) <- variate$I
    }
    if(len$R > 0){## real
        Rmean <- array(grep('Rmean\\[', rownames(mcsamples)),
                       dim=c(len$R,nclusters), dimnames=list(variate$R,NULL))
        Rvar <- array(grep('Rvar\\[', rownames(mcsamples)),
                      dim=c(len$R,nclusters), dimnames=list(variate$R,NULL))
    }
    if(len$D > 0){## two-bounded
        Dmean <- array(grep('Dmean\\[', rownames(mcsamples)),
                       dim=c(len$D,nclusters), dimnames=list(variate$D,NULL))
        Dvar <- array(grep('Dvar\\[', rownames(mcsamples)),
                      dim=c(len$D,nclusters), dimnames=list(variate$D,NULL))
        Dbounds <- -qnorm(varinfo[['n']][variate$D])
    }
    ##
    ## ##
    ## if(Yn$C > 0){## categorical
    ## totake <- sapply(Yv$C,function(x)which(variate$C == x))
    ## YCprob <- sapply(paste0('Cprob\\[',totake,','),
    ##                         grep,rownames(mcsamples))
    ## ncategories <- length(YCprob)/length(totake)/nclusters
    ## seqcategories <- seq_len(ncategories)
    ## YCprob <- aperm(array(YCprob,
    ##                        dim=c(nclusters,ncategories,length(totake)),
    ##                        dimnames=list(NULL,NULL,Yv$C)),
    ##                 c(2,1,3))
    ## }
    ## if(Yn$B > 0){## binary
    ## totake <- sapply(Yv$B,function(x)which(variate$B == x))
    ## YBprob <- array(t(vapply(paste0('Bprob\\[',totake,','), grep,
    ##                        numeric(nclusters),
    ##                        rownames(mcsamples))),
    ##                  dim=c(length(totake),nclusters), dimnames=list(Yv$B,NULL))
    ## }
    ## if(Yn$I > 0){## integer ordinal
    ## totake <- sapply(Yv$I,function(x)which(variate$I == x))
    ## YImean <- array(t(vapply(paste0('Imean\\[',totake,','), grep,
    ##                        numeric(nclusters),
    ##                        rownames(mcsamples))),
    ##                  dim=c(length(totake),nclusters), dimnames=list(Yv$I,NULL))
    ## YIvar <- array(t(vapply(paste0('Ivar\\[',totake,','), grep,
    ##                        numeric(nclusters),
    ##                        rownames(mcsamples))),
    ##                  dim=c(length(totake),nclusters), dimnames=list(Yv$I,NULL))
    ## YIlefts <- lapply(Yv$I,function(v){ qnorm((0:(varinfo[['n']][v]-1L))/varinfo[['n']][v]) })
    ## YIrights <- lapply(Yv$I,function(v){ qnorm((1:(varinfo[['n']][v]))/varinfo[['n']][v]) })
    ## names(YIlefts) <- names(YIrights) <- Yv$I
    ## }
    ## if(Yn$R > 0){## real
    ## totake <- sapply(Yv$R,function(x)which(variate$R == x))
    ## YRmean <- array(t(vapply(paste0('Rmean\\[',totake,','), grep,
    ##                        numeric(nclusters),
    ##                        rownames(mcsamples))),
    ##                  dim=c(length(totake),nclusters), dimnames=list(Yv$R,NULL))
    ## YRvar <- array(t(vapply(paste0('Rvar\\[',totake,','), grep,
    ##                        numeric(nclusters),
    ##                        rownames(mcsamples))),
    ##                  dim=c(length(totake),nclusters), dimnames=list(Yv$R,NULL))
    ## }
    ## if(Yn$O > 0){## one-side censored
    ## totake <- sapply(Yv$O,function(x)which(variate$O == x))
    ## YOmean <- array(t(vapply(paste0('Omean\\[',totake,','), grep,
    ##                        numeric(nclusters),
    ##                        rownames(mcsamples))),
    ##                  dim=c(length(totake),nclusters), dimnames=list(Yv$O,NULL))
    ## YOvar <- array(t(vapply(paste0('Ovar\\[',totake,','), grep,
    ##                        numeric(nclusters),
    ##                        rownames(mcsamples))),
    ##                  dim=c(length(totake),nclusters), dimnames=list(Yv$O,NULL))
    ## YOlefts <- transf(rbind(varinfo[['tmax']][Yv$O]),varinfo,Oout='')
    ## }
    ## if(Yn$D > 0){## two-bounded
    ## totake <- sapply(Yv$D,function(x)which(variate$D == x))
    ## YDmean <- array(t(vapply(paste0('Dmean\\[',totake,','), grep,
    ##                        numeric(nclusters),
    ##                        rownames(mcsamples))),
    ##                  dim=c(length(totake),nclusters), dimnames=list(Yv$D,NULL))
    ## YDvar <- array(t(vapply(paste0('Dvar\\[',totake,','), grep,
    ##                        numeric(nclusters),
    ##                        rownames(mcsamples))),
    ##                dim=c(length(totake),nclusters), dimnames=list(Yv$D,NULL))
    ## YDbounds <- -qnorm(varinfo[['n']][Yv$D])
    ## }
    ## ##
    ## ##
    ## if(Xn$C > 0){## categorical
    ## totake <- sapply(Xv$C,function(x)which(variate$C == x))
    ## XCprob <- sapply(paste0('Cprob\\[',totake,','),
    ##                         grep,rownames(mcsamples))
    ## ncategories <- length(XCprob)/length(totake)/nclusters
    ## seqcategories <- seq_len(ncategories)
    ## XCprob <- aperm(array(XCprob,
    ##                        dim=c(nclusters,ncategories,length(totake)),
    ##                        dimnames=list(NULL,NULL,Xv$C)),
    ##                 c(2,1,3))
    ## }
    ## if(Xn$B > 0){## binary
    ## totake <- sapply(Xv$B,function(x)which(variate$B == x))
    ## XBprob <- array(t(vapply(paste0('Bprob\\[',totake,','), grep,
    ##                        numeric(nclusters),
    ##                        rownames(mcsamples))),
    ##                  dim=c(length(totake),nclusters), dimnames=list(Xv$B,NULL))
    ## }
    ## if(Xn$I > 0){## integer ordinal
    ## totake <- sapply(Xv$I,function(x)which(variate$I == x))
    ## XImean <- array(t(vapply(paste0('Imean\\[',totake,','), grep,
    ##                        numeric(nclusters),
    ##                        rownames(mcsamples))),
    ##                  dim=c(length(totake),nclusters), dimnames=list(Xv$I,NULL))
    ## XIvar <- array(t(vapply(paste0('Ivar\\[',totake,','), grep,
    ##                        numeric(nclusters),
    ##                        rownames(mcsamples))),
    ##                  dim=c(length(totake),nclusters), dimnames=list(Xv$I,NULL))
    ## XIlefts <- lapply(Xv$I,function(v){ qnorm((0:(varinfo[['n']][v]-1L))/varinfo[['n']][v]) })
    ## XIrights <- lapply(Xv$I,function(v){ qnorm((1:(varinfo[['n']][v]))/varinfo[['n']][v]) })
    ## names(XIlefts) <- names(XIrights) <- Xv$I
    ## }
    ## if(Xn$R > 0){## real
    ## totake <- sapply(Xv$R,function(x)which(variate$R == x))
    ## XRmean <- array(t(vapply(paste0('Rmean\\[',totake,','), grep,
    ##                        numeric(nclusters),
    ##                        rownames(mcsamples))),
    ##                  dim=c(length(totake),nclusters), dimnames=list(Xv$R,NULL))
    ## XRvar <- array(t(vapply(paste0('Rvar\\[',totake,','), grep,
    ##                        numeric(nclusters),
    ##                        rownames(mcsamples))),
    ##                  dim=c(length(totake),nclusters), dimnames=list(Xv$R,NULL))
    ## }
    ## if(Xn$O > 0){## one-side censored
    ## totake <- sapply(Xv$O,function(x)which(variate$O == x))
    ## XOmean <- array(t(vapply(paste0('Omean\\[',totake,','), grep,
    ##                        numeric(nclusters),
    ##                        rownames(mcsamples))),
    ##                  dim=c(length(totake),nclusters), dimnames=list(Xv$O,NULL))
    ## XOvar <- array(t(vapply(paste0('Ovar\\[',totake,','), grep,
    ##                        numeric(nclusters),
    ##                        rownames(mcsamples))),
    ##                  dim=c(length(totake),nclusters), dimnames=list(Xv$O,NULL))
    ## XOlefts <- transf(rbind(varinfo[['tmax']][Xv$O]),varinfo,Oout='')
    ## }
    ## if(Xn$D > 0){## two-bounded
    ## totake <- sapply(Xv$D,function(x)which(variate$D == x))
    ## XDmean <- array(t(vapply(paste0('Dmean\\[',totake,','), grep,
    ##                        numeric(nclusters),
    ##                        rownames(mcsamples))),
    ##                  dim=c(length(totake),nclusters), dimnames=list(Xv$D,NULL))
    ## XDvar <- array(t(vapply(paste0('Dvar\\[',totake,','), grep,
    ##                        numeric(nclusters),
    ##                        rownames(mcsamples))),
    ##                dim=c(length(totake),nclusters), dimnames=list(Xv$D,NULL))
    ## XDbounds <- -qnorm(varinfo[['n']][Xv$D])
    ## }
    ##
    Y2 <- transf(Y,varinfo, Iout='index', Dout='index', Oout='index')
    if(!is.null(X)){
        X2 <- transf(X,varinfo,Iout='index', Dout='index', Oout='index')
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
        if(all(is.na(x))){
            probX <- t(log(W))
        }else{
            x <- x[!is.na(x),,drop=F]
            xv <- lapply(variatetypes, function(xx){out <- rownames(x)[varinfo[['type']][rownames(x)]==xx]
                names(out) <- NULL
                out})
            xn <- lapply(xv,length)
            names(xv) <- names(xn) <- variatetypes
            ##
            probX <- t( # rows: MCsamples, cols: clusters
                log(W) + 
                (if(Xn$R > 0){# continuous
                     colSums(
                         array(dnorm(x=x[Xt$R,],
                                     mean=Rmean[Xi$R,,],
                                     sd=sqrt(Rvar[Xi$R,,]),log=T),
                               dim=c(Xn$R, nclusters, nsamples)),
                         na.rm=F)
                 }else{0}) +
                (if(Xn$C > 0){# censored
                     colSums(
                         array(
                             t(sapply(Xt$C, function(v){
                                 if(is.finite(x[v,])){
                                     (dnorm(x=x[v,],
                                            mean=Cmean[v,,],
                                            sd=sqrt(Cvar[v,,]),log=T))
                                 }else{
                                     (pnorm(q=Cbounds[v,1.5+sign(x[v,])/2],
                                            mean=Cmean[v,,],
                                            sd=sqrt(Cvar[v,,]),
                                            lower.tail=(x[v,]<0),
                                            log.p=T))
                                 }
                             })),
                             dim=c(Xn$C, nclusters, length(subsamples))),
                         na.rm=F)
                 }else{0}) +
                (if(Xn$D > 0){# continuous
                     colSums(
                         array(dnorm(x=x[Xt$D,],
                                     mean=Dmean[Xi$D,,],
                                     sd=sqrt(Dvar[Xi$D,,]),log=T),
                               dim=c(Xn$D, nclusters, nsamples)),
                         na.rm=F)
                 }else{0}) +
                (if(Xn$O > 0){
                     vv <- cbind(Oseq,x[Xt$O,])
                     colSums(
                         array(log(
                             pnorm(q=Oright[vv],
                                   mean=Omean[Xi$O,,],
                                   sd=sqrt(Ovar[Xi$O,,])) -
                             pnorm(q=Oleft[vv],
                                   mean=Omean[Xi$O,,],
                                   sd=sqrt(Ovar[Xi$O,,]))
                         ),
                         dim=c(Xn$O, nclusters, nsamples)),
                         na.rm=F)
                 }else{0}) +
                (if(Xn$B > 0){
                     colSums(
                         array(log( x[Xt$B,]*Bprob[Xi$B,,] +
                                    (1-x[Xt$B,])*(1-Bprob[Xi$B,,]) ),
                               dim=c(Xn$B, nclusters, nsamples)),
                         na.rm=F)
                 }else{0})
            ) # end probX
        }
        ##
        y <- y[!is.na(y),,drop=F]
        yv <- lapply(variatetypes, function(xx){
            out <- rownames(y)[varinfo[['type']][rownames(y)]==xx]
            names(out) <- NULL
            out})
        yn <- lapply(yv,length)
        names(yv) <- names(yn) <- variatetypes
        ##
        if(all(is.na(y))){
            probY <- NA
        }else{
            probY <- t( # rows: MCsamples, cols: clusters
            (if(yn$D > 0){
                 colSums(
                     array(
                         t(sapply(yv$D, function(v){
                             if(is.finite(y[v,])){
                                 (dnorm(x=y[v,],
                                        mean=mcsamples[Dmean[v,],],
                                        sd=sqrt(mcsamples[Dvar[v,],]),log=T))
                             }else{
                                 (pnorm(q=Dbounds[v,]*sign(y[v,]),
                                        mean=mcsamples[Dmean[v,],],
                                        sd=sqrt(mcsamples[Dvar[v,],]),
                                        lower.tail=(y[v,]<0),
                                        log.p=T))
                             }
                         })),
                         dim=c(yn$D, nclusters, length(subsamples))),
                     na.rm=F)
             }else{0}) +
            ## (if(yn$O > 0){
            ##      colSums(
            ##          array(
            ##              t(sapply(yv$O, function(v){
            ##                  if(is.finite(y[v,])){
            ##                      (dnorm(x=y[v,],
            ##                              mean=mcsamples[YOmean[v,],],
            ##                              sd=sqrt(mcsamples[YOvar[v,],]),log=T))
            ##                  }else{
            ##                      (pnorm(q=YOlefts[1,v],
            ##                              mean=mcsamples[YOmean[v,],],
            ##                              sd=sqrt(mcsamples[YOvar[v,],]),
            ##                              lower.tail=F,
            ##                              log.p=T))
            ##                  }
            ##              })),
            ##              dim=c(yn$O, nclusters, length(subsamples))),
            ##          na.rm=F)
            ##  }else{0}) +
            (if(yn$R > 0){
                 colSums(
                     array(dnorm(x=y[yv$R,],
                                 mean=mcsamples[Rmean[yv$R,],],
                                 sd=sqrt(mcsamples[Rvar[yv$R,],]),log=T),
                           dim=c(yn$R, nclusters, length(subsamples))),
                     na.rm=F)
             }else{0}) +
            (if(yn$B > 0){
                 colSums(
                     array(log( y[yv$B,]*mcsamples[Bprob[yv$B,],] +
                                (1-y[yv$B,])*(1-mcsamples[Bprob[yv$B,],]) ),
                           dim=c(yn$B, nclusters, length(subsamples))),
                     na.rm=F)
             }else{0}) +
            (if(yn$I > 0){
                 colSums(
                     array(
                         t(sapply(yv$I, function(v){
                             log(pnorm(q=Irights[[v]][y[v,]],
                                       mean=mcsamples[Imean[v,],],
                                       sd=sqrt(mcsamples[Ivar[v,],])) -
                                 pnorm(q=Ilefts[[v]][y[v,]],
                                       mean=mcsamples[Imean[v,],],
                                       sd=sqrt(mcsamples[Ivar[v,],])))
                             ## y2 <- YIrights[[v]][y[v,]]
                             ## (pnorm(q=y2,
                             ##         mean=mcsamples[YImean[v,],],
                             ##         sd=sqrt(mcsamples[YIvar[v,],]), log.p=T) +
                             ##   log1p(-pnorm(q=YIlefts[[v]][y[v,]],
                             ##                mean=mcsamples[YImean[v,],],
                             ##                sd=sqrt(mcsamples[YIvar[v,],]))/
                             ##         pnorm(q=y2,
                             ##               mean=mcsamples[YImean[v,],],
                             ##               sd=sqrt(mcsamples[YIvar[v,],]))))
                         })),
                         dim=c(yn$I, nclusters, length(subsamples))),
                     na.rm=F)
             }else{0})
            ) # end probY
        }
        ##
        ## ## Other approaches tested for roundoff error
        ## testc2 <- rowSums(exp(probX2 + probY - log(rowSums(exp(probX2)))))
        ##
        ## maxp <- apply(probX2 + probY - log(rowSums(exp(probX2))),1,max,na.rm=T)
        ## testc2f <- exp(maxp)*rowSums(exp(probX2 + probY - log(rowSums(exp(probX2)))-maxp))
        ##
        ## testc1d <- rowSums(exp(probX+probY))/rowSums(exp(probX))
        ##        
        ## testc1 <- rowSums(exp(probX + probY - log(rowSums(exp(probX)))))
        ##
        ## maxp <- apply(probX + probY - log(rowSums(exp(probX))),1,max,na.rm=T)
        ## testc1f <- exp(maxp)*rowSums(exp(probX + probY - log(rowSums(exp(probX)))-maxp))
        ##
        ## ## Comparison using Rmpfr library
        ## dprobX <- mpfr(probX, precBits=200)
        ## dprobY <- mpfr(probY, precBits=200)
        ##
        ## dtestc1 <- rowSums(exp(dprobX + dprobY - log(rowSums(exp(dprobX)))))
        ## dtestc1d <- rowSums(exp(dprobX+dprobY))/rowSums(exp(dprobX))
        ##
        ##
        ## This seems to minimize roundoff error
        if(all(is.na(x))){
            out <- rowSums(exp(probX+probY)) 
        }else{
            probX <- probX - apply(probX, 1, max, na.rm=T)
            out <- rowSums(exp(probX+probY))/rowSums(exp(probX))
        }
        fn(out)
    } *
        (if(jacobian){exp(-rowSums(
                               log(invjacobian(Y,varinfo=varinfo))
                             , na.rm=T))}else{1L})
}

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
    vn <- vnames <- Yt <- Xt <- Yi <- Xi <- Yn <- Xn <- Yseq <- Xseq <- vseq <- list()
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
        Yseq[[atype]] <- 1:Yn[[atype]]
        Xseq[[atype]] <- 1:Xn[[atype]]
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
                             t(sapply(Xseq$C, function(v){
                                 v1 <- Xt$C[v]
                                 v2 <- Xi$C[v]
                                 if(is.finite(x[v1,])){
                                     (dnorm(x=x[v1,],
                                            mean=Cmean[v2,,],
                                            sd=sqrt(Cvar[v2,,]),log=T))
                                 }else{
                                     (pnorm(q=Cbounds[v2,1.5+sign(x[v1,])/2],
                                            mean=Cmean[v2,,],
                                            sd=sqrt(Cvar[v2,,]),
                                            lower.tail=(x[v1,]<0),
                                            log.p=T))
                                 }
                             })),
                             dim=c(Xn$C, nclusters, nsamples)),
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
                     v2 <- cbind(Oseq,x[Xt$O,])
                     colSums(
                         array(log(
                             pnorm(q=Oright[v2],
                                   mean=Omean[Xi$O,,],
                                   sd=sqrt(Ovar[Xi$O,,])) -
                             pnorm(q=Oleft[v2],
                                   mean=Omean[Xi$O,,],
                                   sd=sqrt(Ovar[Xi$O,,]))
                         ),
                         dim=c(Xn$O, nclusters, nsamples)),
                         na.rm=F)
                 }else{0}) +
                (if(Xn$N > 0){
                     colSums(
                         array(
                             t(sapply(Xseq$N, function(v){
                                 Nprob[Xi$N[v],,Xt$N[v],]
                                 })),
                             dim=c(Xn$N, nclusters, nsamples)),
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

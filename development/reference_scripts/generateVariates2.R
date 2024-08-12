generateVariates <- function(n, Ynames, X=NULL, mcsamples, varinfo){
    if(!is.null(X) && length(dim(X)) != 2){stop('X must be NULL or have two dimensions')}
    ##
    if(!is.null(X) && ncol(X) == 0){X <- NULL}
    subsamples <- sample(1:nrow(mcsamples), n, replace=(n > nrow(mcsamples)))
    seqn <- 1:n
    mcsamples <- t(mcsamples[subsamples,,drop=FALSE])
    Yv <- Ynames
    Yvn <- length(Yv)
    Xv <- colnames(X)
    Xvn <- length(Xv)
    Vv <- varinfo[['name']]
    if(length(intersect(Yv, Xv)) > 0){warning('overlap in Y and X variates')}
    if(!all(Yv %in% Vv)){warning('unknown Y variates')}
    if(!all(Xv %in% Vv)){warning('unknown X variates')}
    ##
    variate <- lapply(variatetypes, function(x)names(varinfo[['type']])[varinfo[['type']]==x])
    len <- lapply(variate,length)
    names(variate) <- names(len) <- variatetypes
    ##
    Wi <- grep('W',rownames(mcsamples))
    ncomponents <- length(Wi)
    W <- mcsamples[Wi,] # rows: components, cols: MC samples
    ## seqcomponents <- seq_len(ncomponents)
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
    if(Yn$C > 0){## categorical
    totake <- sapply(Yv$C,function(x)which(variate$C == x))
    YCprob <- sapply(paste0('Cprob\\[',totake,','),
                            grep,rownames(mcsamples))
    ncategories <- length(YCprob)/length(totake)/ncomponents
    seqcategories <- seq_len(ncategories)
    YCprob <- aperm(array(YCprob,
                           dim=c(ncomponents,ncategories,length(totake)),
                           dimnames=list(NULL,NULL,Yv$C)),
                    c(2,1,3))
    }
    if(Yn$B > 0){## binary
    totake <- sapply(Yv$B,function(x)which(variate$B == x))
    YBprob <- array(t(vapply(paste0('Bprob\\[',totake,','), grep,
                           numeric(ncomponents),
                           rownames(mcsamples))),
                     dim=c(length(totake),ncomponents), dimnames=list(Yv$B,NULL))
    }
    if(Yn$I > 0){## integer ordinal
    totake <- sapply(Yv$I,function(x)which(variate$I == x))
    YImean <- array(t(vapply(paste0('Imean\\[',totake,','), grep,
                           numeric(ncomponents),
                           rownames(mcsamples))),
                     dim=c(length(totake),ncomponents), dimnames=list(Yv$I,NULL))
    YIvar <- array(t(vapply(paste0('Ivar\\[',totake,','), grep,
                           numeric(ncomponents),
                           rownames(mcsamples))),
                     dim=c(length(totake),ncomponents), dimnames=list(Yv$I,NULL))
    ## YIlefts <- lapply(Yv$I,function(v){ qnorm((0:(varinfo[['n']][v]-1L))/varinfo[['n']][v]) })
    ## YIrights <- lapply(Yv$I,function(v){ qnorm((1:(varinfo[['n']][v]))/varinfo[['n']][v]) })
    ## names(YIlefts) <- names(YIrights) <- Yv$I
    }
    if(Yn$R > 0){## real
    totake <- sapply(Yv$R,function(x)which(variate$R == x))
    YRmean <- array(t(vapply(paste0('Rmean\\[',totake,','), grep,
                           numeric(ncomponents),
                           rownames(mcsamples))),
                     dim=c(length(totake),ncomponents), dimnames=list(Yv$R,NULL))
    YRvar <- array(t(vapply(paste0('Rvar\\[',totake,','), grep,
                           numeric(ncomponents),
                           rownames(mcsamples))),
                     dim=c(length(totake),ncomponents), dimnames=list(Yv$R,NULL))
    }
    if(Yn$O > 0){## one-side censored
    totake <- sapply(Yv$O,function(x)which(variate$O == x))
    YOmean <- array(t(vapply(paste0('Omean\\[',totake,','), grep,
                           numeric(ncomponents),
                           rownames(mcsamples))),
                     dim=c(length(totake),ncomponents), dimnames=list(Yv$O,NULL))
    YOvar <- array(t(vapply(paste0('Ovar\\[',totake,','), grep,
                           numeric(ncomponents),
                           rownames(mcsamples))),
                     dim=c(length(totake),ncomponents), dimnames=list(Yv$O,NULL))
    ## YOlefts <- c(sapply(Yv$O, function(v){
    ##     out <- varinfo[['t']][[v]](varinfo[['max']][v])
    ##     names(out) <- NULL
    ##     out}))
    }
    if(Yn$D > 0){## two-bounded
    totake <- sapply(Yv$D,function(x)which(variate$D == x))
    YDmean <- array(t(vapply(paste0('Dmean\\[',totake,','), grep,
                           numeric(ncomponents),
                           rownames(mcsamples))),
                     dim=c(length(totake),ncomponents), dimnames=list(Yv$D,NULL))
    YDvar <- array(t(vapply(paste0('Dvar\\[',totake,','), grep,
                           numeric(ncomponents),
                           rownames(mcsamples))),
                   dim=c(length(totake),ncomponents), dimnames=list(Yv$D,NULL))
    }
    ##
    ##
    if(Xn$C > 0){## categorical
    totake <- sapply(Xv$C,function(x)which(variate$C == x))
    XCprob <- sapply(paste0('Cprob\\[',totake,','),
                            grep,rownames(mcsamples))
    ncategories <- length(XCprob)/length(totake)/ncomponents
    seqcategories <- seq_len(ncategories)
    XCprob <- aperm(array(XCprob,
                           dim=c(ncomponents,ncategories,length(totake)),
                           dimnames=list(NULL,NULL,Xv$C)),
                    c(2,1,3))
    }
    if(Xn$B > 0){## binary
    totake <- sapply(Xv$B,function(x)which(variate$B == x))
    XBprob <- array(t(vapply(paste0('Bprob\\[',totake,','), grep,
                           numeric(ncomponents),
                           rownames(mcsamples))),
                     dim=c(length(totake),ncomponents), dimnames=list(Xv$B,NULL))
    }
    if(Xn$I > 0){## integer ordinal
    totake <- sapply(Xv$I,function(x)which(variate$I == x))
    XImean <- array(t(vapply(paste0('Imean\\[',totake,','), grep,
                           numeric(ncomponents),
                           rownames(mcsamples))),
                     dim=c(length(totake),ncomponents), dimnames=list(Xv$I,NULL))
    XIvar <- array(t(vapply(paste0('Ivar\\[',totake,','), grep,
                           numeric(ncomponents),
                           rownames(mcsamples))),
                     dim=c(length(totake),ncomponents), dimnames=list(Xv$I,NULL))
    XIlefts <- lapply(Xv$I,function(v){ qnorm((0:(varinfo[['n']][v]-1L))/varinfo[['n']][v]) })
    XIrights <- lapply(Xv$I,function(v){ qnorm((1:(varinfo[['n']][v]))/varinfo[['n']][v]) })
    names(XIlefts) <- names(XIrights) <- Xv$I
    }
    if(Xn$R > 0){## real
    totake <- sapply(Xv$R,function(x)which(variate$R == x))
    XRmean <- array(t(vapply(paste0('Rmean\\[',totake,','), grep,
                           numeric(ncomponents),
                           rownames(mcsamples))),
                     dim=c(length(totake),ncomponents), dimnames=list(Xv$R,NULL))
    XRvar <- array(t(vapply(paste0('Rvar\\[',totake,','), grep,
                           numeric(ncomponents),
                           rownames(mcsamples))),
                     dim=c(length(totake),ncomponents), dimnames=list(Xv$R,NULL))
    }
    if(Xn$O > 0){## one-side censored
    totake <- sapply(Xv$O,function(x)which(variate$O == x))
    XOmean <- array(t(vapply(paste0('Omean\\[',totake,','), grep,
                           numeric(ncomponents),
                           rownames(mcsamples))),
                     dim=c(length(totake),ncomponents), dimnames=list(Xv$O,NULL))
    XOvar <- array(t(vapply(paste0('Ovar\\[',totake,','), grep,
                           numeric(ncomponents),
                           rownames(mcsamples))),
                     dim=c(length(totake),ncomponents), dimnames=list(Xv$O,NULL))
    XOlefts <- transf(rbind(varinfo[['tmax']][Xv$O]),varinfo,Oout='')
    }
    if(Xn$D > 0){## two-bounded
    totake <- sapply(Xv$D,function(x)which(variate$D == x))
    XDmean <- array(t(vapply(paste0('Dmean\\[',totake,','), grep,
                           numeric(ncomponents),
                           rownames(mcsamples))),
                     dim=c(length(totake),ncomponents), dimnames=list(Xv$D,NULL))
    XDvar <- array(t(vapply(paste0('Dvar\\[',totake,','), grep,
                           numeric(ncomponents),
                           rownames(mcsamples))),
                   dim=c(length(totake),ncomponents), dimnames=list(Xv$D,NULL))
    XDbounds <- -qnorm(varinfo[['n']][Xv$D])
    }
    ##
    ##    Y2 <- transf(Y,varinfo, Iout='index', Dout='index', Oout='index')
    if(!is.null(X)){
        X2 <- transf(X,varinfo,Iout='index', Dout='index', Oout='index')
        ## if(nrow(X2) < n){
        ##     warning('*Note: X has fewer entries than n. Recycling*')
        ##     X2 <- t(matrix(rep(t(X2), ceiling(n/nrow(X2))), nrow=ncol(X2), dimnames=list(colnames(X2),NULL)))[1:n,,drop=FALSE]
        ## }
        ## if(nrow(X2) > n){
        ##     warning('*Note: X has more entries than n. Skipping some*')
        ##     Y2 <- t(matrix(rep(t(Y2), ceiling(nrow(X2)/nrow(Y2))), nrow=ncol(Y2), dimnames=list(colnames(Y2),NULL)))[1:nrow(X2),,drop=FALSE]
        ## }
    }else{X2 <- t(NA)}
    ## ndata <- nrow(Y2)
    ##
    ##
    out <- foreach(x=t(X2), .combine=rbind, .inorder=T)%dorng%{
        ## ## for debugging
        ## for(iii in 1:nrow(Y2)){
        ## print(iii)
        ## iii <- iii+1
        ##     y <- t(Y2)[,1,drop=F]
        ##     x <- t(X2)[,1,drop=F]
            ##
        ## if(any(is.na(y))){
        ##     y <- y[!is.na(y),,drop=F]
        ##     yv <- lapply(variatetypes, function(xx){
        ##         out <- rownames(y)[varinfo[['type']][rownames(y)]==xx]
        ##         names(out) <- NULL
        ##         out})
        ##     yn <- lapply(yv,length)
        ##     names(yv) <- names(yn) <- variatetypes
        ## }else{
            yv <- Yv
            yn <- Yn
        ## }
        ##            
        if(all(is.na(x))){
                probX <- t(W)
        }else{
            if(any(is.na(x))){
                x <- x[!is.na(x),,drop=F]
                xv <- lapply(variatetypes, function(xx){out <- rownames(x)[varinfo[['type']][rownames(x)]==xx]
                    names(out) <- NULL
                    out})
                xn <- lapply(xv,length)
                names(xv) <- names(xn) <- variatetypes
            }else{
                xv <- Xv
                xn <- Xn
                }
                ##
        probX <- t(exp( # rows: MCsamples, cols: components
            log(W) + 
                (if(xn$D > 0){
                     colSums(
                         array(
                             t(sapply(xv$D, function(v){
                                 if(is.finite(x[v,])){
                                     (dnorm(x=x[v,],
                                             mean=mcsamples[XDmean[v,],],
                                             sd=sqrt(mcsamples[XDvar[v,],]),log=T))
                                 }else{
                                     (pnorm(q=XDbounds[v,]*sign(x[v,]),
                                             mean=mcsamples[XDmean[v,],],
                                             sd=sqrt(mcsamples[XDvar[v,],]),
                                             lower.tail=(x[v,]<0),
                                             log.p=T))
                                 }
                             })),
                             dim=c(xn$D, ncomponents, length(subsamples))),
                         na.rm=F)
                 }else{0}) +
                (if(xn$O > 0){
                     colSums(
                         array(
                             t(sapply(xv$O, function(v){
                                 if(is.finite(x[v,])){
                                     (dnorm(x=x[v,],
                                             mean=mcsamples[XOmean[v,],],
                                             sd=sqrt(mcsamples[XOvar[v,],]),log=T))
                                 }else{
                                     (pnorm(q=XOlefts[1,v],
                                             mean=mcsamples[XOmean[v,],],
                                             sd=sqrt(mcsamples[XOvar[v,],]),
                                             lower.tail=F,
                                             log.p=T))
                                 }
                             })),
                             dim=c(xn$O, ncomponents, length(subsamples))),
                         na.rm=F)
                 }else{0}) +
                (if(xn$R > 0){
                     colSums(
                         array(dnorm(x=x[xv$R,],
                                     mean=mcsamples[XRmean,],
                                     sd=sqrt(mcsamples[XRvar,]),log=T),
                               dim=c(xn$R, ncomponents, length(subsamples))),
                         na.rm=F)
                 }else{0}) +
                (if(xn$B > 0){
                     colSums(
                         array(log( x[xv$B,]*mcsamples[XBprob,] +
                               (1-x[xv$B,])*(1-mcsamples[XBprob,]) ),
                               dim=c(xn$B, ncomponents, length(subsamples))),
                         na.rm=F)
                 }else{0}) +
                (if(xn$I > 0){
                     colSums(
                         array(
                             t(sapply(xv$I, function(v){
                                 log(pnorm(q=XIrights[[v]][x[v,]],
                                        mean=mcsamples[XImean[v,],],
                                        sd=sqrt(mcsamples[XIvar[v,],])) -
                                  pnorm(q=XIlefts[[v]][x[v,]],
                                        mean=mcsamples[XImean[v,],],
                                        sd=sqrt(mcsamples[XIvar[v,],])))
                                 ## x2 <- XIrights[[v]][x[v,]]
                                 ## (pnorm(q=x2,
                                 ##         mean=mcsamples[XImean[v,],],
                                 ##         sd=sqrt(mcsamples[XIvar[v,],]), log.p=T) +
                                 ##   log1p(-pnorm(q=XIlefts[[v]][x[v,]],
                                 ##                mean=mcsamples[XImean[v,],],
                                 ##                sd=sqrt(mcsamples[XIvar[v,],]))/
                                 ##         pnorm(q=x2,
                                 ##               mean=mcsamples[XImean[v,],],
                                 ##               sd=sqrt(mcsamples[XIvar[v,],]))))
                             })),
                             dim=c(xn$I, ncomponents, length(subsamples))),
                         na.rm=F)
                 }else{0})
        )) # end probX
        }
    ##
    Ws <- extraDistr::rcat(n=n, prob=probX)
        probY <- cbind( # rows: samples, cols: variates
        (if(Yn$D > 0){
             matrix(rnorm(n=n*Yn$D,
                          mean=mcsamples[cbind(c(t(YDmean[,Ws])),seqn)],
                          sd=sqrt(mcsamples[cbind(c(t(YDvar[,Ws])),seqn)])),
                    nrow=n, ncol=Yn$D, dimnames=list(NULL,Yv$D))
                 }else{NULL}),
        (if(Yn$O > 0){
             matrix(rnorm(n=n*Yn$O,
                          mean=mcsamples[cbind(c(t(YOmean[,Ws])),seqn)],
                          sd=sqrt(mcsamples[cbind(c(t(YOvar[,Ws])),seqn)])),
                    nrow=n, ncol=Yn$O, dimnames=list(NULL,Yv$O))
                 }else{NULL}),
        (if(Yn$R > 0){
             matrix(rnorm(n=n*Yn$R,
                          mean=mcsamples[cbind(c(t(YRmean[,Ws])),seqn)],
                          sd=sqrt(mcsamples[cbind(c(t(YRvar[,Ws])),seqn)])),
                    nrow=n, ncol=Yn$R, dimnames=list(NULL,Yv$R))
                 }else{NULL}),
        (if(Yn$B > 0){
             matrix(extraDistr::rbern(n=n*Yn$B,
                                      prob=mcsamples[cbind(c(t(YBprob[,Ws])),seqn)]),
                    nrow=n, ncol=Yn$B, dimnames=list(NULL,Yv$B))
         }else{NULL}),
        (if(Yn$I > 0){
             matrix(rnorm(n=n*Yn$I,
                          mean=mcsamples[cbind(c(t(YImean[,Ws])),seqn)],
                          sd=sqrt(mcsamples[cbind(c(t(YIvar[,Ws])),seqn)])),
                    nrow=n, ncol=Yn$I, dimnames=list(NULL,Yv$I))
                 }else{NULL})
        ) # end probY
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
    probY
    }
    attr(out,'rng') <- NULL
    attr(out,'doRNG_version') <- NULL
    out <- t(invtransf(z=out[,Ynames,drop=F], varinfo, Oout='censored'))
    dim(out) <- c(length(Ynames),n,nrow(X2))
    dimnames(out) <- list(Ynames,NULL,rownames(X2))
    ##
    out
}

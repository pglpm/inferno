## Author: PGL  Porta Mana
## Created: 2022-09-08T17:03:24+0200
## Last-Updated: 2023-02-21T14:29:05+0100
#########################################
## Inference of exchangeable variates (nonparametric density regression)
## using effectively-infinite mixture of product kernels
## Auxiliary functions
#########################################

sd2iqr <- 0.5/qnorm(0.75)

variatetypes <- c( 'categorical'='C', 'binary'='B', 'ordinal'='I', 'continuous'='R', 'one_censored'='O', 'two_censored'='D' )

## Create interval bounds for transformed integer variates; pad with +Inf
## createbounds <- function(n, nmax=n){
##     if(n>1){ c( qnorm((1:(n-1))/n), rep(+Inf, nmax-n) ) } else {qnorm(c(n, 1-n))}
## }


## m <- 0
## M <- 100
## n <- 7
## sca <- (M-m)/(n-1)
## loc <- m
## tests <- cbind('test'=c(NA,seq(m,M,length.out=n)))
## testvarinfo <- rbind('test'=c('min'=m, 'max'=M, 'n'=n, 'location'=loc, 'scale'=sca, 'type'=1))


## m <- 0
## M <- 100
## n <- 2^-6
## sca <- (M-m)/(1-2*n)
## loc <- m-n*sca
## tests <- cbind('test'=c(NA,seq(m,M,length.out=20)))
## testvarinfo <- rbind('test'=c('min'=m, 'max'=M, 'n'=n, 'location'=loc, 'scale'=sca, 'type'=-2))
## cbind(tests,transf(tests,testvarinfo,Tout='left'),transf(tests,testvarinfo,Tout='init'),transf(tests,testvarinfo,Tout='data'),transf(tests,testvarinfo,Tout='right'))

#### Transformation types:
## 'identity':
## direct: (x-mu)/sigma
## inverse: y*sigma+mu
## invjacobian: sigma
##
## 'integer':
## direct: round((x-mu)/sigma)
## inverse: y*sigma+mu
## invjacobian: rep(1L, length(x))
##
## 'binary':
## direct: round((x-mu)/sigma)
## inverse: y*sigma+mu
## invjacobian: rep(1L, length(x))
##
## 'log':
## direct: (log(x-min)-mu)/sigma
## inverse: exp(y*sigma+mu)
## invjacobian: x*sigma
##
## 'probit':
## direct: (qnorm((x-min)/(max-min))-mu)/sigma
## inverse: exp(y*sigma+mu)
## invjacobian: x*sigma


## Transformation from variate to internal variable
transf <- function(x, varinfo, Iout='init', Oout='data', Dout='data', Bout='data', variates=NULL, Ifunction='qnorm'){ # 'in' 'data' 'aux'
    x <- cbind(data.matrix(x))
    if(!is.null(variates)){colnames(x) <- variates}
    matrix(sapply(colnames(x), function(v){
        ##
        datum <- x[,v,drop=F]
        info <- lapply(varinfo,function(xxx)xxx[v])
        ##
        if(info$type == 'R'){
            if (info$transform == 'log'){
                datum <- log(datum-info$min)
            }else if (info$transform == 'probit'){
                datum <- qnorm((datum-info$min)/(info$max-info$min))
            }
            datum <- (datum-info$location)/info$scale
            ##
        } else if(info$type == 'I'){ # integer, discrete ordinal
            if(Ifunction == 'qnorm'){Ifunction <- qnorm
            } else if(Ifunction == 'qt'){Ifunction <- function(xxx){2.25*qt(xxx,df=2.25)}}
            datum <- round((datum-info$location)/info$scale) # output is in range 0 to n-1
            if(Iout == 'init'){ # in sampling functions or init MCMC
                datum[is.na(datum)] <- round(info$n/2-0.5)
                datum <- Ifunction((datum+0.5)/info$n)
            } else if(Iout == 'left'){ # as left for MCMC
                datum <- Ifunction(pmax(0,datum)/info$n)
                datum[is.na(datum)] <- -Inf
            } else if(Iout == 'right'){ # as right for MCMC
                datum <- Ifunction(pmin(info$n,datum+1)/info$n)
                datum[is.na(datum)] <- +Inf
            } else if(Iout == 'index'){ # in output functions
                datum <- datum+1L
            } else if(Iout == 'correct'){ # in output functions
                datum <- datum * info$scale + info$location
            }
            ##
        } else if(info$type == 'B'){ # binary
            datum <- round((datum-info$location)/info$scale) # output in range 0 to 1
            if(Bout == 'correct'){ # in output functions
                datum <- datum * info$scale + info$location
            }
            ##
        } else if(info$type == 'D'){ # censored
            if (info$transform == 'log'){
                datum <- log(datum-info$min)
                rightbound <- log(info$tmax-info$min)
                leftbound <- log(info$tmin-info$min)
            }else if (info$transform == 'probit'){
                datum <- qnorm((datum-info$min)/(info$max-info$min))
                rightbound <- qnorm((info$tmax-info$min)/(info$max-info$min))
                leftbound <- qnorm((info$tmin-info$min)/(info$max-info$min))
            }
            datum <- (datum-info$location)/info$scale
            rightbound <- (rightbound-info$location)/info$scale
            leftbound <- (leftbound-info$location)/info$scale
            if(Dout == 'left'){ # in MCMC
                sel <- is.na(datum) | (x[,v] < info$tmax)
                datum[sel] <- -Inf
                datum[!sel] <- rightbound
            } else if(Dout == 'right'){ # in MCMC
                sel <- is.na(datum) | (x[,v] > info$tmin)
                datum[sel] <- +Inf
                datum[!sel] <- leftbound
            } else if(Dout == 'data'){ # data in MCMC
                sel <- is.na(datum) | (x[,v] >= info$tmax) | (x[,v] <= info$tmin)
                datum[sel] <- NA
            } else if(Dout == 'init'){ #init in MCMC
                datum[is.na(datum)] <- 0L
                datum[x[,v] >= info$tmax] <- rightbound + 0.125
                datum[x[,v] <= info$tmin] <- leftbound - 0.125
            } else if(Dout == 'index'){ #in sampling functions
                datum[x[,v] >= info$tmax] <- +Inf
                datum[x[,v] <= info$tmin] <- -Inf
            }
            ##
        ## } else if(info$type == 'O'){ # one-censored
        ##     if (info$transform == 'log'){
        ##         datum <- log(datum-info$min)
        ##         rightbound <- log(info$tmax-info$min)
        ##     }else if (info$transform == 'probit'){
        ##         datum <- qnorm((datum-info$min)/(info$max-info$min))
        ##         rightbound <- qnorm((info$tmax-info$min)/(info$max-info$min))
        ##     }
        ##     datum <- (datum-info$location)/info$scale
        ##     rightbound <- (rightbound-info$location)/info$scale
        ##     if(Oout == 'left'){ # in MCMC
        ##         sel <- is.na(datum) | (x[,v] < info$tmax)
        ##         datum[sel] <- -Inf
        ##         datum[!sel] <- rightbound
        ##     } else if(Oout == 'data'){ # data in MCMC
        ##         sel <- is.na(datum) | (x[,v] >= info$tmax)
        ##         datum[sel] <- NA
        ##     } else if(Oout == 'init'){ #init in MCMC
        ##         datum[is.na(datum)] <- 0L
        ##         datum[x[,v] >= info$tmax] <- rightbound+0.125
        ##     } else if(Oout == 'index'){ #in sampling functions
        ##         datum[x[,v] >= info$tmax] <- +Inf
        ##     }
        ##     ##
        }
        datum
    }),ncol=ncol(x),dimnames=list(NULL,colnames(x)))
}

## Inverse of Jacobian, in terms of original variate
invjacobian <- function(x, varinfo, variates=NULL){ # 'in' 'data' 'aux'
    x <- cbind(data.matrix(x))
    if(!is.null(variates)){colnames(x) <- variates}
    matrix(sapply(colnames(x), function(v){
        ##
        datum <- x[,v,drop=F]
        info <- lapply(varinfo,function(xxx)xxx[v])
        ##
        if(info$type == 'R'){
            if(info$transform == 'identity'){
                datum <- 1L + 0L *datum
            }else if (info$transform == 'probit'){
                datum <- dnorm(qnorm((datum-info$min)/(info$max-info$min))) * (info$max-info$min)
            }
            datum <- datum * info$scale
            ##
        }else if(info$type == 'O'){ # one-bounded
            if(info$transform == 'identity'){
                datum <- rep(1L, length(datum))
            }else if (info$transform == 'probit'){
                datum <- dnorm(qnorm((datum-info$min)/(info$max-info$min))) * (info$max-info$min)
            }
            datum <- datum * info$scale
            datum[x[,v] >= info$tmax] <- 1L
            ##
        } else {
            datum <- 1L + 0L *datum
        }
        datum
    }),ncol=ncol(x),dimnames=list(NULL,colnames(x)))
}

## Inverse transformation from internal variable to variate
invtransf <- function(z, varinfo, Oout='censored', variates=NULL){ # 
    z <- cbind(data.matrix(z))
    if(!is.null(variates)){colnames(z) <- variates}
    matrix(sapply(colnames(z), function(v){
        ##
        datum <- z[,v,drop=F]
        info <- lapply(varinfo,function(xxx)xxx[v])
        ##
        if(info$type == 'R'){
            datum <- datum*info$scale + info$location
            if (info$transform == 'log'){
                datum <- exp(datum) + info$min
            } else if (info$transform == 'probit'){
                datum <- pnorm(datum) * (info$max-info$min) + info$min
            }
            ##
        } else if(info$type == 'I'){ # integer, discrete ordinal
            datum <- nimble::rinterval(n=length(datum), datum,
                                       qnorm((1:(info$n-1))/info$n) )
            datum <- datum * info$scale + info$location
            ##
        } else if(info$type == 'B'){ # binary
            datum <- datum * info$scale + info$location
            ##
        } else if(info$type == 'O'){ # one-censored
            datum <- datum*info$scale + info$location
            if (info$transform == 'log'){
                datum <- exp(datum) + info$min
            } else if (info$transform == 'probit'){
                datum <- pnorm(datum) * (info$max-info$min) + info$min
            }
            if(Oout == 'censored'){
                datum[datum > info$tmax] <- info$tmax
            }
        }
        datum
    }),ncol=ncol(z),dimnames=list(NULL,colnames(z)))
}


## Extend range
extendrange <- function(vec,by=0.125){
    (1+by) * range(vec, na.rm=TRUE)- by * mean(range(vec, na.rm=TRUE))
}
## Normalize vector
normalize <- function(freqs){freqs/sum(freqs)}
## Normalize rows of matrix
normalizerows <- function(freqs){freqs/rowSums(freqs)}
##
## Assess stationarity (from LaplacesDemon::LaplacesDemon)
proposeburnin <- function(x, batches=16){
    if(is.null(dim(x))){x <- cbind(x) }
    lx <- nrow(x)
    if(lx%%batches != 0){
        x <- x[1:(batches * trunc(lx/batches)), ]
    }
    lx2 <- nrow(x)
    HD <- LaplacesDemon::BMK.Diagnostic(x, batches = batches)
    Ind <- 1 * (HD > 0.5)
    BurnIn <- lx
    batch.list <- seq(from = 1, to = lx2, by = floor(lx2/batches))
    for (i in 1:(batches-1)) {
        if (sum(Ind[, i:(batches-1)]) == 0) {
            BurnIn <- batch.list[i] - 1
            break
        }
    }
BurnIn
}
##
## Assess thinning (from LaplacesDemon::LaplacesDemon)
proposethinning <- function(x){
    if(is.null(dim(x))){ x <- cbind(x) }
    lx <- nrow(x)
    LIV <- ncol(x)
    acf.rows <- trunc(10 * log10(lx))
    acf.temp <- matrix(1, acf.rows, LIV)
    Rec.Thin <- rep(1, LIV)
    names(Rec.Thin) <- colnames(x)
    for (j in 1:LIV) {
        temp0 <- acf(x[, j], lag.max = acf.rows, plot = FALSE)
        if (length(temp0$acf[-1, 1, 1]) == acf.rows) 
            acf.temp[, j] <- abs(temp0$acf[-1, 1, 1])
        ##ESS1[j] <- LaplacesDemon::ESS(x[, j])
        Rec.Thin[j] <- which(acf.temp[, j] <= 0.1)[1] * 1
    }
    Rec.Thin[which(is.na(Rec.Thin))] <- nrow(acf.temp)
    Rec.Thin
}

## samples of marginal and conditional full-population freq. distributions
samplesFDistribution2 <- function(Y, X=NULL, mcsamples, varinfo, subsamples=1:nrow(mcsamples), jacobian=TRUE, fn=identity){
    if(length(dim(Y)) != 2){stop('Y must have two dimensions')}
    if(!is.null(X) && length(dim(X)) != 2){stop('X must be NULL or have two dimensions')}
    ##
    if(!is.null(X) && ncol(X) == 0){X <- NULL}
    if(length(subsamples) == 1 && !is.numeric(subsamples)){
        subsamples <- seq(1, nrow(mcsamples), length.out=round(abs(as.complex(subsamples))))
    }
    mcsamples <- t(mcsamples[subsamples,,drop=FALSE])
    Yv <- colnames(Y)
    Yvn <- length(Yv)
    Xv <- colnames(X)
    Xvn <- length(Xv)
    Vv <- varinfo[['name']]
    if(length(intersect(Yv, Xv)) > 0){cat('WARNING: overlap in Y and X variates\n')}
    if(!all(Yv %in% Vv)){cat('Warning: unknown Y variates\n')}
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
    ndata <- nrow(Y2)
    ##
    ##
##foreach(y=t(Y2), x=t(X2), .combine=rbind, .inorder=T)%dopar%{
        ## ## for debugging
        ## for(iii in 1:nrow(Y2)){
        ## print(iii)
        ## iii <- iii+1
        ##     y <- t(Y2)[,1,drop=F]
        ##     x <- t(X2)[,1,drop=F]
            ##
    ##
        if(is.null(X)){
                probX <- log(c(W))
        }else{
            ## x <- x[!is.na(x),,drop=F]
            x <- t(X2)
            xv <- lapply(variatetypes, function(xx){
                out <- rownames(x)[varinfo[['type']][rownames(x)]==xx]
                names(out) <- NULL
                out
            })
            xn <- lapply(xv,length)
            names(xv) <- names(xn) <- variatetypes
            ##
            probX <- ( # rows: MCsamples, cols: clusters
                log(c(W)) + aperm(
                (if(xn$D > 0){
                     colSums(
                         array(
                             t(sapply(xv$D, function(v){
                                 if(is.finite(x[v,])){
                                     (dnorm(x=x[v,],
                                             mean=mcsamples[Dmean[rep(v,ndata),],],
                                             sd=sqrt(mcsamples[Dvar[rep(v,ndata),],]),log=T))
                                 }else{
                                     (pnorm(q=Dbounds[rep(v,ndata),]*sign(x[v,]),
                                             mean=mcsamples[Dmean[rep(v,ndata),],],
                                             sd=sqrt(mcsamples[Dvar[rep(v,ndata),],]),
                                             lower.tail=(x[v,]<0),
                                             log.p=T))
                                 }
                             })),
                             dim=c(xn$D, ndata, nclusters, length(subsamples))),
                         na.rm=T)
                 }else{0}) +
                ## (if(xn$O > 0){
                ##      colSums(
                ##          array(
                ##              t(sapply(xv$O, function(v){
                ##                  if(is.finite(x[v,])){
                ##                      (dnorm(x=x[v,],
                ##                              mean=mcsamples[XOmean[v,],],
                ##                              sd=sqrt(mcsamples[XOvar[v,],]),log=T))
                ##                  }else{
                ##                      (pnorm(q=XOlefts[1,v],
                ##                              mean=mcsamples[XOmean[v,],],
                ##                              sd=sqrt(mcsamples[XOvar[v,],]),
                ##                              lower.tail=F,
                ##                              log.p=T))
                ##                  }
                ##              })),
                ##              dim=c(xn$O, nclusters, length(subsamples))),
                ##          na.rm=F)
                ##  }else{0}) +
                (if(xn$R > 0){
                     colSums(
                         array(dnorm(x=x[xv$R,],
                                     mean=mcsamples[Rmean[rep(xv$R,ndata),],],
                                     sd=sqrt(mcsamples[Rvar[rep(xv$R,ndata),],]),log=T),
                               dim=c(xn$R, ndata, nclusters, length(subsamples))),
                         na.rm=T)
                 }else{0}) +
            (if(xn$B > 0){
                 colSums(
                     array(log( c(x[xv$B,])*mcsamples[Bprob[rep(xv$B,ndata),],] +
                                c(1-x[xv$B,])*(1-mcsamples[Bprob[rep(xv$B,ndata),],]) ),
                           dim=c(xn$B, ndata, nclusters, length(subsamples))),
                     na.rm=T)
                 }else{0}) +
                (if(xn$I > 0){
                     colSums(
                         array(
                             t(sapply(xv$I, function(v){
                                 log(pnorm(q=Irights[[v]][x[v,]],
                                           mean=mcsamples[Imean[rep(v,ndata),],],
                                           sd=sqrt(mcsamples[Ivar[rep(v,ndata),],])) -
                                  pnorm(q=Ilefts[[v]][x[v,]],
                                        mean=mcsamples[Imean[rep(v,ndata),],],
                                        sd=sqrt(mcsamples[Ivar[rep(v,ndata),],])))
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
                             dim=c(xn$I, ndata, nclusters, length(subsamples))),
                         na.rm=T)
                 }else{0})
            , c(2,3,1))
            )# end probX
        }
    ##
    y <- t(Y2)
    yv <- lapply(variatetypes, function(xx){
        out <- rownames(y)[varinfo[['type']][rownames(y)]==xx]
        names(out) <- NULL
        out
    })
    yn <- lapply(yv,length)
    names(yv) <- names(yn) <- variatetypes
    ##
    probY <- aperm( # rows: MCsamples, cols: clusters
    (if(yn$D > 0){
         colSums(
             array(
                 t(sapply(yv$D, function(v){
                     if(is.finite(y[v,])){
                         (dnorm(x=y[v,],
                                mean=mcsamples[Dmean[rep(v,ndata),],],
                                sd=sqrt(mcsamples[Dvar[rep(v,ndata),],]),log=T))
                     }else{
                         (pnorm(q=Dbounds[rep(v,ndata),]*sign(y[v,]),
                                mean=mcsamples[Dmean[rep(v,ndata),],],
                                sd=sqrt(mcsamples[Dvar[rep(v,ndata),],]),
                                lower.tail=(y[v,]<0),
                                log.p=T))
                     }
                 })),
                 dim=c(yn$D, ndata, nclusters, length(subsamples))),
             na.rm=T)
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
                                     mean=mcsamples[Rmean[rep(yv$R,ndata),],],
                                     sd=sqrt(mcsamples[Rvar[rep(yv$R,ndata),],]),log=T),
                               dim=c(yn$R, ndata, nclusters, length(subsamples))),
                         na.rm=T)
                 }else{0}) +
                (if(yn$B > 0){
                     colSums(
                         array(log( c(y[yv$B,])*mcsamples[Bprob[rep(yv$B,ndata),],] +
                               c(1-y[yv$B,])*(1-mcsamples[Bprob[rep(yv$B,ndata),],]) ),
                               dim=c(yn$B, ndata, nclusters, length(subsamples))),
                         na.rm=T)
                 }else{0}) +
                (if(yn$I > 0){
                     colSums(
                         array(
                             t(sapply(yv$I, function(v){
                                 log(pnorm(q=Irights[[v]][y[v,]],
                                        mean=mcsamples[Imean[rep(v,ndata),],],
                                        sd=sqrt(mcsamples[Ivar[rep(v,ndata),],])) -
                                  pnorm(q=Ilefts[[v]][y[v,]],
                                        mean=mcsamples[Imean[rep(v,ndata),],],
                                        sd=sqrt(mcsamples[Ivar[rep(v,ndata),],])))
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
                             dim=c(yn$I, ndata, nclusters, length(subsamples))),
                         na.rm=T)
                 }else{0})
        , c(2,3,1)) # end probY
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
        if(is.null(X)){
            out <- colSums(exp(probX+probY)) 
        }else{
            probX <- apply(probX, c(2,3), function(xxx){xxx-max(xxx, na.rm=T)})
            out <- colSums(exp(probX+probY))/colSums(exp(probX))
        }
    fn(t(out)) *
        (if(jacobian){exp(-rowSums(
                           log(invjacobian(Y,varinfo=varinfo))
                           , na.rm=T))}else{1L})
}




## samples of marginal and conditional full-population freq. distributions
samplesFDistribution <- function(Y, X=NULL, mcsamples, varinfo, subsamples=1:nrow(mcsamples), jacobian=TRUE, fn=identity){
    if(length(dim(Y)) != 2){stop('Y must have two dimensions')}
    if(!is.null(X) && length(dim(X)) != 2){stop('X must be NULL or have two dimensions')}
    ##
    if(!is.null(X) && ncol(X) == 0){X <- NULL}
    if(length(subsamples) == 1 && !is.numeric(subsamples)){
        subsamples <- seq(1, nrow(mcsamples), length.out=round(abs(as.complex(subsamples))))
    }
    mcsamples <- t(mcsamples[subsamples,,drop=FALSE])
    Yv <- colnames(Y)
    Yvn <- length(Yv)
    Xv <- colnames(X)
    Xvn <- length(Xv)
    Vv <- varinfo[['name']]
    if(length(intersect(Yv, Xv)) > 0){cat('WARNING: overlap in Y and X variates\n')}
    if(!all(Yv %in% Vv)){cat('Warning: unknown Y variates\n')}
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
                (if(xn$D > 0){
                     colSums(
                         array(
                             t(sapply(xv$D, function(v){
                                 if(is.finite(x[v,])){
                                     (dnorm(x=x[v,],
                                             mean=mcsamples[Dmean[v,],],
                                             sd=sqrt(mcsamples[Dvar[v,],]),log=T))
                                 }else{
                                     (pnorm(q=Dbounds[v,]*sign(x[v,]),
                                             mean=mcsamples[Dmean[v,],],
                                             sd=sqrt(mcsamples[Dvar[v,],]),
                                             lower.tail=(x[v,]<0),
                                             log.p=T))
                                 }
                             })),
                             dim=c(xn$D, nclusters, length(subsamples))),
                         na.rm=F)
                 }else{0}) +
                ## (if(xn$O > 0){
                ##      colSums(
                ##          array(
                ##              t(sapply(xv$O, function(v){
                ##                  if(is.finite(x[v,])){
                ##                      (dnorm(x=x[v,],
                ##                              mean=mcsamples[XOmean[v,],],
                ##                              sd=sqrt(mcsamples[XOvar[v,],]),log=T))
                ##                  }else{
                ##                      (pnorm(q=XOlefts[1,v],
                ##                              mean=mcsamples[XOmean[v,],],
                ##                              sd=sqrt(mcsamples[XOvar[v,],]),
                ##                              lower.tail=F,
                ##                              log.p=T))
                ##                  }
                ##              })),
                ##              dim=c(xn$O, nclusters, length(subsamples))),
                ##          na.rm=F)
                ##  }else{0}) +
                (if(xn$R > 0){
                     colSums(
                         array(dnorm(x=x[xv$R,],
                                     mean=mcsamples[Rmean[xv$R,],],
                                     sd=sqrt(mcsamples[Rvar[xv$R,],]),log=T),
                               dim=c(xn$R, nclusters, length(subsamples))),
                         na.rm=F)
                 }else{0}) +
                (if(xn$B > 0){
                     colSums(
                         array(log( x[xv$B,]*mcsamples[Bprob[xv$B,],] +
                               (1-x[xv$B,])*(1-mcsamples[Bprob[xv$B,],]) ),
                               dim=c(xn$B, nclusters, length(subsamples))),
                         na.rm=F)
                 }else{0}) +
                (if(xn$I > 0){
                     colSums(
                         array(
                             t(sapply(xv$I, function(v){
                                 log(pnorm(q=Irights[[v]][x[v,]],
                                        mean=mcsamples[Imean[v,],],
                                        sd=sqrt(mcsamples[Ivar[v,],])) -
                                  pnorm(q=Ilefts[[v]][x[v,]],
                                        mean=mcsamples[Imean[v,],],
                                        sd=sqrt(mcsamples[Ivar[v,],])))
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
                             dim=c(xn$I, nclusters, length(subsamples))),
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



## Samples of variates
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
    nclusters <- length(Wi)
    W <- mcsamples[Wi,] # rows: clusters, cols: MC samples
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
    if(Yn$C > 0){## categorical
    totake <- sapply(Yv$C,function(x)which(variate$C == x))
    YCprob <- sapply(paste0('Cprob\\[',totake,','),
                            grep,rownames(mcsamples))
    ncategories <- length(YCprob)/length(totake)/nclusters
    seqcategories <- seq_len(ncategories)
    YCprob <- aperm(array(YCprob,
                           dim=c(nclusters,ncategories,length(totake)),
                           dimnames=list(NULL,NULL,Yv$C)),
                    c(2,1,3))
    }
    if(Yn$B > 0){## binary
    totake <- sapply(Yv$B,function(x)which(variate$B == x))
    YBprob <- array(t(vapply(paste0('Bprob\\[',totake,','), grep,
                           numeric(nclusters),
                           rownames(mcsamples))),
                     dim=c(length(totake),nclusters), dimnames=list(Yv$B,NULL))
    }
    if(Yn$I > 0){## integer ordinal
    totake <- sapply(Yv$I,function(x)which(variate$I == x))
    YImean <- array(t(vapply(paste0('Imean\\[',totake,','), grep,
                           numeric(nclusters),
                           rownames(mcsamples))),
                     dim=c(length(totake),nclusters), dimnames=list(Yv$I,NULL))
    YIvar <- array(t(vapply(paste0('Ivar\\[',totake,','), grep,
                           numeric(nclusters),
                           rownames(mcsamples))),
                     dim=c(length(totake),nclusters), dimnames=list(Yv$I,NULL))
    ## YIlefts <- lapply(Yv$I,function(v){ qnorm((0:(varinfo[['n']][v]-1L))/varinfo[['n']][v]) })
    ## YIrights <- lapply(Yv$I,function(v){ qnorm((1:(varinfo[['n']][v]))/varinfo[['n']][v]) })
    ## names(YIlefts) <- names(YIrights) <- Yv$I
    }
    if(Yn$R > 0){## real
    totake <- sapply(Yv$R,function(x)which(variate$R == x))
    YRmean <- array(t(vapply(paste0('Rmean\\[',totake,','), grep,
                           numeric(nclusters),
                           rownames(mcsamples))),
                     dim=c(length(totake),nclusters), dimnames=list(Yv$R,NULL))
    YRvar <- array(t(vapply(paste0('Rvar\\[',totake,','), grep,
                           numeric(nclusters),
                           rownames(mcsamples))),
                     dim=c(length(totake),nclusters), dimnames=list(Yv$R,NULL))
    }
    if(Yn$O > 0){## one-side censored
    totake <- sapply(Yv$O,function(x)which(variate$O == x))
    YOmean <- array(t(vapply(paste0('Omean\\[',totake,','), grep,
                           numeric(nclusters),
                           rownames(mcsamples))),
                     dim=c(length(totake),nclusters), dimnames=list(Yv$O,NULL))
    YOvar <- array(t(vapply(paste0('Ovar\\[',totake,','), grep,
                           numeric(nclusters),
                           rownames(mcsamples))),
                     dim=c(length(totake),nclusters), dimnames=list(Yv$O,NULL))
    ## YOlefts <- c(sapply(Yv$O, function(v){
    ##     out <- varinfo[['t']][[v]](varinfo[['max']][v])
    ##     names(out) <- NULL
    ##     out}))
    }
    if(Yn$D > 0){## two-bounded
    totake <- sapply(Yv$D,function(x)which(variate$D == x))
    YDmean <- array(t(vapply(paste0('Dmean\\[',totake,','), grep,
                           numeric(nclusters),
                           rownames(mcsamples))),
                     dim=c(length(totake),nclusters), dimnames=list(Yv$D,NULL))
    YDvar <- array(t(vapply(paste0('Dvar\\[',totake,','), grep,
                           numeric(nclusters),
                           rownames(mcsamples))),
                   dim=c(length(totake),nclusters), dimnames=list(Yv$D,NULL))
    }
    ##
    ##
    if(Xn$C > 0){## categorical
    totake <- sapply(Xv$C,function(x)which(variate$C == x))
    XCprob <- sapply(paste0('Cprob\\[',totake,','),
                            grep,rownames(mcsamples))
    ncategories <- length(XCprob)/length(totake)/nclusters
    seqcategories <- seq_len(ncategories)
    XCprob <- aperm(array(XCprob,
                           dim=c(nclusters,ncategories,length(totake)),
                           dimnames=list(NULL,NULL,Xv$C)),
                    c(2,1,3))
    }
    if(Xn$B > 0){## binary
    totake <- sapply(Xv$B,function(x)which(variate$B == x))
    XBprob <- array(t(vapply(paste0('Bprob\\[',totake,','), grep,
                           numeric(nclusters),
                           rownames(mcsamples))),
                     dim=c(length(totake),nclusters), dimnames=list(Xv$B,NULL))
    }
    if(Xn$I > 0){## integer ordinal
    totake <- sapply(Xv$I,function(x)which(variate$I == x))
    XImean <- array(t(vapply(paste0('Imean\\[',totake,','), grep,
                           numeric(nclusters),
                           rownames(mcsamples))),
                     dim=c(length(totake),nclusters), dimnames=list(Xv$I,NULL))
    XIvar <- array(t(vapply(paste0('Ivar\\[',totake,','), grep,
                           numeric(nclusters),
                           rownames(mcsamples))),
                     dim=c(length(totake),nclusters), dimnames=list(Xv$I,NULL))
    XIlefts <- lapply(Xv$I,function(v){ qnorm((0:(varinfo[['n']][v]-1L))/varinfo[['n']][v]) })
    XIrights <- lapply(Xv$I,function(v){ qnorm((1:(varinfo[['n']][v]))/varinfo[['n']][v]) })
    names(XIlefts) <- names(XIrights) <- Xv$I
    }
    if(Xn$R > 0){## real
    totake <- sapply(Xv$R,function(x)which(variate$R == x))
    XRmean <- array(t(vapply(paste0('Rmean\\[',totake,','), grep,
                           numeric(nclusters),
                           rownames(mcsamples))),
                     dim=c(length(totake),nclusters), dimnames=list(Xv$R,NULL))
    XRvar <- array(t(vapply(paste0('Rvar\\[',totake,','), grep,
                           numeric(nclusters),
                           rownames(mcsamples))),
                     dim=c(length(totake),nclusters), dimnames=list(Xv$R,NULL))
    }
    if(Xn$O > 0){## one-side censored
    totake <- sapply(Xv$O,function(x)which(variate$O == x))
    XOmean <- array(t(vapply(paste0('Omean\\[',totake,','), grep,
                           numeric(nclusters),
                           rownames(mcsamples))),
                     dim=c(length(totake),nclusters), dimnames=list(Xv$O,NULL))
    XOvar <- array(t(vapply(paste0('Ovar\\[',totake,','), grep,
                           numeric(nclusters),
                           rownames(mcsamples))),
                     dim=c(length(totake),nclusters), dimnames=list(Xv$O,NULL))
    XOlefts <- transf(rbind(varinfo[['tmax']][Xv$O]),varinfo,Oout='')
    }
    if(Xn$D > 0){## two-bounded
    totake <- sapply(Xv$D,function(x)which(variate$D == x))
    XDmean <- array(t(vapply(paste0('Dmean\\[',totake,','), grep,
                           numeric(nclusters),
                           rownames(mcsamples))),
                     dim=c(length(totake),nclusters), dimnames=list(Xv$D,NULL))
    XDvar <- array(t(vapply(paste0('Dvar\\[',totake,','), grep,
                           numeric(nclusters),
                           rownames(mcsamples))),
                   dim=c(length(totake),nclusters), dimnames=list(Xv$D,NULL))
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
        probX <- t(exp( # rows: MCsamples, cols: clusters
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
                             dim=c(xn$D, nclusters, length(subsamples))),
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
                             dim=c(xn$O, nclusters, length(subsamples))),
                         na.rm=F)
                 }else{0}) +
                (if(xn$R > 0){
                     colSums(
                         array(dnorm(x=x[xv$R,],
                                     mean=mcsamples[XRmean,],
                                     sd=sqrt(mcsamples[XRvar,]),log=T),
                               dim=c(xn$R, nclusters, length(subsamples))),
                         na.rm=F)
                 }else{0}) +
                (if(xn$B > 0){
                     colSums(
                         array(log( x[xv$B,]*mcsamples[XBprob,] +
                               (1-x[xv$B,])*(1-mcsamples[XBprob,]) ),
                               dim=c(xn$B, nclusters, length(subsamples))),
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
                             dim=c(xn$I, nclusters, length(subsamples))),
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








#### OLD FUNCTIONS ####
## if(FALSE){

## ##
## ## Construct a list of parameter samples from the raw MCMC samples for the second monitored set
## finalstate2list <- function(finalstate, realVars, integerVars, categoryVars, binaryVars, compoundgamma){
##     if(!is.vector(finalstate)){print('ERROR!')}
##     nclusters <- sum(grepl('^q\\[', names(finalstate)))
##     ncategories <- sum(grepl('^probC\\[1, 1, ', colnames(mcsamples)))
##     nrcovs <- length(realVars)
##     nicovs <- length(integerVars)
##     nccovs <- length(categoryVars)
##     nbcovs <- length(binaryVars)
##     ##
##     parmNames <- c('q',
##     (if(nrcovs > 0){c('meanR', 'varR', (if(compoundgamma){'varRrate'}))}),
##     (if(nicovs > 0){c('probI', 'sizeI')}),
##     (if(nccovs > 0){'probC'}),
##     (if(nbcovs > 0){'probB'}),
##     'C'
##     )
##     parmList <- foreach(var=parmNames)%do%{
##         out <- finalstate[grepl(paste0('^',var,'\\['), names(finalstate))]
##         if(var == 'meanR'||var == 'varR'){
##             dim(out) <- c(nrcovs, nclusters)
##             dimnames(out) <- list(realVars, NULL)
##         }else if(var == 'probI'||var == 'sizeI'){
##             dim(out) <- c(nicovs, nclusters)
##             dimnames(out) <- list(integerVars, NULL)
##         }else if(var == 'probC'){
##             dim(out) <- c(nccovs, nclusters, ncategories)
##             dimnames(out) <- list(categoryVars, NULL, NULL)
##         }else if(var == 'probB'){
##             dim(out) <- c(nbcovs, nclusters)
##             dimnames(out) <- list(binaryVars, NULL)
##         } # 'q' and 'C' and 'varRrate' are vectors with no names
##         out
##     }
##     names(parmList) <- parmNames
##     parmList
## }
## ##
## ## Calculates the probability of the data (likelihood of parameters) for several MCMC samples
## llSamples <- function(dat, parmList){
##     realCovs <- dimnames(parmList$meanR)[[2]]
##     integerCovs <- dimnames(parmList$probI)[[2]]
##     categoryCovs <- dimnames(parmList$probC)[[2]]
##     binaryCovs <- dimnames(parmList$probB)[[2]]
##     covNames <- c(realCovs, integerCovs, categoryCovs, binaryCovs)
##     nrcovs <- length(realCovs)
##     nicovs <- length(integerCovs)
##     nccovs <- length(categoryCovs)
##     nbcovs <- length(binaryCovs)
##     ndata <- nrow(dat[[which.max(c(nrcovs,nicovs,nccovs,nbcovs))]])
##     q <- parmList$q
##     ##
##     foreach(asample=seq_len(nrow(q)), .combine=c, .inorder=TRUE)%dopar%{
##         sum( log( colSums(
##             exp(
##                 log(q[asample,]) +
##                 t(vapply(seq_len(ncol(q)), function(acluster){
##                     #### rows: variates, columns: data
##                     ## real variates
##                     (if(nrcovs > 0){
##                         colSums(dnorm(t(dat$Real), mean=parmList$meanR[asample,,acluster], sd=sqrt(parmList$varR[asample,,acluster]), log=TRUE), na.rm=T)
##                         }else{0}) +
##                         ## integer variates
##                         (if(nicovs > 0){
##                             colSums(dbinom(t(dat$Integer), prob=parmList$probI[asample,,acluster], size=parmList$sizeI[asample,,acluster], log=TRUE), na.rm=T)
##                             }else{0}) +
##                         ## category variates
##                         (if(nccovs > 0){
##                             rowSums(log(sapply(seq_len(nccovs), function(acov){
##                                 parmList$probC[asample,acov,acluster,dat$Category[,acov]]})), na.rm=T)
##                             }else{0}) +
##                     ## binary variates
##                     (if(nbcovs > 0){
##                             colSums(log(
##                                 parmList$probB[asample,,acluster] * t(dat$Binary) +
##                                 (1-parmList$probB[asample,,acluster]) * (1-t(dat$Binary))
##                                 ), na.rm=T)
##                     }else{0})
##     }, numeric(ndata)))
##             )
##         ) ) )
##     }
## }
## ##
## ## Calculates the probability of the data (likelihood of parameters) for several MCMC samples
## ## using mcsamples directly
## llSamplesmc <- function(dat, mcsamples){
##     nrcovs <- sum(ncol(dat$Real))
##     nicovs <- sum(ncol(dat$Integer))
##     nccovs <- sum(ncol(dat$Category))
##     nbcovs <- sum(ncol(dat$Binary))
##     ndata <- nrow(dat[[which.max(c(nrcovs,nicovs,nccovs,nbcovs))]])
##     ##
##     Qi <- grep('q',colnames(mcsamples))
##     nclusters <- length(Qi)
##     sclusters <- seq_len(nclusters)
##     if(nrcovs > 0){
##         meanRi <- grep('meanR',colnames(mcsamples))
##         varRi <- grep('varR',colnames(mcsamples))
##         dim(meanRi) <- dim(varRi) <- c(nrcovs, nclusters)
##         }
##     if(nicovs > 0){
##         probIi <- grep('probI',colnames(mcsamples))
##         sizeIi <- grep('sizeI',colnames(mcsamples))
##         dim(probIi) <- dim(sizeIi) <- c(nicovs,nclusters)
##         }
##     if(nccovs > 0){
##         probCi <- grep('probC',colnames(mcsamples))
##         dim(probCi) <- c(nccovs,nclusters,length(probCi)/nccovs/nclusters)
##         }
##     if(nbcovs > 0){
##         probBi <- grep('probB',colnames(mcsamples))
##         dim(probBi) <- c(nbcovs,nclusters)
##         }
##     ##
##     foreach(asample=t(mcsamples), .combine=c, .inorder=TRUE)%dopar%{
##         sum( log( colSums(
##             exp(
##                 log(asample[Qi]) +
##                 t(vapply(sclusters, function(acluster){
##                     #### rows: variates, columns: data
##                     ## real variates
##                     (if(nrcovs > 0){
##                          colSums(dnorm(t(dat$Real), mean=asample[meanRi[,acluster]], sd=sqrt(asample[varRi[,acluster]]), log=TRUE), na.rm=T)
##                         }else{0}) +
##                         ## integer variates
##                         (if(nicovs > 0){
##                             colSums(dbinom(t(dat$Integer), prob=asample[probIi[,acluster]], size=asample[sizeIi[,acluster]], log=TRUE), na.rm=T)
##                             }else{0}) +
##                         ## category variates
##                         (if(nccovs > 0){
##                              rowSums(log(sapply(seq_len(nccovs), function(acov){
##                                  asample[probCi[acov,acluster,dat$Category[,acov]]]})), na.rm=T)
##                             }else{0}) +
##                     ## binary variates
##                     (if(nbcovs > 0){
##                          colSums(log(
##                              asample[probBi[,acluster]] * t(dat$Binary) +
##                              (1-asample[probBi[,acluster]]) * (1-t(dat$Binary))
##                                 ), na.rm=T)
##                     }else{0})
##     }, numeric(ndata)))
##             )
##         ) ) )
##     }
## }
## ##
## ## Calculate the distribution of means and variances
## samplesMVmc <- function(Ynames=NULL, X=NULL, mcsamples, variateparameters, fromsamples=nrow(mcsamples), inorder=FALSE){
##     allvariates <- rownames(variateparameters)
##     if(!is.null(X)){
##         X <- data.matrix(rbind(X))
##     }
##     Xnames <- colnames(X)
##     if(is.null(Ynames)){
##         Ynames <- allvariates
##     }
##     if(!all(Ynames %in% allvariates)){
##         warning('*Some requested variates missing from variateparameters argument*')
##         warning(paste0('*Discarding: ',paste0(setdiff(Ynames,allvariates),collapse=' ')))
##         Ynames <- intersect(Ynames,allvariates)
##     }
##     if(!all(Xnames %in% allvariates)){
##         warning('*Some conditional variates missing from variateparameters argument*')
##         warning(paste0('*Discarding: ',paste0(setdiff(Xnames,allvariates),collapse=' ')))
##         Xnames <- intersect(Xnames,allvariates)
##         X <- X[,Xnames,drop=F]
##     }
##     if(length(intersect(Ynames,Xnames)) > 0){
##         warning('*Predictor and predictand have variates in common.*')
##         warning(paste0('*Discarding: ',paste0(intersect(Ynames,Xnames),collapse=' '),' from predictor'))
##         Xnames <- setdiff(Xnames,Ynames)
##     }
##     if(length(Xnames) == 0){
##         X <- Xnames <- NULL
##     }
##     ##
##     YXnames <- c(Ynames, Xnames)
##     ## variateparameters <- variateparameters[YXnames,,drop=F]
##     ##
##     rY <- Ynames[variateparameters[Ynames,'type'] == 0]
##     iY <- Ynames[variateparameters[Ynames,'type'] == 3]
##     cY <- Ynames[variateparameters[Ynames,'type'] == 1]
##     bY <- Ynames[variateparameters[Ynames,'type'] == 2]
##     cNamesY <- c(rY, iY, cY, bY)
##     nrY <- length(rY)
##     niY <- length(iY)
##     ncY <- length(cY)
##     nbY <- length(bY)
##     ##
##     rX <- Xnames[variateparameters[Xnames,'type'] == 0]
##     iX <- Xnames[variateparameters[Xnames,'type'] == 3]
##     cX <- Xnames[variateparameters[Xnames,'type'] == 1]
##     bX <- Xnames[variateparameters[Xnames,'type'] == 2]
##     ## cNamesX <- c(rX, iX, cX, bX)
##     nrX <- length(rX)
##     niX <- length(iX)
##     ncX <- length(cX)
##     nbX <- length(bX)
##     ##
##     ##
##     if(!('location' %in% colnames(variateparameters))){
##         locations <- integer(length(YXnames))
##         names(locations) <- YXnames
##     }else{locations <- variateparameters[,'location']
##         names(locations) <- rownames(variateparameters)
##     }
##     if(!('scale' %in% colnames(variateparameters))){
##         scales <- locations * 0L + 1L
##     }else{scales <- variateparameters[,'scale']
##             names(scales) <- rownames(variateparameters)
##     }
##     ##
##     if(!is.null(X)){
##         X <- (t(X)-locations[Xnames])/scales[Xnames]
##     }
##     ##
##     Qi <- grep('q',colnames(mcsamples))
##     nclusters <- length(Qi)
##     sclusters <- seq_len(nclusters)
##     if(nrY > 0){
##         totake <- variateparameters[rY,'index']
##         YmeanRi <- sapply(paste0('meanR\\[',totake,','),grep,colnames(mcsamples))
##         YvarRi <- sapply(paste0('varR\\[',totake,','),grep,colnames(mcsamples))
##         dim(YmeanRi) <- dim(YvarRi) <- c(nclusters,nrY)
##         colnames(YmeanRi) <- colnames(YvarRi) <- rY
##     }
##     if(niY > 0){
##         totake <- variateparameters[iY,'index']
##         YprobIi <- sapply(paste0('probI\\[',totake,','),grep,colnames(mcsamples))
##         YsizeIi <- sapply(paste0('sizeI\\[',totake,','),grep,colnames(mcsamples))
##         dim(YprobIi) <- dim(YsizeIi) <- c(nclusters,niY)
##         colnames(YprobIi) <- colnames(YsizeIi) <- iY
##     }
##     if(ncY > 0){
##         totake <- variateparameters[cY,'index']
##         YprobCi <- t(sapply(paste0('probC\\[',totake,','),grep,colnames(mcsamples)))
##         ncategories <- length(YprobCi)/ncY/nclusters
##         dim(YprobCi) <- c(ncY,nclusters,ncategories)
##         YprobCi <- aperm(YprobCi)
##         dim(YprobCi) <- c(ncategories*nclusters,ncY)
##         colnames(YprobCi) <- cY
##         scategories <- seq_len(ncategories)
##     }
##     if(nbY > 0){
##         totake <- variateparameters[bY,'index']
##         YprobBi <- sapply(paste0('probB\\[',totake,','),grep,colnames(mcsamples))
##         dim(YprobBi) <- c(nclusters,nbY)
##         colnames(YprobBi) <- bY
##         ## sbins <- 0:1
##     }
##     ##
##     if(nrX > 0){
##         totake <- variateparameters[rX,'index']
##         XmeanRi <- t(sapply(paste0('meanR\\[',totake,','),grep,colnames(mcsamples)))
##         XvarRi <- t(sapply(paste0('varR\\[',totake,','),grep,colnames(mcsamples)))
##         ## dim(meanRi) <- dim(varRi) <- c(nrcovs, nclusters)
##         ## rownames(meanRi) <- rownames(varRi) <- rCovs
##     }
##     if(niX > 0){
##         totake <- variateparameters[iX,'index']
##         XprobIi <- t(sapply(paste0('probI\\[',totake,','),grep,colnames(mcsamples)))
##         XsizeIi <- t(sapply(paste0('sizeI\\[',totake,','),grep,colnames(mcsamples)))
##         dim(XprobIi) <- dim(XsizeIi) <- c(niX,nclusters)
##         rownames(XprobIi) <- rownames(XsizeIi) <- iX
##     }
##     if(ncX > 0){
##         totake <- variateparameters[cX,'index']
##         XprobCi <- t(sapply(paste0('probC\\[',totake,','),grep,colnames(mcsamples)))
##         if(!exists('ncategories')){
##             ncategories <- length(XprobCi)/ncX/nclusters
##             dim(XprobCi) <- c(ncX,nclusters,ncategories)
##             dimnames(XprobCi) <- list(cX, NULL, NULL)
##             scategories <- seq_len(ncategories)
##         }
##     }
##     if(nbX > 0){
##         totake <- variateparameters[bX,'index']
##         XprobBi <- t(sapply(paste0('probB\\[',totake,','),grep,colnames(mcsamples)))
##         ## dim(probBi) <- c(nbcovs,nclusters)
##         ## rownames(probBi) <- bCovs
##         ## sbins <- 0:1
##     }
##     ##
##     if(length(fromsamples) == 1){
##         if(inorder){
##             fromsamples <- round(seq(1, nrow(mcsamples), length.out=fromsamples))
##         }else{
##             fromsamples <- sample(1:nrow(mcsamples),fromsamples,replace=(fromsamples > nrow(mcsamples)))
##         }
##     }
##     mcsamples <- t(mcsamples[fromsamples,,drop=F])
##     nsamples <- length(fromsamples)
##     ndata <- max(ncol(X),nsamples)
##     q <- mcsamples[Qi,] # rows=clusters
##     ##
##     if(!is.null(X)){
##             ## pX: cols=clusters, rows=datapoints
##             pX <- t(exp(
##                 log(q) + t(#rows=clusters
##                 vapply(sclusters, function(acluster){#cols=clusters
##                     ## real covariates (rows)
##                     (if(nrX > 0){
##                         colSums(dnorm(x=X[rX,,drop=F], mean=mcsamples[XmeanRi[,acluster],,drop=F], sd=sqrt(mcsamples[XvarRi[,acluster],,drop=F]), log=TRUE), na.rm=TRUE)
##                     }else{0}) +
##                         ## integer covariates (rows)
##                         (if(niX > 0){
##                             colSums(dbinom(x=X[iX,,drop=FALSE], prob=mcsamples[XprobIi[,acluster],,drop=F], size=mcsamples[XsizeIi[,acluster],,drop=F], log=TRUE), na.rm=TRUE)
##                         }else{0}) +
##                         ## category variates (cols)
##                         (if(ncX > 0){
##                              rowSums(sapply(cX,function(acov){
##                                  extraDistr::dcat(x=X[acov,], prob=t(mcsamples[XprobCi[acov,acluster,],,drop=F]), log=TRUE)}), na.rm=TRUE)
##                             }else{0}) +
##                         ## binary covariates (rows)
##                         (if(nbX > 0){
##                              colSums(matrix(log(
##                                  c(mcsamples[XprobBi[,acluster],]) * c(X[bX,]) +
##                                  c(1-mcsamples[XprobBi[,acluster],])* c(1-X[bX,])
##                              ), nrow=nbX),
##                              na.rm=TRUE)
##                          }else{0})
##                 }, numeric(ndata)))
##             ))
##         pX <- t(pX/rowSums(pX)) # rows=clusters, cols=datapoints
##     }else{
##         pX <- q
##     }
##     ## 
##     ## mcsamples <- t(mcsamples)
##     ##
##     if(nrY > 0){
##         rM <- t(sapply(rY, function(acov){
##             out <- mcsamples[YmeanRi[,acov],]
##             rbind(out2 <- colSums(pX*out),
##               colSums(pX * (mcsamples[YvarRi[,acov],] + out*out)) -out2*out2)
##         })) # Y, moments, samples
##     }else{rM <- NULL}
##     ##
##     if(niY > 0){
##         iM <- t(sapply(iY, function(acov){
##             out <- mcsamples[YprobIi[,acov],]*mcsamples[YsizeIi[acov,],]
##             rbind(out2 <- colSums(pX*out),
##                   colSums(pX*out*
##                           (out-mcsamples[YprobIi[,acov],]+1)) -out2*out2)
##         })) # Y, moments, samples
##     }else{iM <- NULL}
##     ##
##     if(ncY > 0){
##         cM <- t(sapply(cY, function(acov){
##             out <- matrix(mcsamples[YprobCi[,acov],],nrow=ncategories)*scategories
##             rbind(colSums(pX*colSums(out)),
##                   colSums(pX*colSums(out*scategories)))
##         })) # Y, moments, samples
##     }else{cM <- NULL}
##     if(exists('out')){rm(out)}
##     ##
##     if(nbY > 0){
##         bM <- t(sapply(bY, function(acov){
##             rbind(out2 <- colSums(pX*mcsamples[YprobBi[,acov],]),
##                   out2*(1-out2))
##         })) # Y, moments, samples
##     }else{bM <- NULL}
##     if(exists('out2')){rm(out2)}
##     ##
##     array(rbind(rM,iM,cM,bM),
##                 dim=c(length(cNamesY), 2, ndata),
##                 dimnames=list(cNamesY,c('mean','var'),NULL))[order(match(cNamesY,Ynames)),,,drop=F] * c(scales[Ynames],scales[Ynames]*scales[Ynames])+c(locations[Ynames],locations[Ynames]*0L)
## }
## ##
## ## tottime <- Sys.time()
## ## for(i in 1:1e6){extraDistr::dbern(x=testx, prob=testp, log=TRUE)}
## ## Sys.time()-tottime
## ## ## Time difference of 56.0237 secs
## ## tottime <- Sys.time()
## ## for(i in 1:1e6){log(testp*testx + (1-testp)*(1L-testx))}
## ## Sys.time()-tottime
## ## ## Time difference of 31.73911 secs
## ## tottime <- Sys.time()
## ## for(i in 1:1e6){log(testx + (-1L)^testx*(1-testp))}
## ## Sys.time()-tottime
## ## ## Time difference of 51.8646 secs
## ##
## ## Calculates means and covariances of the sampled frequency distributions
## moments12Samples <- function(parmList){
##     realCovs <- dimnames(parmList$meanR)[[2]]
##     integerCovs <- dimnames(parmList$probI)[[2]]
##     binaryCovs <- dimnames(parmList$probB)[[2]]
##     covNames <- c(realCovs, integerCovs, binaryCovs)
##     nrcovs <- length(realCovs)
##     nicovs <- length(integerCovs)
##     nbcovs <- length(binaryCovs)
##     ncovs <- length(covNames)
##     q <- t(parmList$q)
##     ##
## if(length(realCovs) > 0){
##     meansr <- aperm(parmList$meanR, c(3, 1, 2))
##     quadr <- aperm(1/parmList$tauR + parmList$meanR * parmList$meanR, c(3, 1, 2))
## }else{
##     meansr <- NULL
##     quadr <- NULL
## }
##     if(length(integerCovs) > 0){
##     meansi <- aperm(parmList$probI * parmList$sizeI, c(3, 1, 2))
##     quadi <- aperm(parmList$probI * parmList$sizeI *
##                     (1 + parmList$probI * (parmList$sizeI - 1)), c(3, 1, 2))
##     }else{
##         meansi <- NULL
##         quadi <- NULL
## }
##     if(length(binaryCovs) > 0){
##     meansb <- aperm(parmList$probB, c(3, 1, 2))
##     quadb <- aperm(parmList$probB * parmList$probB, c(3, 1, 2))
##     }else{
##         meansb <- NULL
##         quadb <- NULL
## }
##     clustermeans <- c(meansr, meansi, meansb)
##     dim(clustermeans) <- c(dim((q)), ncovs)
##     mixmeans <- colSums(c(q) * clustermeans)
##     dimnames(mixmeans) <- list(NULL, paste0('MEAN_', covNames))
##     ##
##     mixvars <- c(quadr, quadi, quadb)
##     dim(mixvars) <- c(dim((q)), ncovs)
##     mixvars <- colSums(c(q) * mixvars) - mixmeans*mixmeans
##     dimnames(mixvars) <- list(NULL, paste0('VAR_', covNames))
##     ##
##     mixcovars <- foreach(cov1=seq_len(ncovs-1), .combine=cbind)%:%foreach(cov2=(cov1+1):ncovs, .combine=cbind)%do%{
##        colSums(c(q)*clustermeans[,,cov1,drop=FALSE]*clustermeans[,,cov2,drop=FALSE]) - mixmeans[,cov1,drop=FALSE]*mixmeans[,cov2,drop=FALSE]
##     }
##     colnames(mixcovars) <- foreach(cov1=seq_len(ncovs-1), .combine=c)%:%foreach(cov2=(cov1+1):ncovs, .combine=c)%do%{paste0('COV_',covNames[cov1],'|',covNames[cov2])}
##     ##
##     list(
##         Dcovars=((matrix(colSums(c(q) * apply(
##         (clustermeans -
##          array(rep(c(mixmeans), each=nrow(q)), dim=dim(clustermeans)))/array(rep(sqrt(c(mixvars)), each=nrow(q)), dim=dim(clustermeans)),
##         c(1,2), prod)),
##         ncol=1, dimnames=list(NULL, 'Dcov')))),
##         means=mixmeans,
##         vars=(mixvars),
##         covars=mixcovars
##         )
## }
## ##
## ## Gives samples of frequency distributions of any set of Y conditional on any set of X
## ## uses mcsamples
## samplesFmc <- function(Y, X=NULL, mcsamples, variateparameters, fromsamples=nrow(mcsamples), inorder=FALSE){
##     ##
##     rCovs <- rownames(variateparameters)[variateparameters[,'type'] == 0]
##     iCovs <- rownames(variateparameters)[variateparameters[,'type'] == 3]
##     cCovs <- rownames(variateparameters)[variateparameters[,'type'] == 1]
##     bCovs <- rownames(variateparameters)[variateparameters[,'type'] == 2]
##     cNames <- c(rCovs, iCovs, cCovs, bCovs)
##     nrcovs <- length(rCovs)
##     nicovs <- length(iCovs)
##     nccovs <- length(cCovs)
##     nbcovs <- length(bCovs)
##     if(!('location' %in% colnames(variateparameters))){
##         locations <- integer(length(cNames))
##         names(locations) <- cNames
##     }else{locations <- variateparameters[,'location']
##         names(locations) <- rownames(variateparameters)}
##     if(!('scale' %in% colnames(variateparameters))){
##         scales <- locations * 0L + 1L
##     }else{scales <- variateparameters[,'scale']
##         names(scales) <- rownames(variateparameters)}
##     ##
##     Y <- data.matrix(rbind(Y))
##     cnY <- colnames(Y)
##     if(!all(cnY %in% rownames(variateparameters))){
##         warning('*WARNING: some requested predictands missing from variateparameters argument*')
##         warning(paste0('*Discarding: ',paste0(setdiff(cnY,rownames(variateparameters)),collapse=' ')))
##         cnY <- intersect(cnY,rownames(variateparameters))
##     }
##     Y <- t((t(Y[,cnY,drop=F])-locations[cnY])/scales[cnY])
##     rY <- cnY[cnY %in% rCovs]
##     iY <- cnY[cnY %in% iCovs]
##     cY <- cnY[cnY %in% cCovs]
##     bY <- cnY[cnY %in% bCovs]
##     ##
##     if(!is.null(X)){
##         X <- data.matrix(rbind(X))
##         cnX <- colnames(X)
##         if(!all(cnX %in% rownames(variateparameters))){
##             warning('*WARNING: some requested predictands missing from variateparameters argument*')
##             warning(paste0('*Discarding: ',paste0(setdiff(cnX,rownames(variateparameters)),collapse=' ')))
##             cnX <- intersect(cnX,rownames(variateparameters))
##         }
##         X <- t((t(X[,cnX,drop=F])-locations[cnX])/scales[cnX])
##         rX <- cnX[cnX %in% rCovs]
##         if(length(intersect(rX,rY)) > 0){
##             warning('*WARNING: predictor and predictand have real variates in common. Removing from predictor*')
##             rX <- setdiff(rX,rY)
##         }
##         iX <- cnX[cnX %in% iCovs]
##         if(length(intersect(iX,iY)) > 0){
##             warning('*WARNING: predictor and predictand have integer variates in common. Removing from predictor*')
##             iX <- setdiff(iX,iY)
##         }
##         cX <- cnX[cnX %in% cCovs]
##         if(length(intersect(cX,cY)) > 0){
##             warning('*WARNING: predictor and predictand have category variates in common. Removing from predictor*')
##             cX <- setdiff(cX,cY)
##         }
##         bX <- cnX[cnX %in% bCovs]
##         if(length(intersect(bX,bY)) > 0){
##             warning('*WARNING: predictor and predictand have binary variates in common. Removing from predictor*')
##             bX <- setdiff(bX,bY)
##         }
##         if(length(c(rX,iX,cX,bX)) == 0){X <- NULL}
##     }
##     ##
##     if(!is.null(X)){
##         if(nrow(X) < nrow(Y)){
##         warning('*Note: X has fewer data than Y. Recycling*')
##         X <- t(matrix(rep(t(X), ceiling(nrow(Y)/nrow(X))), nrow=ncol(X), dimnames=list(cnX,NULL)))[1:nrow(Y),,drop=FALSE]
##     }
##     if(nrow(X) > nrow(Y)){
##         warning('*Note: X has more data than Y. Recycling*')
##         Y <- t(matrix(rep(t(Y), ceiling(nrow(X)/nrow(Y))), nrow=ncol(Y), dimnames=list(cnY,NULL)))[1:nrow(X),,drop=FALSE]
##     }
##     }
##     ndata <- nrow(Y)
##     ##
##     Qi <- grep('q',colnames(mcsamples))
##     nclusters <- length(Qi)
##     sclusters <- seq_len(nclusters)
##     if(nrcovs > 0){
##         meanRi <- grep('meanR',colnames(mcsamples))
##         varRi <- grep('varR',colnames(mcsamples))
##         dim(meanRi) <- dim(varRi) <- c(nrcovs, nclusters)
##         rownames(meanRi) <- rownames(varRi) <- rCovs
##         }
##     if(nicovs > 0){
##         probIi <- grep('probI',colnames(mcsamples))
##         sizeIi <- grep('sizeI',colnames(mcsamples))
##         dim(probIi) <- dim(sizeIi) <- c(nicovs,nclusters)
##         rownames(probIi) <- rownames(sizeIi) <- iCovs
##         }
##     if(nccovs > 0){
##         probCi <- grep('probC',colnames(mcsamples))
##         dim(probCi) <- c(nccovs,nclusters,length(probCi)/nccovs/nclusters)
##         dimnames(probCi) <- list(cCovs, NULL, NULL)
##         }
##     if(nbcovs > 0){
##         probBi <- grep('probB',colnames(mcsamples))
##         dim(probBi) <- c(nbcovs,nclusters)
##         rownames(probBi) <- bCovs
##         }
##     ##
##     if(length(fromsamples) == 1){
##         if(inorder){
##             fromsamples <- round(seq(1, nrow(mcsamples), length.out=fromsamples))
##         }else{
##             fromsamples <- sample(1:nrow(mcsamples),fromsamples,replace=(fromsamples > nrow(mcsamples)))
##         }
##     }
##     ##
##     if(!is.null(X)){
##         freqs <- foreach(asample=t(mcsamples[fromsamples,,drop=F]), .combine=cbind, .inorder=inorder)%dopar%{
##             ## pX: rows=clusters, cols=datapoints
##             pX <- exp(
##                 log(asample[Qi]) +
##                 t(rbind(vapply(sclusters, function(acluster){
##                     ## real covariates
##                     (if(length(rX)>0){
##                         colSums(dnorm(x=t(X[,rX,drop=FALSE]), mean=asample[meanRi[rX,acluster]], sd=sqrt(asample[varRi[rX,acluster]]), log=TRUE), na.rm=TRUE)
##                     }else{0}) +
##                         ## integer covariates
##                         (if(length(iX)>0){
##                             colSums(dbinom(x=t(X[,iX,drop=FALSE]), prob=asample[probIi[iX,acluster]], size=asample[sizeIi[iX,acluster]], log=TRUE), na.rm=TRUE)
##                         }else{0}) +
##                         ## category variates
##                         (if(length(cX)>0){
##                              rowSums(log(sapply(cX, function(acov){
##                                  asample[probCi[acov,acluster,X[,acov]]]})), na.rm=T)
##                             }else{0}) +
##                         ## binary covariates
##                         (if(length(bX)>0){
##                              colSums(log(
##                                  asample[probBi[bX,acluster]] * t(X[,bX,drop=FALSE]) +
##                                  (1-asample[probBi[bX,acluster]]) * (1-t(X[,bX,drop=FALSE]))
##                              ), na.rm=TRUE)
##                         }else{0})
##                 }, numeric(ndata))))
##             )
##             ## pY: rows=clusters, cols=datapoints
##             pY <- exp(
##                 t(rbind(vapply(sclusters, function(acluster){
##                     ## real covariates
##                     (if(length(rY)>0){
##                         colSums(dnorm(x=t(Y[,rY,drop=FALSE]), mean=asample[meanRi[rY,acluster]], sd=sqrt(asample[varRi[rY,acluster]]), log=TRUE), na.rm=TRUE)
##                     }else{0}) +
##                         ## integer covariates
##                         (if(length(iY)>0){
##                             colSums(dbinom(x=t(Y[,iY,drop=FALSE]), prob=asample[probIi[iY,acluster]], size=asample[sizeIi[iY,acluster]], log=TRUE), na.rm=TRUE)
##                         }else{0}) +
##                         ## category variates
##                         (if(length(cY)>0){
##                             rowSums(log(sapply(cY, function(acov){
##                                  asample[probCi[acov,acluster,Y[,acov]]]})), na.rm=T)
##                             }else{0}) +
##                         ## binary covariates
##                         (if(length(bY)>0){
##                             colSums(log(
##                                 asample[probBi[bY,acluster]] * t(Y[,bY,drop=FALSE]) +
##                                 (1-asample[probBi[bY,acluster]]) * (1-t(Y[,bY,drop=FALSE]))
##                             ), na.rm=TRUE)
##                         }else{0})
##                 }, numeric(ndata))))
##             )
##             ##
##             colSums(pX * pY)/colSums(pX)
##         }
##     }else{
##         freqs <- foreach(asample=t(mcsamples[fromsamples,,drop=F]), .combine=cbind, .inorder=inorder)%dopar%{
##             ## pY: rows=clusters, cols=datapoints
##             pY <- exp(
##                 log(asample[Qi]) +
##                 t(rbind(vapply(sclusters, function(acluster){
##                     ## real covariates
##                     (if(length(rY)>0){
##                         colSums(dnorm(x=t(Y[,rY,drop=FALSE]), mean=asample[meanRi[rY,acluster]], sd=sqrt(asample[varRi[rY,acluster]]), log=TRUE), na.rm=TRUE)
##                     }else{0}) +
##                         ## integer covariates
##                         (if(length(iY)>0){
##                             colSums(dbinom(x=t(Y[,iY,drop=FALSE]), prob=asample[probIi[iY,acluster]], size=asample[sizeIi[iY,acluster]], log=TRUE), na.rm=TRUE)
##                         }else{0}) +
##                         ## category variates
##                         (if(length(cY)>0){
##                             rowSums(log(sapply(cY, function(acov){
##                                  asample[probCi[acov,acluster,Y[,acov]]]})), na.rm=T)
##                             }else{0}) +
##                         ## binary covariates
##                         (if(length(bY)>0){
##                             colSums(log(
##                                 asample[probBi[bY,acluster]] * t(Y[,bY,drop=FALSE]) +
##                                 (1-asample[probBi[bY,acluster]]) * (1-t(Y[,bY,drop=FALSE]))
##                             ), na.rm=TRUE)
##                         }else{0})
##                 }, numeric(ndata))))
##             )
##             ##
##             colSums(pY)
##         }
##     }
##     freqs/prod(scales[cnY])
## }
## ##
## ## Gives samples of frequency distributions of any set of Y conditional on any set of X
## samplesF <- function(Y, X=NULL, parmList, nfsamples=NULL, inorder=FALSE, transform=NULL){
##     rCovs <- dimnames(parmList$meanR)[[2]]
##     iCovs <- dimnames(parmList$probI)[[2]]
##     cCovs <- dimnames(parmList$probC)[[2]]
##     bCovs <- dimnames(parmList$probB)[[2]]
##     cNames <- c(rCovs, iCovs, cCovs, bCovs)
##     if(is.null(transform)){
##         transform <- matrix(c(0,1),nrow=length(cNames),ncol=2,byrow=TRUE,dimnames=list(cNames, c('location','scale')))
##     }
##     nrcovs <- length(rCovs)
##     nicovs <- length(iCovs)
##     nccovs <- length(cCovs)
##     nbcovs <- length(bCovs)
##     ##
##     Y <- data.matrix(rbind(Y))
##     Y <- t((t(Y)-transform[colnames(Y),'location'])/transform[colnames(Y),'scale'])
##     rY <- colnames(Y)[colnames(Y) %in% rCovs]
##     iY <- colnames(Y)[colnames(Y) %in% iCovs]
##     cY <- colnames(Y)[colnames(Y) %in% cCovs]
##     bY <- colnames(Y)[colnames(Y) %in% bCovs]
##     ##
##     if(!is.null(X)){
##         X <- data.matrix(rbind(X))
##         X <- t((t(X)-transform[colnames(X),'location'])/transform[colnames(X),'scale'])
##         rX <- colnames(X)[colnames(X) %in% rCovs]
##         if(length(intersect(rX,rY))>0){
##             warning('*WARNING: predictor and predictand have real variates in common. Removing from predictor*')
##             rX <- setdiff(rX,rY)
##         }
##         iX <- colnames(X)[colnames(X) %in% iCovs]
##         if(length(intersect(iX,iY))>0){
##             warning('*WARNING: predictor and predictand have integer variates in common. Removing from predictor*')
##             iX <- setdiff(iX,iY)
##         }
##         cX <- colnames(X)[colnames(X) %in% cCovs]
##         if(length(intersect(cX,cY))>0){
##             warning('*WARNING: predictor and predictand have category variates in common. Removing from predictor*')
##             cX <- setdiff(cX,cY)
##         }
##         bX <- colnames(X)[colnames(X) %in% bCovs]
##         if(length(intersect(bX,bY))>0){
##             warning('*WARNING: predictor and predictand have binary variates in common. Removing from predictor*')
##             bX <- setdiff(bX,bY)
##         }
##         if(length(c(rX,iX,cX,bX))==0){X <- NULL}
##     }
##     ##
##     if(!is.null(X)){
##         if(nrow(X) < nrow(Y)){
##         warning('*Note: X has fewer data than Y. Recycling*')
##         X <- t(matrix(rep(t(X), ceiling(nrow(Y)/nrow(X))), nrow=ncol(X), dimnames=list(colnames(X),NULL)))[1:nrow(Y),,drop=FALSE]
##     }
##     if(nrow(X) > nrow(Y)){
##         warning('*Note: X has more data than Y. Recycling*')
##         Y <- t(matrix(rep(t(Y), ceiling(nrow(X)/nrow(Y))), nrow=ncol(Y), dimnames=list(colnames(Y),NULL)))[1:nrow(X),,drop=FALSE]
##     }
##     }
##     nydata <- nrow(Y)
##     ##
##     q <- parmList$q
##     nclusters <- ncol(q)
##     sclusters <- seq_len(nclusters)
##     if(is.numeric(nfsamples)){
##         fsubsamples <- round(seq(1, nrow(q), length.out=nfsamples))
##     }else{
##         nfsamples <- nrow(q)
##         fsubsamples <- seq_len(nfsamples)
##     }
##     ##
##     if(!is.null(X)){
##         freqs <- foreach(asample=fsubsamples, .combine=cbind, .inorder=inorder)%dopar%{
##             ## pX: rows=clusters, cols=datapoints
##             pX <- exp(
##                 log(q[asample,]) +
##                 t(rbind(vapply(sclusters, function(acluster){
##                     ## real covariates
##                     (if(length(rX)>0){
##                         colSums(dnorm(x=t(X[,rX,drop=FALSE]), mean=parmList$meanR[asample,rX,acluster], sd=sqrt(parmList$varR[asample,rX,acluster]), log=TRUE), na.rm=TRUE)
##                     }else{0}) +
##                         ## integer covariates
##                         (if(length(iX)>0){
##                             colSums(dbinom(x=t(X[,iX,drop=FALSE]), prob=parmList$probI[asample,iX,acluster], size=parmList$sizeI[asample,iX,acluster], log=TRUE), na.rm=TRUE)
##                         }else{0}) +
##                         ## category variates
##                         (if(length(cX)>0){
##                             rowSums(log(sapply(cX, function(acov){
##                                 parmList$probC[asample,acov,acluster,X[,acov]]})), na.rm=T)
##                             }else{0}) +
##                         ## binary covariates
##                         (if(length(bX)>0){
##                             colSums(log(
##                                 parmList$probB[asample,bX,acluster] * t(X[,bX,drop=FALSE]) +
##                                 (1-parmList$probB[asample,bX,acluster]) * (1-t(X[,bX,drop=FALSE]))
##                             ), na.rm=TRUE)
##                         }else{0})
##                 }, numeric(nydata))))
##             )
##             ## pY: rows=clusters, cols=datapoints
##             pY <- exp(
##                 t(rbind(vapply(sclusters, function(acluster){
##                     ## real covariates
##                     (if(length(rY)>0){
##                         colSums(dnorm(x=t(Y[,rY,drop=FALSE]), mean=parmList$meanR[asample,rY,acluster], sd=sqrt(parmList$varR[asample,rY,acluster]), log=TRUE), na.rm=TRUE)
##                     }else{0}) +
##                         ## integer covariates
##                         (if(length(iY)>0){
##                             colSums(dbinom(x=t(Y[,iY,drop=FALSE]), prob=parmList$probI[asample,iY,acluster], size=parmList$sizeI[asample,iY,acluster], log=TRUE), na.rm=TRUE)
##                         }else{0}) +
##                         ## category variates
##                         (if(length(cY)>0){
##                             rowSums(log(sapply(cY, function(acov){
##                                 parmList$probC[asample,acov,acluster,Y[,acov]]})), na.rm=T)
##                             }else{0}) +
##                         ## binary covariates
##                         (if(length(bY)>0){
##                             colSums(log(
##                                 parmList$probB[asample,bY,acluster] * t(Y[,bY,drop=FALSE]) +
##                                 (1-parmList$probB[asample,bY,acluster]) * (1-t(Y[,bY,drop=FALSE]))
##                             ), na.rm=TRUE)
##                         }else{0})
##                 }, numeric(nydata))))
##             )
##             ##
##             colSums(pX * pY)/colSums(pX)
##         }
##     }else{
##         freqs <- foreach(asample=fsubsamples, .combine=cbind, .inorder=inorder)%dopar%{
##             ## pY: rows=clusters, cols=datapoints
##             pY <- exp(
##                 log(q[asample,]) +
##                 t(rbind(vapply(sclusters, function(acluster){
##                     ## real covariates
##                     (if(length(rY)>0){
##                         colSums(dnorm(x=t(Y[,rY,drop=FALSE]), mean=parmList$meanR[asample,rY,acluster], sd=sqrt(parmList$varR[asample,rY,acluster]), log=TRUE), na.rm=TRUE)
##                     }else{0}) +
##                         ## integer covariates
##                         (if(length(iY)>0){
##                             colSums(dbinom(x=t(Y[,iY,drop=FALSE]), prob=parmList$probI[asample,iY,acluster], size=parmList$sizeI[asample,iY,acluster], log=TRUE), na.rm=TRUE)
##                         }else{0}) +
##                         ## category variates
##                         (if(length(cY)>0){
##                             rowSums(log(sapply(cY, function(acov){
##                                 parmList$probC[asample,acov,acluster,Y[,acov]]})), na.rm=T)
##                             }else{0}) +
##                         ## binary covariates
##                         (if(length(bY)>0){
##                             colSums(log(
##                                 parmList$probB[asample,bY,acluster] * t(Y[,bY,drop=FALSE]) +
##                                 (1-parmList$probB[asample,bY,acluster]) * (1-t(Y[,bY,drop=FALSE]))
##                             ), na.rm=TRUE)
##                         }else{0})
##                 }, numeric(nydata))))
##             )
##             ##
##             colSums(pY)
##         }
##     }
##     freqs/prod(transform[colnames(Y),'scale'])
## }
## ##
## ## Gives samples of sums of log-frequency distributions of any set of Y conditional on any set of X
## ## uses mcsamples
## logsumsamplesFmc <- function(Y, X=NULL, mcsamples, variateparameters, fromsamples=nrow(mcsamples), inorder=FALSE){
##     ##
##     rCovs <- rownames(variateparameters)[variateparameters[,'type']==0]
##     iCovs <- rownames(variateparameters)[variateparameters[,'type']==3]
##     cCovs <- rownames(variateparameters)[variateparameters[,'type']==1]
##     bCovs <- rownames(variateparameters)[variateparameters[,'type']==2]
##     cNames <- c(rCovs, iCovs, cCovs, bCovs)
##     nrcovs <- length(rCovs)
##     nicovs <- length(iCovs)
##     nccovs <- length(cCovs)
##     nbcovs <- length(bCovs)
##     ##
##     Y <- data.matrix(rbind(Y))
##     rY <- colnames(Y)[colnames(Y) %in% rCovs]
##     iY <- colnames(Y)[colnames(Y) %in% iCovs]
##     cY <- colnames(Y)[colnames(Y) %in% cCovs]
##     bY <- colnames(Y)[colnames(Y) %in% bCovs]
##     ##
##     if(!is.null(X)){
##         X <- data.matrix(rbind(X))
##         rX <- colnames(X)[colnames(X) %in% rCovs]
##         if(length(intersect(rX,rY))>0){
##             warning('*WARNING: predictor and predictand have real variates in common. Removing from predictor*')
##             rX <- setdiff(rX,rY)
##         }
##         iX <- colnames(X)[colnames(X) %in% iCovs]
##         if(length(intersect(iX,iY))>0){
##             warning('*WARNING: predictor and predictand have integer variates in common. Removing from predictor*')
##             iX <- setdiff(iX,iY)
##         }
##         cX <- colnames(X)[colnames(X) %in% cCovs]
##         if(length(intersect(cX,cY))>0){
##             warning('*WARNING: predictor and predictand have category variates in common. Removing from predictor*')
##             cX <- setdiff(cX,cY)
##         }
##         bX <- colnames(X)[colnames(X) %in% bCovs]
##         if(length(intersect(bX,bY))>0){
##             warning('*WARNING: predictor and predictand have binary variates in common. Removing from predictor*')
##             bX <- setdiff(bX,bY)
##         }
##         if(length(c(rX,iX,cX,bX))==0){X <- NULL}
##     }
##     ##
##     if(!is.null(X)){
##         if(nrow(X) < nrow(Y)){
##         warning('*Note: X has fewer data than Y. Recycling*')
##         X <- t(matrix(rep(t(X), ceiling(nrow(Y)/nrow(X))), nrow=ncol(X), dimnames=list(colnames(X),NULL)))[1:nrow(Y),,drop=FALSE]
##     }
##     if(nrow(X) > nrow(Y)){
##         warning('*Note: X has more data than Y. Recycling*')
##         Y <- t(matrix(rep(t(Y), ceiling(nrow(X)/nrow(Y))), nrow=ncol(Y), dimnames=list(colnames(Y),NULL)))[1:nrow(X),,drop=FALSE]
##     }
##     }
##     ndata <- nrow(Y)
##     ##
##     Qi <- grep('q',colnames(mcsamples))
##     nclusters <- length(Qi)
##     sclusters <- seq_len(nclusters)
##     if(nrcovs>0){
##         meanRi <- grep('meanR',colnames(mcsamples))
##         varRi <- grep('varR',colnames(mcsamples))
##         dim(meanRi) <- dim(varRi) <- c(nrcovs, nclusters)
##         rownames(meanRi) <- rownames(varRi) <- rCovs
##         }
##     if(nicovs>0){
##         probIi <- grep('probI',colnames(mcsamples))
##         sizeIi <- grep('sizeI',colnames(mcsamples))
##         dim(probIi) <- dim(sizeIi) <- c(nicovs,nclusters)
##         rownames(probIi) <- rownames(sizeIi) <- iCovs
##         }
##     if(nccovs>0){
##         probCi <- grep('probC',colnames(mcsamples))
##         dim(probCi) <- c(nccovs,nclusters,length(probCi)/nccovs/nclusters)
##         dimnames(probCi) <- list(cCovs, NULL, NULL)
##         }
##     if(nbcovs>0){
##         probBi <- grep('probB',colnames(mcsamples))
##         dim(probBi) <- c(nbcovs,nclusters)
##         rownames(probBi) <- bCovs
##         }
##     ##
##     if(length(fromsamples) == 1){
##         if(inorder){
##             fromsamples <- round(seq(1, nrow(mcsamples), length.out=fromsamples))
##         }else{
##             fromsamples <- sample(1:nrow(mcsamples),fromsamples,replace=(fromsamples>nrow(mcsamples)))
##         }
##     }
##     ##
##     if(!is.null(X)){
##         freqs <- foreach(asample=t(mcsamples[fromsamples,,drop=F]), .combine=cbind, .inorder=inorder)%dopar%{
##             ## pX: rows=clusters, cols=datapoints
##             pX <- exp(
##                 log(asample[Qi]) +
##                 t(rbind(vapply(sclusters, function(acluster){
##                     ## real covariates
##                     (if(length(rX)>0){
##                         colSums(dnorm(x=t(X[,rX,drop=FALSE]), mean=asample[meanRi[rX,acluster]], sd=sqrt(asample[varRi[rX,acluster]]), log=TRUE), na.rm=TRUE)
##                     }else{0}) +
##                         ## integer covariates
##                         (if(length(iX)>0){
##                             colSums(dbinom(x=t(X[,iX,drop=FALSE]), prob=asample[probIi[iX,acluster]], size=asample[sizeIi[iX,acluster]], log=TRUE), na.rm=TRUE)
##                         }else{0}) +
##                         ## category variates
##                         (if(length(cX)>0){
##                              rowSums(log(sapply(cX, function(acov){
##                                  asample[probCi[acov,acluster,X[,acov]]]})), na.rm=T)
##                             }else{0}) +
##                         ## binary covariates
##                         (if(length(bX)>0){
##                              colSums(log(
##                                  asample[probBi[bX,acluster]] * t(X[,bX,drop=FALSE]) +
##                                  (1-asample[probBi[bX,acluster]]) * (1-t(X[,bX,drop=FALSE]))
##                              ), na.rm=TRUE)
##                         }else{0})
##                 }, numeric(ndata))))
##             )
##             ## pY: rows=clusters, cols=datapoints
##             pY <- exp(
##                 t(rbind(vapply(sclusters, function(acluster){
##                     ## real covariates
##                     (if(length(rY)>0){
##                         colSums(dnorm(x=t(Y[,rY,drop=FALSE]), mean=asample[meanRi[rY,acluster]], sd=sqrt(asample[varRi[rY,acluster]]), log=TRUE), na.rm=TRUE)
##                     }else{0}) +
##                         ## integer covariates
##                         (if(length(iY)>0){
##                             colSums(dbinom(x=t(Y[,iY,drop=FALSE]), prob=asample[probIi[iY,acluster]], size=asample[sizeIi[iY,acluster]], log=TRUE), na.rm=TRUE)
##                         }else{0}) +
##                         ## category variates
##                         (if(length(cY)>0){
##                             rowSums(log(sapply(cY, function(acov){
##                                  asample[probCi[acov,acluster,Y[,acov]]]})), na.rm=T)
##                             }else{0}) +
##                         ## binary covariates
##                         (if(length(bY)>0){
##                             colSums(log(
##                                 asample[probBi[bY,acluster]] * t(Y[,bY,drop=FALSE]) +
##                                 (1-asample[probBi[bY,acluster]]) * (1-t(Y[,bY,drop=FALSE]))
##                             ), na.rm=TRUE)
##                         }else{0})
##                 }, numeric(ndata))))
##             )
##             ##
##             sum(log(colSums(pX * pY)))-sum(log(colSums(pX)))
##         }
##     }else{
##         freqs <- foreach(asample=t(mcsamples[fromsamples,,drop=F]), .combine=cbind, .inorder=inorder)%dopar%{
##             ## pY: rows=clusters, cols=datapoints
##             pY <- exp(
##                 log(asample[Qi]) +
##                 t(rbind(vapply(sclusters, function(acluster){
##                     ## real covariates
##                     (if(length(rY)>0){
##                         colSums(dnorm(x=t(Y[,rY,drop=FALSE]), mean=asample[meanRi[rY,acluster]], sd=sqrt(asample[varRi[rY,acluster]]), log=TRUE), na.rm=TRUE)
##                     }else{0}) +
##                         ## integer covariates
##                         (if(length(iY)>0){
##                             colSums(dbinom(x=t(Y[,iY,drop=FALSE]), prob=asample[probIi[iY,acluster]], size=asample[sizeIi[iY,acluster]], log=TRUE), na.rm=TRUE)
##                         }else{0}) +
##                         ## category variates
##                         (if(length(cY)>0){
##                             rowSums(log(sapply(cY, function(acov){
##                                  asample[probCi[acov,acluster,Y[,acov]]]})), na.rm=T)
##                             }else{0}) +
##                         ## binary covariates
##                         (if(length(bY)>0){
##                             colSums(log(
##                                 asample[probBi[by,acluster]] * t(Y[,bY,drop=FALSE]) +
##                                 (1-asample[probBi[by,acluster]]) * (1-t(Y[,bY,drop=FALSE]))
##                             ), na.rm=TRUE)
##                         }else{0})
##                 }, numeric(ndata))))
##             )
##             ##
##             sum(log(colSums(pY)))
##         }
##     }
##     freqs
## }
## ##
## ## Gives samples of sums of log-frequency distributions of any set of Y conditional on any set of X
## logsumsamplesF <- function(Y, X=NULL, parmList, nfsamples=NULL, inorder=FALSE){
##     rCovs <- dimnames(parmList$meanR)[[2]]
##     iCovs <- dimnames(parmList$probI)[[2]]
##     cCovs <- dimnames(parmList$probC)[[2]]
##     bCovs <- dimnames(parmList$probB)[[2]]
##     covNames <- c(rCovs, iCovs, cCovs, bCovs)
##     nrcovs <- length(rCovs)
##     nicovs <- length(iCovs)
##     nccovs <- length(cCovs)
##     nbcovs <- length(bCovs)
##     ##
##     Y <- data.matrix(rbind(Y))
##     rY <- colnames(Y)[colnames(Y) %in% rCovs]
##     iY <- colnames(Y)[colnames(Y) %in% iCovs]
##     cY <- colnames(Y)[colnames(Y) %in% cCovs]
##     bY <- colnames(Y)[colnames(Y) %in% bCovs]
##     ##
##     if(!is.null(X)){
##         X <- data.matrix(rbind(X))
##         rX <- colnames(X)[colnames(X) %in% rCovs]
##         if(length(intersect(rX,rY))>0){
##             warning('*WARNING: predictor and predictand have real variates in common. Removing from predictor*')
##             rX <- setdiff(rX,rY)
##         }
##         iX <- colnames(X)[colnames(X) %in% iCovs]
##         if(length(intersect(iX,iY))>0){
##             warning('*WARNING: predictor and predictand have integer variates in common. Removing from predictor*')
##             iX <- setdiff(iX,iY)
##         }
##         cX <- colnames(X)[colnames(X) %in% cCovs]
##         if(length(intersect(cX,cY))>0){
##             warning('*WARNING: predictor and predictand have category variates in common. Removing from predictor*')
##             cX <- setdiff(cX,cY)
##         }
##         bX <- colnames(X)[colnames(X) %in% bCovs]
##         if(length(intersect(bX,bY))>0){
##             warning('*WARNING: predictor and predictand have binary variates in common. Removing from predictor*')
##             bX <- setdiff(bX,bY)
##         }
##         if(length(c(rX,iX,cX,bX))==0){X <- NULL}
##     }
##     ##
##     if(!is.null(X)){
##         if(nrow(X) < nrow(Y)){
##         warning('*Note: X has fewer data than Y. Recycling*')
##         X <- t(matrix(rep(t(X), ceiling(nrow(Y)/nrow(X))), nrow=ncol(X), dimnames=list(colnames(X),NULL)))[1:nrow(Y),,drop=FALSE]
##     }
##     if(nrow(X) > nrow(Y)){
##         warning('*Note: X has more data than Y. Recycling*')
##         Y <- t(matrix(rep(t(Y), ceiling(nrow(X)/nrow(Y))), nrow=ncol(Y), dimnames=list(colnames(Y),NULL)))[1:nrow(X),,drop=FALSE]
##     }
##     }
##     nydata <- nrow(Y)
##     ##
##     q <- parmList$q
##     nclusters <- ncol(q)
##     sclusters <- seq_len(nclusters)
##     if(is.numeric(nfsamples)){
##         fsubsamples <- round(seq(1, nrow(q), length.out=nfsamples))
##     }else{
##         nfsamples <- nrow(q)
##         fsubsamples <- seq_len(nfsamples)
##     }
##     ##
##     if(!is.null(X)){
##         freqs <- foreach(asample=fsubsamples, .combine=c, .inorder=inorder)%dopar%{
##             ## pX: rows=clusters, cols=datapoints
##             pX <- exp(
##                 log(q[asample,]) +
##                 t(rbind(vapply(sclusters, function(acluster){
##                     ## real covariates
##                     (if(length(rX)>0){
##                         colSums(dnorm(x=t(X[,rX,drop=FALSE]), mean=parmList$meanR[asample,rX,acluster], sd=sqrt(parmList$varR[asample,rX,acluster]), log=TRUE), na.rm=TRUE)
##                     }else{0}) +
##                         ## integer covariates
##                         (if(length(iX)>0){
##                             colSums(dbinom(x=t(X[,iX,drop=FALSE]), prob=parmList$probI[asample,iX,acluster], size=parmList$sizeI[asample,iX,acluster], log=TRUE), na.rm=TRUE)
##                         }else{0}) +
##                         ## category variates
##                         (if(length(cX)>0){
##                             rowSums(log(sapply(cX, function(acov){
##                                 parmList$probC[asample,acov,acluster,X[,acov]]})), na.rm=T)
##                             }else{0}) +
##                         ## binary covariates
##                         (if(length(bX)>0){
##                             colSums(log(
##                                 parmList$probB[asample,bX,acluster] * t(X[,bX,drop=FALSE]) +
##                                 (1-parmList$probB[asample,bX,acluster]) * (1-t(X[,bX,drop=FALSE]))
##                             ), na.rm=TRUE)
##                         }else{0})
##                 }, numeric(nydata))))
##             )
##             ## pY: rows=clusters, cols=datapoints
##             pY <- exp(
##                 t(rbind(vapply(sclusters, function(acluster){
##                     ## real covariates
##                     (if(length(rY)>0){
##                         colSums(dnorm(x=t(Y[,rY,drop=FALSE]), mean=parmList$meanR[asample,rY,acluster], sd=sqrt(parmList$varR[asample,rY,acluster]), log=TRUE), na.rm=TRUE)
##                     }else{0}) +
##                         ## integer covariates
##                         (if(length(iY)>0){
##                             colSums(dbinom(x=t(Y[,iY,drop=FALSE]), prob=parmList$probI[asample,iY,acluster], size=parmList$sizeI[asample,iY,acluster], log=TRUE), na.rm=TRUE)
##                         }else{0}) +
##                         ## category variates
##                         (if(length(cY)>0){
##                             rowSums(log(sapply(cY, function(acov){
##                                 parmList$probC[asample,acov,acluster,Y[,acov]]})), na.rm=T)
##                             }else{0}) +
##                         ## binary covariates
##                         (if(length(bY)>0){
##                             colSums(log(
##                                 parmList$probB[asample,bY,acluster] * t(Y[,bY,drop=FALSE]) +
##                                 (1-parmList$probB[asample,bY,acluster]) * (1-t(Y[,bY,drop=FALSE]))
##                             ), na.rm=TRUE)
##                         }else{0})
##                 }, numeric(nydata))))
##             )
##             ##
##             sum(log(colSums(pX * pY)))-sum(log(colSums(pX)))
##         }
##     }else{
##         freqs <- foreach(asample=fsubsamples, .combine=c, .inorder=inorder)%dopar%{
##             ## pY: rows=clusters, cols=datapoints
##             pY <- exp(
##                 log(q[asample,]) +
##                 t(rbind(vapply(sclusters, function(acluster){
##                     ## real covariates
##                     (if(length(rY)>0){
##                         colSums(dnorm(x=t(Y[,rY,drop=FALSE]), mean=parmList$meanR[asample,rY,acluster], sd=sqrt(parmList$varR[asample,rY,acluster]), log=TRUE), na.rm=TRUE)
##                     }else{0}) +
##                         ## integer covariates
##                         (if(length(iY)>0){
##                             colSums(dbinom(x=t(Y[,iY,drop=FALSE]), prob=parmList$probI[asample,iY,acluster], size=parmList$sizeI[asample,iY,acluster], log=TRUE), na.rm=TRUE)
##                         }else{0}) +
##                         ## category variates
##                         (if(length(cY)>0){
##                             rowSums(log(sapply(cY, function(acov){
##                                 parmList$probC[asample,acov,acluster,Y[,acov]]})), na.rm=T)
##                             }else{0}) +
##                         ## binary covariates
##                         (if(length(bY)>0){
##                             colSums(log(
##                                 parmList$probB[asample,bY,acluster] * t(Y[,bY,drop=FALSE]) +
##                                 (1-parmList$probB[asample,bY,acluster]) * (1-t(Y[,bY,drop=FALSE]))
##                             ), na.rm=TRUE)
##                         }else{0})
##                 }, numeric(nydata))))
##             )
##             ##
##             sum(log(colSums(pY)))
##         }
##     }
##     freqs
## }
## ##
## ## Gives samples of means of distributions of any set of Y conditional on any set of X
## samplesMeans <- function(Y, X=NULL, parmList, nfsamples=NULL, inorder=FALSE){
##     rCovs <- dimnames(parmList$meanR)[[2]]
##     iCovs <- dimnames(parmList$probI)[[2]]
##     bCovs <- dimnames(parmList$probB)[[2]]
##     covNames <- c(rCovs, iCovs, bCovs)
##     nrcovs <- length(rCovs)
##     nicovs <- length(iCovs)
##     nbcovs <- length(bCovs)
##     ##
##     if(length(Y)==0 | length(intersect(Y, covNames))==0){stop('invalid variates in Y')}
##     Y <- intersect(Y, covNames)
##     rY <- Y[Y %in% rCovs]
##     iY <- Y[Y %in% iCovs]
##     bY <- Y[Y %in% bCovs]
##     ##
##     if(!is.null(X)){
##         X <- data.matrix(rbind(X))
##         rX <- colnames(X)[colnames(X) %in% rCovs]
##         if(length(intersect(rX,rY))>0){
##             warning('*WARNING: predictor and predictand have real variates in common. Removing from predictor*')
##             rX <- setdiff(rX,rY)
##         }
##         iX <- colnames(X)[colnames(X) %in% iCovs]
##         if(length(intersect(iX,iY))>0){
##             warning('*WARNING: predictor and predictand have integer variates in common. Removing from predictor*')
##             iX <- setdiff(iX,iY)
##         }
##         bX <- colnames(X)[colnames(X) %in% bCovs]
##         if(length(intersect(bX,bY))>0){
##             warning('*WARNING: predictor and predictand have binary variates in common. Removing from predictor*')
##             bX <- setdiff(bX,bY)
##         }
##         if(length(c(rX,iX,bX))==0){X <- NULL}
##     }
##     ##
##     ## if(!is.null(X)){
##     ##     if(nrow(X) < nrow(Y)){
##     ##     warning('*Note: X has fewer data than Y. Recycling*')
##     ##     X <- t(matrix(rep(t(X), ceiling(nrow(Y)/nrow(X))), nrow=ncol(X), dimnames=list(colnames(X),NULL)))[1:nrow(Y),,drop=FALSE]
##     ## }
##     ## if(nrow(X) > nrow(Y)){
##     ##     warning('*Note: X has more data than Y. Recycling*')
##     ##     Y <- t(matrix(rep(t(Y), ceiling(nrow(X)/nrow(Y))), nrow=ncol(Y), dimnames=list(Y,NULL)))[1:nrow(X),,drop=FALSE]
##     ## }
##     ## }
##     ndata <- nrow(X)
##     ##
##     q <- parmList$q
##     nclusters <- ncol(q)
##     sclusters <- seq_len(nclusters)
##     if(is.numeric(nfsamples)){
##         fsubsamples <- round(seq(1, nrow(q), length.out=nfsamples))
##     }else{
##         nfsamples <- nrow(q)
##         fsubsamples <- seq_len(nfsamples)
##     }
##     ##
##     if(!is.null(X)){
##         freqs <- foreach(asample=fsubsamples, .combine=cbind, .inorder=inorder)%dopar%{
##             ## pX: rows=clusters, cols=datapoints
##             pX <- exp(
##                 log(q[asample,]) +
##                 t(rbind(vapply(sclusters, function(acluster){
##                     ## real covariates
##                     (if(length(rX)>0){
##                         colSums(dnorm(x=t(X[,rX,drop=FALSE]), mean=parmList$meanR[asample,rX,acluster], sd=sqrt(parmList$varR[asample,rX,acluster]), log=TRUE))
##                     }else{0}) +
##                         ## integer covariates
##                         (if(length(iX)>0){
##                             colSums(dbinom(x=t(X[,iX,drop=FALSE]), prob=parmList$probI[asample,iX,acluster], size=parmList$sizeI[asample,iX,acluster], log=TRUE))
##                         }else{0}) +
##                         ## binary covariates
##                         (if(length(bX)>0){
##                             colSums(log(
##                                 parmList$probB[asample,bX,acluster] * t(X[,bX,drop=FALSE]) +
##                                 (1-parmList$probB[asample,bX,acluster]) * (1-t(X[,bX,drop=FALSE]))
##                             ))
##                         }else{0})
##                 }, numeric(ndata))))
##             )
##             ## pY: rows=means, cols=clusters
##             pY <- rbind(
##                     ## real covariates
##                 if(length(rY)>0){
##                     parmList$meanR[asample,rY,]
##                     },
##                         ## integer covariates
##                 if(length(iY)>0){
##                     parmList$probI[asample,iY,] * parmList$sizeI[asample,iY,]
##                     },
##                 ## binary covariates
##                 if(length(bY)>0){
##                     parmList$probB[asample,bY,]
##                 }
##                 )
##             ##
##             t(sapply(1:nrow(pY),function(amean){
##                 colSums(pY[amean,] * pX)/colSums(pX)
##             }))
##         }
##     }else{
##         freqs <- foreach(asample=fsubsamples, .combine=cbind, .inorder=inorder)%dopar%{
##             ## pY: rows=means, cols=clusters
##             pY <- rbind(
##                     ## real covariates
##                 if(length(rY)>0){
##                     parmList$meanR[asample,rY,]
##                     },
##                         ## integer covariates
##                 if(length(iY)>0){
##                     parmList$probI[asample,iY,] * parmList$sizeI[asample,iY,]
##                     },
##                 ## binary covariates
##                 if(length(bY)>0){
##                     parmList$probB[asample,bY,]
##                 }
##                 )
##             ##
##             colSums(q[asample,]*t(pY))
##         }
##     }
##     dim(freqs) <- c(length(Y), ndata, nfsamples)
##     dimnames(freqs) <- list(Y,NULL,NULL)
##     freqs
## }
## ##
## ## Gives samples of variances of distributions of any set of Y conditional on any set of X
## samplesVars <- function(Y, X=NULL, parmList, nfsamples=NULL, inorder=FALSE){
##     rCovs <- dimnames(parmList$meanR)[[2]]
##     iCovs <- dimnames(parmList$probI)[[2]]
##     bCovs <- dimnames(parmList$probB)[[2]]
##     covNames <- c(rCovs, iCovs, bCovs)
##     nrcovs <- length(rCovs)
##     nicovs <- length(iCovs)
##     nbcovs <- length(bCovs)
##     ##
##     if(length(Y)==0 | length(intersect(Y, covNames))==0){stop('invalid variates in Y')}
##     Y <- intersect(Y, covNames)
##     rY <- Y[Y %in% rCovs]
##     iY <- Y[Y %in% iCovs]
##     bY <- Y[Y %in% bCovs]
##     ##
##     if(!is.null(X)){
##         X <- data.matrix(rbind(X))
##         rX <- colnames(X)[colnames(X) %in% rCovs]
##         if(length(intersect(rX,rY))>0){
##             warning('*WARNING: predictor and predictand have real variates in common. Removing from predictor*')
##             rX <- setdiff(rX,rY)
##         }
##         iX <- colnames(X)[colnames(X) %in% iCovs]
##         if(length(intersect(iX,iY))>0){
##             warning('*WARNING: predictor and predictand have integer variates in common. Removing from predictor*')
##             iX <- setdiff(iX,iY)
##         }
##         bX <- colnames(X)[colnames(X) %in% bCovs]
##         if(length(intersect(bX,bY))>0){
##             warning('*WARNING: predictor and predictand have binary variates in common. Removing from predictor*')
##             bX <- setdiff(bX,bY)
##         }
##         if(length(c(rX,iX,bX))==0){X <- NULL}
##     }
##     ##
##     ## if(!is.null(X)){
##     ##     if(nrow(X) < nrow(Y)){
##     ##     warning('*Note: X has fewer data than Y. Recycling*')
##     ##     X <- t(matrix(rep(t(X), ceiling(nrow(Y)/nrow(X))), nrow=ncol(X), dimnames=list(colnames(X),NULL)))[1:nrow(Y),,drop=FALSE]
##     ## }
##     ## if(nrow(X) > nrow(Y)){
##     ##     warning('*Note: X has more data than Y. Recycling*')
##     ##     Y <- t(matrix(rep(t(Y), ceiling(nrow(X)/nrow(Y))), nrow=ncol(Y), dimnames=list(Y,NULL)))[1:nrow(X),,drop=FALSE]
##     ## }
##     ## }
##     ndata <- nrow(X)
##     ##
##     q <- parmList$q
##     nclusters <- ncol(q)
##     sclusters <- seq_len(nclusters)
##     if(is.numeric(nfsamples)){
##         fsubsamples <- round(seq(1, nrow(q), length.out=nfsamples))
##     }else{
##         nfsamples <- nrow(q)
##         fsubsamples <- seq_len(nfsamples)
##     }
##     ##
##     if(!is.null(X)){
##         freqs <- foreach(asample=fsubsamples, .combine=cbind, .inorder=inorder)%dopar%{
##             ## pX: rows=clusters, cols=datapoints
##             pX <- exp(
##                 log(q[asample,]) +
##                 t(rbind(vapply(sclusters, function(acluster){
##                     ## real covariates
##                     (if(length(rX)>0){
##                         colSums(dnorm(x=t(X[,rX,drop=FALSE]), mean=parmList$meanR[asample,rX,acluster], sd=sqrt(parmList$varR[asample,rX,acluster]), log=TRUE))
##                     }else{0}) +
##                         ## integer covariates
##                         (if(length(iX)>0){
##                             colSums(dbinom(x=t(X[,iX,drop=FALSE]), prob=parmList$probI[asample,iX,acluster], size=parmList$sizeI[asample,iX,acluster], log=TRUE))
##                         }else{0}) +
##                         ## binary covariates
##                         (if(length(bX)>0){
##                             colSums(log(
##                                 parmList$probB[asample,bX,acluster] * t(X[,bX,drop=FALSE]) +
##                                 (1-parmList$probB[asample,bX,acluster]) * (1-t(X[,bX,drop=FALSE]))
##                             ))
##                         }else{0})
##                 }, numeric(ndata))))
##             )
##             ## pY: rows=means, cols=clusters
##             pY <- rbind(
##                     ## real covariates
##                 if(length(rY)>0){
##                     parmList$varR[asample,rY,]
##                     },
##                         ## integer covariates
##                 if(length(iY)>0){
##                     (1-parmList$probI[asample,iY,]) * parmList$probI[asample,iY,] * parmList$sizeI[asample,iY,]
##                     },
##                 ## binary covariates
##                 if(length(bY)>0){
##                     (1-parmList$probB[asample,bY,]) * parmList$probB[asample,bY,]
##                 }
##                 )
##             ##
##             t(sapply(1:nrow(pY),function(amean){
##                 colSums(pY[amean,] * pX)/colSums(pX)
##             }))
##         }
##     }else{
##         freqs <- foreach(asample=fsubsamples, .combine=cbind, .inorder=inorder)%dopar%{
##             ## pY: rows=means, cols=clusters
##             pY <- rbind(
##                     ## real covariates
##                 if(length(rY)>0){
##                     parmList$varR[asample,rY,]
##                     },
##                         ## integer covariates
##                 if(length(iY)>0){
##                     (1-parmList$probI[asample,iY,]) * parmList$probI[asample,iY,] * parmList$sizeI[asample,iY,]
##                     },
##                 ## binary covariates
##                 if(length(bY)>0){
##                     (1-parmList$probB[asample,bY,]) * parmList$probB[asample,bY,]
##                 }
##                 )
##             ##
##             colSums(q[asample,]*t(pY))
##         }
##     }
##     dim(freqs) <- c(length(Y), ndata, nfsamples)
##     dimnames(freqs) <- list(Y,NULL,NULL)
##     freqs
## }
## ##
## ## Gives samples of variate values. Uses mcsamples
## samplesXmc <- function(mcsamples, variateparameters, Xnames=NULL, pointspermcsample=2, fromsamples=nrow(mcsamples), inorder=FALSE, seed=NULL){
##     allvariates <- rownames(variateparameters)
##     if(!is.null(Xnames)){
##         if(!all(Xnames %in% allvariates)){
##             warning('*WARNING: some requested variates missing from variateparameters argument*')
##             warning(paste0('*Discarding: ',paste0(setdiff(Xnames,allvariates),collapse=' ')))
##             Xnames <- intersect(Xnames,allvariates)
##         }
##     }else{
##         Xnames <- allvariates
##     }
##     ##
##     rCovs <- Xnames[variateparameters[Xnames,'type']==0]
##     iCovs <- Xnames[variateparameters[Xnames,'type']==3]
##     cCovs <- Xnames[variateparameters[Xnames,'type']==1]
##     bCovs <- Xnames[variateparameters[Xnames,'type']==2]
##     cNames <- c(rCovs, iCovs, cCovs, bCovs)
##     nrcovs <- length(rCovs)
##     nicovs <- length(iCovs)
##     nccovs <- length(cCovs)
##     nbcovs <- length(bCovs)
##     ##
##     if(!('location' %in% colnames(variateparameters))){
##         locations <- integer(length(Xnames))
##         names(locations) <- Xnames
##     }else{locations <- variateparameters[Xnames,'location']
##     names(locations) <- Xnames}
##     if(!('scale' %in% colnames(variateparameters))){
##         scales <- locations * 0L + 1L
##     }else{scales <- variateparameters[Xnames,'scale']
##     names(scales) <- Xnames}
##     ##
##     Qi <- grep('q',colnames(mcsamples))
##     nclusters <- length(Qi)
##     sclusters <- seq_len(nclusters)
##     if(nrcovs>0){
##         totake <- variateparameters[rCovs,'index']
##         meanRi <- t(sapply(paste0('meanR\\[',totake,','),grep,colnames(mcsamples)))
##         varRi <- t(sapply(paste0('varR\\[',totake,','),grep,colnames(mcsamples)))
##         dim(meanRi) <- dim(varRi) <- c(nrcovs, nclusters)
##         rownames(meanRi) <- rownames(varRi) <- rCovs
##         }
##     if(nicovs>0){
##         totake <- variateparameters[iCovs,'index']
##         probIi <- t(sapply(paste0('probI\\[',totake,','),grep,colnames(mcsamples)))
##         sizeIi <- t(sapply(paste0('sizeI\\[',totake,','),grep,colnames(mcsamples)))
##         dim(probIi) <- dim(sizeIi) <- c(nicovs,nclusters)
##         rownames(probIi) <- rownames(sizeIi) <- iCovs
##         }
##     if(nccovs>0){
##         totake <- variateparameters[cCovs,'index']
##         probCi <- t(sapply(paste0('probC\\[',totake,','),grep,colnames(mcsamples)))
##         ncategories <- length(probCi)/nccovs/nclusters
##         dim(probCi) <- c(nccovs,nclusters,ncategories)
##         dimnames(probCi) <- list(cCovs, NULL, NULL)
##         scategories <- seq_len(ncategories)
##         }
##     if(nbcovs>0){
##         totake <- variateparameters[bCovs,'index']
##         probBi <- t(sapply(paste0('probB\\[',totake,','),grep,colnames(mcsamples)))
##         dim(probBi) <- c(nbcovs,nclusters)
##         rownames(probBi) <- bCovs
##         sbins <- 0:1
##         }
##     ##
##     if(length(fromsamples) == 1){
##         if(inorder){
##             fromsamples <- round(seq(1, nrow(mcsamples), length.out=fromsamples))
##         }else{
##             fromsamples <- sample(1:nrow(mcsamples),fromsamples,replace=(fromsamples>nrow(mcsamples)))
##         }
##     }
##     ##
##     if(!is.null(seed)){set.seed(seed)}
##     XX <- foreach(asample=t(mcsamples[fromsamples,,drop=F]), .combine=cbind, .inorder=inorder)%dorng%{
##         asample <- asample[,1]
##         acluster <- sample(x=sclusters, size=pointspermcsample, prob=asample[Qi], replace=TRUE)
        
##         ##
##         (if(nrcovs>0){
##              rX <- rnorm(n=pointspermcsample*nrcovs, mean=asample[meanRi[,acluster]], sd=sqrt(asample[varRi[,acluster]]))
##              dim(rX) <- c(nrcovs, pointspermcsample)
##          }else{rX <- NULL})
##         ##        
##         (if(nicovs>0){
##              iX <- rbinom(n=pointspermcsample*nicovs, prob=asample[probIi[,acluster]], size=asample[sizeIi[,acluster]])
##              dim(iX) <- c(nicovs, pointspermcsample)#, dimnames=list(NULL, iCovs))
##          }else{iX <- NULL})
##         ##        
##         (if(nccovs>0){
##              cX <- sapply(acluster, function(clus){sapply(cCovs,function(acov){
##                  sample(x=scategories, size=1, prob=asample[probCi[acov,clus,]])
##              },USE.NAMES=F)})
##          }else{cX <- NULL})
##         ##
##         (if(nbcovs>0){
##         bX <- sapply(acluster, function(clus){sapply(bCovs,function(acov){
##             sample(x=sbins, size=1, prob=c(1-asample[probBi[acov,clus]], asample[probBi[acov,clus]]))
##         },USE.NAMES=F)})
##          }else{bX <- NULL})
##         ##
##         rbind(rX, iX, cX, bX)[order(match(cNames,Xnames)),,drop=F]
##     }
##     attr(XX, 'rng') <- attr(XX, 'doRNG_version') <- NULL
##     dim(XX) <- c(length(Xnames), pointspermcsample, length(fromsamples))
##     rownames(XX) <- Xnames
##     ## rows: mcsamples, cols: samples from one mcsample, 3rd dim: variates
##     aperm(XX*scales+locations)
## }
## ##
## ## Gives samples of variate values
## samplesX <- function(parmListt, nperf=2, nfsamples=NULL, inorder=FALSE, seed=NULL){
##     rCovs <- dimnames(parmList$meanR)[[2]]
##     iCovs <- dimnames(parmList$probI)[[2]]
##     bCovs <- dimnames(parmList$probB)[[2]]
##     covNames <- c(rCovs, iCovs, bCovs)
##     nrcovs <- length(rCovs)
##     nicovs <- length(iCovs)
##     nbcovs <- length(bCovs)
##     ##
##     q <- parmList$q
##     nclusters <- ncol(q)
##     sclusters <- seq_len(nclusters)
##     if(is.numeric(nfsamples)){
##         fsubsamples <- round(seq(1, nrow(q), length.out=nfsamples))
##     }else{
##         nfsamples <- nrow(q)
##         fsubsamples <- seq_len(nfsamples)
##     }
##     ##
##     if(!is.null(seed)){set.seed(seed)}
##     XX <- foreach(asample=fsubsamples, .combine=rbind, .inorder=inorder)%dorng%{
##         acluster <- sample(x=sclusters, size=nperf, prob=q[asample,], replace=TRUE)
##         ##
##         rX <- matrix(rnorm(n=nperf*nrcovs, mean=t(parmList$meanR[asample,,acluster]), sd=sqrt(t(parmList$varR[asample,,acluster]))),
##                      nrow=nperf, dimnames=list(NULL,rCovs))
##         ##        
##         iX <- matrix(rbinom(n=nperf*nicovs, prob=t(parmList$probI[asample,,acluster]), size=t(parmList$sizeI[asample,,acluster])),
##                      nrow=nperf, dimnames=list(NULL, iCovs))
##         ##
##         (if(nbcovs>0){
##         bX <- sapply(bCovs, function(acov){sapply(acluster,function(clus){
##             sample(x=0:1, size=1, prob=c(1-parmList$probB[asample,acov,clus], parmList$probB[asample,acov,clus]))
##         })})}else{bX <- NULL})
##         ##
##         cbind(rX, iX, bX)
##     }
##     attr(XX, 'rng') <- NULL
##     XX
## }
## ## ##
## ## ## Gives samples of long-run mutual information for two sets of variates
## samplesMI2 <- function(Y, X=Y, X2=NULL, parmList, base=2L, nperf=1024, nfsamples=NULL, inorder=FALSE, seed=NULL){
##     rCovs <- dimnames(parmList$meanR)[[2]]
##     iCovs <- dimnames(parmList$probI)[[2]]
##     bCovs <- dimnames(parmList$probB)[[2]]
##     covNames <- c(rCovs, iCovs, bCovs)
##     ##
##     rY <- Y[Y %in% rCovs]
##     iY <- Y[Y %in% iCovs]
##     bY <- Y[Y %in% bCovs]
##     rX <- X[X %in% rCovs]
##     iX <- X[X %in% iCovs]
##     bX <- X[X %in% bCovs]
##     rX2 <- X2[X2 %in% rCovs]
##     iX2 <- X2[X2 %in% iCovs]
##     bX2 <- X2[X2 %in% bCovs]
##     ##
##     rZ <- union(rY,union(rX,rX2))
##     iZ <- union(iY,union(iX,iX2))
##     bZ <- union(bY,union(bX,bX2))
##     nrcovs <- length(rZ)
##     nicovs <- length(iZ)
##     nbcovs <- length(bZ)
##     ##
##     q <- parmList$q
##     nclusters <- ncol(q)
##     sclusters <- seq_len(nclusters)
##     if(is.numeric(nfsamples)){
##         fsubsamples <- round(seq(1, nrow(q), length.out=nfsamples))
##     }else{
##         nfsamples <- nrow(q)
##         fsubsamples <- seq_len(nfsamples)
##     }
##     ##
##     if(!is.null(seed)){set.seed(seed)}
##     MI <- foreach(asample=fsubsamples, .combine=cbind, .inorder=inorder)%dorng%{
##         someclusters <- sample(x=sclusters, size=nperf, prob=q[asample,], replace=TRUE)
##         ##
##         Zs <- cbind(
##             (if(nrcovs>0){matrix(rnorm(n=nperf*nrcovs, mean=t(parmList$meanR[asample,rZ,someclusters]), sd=sqrt(t(parmList$varR[asample,rZ,someclusters]))),
##                      nrow=nperf, dimnames=list(NULL,rZ))}else{NULL}),
##         ##        
##         (if(nicovs>0){matrix(rbinom(n=nperf*nicovs, prob=t(parmList$probI[asample,iZ,someclusters]), size=t(parmList$sizeI[asample,iZ,someclusters])),
##                      nrow=nperf, dimnames=list(NULL, iZ))}else{NULL}),
##         ##
##         (if(nbcovs>0){sapply(bZ, function(acov){sapply(someclusters,function(acluster){
##             sample(x=0:1, size=1, prob=c(1-parmList$probB[asample,acov,acluster], parmList$probB[asample,acov,acluster]))
##         })})}else{NULL})
##         )
##         ##
##         ##
##         logprobs <- aperm(sapply(sclusters, simplify='array', function(acluster){rbind(
##                     ## real covariates
##                     (if(nrcovs>0){
##                         dnorm(x=t(Zs[,rZ,drop=FALSE]), mean=parmList$meanR[asample,rZ,acluster], sd=sqrt(parmList$varR[asample,rZ,acluster]), log=TRUE)
##                     }else{NULL}),
##                         ## integer covariates
##                         (if(nicovs>0){
##                             dbinom(x=t(Zs[,iZ,drop=FALSE]), prob=parmList$probI[asample,iZ,acluster], size=parmList$sizeI[asample,iZ,acluster], log=TRUE)
##                         }else{NULL}),
##                         ## binary covariates
##                         (if(nbcovs>0){
##                             log(
##                                 parmList$probB[asample,bZ,acluster] * t(Zs[,bZ,drop=FALSE]) +
##                                 (1-parmList$probB[asample,bZ,acluster]) * (1-t(Zs[,bZ,drop=FALSE]))
##                             )
##                         }else{NULL})
##         )}), c(1,3,2))
##         ##
##         pY <- colSums(exp(log(q[asample,]) + colSums(logprobs[Y,,,drop=FALSE], na.rm=TRUE)), na.rm=TRUE)
##         pX <- colSums(exp(log(q[asample,]) + colSums(logprobs[X,,,drop=FALSE], na.rm=TRUE)), na.rm=TRUE)
##         pYX <- colSums(exp(log(q[asample,]) + colSums(logprobs[union(Y,X),,,drop=FALSE], na.rm=TRUE)), na.rm=TRUE)
##         if(!is.null(X2)){
##             pX2 <- colSums(exp(log(q[asample,]) + colSums(logprobs[X2,,,drop=FALSE], na.rm=TRUE)), na.rm=TRUE)
##             pYX2 <- colSums(exp(log(q[asample,]) + colSums(logprobs[union(Y,X2),,,drop=FALSE], na.rm=TRUE)), na.rm=TRUE)
##             }
##         ##
##         ##mean(log2(pYX/(pX*pY)), na.rm=TRUE)
##         c(mean(log2(pYX)-log2(pX)-log2(pY), na.rm=TRUE),
##           if(!is.null(X2)){ mean(log2(pYX2)-log2(pX2)-log2(pY), na.rm=TRUE) })
##     }/log2(base)
##     attr(MI, 'rng') <- NULL
##     MI
## }
## ## ##
## ## ## Gives samples of long-run mutual information for several sets of variates
## samplesMI <- function(Y, X, parmList, base=2L, nperf=1024, nfsamples=NULL, inorder=FALSE, seed=NULL){
##     rCovs <- dimnames(parmList$meanR)[[2]]
##     iCovs <- dimnames(parmList$probI)[[2]]
##     bCovs <- dimnames(parmList$probB)[[2]]
##     covNames <- c(rCovs, iCovs, bCovs)
##     ##
##     rY <- Y[Y %in% rCovs]
##     iY <- Y[Y %in% iCovs]
##     bY <- Y[Y %in% bCovs]
##     if(!is.list(X)){X <- list(X)}
##     X2 <- unlist(X)
##     rX2 <- X2[X2 %in% rCovs]
##     iX2 <- X2[X2 %in% iCovs]
##     bX2 <- X2[X2 %in% bCovs]
##     ##
##     rZ <- union(rY,rX2)
##     iZ <- union(iY,iX2)
##     bZ <- union(bY,bX2)
##     nrcovs <- length(rZ)
##     nicovs <- length(iZ)
##     nbcovs <- length(bZ)
##     ##
##     q <- parmList$q
##     nclusters <- ncol(q)
##     sclusters <- seq_len(nclusters)
##     if(is.numeric(nfsamples)){
##         fsubsamples <- round(seq(1, nrow(q), length.out=nfsamples))
##     }else{
##         nfsamples <- nrow(q)
##         fsubsamples <- seq_len(nfsamples)
##     }
##     ##
##     if(!is.null(seed)){set.seed(seed)}
##     MI <- foreach(asample=fsubsamples, .combine=cbind, .inorder=inorder)%dorng%{
##         someclusters <- sample(x=sclusters, size=nperf, prob=q[asample,], replace=TRUE)
##         ##
##         Zs <- cbind(
##             (if(nrcovs>0){matrix(rnorm(n=nperf*nrcovs, mean=t(parmList$meanR[asample,rZ,someclusters]), sd=sqrt(t(parmList$varR[asample,rZ,someclusters]))),
##                      nrow=nperf, dimnames=list(NULL,rZ))}else{NULL}),
##         ##        
##         (if(nicovs>0){matrix(rbinom(n=nperf*nicovs, prob=t(parmList$probI[asample,iZ,someclusters]), size=t(parmList$sizeI[asample,iZ,someclusters])),
##                      nrow=nperf, dimnames=list(NULL, iZ))}else{NULL}),
##         ##
##         (if(nbcovs>0){sapply(bZ, function(acov){sapply(someclusters,function(acluster){
##             sample(x=0:1, size=1, prob=c(1-parmList$probB[asample,acov,acluster], parmList$probB[asample,acov,acluster]))
##         })})}else{NULL})
##         )
##         ##
##         ##
##         logprobs <- aperm(sapply(sclusters, simplify='array', function(acluster){rbind(
##                     ## real covariates
##                     (if(nrcovs>0){
##                         dnorm(x=t(Zs[,rZ,drop=FALSE]), mean=parmList$meanR[asample,rZ,acluster], sd=sqrt(parmList$varR[asample,rZ,acluster]), log=TRUE)
##                     }else{NULL}),
##                         ## integer covariates
##                         (if(nicovs>0){
##                             dbinom(x=t(Zs[,iZ,drop=FALSE]), prob=parmList$probI[asample,iZ,acluster], size=parmList$sizeI[asample,iZ,acluster], log=TRUE)
##                         }else{NULL}),
##                         ## binary covariates
##                         (if(nbcovs>0){
##                             log(
##                                 parmList$probB[asample,bZ,acluster] * t(Zs[,bZ,drop=FALSE]) +
##                                 (1-parmList$probB[asample,bZ,acluster]) * (1-t(Zs[,bZ,drop=FALSE]))
##                             )
##                         }else{NULL})
##         )}), c(1,3,2))
##         ##
##         pY <- colSums(exp(log(q[asample,]) + colSums(logprobs[Y,,,drop=FALSE], na.rm=TRUE)), na.rm=TRUE)
##         ##
##         sapply(X, function(anX){
##             pX <- colSums(exp(log(q[asample,]) + colSums(logprobs[anX,,,drop=FALSE], na.rm=TRUE)), na.rm=TRUE)
##             pYX <- colSums(exp(log(q[asample,]) + colSums(logprobs[union(Y,anX),,,drop=FALSE], na.rm=TRUE)), na.rm=TRUE)
##             mean(log2(pYX)-log2(pX), na.rm=TRUE)
##         }) - mean(log2(pY), na.rm=TRUE)
##     }/log2(base)
##     attr(MI, 'rng') <- NULL
##     ## rownames(MI) <- paste0(Y,'|',sapply(X,function(anX){paste0(anX,collapse = ',')}))
##     rownames(MI) <- names(X)
##     MI
## }

## ## samplesMI2 <- function(Y, X=Y, X2=NULL, parmList, base=2L, nperf=1024, nfsamples=NULL, inorder=FALSE){
## ##     rCovs <- dimnames(parmList$meanR)[[2]]
## ##     iCovs <- dimnames(parmList$probI)[[2]]
## ##     bCovs <- dimnames(parmList$probB)[[2]]
## ##     covNames <- c(rCovs, iCovs, bCovs)
## ##     ##
## ##     rY <- Y[Y %in% rCovs]
## ##     iY <- Y[Y %in% iCovs]
## ##     bY <- Y[Y %in% bCovs]
## ##     rX <- X[X %in% rCovs]
## ##     iX <- X[X %in% iCovs]
## ##     bX <- X[X %in% bCovs]
## ##     rX2 <- X2[X2 %in% rCovs]
## ##     iX2 <- X2[X2 %in% iCovs]
## ##     bX2 <- X2[X2 %in% bCovs]
## ##     ##
## ##     ## rX <- X[X %in% rCovs]
## ##     ##     if(length(intersect(rX,rY))>0){
## ##     ##         warning('*WARNING: predictor and predictand have real variates in common. Removing from predictor*')
## ##     ##         rX <- setdiff(rX,rY)
## ##     ##     }
## ##     ##     iX <- X[X %in% iCovs]
## ##     ##     if(length(intersect(iX,iY))>0){
## ##     ##         warning('*WARNING: predictor and predictand have integer variates in common. Removing from predictor*')
## ##     ##         iX <- setdiff(iX,iY)
## ##     ##     }
## ##     ##     bX <- X[X %in% bCovs]
## ##     ##     if(length(intersect(bX,bY))>0){
## ##     ##         warning('*WARNING: predictor and predictand have binary variates in common. Removing from predictor*')
## ##     ##         bX <- setdiff(bX,bY)
## ##     ##     }
## ##     ##     if(length(c(rX,iX,bX))==0){stop('X is Null')}
## ##     ##
## ##     rZ <- union(rY,union(rX,rX2))
## ##     iZ <- union(iY,union(iX,iX2))
## ##     bZ <- union(bY,union(bX,bX2))
## ##     nrcovs <- length(rZ)
## ##     nicovs <- length(iZ)
## ##     nbcovs <- length(bZ)
## ##     ##
## ##     q <- parmList$q
## ##     nclusters <- ncol(q)
## ##     sclusters <- seq_len(nclusters)
## ##     if(is.numeric(nfsamples)){
## ##         fsubsamples <- round(seq(1, nrow(q), length.out=nfsamples))
## ##     }else{
## ##         nfsamples <- nrow(q)
## ##         fsubsamples <- seq_len(nfsamples)
## ##     }
## ##     ##
## ##     MI <- foreach(asample=fsubsamples, .combine=cbind, .inorder=inorder)%dorng%{
## ##         someclusters <- sample(x=sclusters, size=nperf, prob=q[asample,], replace=TRUE)
## ##         ##
## ##         Zs <- cbind(
## ##             (if(nrcovs>0){matrix(rnorm(n=nperf*nrcovs, mean=t(parmList$meanR[asample,rZ,someclusters]), sd=1/sqrt(t(parmList$tauR[asample,rZ,someclusters]))),
## ##                      nrow=nperf, dimnames=list(NULL,rZ))}else{NULL}),
## ##         ##        
## ##         (if(nicovs>0){matrix(rbinom(n=nperf*nicovs, prob=t(parmList$probI[asample,iZ,someclusters]), size=t(parmList$sizeI[asample,iZ,someclusters])),
## ##                      nrow=nperf, dimnames=list(NULL, iZ))}else{NULL}),
## ##         ##
## ##         (if(nbcovs>0){sapply(bZ, function(acov){sapply(someclusters,function(acluster){
## ##             sample(x=0:1, size=1, prob=c(1-parmList$probB[asample,acov,acluster], parmList$probB[asample,acov,acluster]))
## ##         })})}else{NULL})
## ##         )
## ##         ##
## ##         ##
## ##         logprobs <- aperm(sapply(sclusters, simplify='array', function(acluster){rbind(
## ##                     ## real covariates
## ##                     (if(nrcovs>0){
## ##                         dnorm(x=t(Zs[,rZ,drop=FALSE]), mean=parmList$meanR[asample,rZ,acluster], sd=1/sqrt(parmList$tauR[asample,rZ,acluster]), log=TRUE)
## ##                     }else{NULL}),
## ##                         ## integer covariates
## ##                         (if(nicovs>0){
## ##                             dbinom(x=t(Zs[,iZ,drop=FALSE]), prob=parmList$probI[asample,iZ,acluster], size=parmList$sizeI[asample,iZ,acluster], log=TRUE)
## ##                         }else{NULL}),
## ##                         ## binary covariates
## ##                         (if(nbcovs>0){
## ##                              log(
## ##                              (1-t(Zs[,bZ,drop=FALSE])) +
## ##                              (2*t(Zs[,bZ,drop=FALSE])-1) *
## ##                              parmList$probB[asample,bZ,acluster]
## ##                              )
## ##                         }else{NULL})
## ##         )}), c(1,3,2))
## ##         ##
## ##         pY <- colSums(exp(log(q[asample,]) + colSums(logprobs[Y,,,drop=FALSE], na.rm=TRUE)), na.rm=TRUE)
## ##         pX <- colSums(exp(log(q[asample,]) + colSums(logprobs[X,,,drop=FALSE], na.rm=TRUE)), na.rm=TRUE)
## ##         pYX <- colSums(exp(log(q[asample,]) + colSums(logprobs[union(Y,X),,,drop=FALSE], na.rm=TRUE)), na.rm=TRUE)
## ##         if(!is.null(X2)){
## ##             pX2 <- colSums(exp(log(q[asample,]) + colSums(logprobs[X2,,,drop=FALSE], na.rm=TRUE)), na.rm=TRUE)
## ##             pYX2 <- colSums(exp(log(q[asample,]) + colSums(logprobs[union(Y,X2),,,drop=FALSE], na.rm=TRUE)), na.rm=TRUE)
## ##             }
## ##         ##
## ##         ##mean(log2(pYX/(pX*pY)), na.rm=TRUE)
## ##         c(mean(log2(pYX)-log2(pX)-log2(pY), na.rm=TRUE),
## ##           if(!is.null(X2)){ mean(log2(pYX2)-log2(pX2)-log2(pY), na.rm=TRUE) })
## ##     }/log2(base)
## ##     attr(MI, 'rng') <- NULL
## ##     MI
## ## }

## ## ## Gives joint samples of Y or conditional samples of Y given X
## ## samplesYX <- function(Y, X=NULL, parmList, nperf=2, nfsamples=NULL, inorder=FALSE){
## ##     rCovs <- dimnames(parmList$meanR)[[2]]
## ##     iCovs <- dimnames(parmList$probI)[[2]]
## ##     bCovs <- dimnames(parmList$probB)[[2]]
## ##     covNames <- c(rCovs, iCovs, bCovs)
## ##     nrcovs <- length(rCovs)
## ##     nicovs <- length(iCovs)
## ##     nbcovs <- length(bCovs)
## ##     ##
## ##     Y <- data.matrix(rbind(Y))
## ##     rY <- colnames(Y)[colnames(Y) %in% rCovs]
## ##     iY <- colnames(Y)[colnames(Y) %in% iCovs]
## ##     bY <- colnames(Y)[colnames(Y) %in% bCovs]
## ##     ##
## ##     if(!is.null(X)){
## ##         X <- data.matrix(rbind(X))
## ##         rX <- colnames(X)[colnames(X) %in% rCovs]
## ##         if(length(intersect(rX,rY))>0){
## ##             warning('*WARNING: predictor and predictand have real variates in common. Removing from predictor*')
## ##             rX <- setdiff(rX,rY)
## ##         }
## ##         iX <- colnames(X)[colnames(X) %in% iCovs]
## ##         if(length(intersect(iX,iY))>0){
## ##             warning('*WARNING: predictor and predictand have integer variates in common. Removing from predictor*')
## ##             iX <- setdiff(iX,iY)
## ##         }
## ##         bX <- colnames(X)[colnames(X) %in% bCovs]
## ##         if(length(intersect(bX,bY))>0){
## ##             warning('*WARNING: predictor and predictand have binary variates in common. Removing from predictor*')
## ##             bX <- setdiff(bX,bY)
## ##         }
## ##         if(length(c(rX,iX,bX))==0){X <- NULL}
## ##     }
## ##     ##
## ##     if(!is.null(X)){
## ##         if(nrow(X) < nrow(Y)){
## ##         warning('*Note: X has fewer data than Y. Recycling*')
## ##         X <- t(matrix(rep(t(X), ceiling(nrow(Y)/nrow(X))), nrow=ncol(X), dimnames=list(colnames(X),NULL)))[1:nrow(Y),,drop=FALSE]
## ##     }
## ##     if(nrow(X) > nrow(Y)){
## ##         warning('*Note: X has more data than Y. Recycling*')
## ##         Y <- t(matrix(rep(t(Y), ceiling(nrow(X)/nrow(Y))), nrow=ncol(Y), dimnames=list(colnames(Y),NULL)))[1:nrow(X),,drop=FALSE]
## ##     }
## ##     }
## ##     nydata <- nrow(Y)
## ##     ##
## ##     q <- parmList$q
## ##     nclusters <- ncol(q)
## ##     if(is.numeric(nfsamples)){
## ##         fsubsamples <- round(seq(1, nrow(q), length.out=nfsamples))
## ##     }else{
## ##         nfsamples <- nrow(q)
## ##         fsubsamples <- seq_len(nfsamples)
## ##     }
## ##     ##
## ##     if(!is.null(X)){
## ##         freqs <- foreach(asample=fsubsamples, .combine=cbind, .inorder=inorder)%dopar%{
## ##             ## pX: rows=clusters, cols=datapoints
## ##             pX <- exp(
## ##                 log(q[asample,]) +
## ##                 t(rbind(vapply(seq_len(nclusters), function(acluster){
## ##                     ## real covariates
## ##                     (if(length(rX)>0){
## ##                         colSums(dnorm(x=t(X[,rX,drop=FALSE]), mean=parmList$meanR[asample,rX,acluster], sd=1/sqrt(parmList$tauR[asample,rX,acluster]), log=TRUE), na.rm=TRUE)
## ##                     }else{0}) +
## ##                         ## integer covariates
## ##                         (if(length(iX)>0){
## ##                             colSums(dbinom(x=t(X[,iX,drop=FALSE]), prob=parmList$probI[asample,iX,acluster], size=parmList$sizeI[asample,iX,acluster], log=TRUE), na.rm=TRUE)
## ##                         }else{0}) +
## ##                         ## binary covariates
## ##                         (if(length(bX)>0){
## ##                             colSums(log(
## ##                                 parmList$probB[asample,bX,acluster] * t(X[,bX,drop=FALSE]) +
## ##                                 (1-parmList$probB[asample,bX,acluster]) * (1-t(X[,bX,drop=FALSE]))
## ##                             ), na.rm=TRUE)
## ##                         }else{0})
## ##                 }, numeric(nydata))))
## ##             )
## ##             ## pY: rows=clusters, cols=datapoints
## ##             pY <- exp(
## ##                 t(rbind(vapply(seq_len(nclusters), function(acluster){
## ##                     ## real covariates
## ##                     (if(length(rY)>0){
## ##                         colSums(dnorm(x=t(Y[,rY,drop=FALSE]), mean=parmList$meanR[asample,rY,acluster], sd=1/sqrt(parmList$tauR[asample,rY,acluster]), log=TRUE), na.rm=TRUE)
## ##                     }else{0}) +
## ##                         ## integer covariates
## ##                         (if(length(iY)>0){
## ##                             colSums(dbinom(x=t(Y[,iY,drop=FALSE]), prob=parmList$probI[asample,iY,acluster], size=parmList$sizeI[asample,iY,acluster], log=TRUE), na.rm=TRUE)
## ##                         }else{0}) +
## ##                         ## binary covariates
## ##                         (if(length(bY)>0){
## ##                             colSums(log(
## ##                                 parmList$probB[asample,bY,acluster] * t(Y[,bY,drop=FALSE]) +
## ##                                 (1-parmList$probB[asample,bY,acluster]) * (1-t(Y[,bY,drop=FALSE]))
## ##                             ), na.rm=TRUE)
## ##                         }else{0})
## ##                 }, numeric(nydata))))
## ##             )
## ##             ##
## ##             colSums(pX * pY)/colSums(pX)
## ##         }
## ##     }else{
## ##         freqs <- foreach(asample=fsubsamples, .combine=cbind, .inorder=inorder)%dopar%{
## ##             ## pY: rows=clusters, cols=datapoints
## ##             pY <- exp(
## ##                 log(q[asample,]) +
## ##                 t(rbind(vapply(seq_len(nclusters), function(acluster){
## ##                     ## real covariates
## ##                     (if(length(rY)>0){
## ##                         colSums(dnorm(x=t(Y[,rY,drop=FALSE]), mean=parmList$meanR[asample,rY,acluster], sd=1/sqrt(parmList$tauR[asample,rY,acluster]), log=TRUE), na.rm=TRUE)
## ##                     }else{0}) +
## ##                         ## integer covariates
## ##                         (if(length(iY)>0){
## ##                             colSums(dbinom(x=t(Y[,iY,drop=FALSE]), prob=parmList$probI[asample,iY,acluster], size=parmList$sizeI[asample,iY,acluster], log=TRUE), na.rm=TRUE)
## ##                         }else{0}) +
## ##                         ## binary covariates
## ##                         (if(length(bY)>0){
## ##                             colSums(log(
## ##                                 parmList$probB[asample,bY,acluster] * t(Y[,bY,drop=FALSE]) +
## ##                                 (1-parmList$probB[asample,bY,acluster]) * (1-t(Y[,bY,drop=FALSE]))
## ##                             ), na.rm=TRUE)
## ##                         }else{0})
## ##                 }, numeric(nydata))))
## ##             )
## ##             ##
## ##             colSums(pY)
## ##         }
## ##     }
## ##     freqs
## ## }
## ##
## ## Function to draw 2D plot of two variates
## plot2dF <- function(xygrid, fsamples, grid=FALSE, labs=TRUE, ticks=TRUE, mar=NULL){
##     xcov <- colnames(xygrid)[1]
##     ycov <- colnames(xygrid)[2]
##     xgrid <- sort(unique(xygrid[,1]))
##     ygrid <- sort(unique(xygrid[,2]))
##     ##
##     ax <- diff(xgrid)[1]/2
##     ay <- diff(ygrid)[1]/2
##     if(ticks){
##     xticks <- if(xcov %in% realCovs){NULL}else{xgrid}
##     yticks <- if(ycov %in% realCovs){NULL}else{ygrid}
##     }else{
##         xticks <- yticks <- FALSE
##     }
##     xlim <- if(xcov %in% realCovs){extendrange(xgrid)}else{range(xgrid)+c(-1,1)/2}
##     ylim <- if(ycov %in% realCovs){extendrange(ygrid)}else{range(ygrid)+c(-1,1)/2}
##     pmax <- max(fsamples)
##     ##
##     tplot(x=NA, y=NA, xlim=xlim, ylim=ylim, xlab=(if(labs){xcov}else{NA}), ylab=(if(labs){ycov}else{NA}), xticks=xticks, yticks=yticks, mar=mar)
##     for(i in 1:nrow(xygrid)){
##         rat <- fsamples[i]/pmax
##         polygon(x=xygrid[i,1]+c(-1,1,1,-1)*ax,
##                 y=xygrid[i,2]+c(-1,-1,1,1)*ay,
##                 border=gray(1-rat), col=gray(1-rat))
##     }
##     if(grid){
##         tplot(x=NA, y=NA, xlim=xlim, ylim=ylim, xlab=xcov, ylab=ycov, xticks=xticks, yticks=yticks, mar=mar, add=T)
##     }
## }


## ####################################################
## #### Functions below are not used at the moment ####
## ####################################################
## ##
## ## Calculates the median and IQR of each sample frequency
## calcSampleMQ <- function(parmList, maxD=1000){
##     realCovs <- dimnames(parmList$meanR)[[2]]
##     integerCovs <- dimnames(parmList$probI)[[2]]
##     covNames <- c(realCovs, integerCovs)
##     ncovs <- length(covNames)
##     q <- parmList$q
##     nsamples <- nrow(q)
##     ##
##     quants <- foreach(asample=seq_len(nsamples), .combine=c)%:%foreach(acov=covNames, .combine=c)%dopar%{
##         if(acov %in% realCovs){
##                 mixq <- function(x){sum(q[asample,] * pnorm(x, mean=parmList$meanR[asample,acov,], sd=sqrt(parmList$varR[asample,acov,])))}
##                 fn <- function(par){
##                     (mixq(par[1]) - 0.5)^2 +
##                         (mixq(par[2]) - 0.25)^2 +
##                         (mixq(par[3]) - 0.75)^2
##                 }
##                 out <- myoptim(par=rep(q[asample,]%*%parmList$meanR[asample,acov,],3), fn=fn)$par
##             ##     out <- sapply(c(0.005, 0.995), function(border){
##             ## optim(0, #rep(q[asample,] %*% parmList$meanR[asample,acov,], 2),
##             ##              fn=fn,
##             ##              gr = function(x) pracma::grad(fn, x), 
##             ##              method = "L-BFGS-B",
##             ##              lower = -Inf, upper = Inf,
##             ##              control = list(factr = 1e-10, pgtol = 0, maxit = 100))
##         }else{
##             searchgrid <- 0:maxD
##             dq <- colSums(c(q[asample,]) * pbinom(matrix(searchgrid, ncol=length(searchgrid), nrow=ncol(q), byrow=TRUE), prob=parmList$probI[asample,acov,], size=parmList$sizeI[asample,acov,]))
##             out <- c(
##                 which.min(abs(dq - 0.5))-1,
##                 which.min(abs(dq - 0.25))-1,
##                 which.min(abs(dq - 0.75))-1
##             )
##         }
##         out
##     }
##     ##
##     dim(quants) <- c(3, ncovs, nsamples)
##     quants <- aperm(quants)
##     dimnames(quants) <- list(NULL, covNames, c('50%', '25%', '75%'))
##     quants
## }
## ##
## ## Calculates the probability of some datapoints for the MCMC samples
## probValuesSamples <- function(X, parmList){
##     realCovs <- dimnames(parmList$meanR)[[2]]
##     integerCovs <- dimnames(parmList$probI)[[2]]
##     covNames <- c(realCovs, integerCovs)
##     ncovs <- length(covNames)
##     q <- parmList$q
##     ndataz <- nrow(X)
##     ##
##     (foreach(asample=seq_len(nrow(q)), .combine=cbind, .inorder=TRUE)%dopar%{
##         colSums(
##             exp(
##                 log(q[asample,]) +
##                 t(vapply(seq_len(ncol(q)), function(acluster){
##                     ## real covariates
##                     colSums(dnorm(t(X[,realCovs]), mean=parmList$meanR[asample,,acluster,drop=FALSE], sd=sqrt(parmList$varR[asample,,acluster,drop=FALSE]), log=TRUE)) +
##                         ## integer covariates
##                     colSums(dbinom(t(X[,integerCovs]), prob=parmList$probI[asample,,acluster,drop=FALSE], size=parmList$sizeI[asample,,acluster,drop=FALSE], log=TRUE))
##     }, numeric(ndataz)))
##             )
##         )
##     })
## }
## ##
## ## Improved optimization functions
## myoptim <- function(par, fn){
##     resu0 <- list(par=par)
##     resu <- optim(par=resu0$par, fn=fn, control=list(factr = 1e-10, maxit=10000))
##     while(any(resu$par!=resu0$par)){
##         resu0 <- resu
##         resu <- optim(par=resu0$par, fn=fn, control=list(factr = 1e-10, maxit=10000))
##     }
##     resu}
## myoptimbounds <- function(par, fn, lower, upper, maxit=100){
##     optim(par, fn=fn,
##           gr = function(x) pracma::grad(fn, x), 
##           method = "L-BFGS-B",
##           lower = lower, upper = upper,
##           control = list(factr = 1e-10, maxit = maxit))
## }
## ##
## ## Calculates the 0.5% and 99.5% quantiles of each sample frequency
## calcSampleQuantiles <- function(parmList){
##     realCovs <- dimnames(parmList$meanR)[[2]]
##     integerCovs <- dimnames(parmList$probI)[[2]]
##     covNames <- c(realCovs, integerCovs)
##     ncovs <- length(covNames)
##     q <- parmList$q
##     nsamples <- nrow(q)
##     ##
##     quants <- foreach(asample=seq_len(nsamples), .combine=c)%:%foreach(acov=covNames, .combine=c)%dopar%{
##         if(acov %in% realCovs){
##                 mixq <- function(x){sum(q[asample,] * pnorm(x, mean=parmList$meanR[asample,acov,], sd=sqrt(parmList$varR[asample,acov,])))}
##                 fn <- function(par){(mixq(par[1]) - 0.005)^2 +
##                                         (mixq(par[2]) - 0.995)^2}
##                 out <- myoptim(par=rep(q[asample,]%*%parmList$meanR[asample,acov,],2), fn=fn)$par
##             ##     out <- sapply(c(0.005, 0.995), function(border){
##             ## optim(0, #rep(q[asample,] %*% parmList$meanR[asample,acov,], 2),
##             ##              fn=fn,
##             ##              gr = function(x) pracma::grad(fn, x), 
##             ##              method = "L-BFGS-B",
##             ##              lower = -Inf, upper = Inf,
##             ##              control = list(factr = 1e-10, pgtol = 0, maxit = 100))
##         }else{
##             searchgrid <- 0:max(parmList$sizeI[asample,acov,])
##             dq <- colSums(c(q[asample,]) * pbinom(matrix(searchgrid, ncol=length(searchgrid), nrow=ncol(q), byrow=TRUE), prob=parmList$probI[asample,acov,], size=parmList$sizeI[asample,acov,]))
##             out <- c(which.min(abs(dq - 0.005))-1, which.min(abs(dq - 0.995))-1)
##         }
##         out
##     }
##     ##
##     dim(quants) <- c(2, ncovs, nsamples)
##     quants <- aperm(quants)
##     dimnames(quants) <- list(NULL, covNames, c('0.5%', '99.5%'))
##     quants
## }
## ##
## ## Calculates the probability of the data (likelihood of parameters) for one MCMC samples
## ## if(!exists('calcLL')){
## ## calcLL <- nimbleFunction( run=function(
## ##     X=double(2), Y=double(2), Q=double(1),
## ##     MeanC=double(2), TauC=double(2),
## ##     ProbD=double(2), SizeD=double(2)
## ##     ){
## ##     returnType(double(0))
## ##     Nclusters <- length(Q)
## ##     Ndata <- dim(X)[1]
## ##     Nrcovs <- dim(X)[2]
## ##     Nicovs <- dim(Y)[2]
## ##     LL <- 0
## ##     for(adatum in 1:Ndata){
## ##         clustersum <- log(Q)
## ##         for(acov in 1:Nrcovs){
## ##             clustersum <- clustersum +
## ##                 dnorm(x=X[adatum,acov], mean=MeanC[acov,], sd=1/sqrt(TauC[acov,]), log=TRUE)
## ##         }
## ##         for(acov in 1:Nicovs){
## ##             clustersum <- clustersum +
## ##                 dbinom(x=Y[adatum,acov], prob=ProbD[acov,], size=SizeD[acov,], log=TRUE)
## ##         }
## ##         LL <- LL + log(sum(exp(clustersum)))
## ##     }
## ##     return(LL)
## ## } )
## ## CcalcLL <- compileNimble(calcLL)
## ## assign('CcalcLL', CcalcLL, envir = .GlobalEnv)
## ## assign('calcLL', calcLL, envir = .GlobalEnv)
## ## }
## ##
## ## Calculates the probability of several datapoints for several MCMC samples
## probJointSamples <- function(dat, parmList, log=FALSE, inorder=FALSE){
##     ndataz <- nrow(dat$Real)
##     q <- parmList$q
##     ##
##     freqs <- foreach(asample=seq_len(nrow(q)), .combine=cbind, .inorder=inorder)%dopar%{
##         colSums(
##             exp(
##                 log(q[asample,]) +
##                 t(vapply(seq_len(ncol(q)), function(acluster){
##                     ## real covariates
##                     colSums(dnorm(t(dat$Real), mean=parmList$meanR[asample,,acluster], sd=sqrt(parmList$varR[asample,,acluster]), log=TRUE)) +
##                         ## integer covariates
##                     colSums(dbinom(t(dat$Integer), prob=parmList$probI[asample,,acluster], size=parmList$sizeI[asample,,acluster], log=TRUE))
##     }, numeric(ndataz)))
##             )
##         )
##     }
##     if(!log){freqs} else {log(freqs)}
## }
## ## Calculates the MCMC posterior probability of several datapoints
## probJointMean <- function(X, parmList){
##     realCovs <- dimnames(parmList$meanR)[[2]]
##     integerCovs <- dimnames(parmList$probI)[[2]]
##     if(is.list(X)){ X <- cbind(X$Real, X$Integer) }
##     ndataz <- nrow(X)
##     q <- parmList$q
##     nsamples <- nrow(q)
##     ##
##     foreach(asample=seq_len(nsamples), .combine='+', .inorder=FALSE)%dopar%{
##         colSums(
##             exp(
##                 log(q[asample,]) +
##                 t(vapply(seq_len(ncol(q)), function(acluster){
##                     ## real covariates
##                     colSums(dnorm(t(X[,realCovs]), mean=parmList$meanR[asample,,acluster], sd=sqrt(parmList$varR[asample,,acluster]), log=TRUE)) +
##                         ## integer covariates
##                     colSums(dbinom(t(X[,integerCovs]), prob=parmList$probI[asample,,acluster], size=parmList$sizeI[asample,,acluster], log=TRUE))
##     }, numeric(ndataz)))
##             )
##         )
##     }/nsamples
## }
## ##
## ## Produces multidimensional samples from several MCMC distribution-samples
## options(doFuture.rng.onMisuse = "ignore")
## samplesFsamples <- function(parmList, nxsamples=1000, nfsamples=NULL, seed=149){
##     realCovs <- dimnames(parmList$meanR)[[2]]
##     integerCovs <- dimnames(parmList$probI)[[2]]
##     ncC <- length(realCovs)
##     ndC <- length(integerCovs)
##     q <- parmList$q
##         if(is.numeric(nfsamples)){
##         fsubsamples <- seq(1, nrow(q), length.out=nfsamples)
##     }else{
##         nfsamples <- nrow(q)
##         fsubsamples <- seq_len(nfsamples)
##     }
##     ##
##     rng <- RNGseq( nfsamples * nxsamples, seed)
##     allsamples <- foreach(afsample=fsubsamples, .combine=cbind, .inorder=FALSE)%:%foreach(axsample=seq_len(nxsamples), r=rng[(afsample-1)*nxsamples + 1:nxsamples], .combine=c, .inorder=FALSE)%dopar%{
##         rngtools::setRNG(r)
##         acluster <- rcat(n=1, prob=q[afsample,])
##         c(
##             ## real covariates
##             rnorm(n=ncC, mean=parmList$meanR[afsample,realCovs,acluster], sd=sqrt(parmList$varR[afsample,realCovs,acluster])),
##             ## integer covariates
##             rbinom(n=ndC, prob=parmList$probI[afsample,integerCovs,acluster], size=parmList$sizeI[afsample,integerCovs,acluster])
##         )
##     }
##     dim(allsamples) <- c(ncC+ndC, nxsamples, nfsamples)
##     dimnames(allsamples) <- list(c(realCovs,integerCovs), NULL, NULL)
##     allsamples
## }
## ##
## ## Calculates the predictive probability for several log_RMSD values conditional on several feature values
## expeRgivenX <- function(maincov, X, parmList){
##     realCovs <- dimnames(parmList$meanR)[[2]]
##     integerCovs <- dimnames(parmList$probI)[[2]]
##     cC <- setdiff(realCovs, maincov)
##     dC <- setdiff(integerCovs, maincov)
##     ndataz <- nrow(X)
##     q <- parmList$q
##     nsamples <- nrow(q)
##     ##
##     foreach(asample=seq_len(nsamples), .combine='+', .inorder=FALSE)%dopar%{
##         ## W: rows=clusters, cols=datapoints
##         W <- exp(
##             log(q[asample,]) +
##             t(vapply(seq_len(ncol(q)), function(acluster){
##                 ## real covariates
##                 if(length(cC)>0){
##                     colSums(dnorm(t(X[,cC]), mean=parmList$meanR[asample,cC,acluster], sd=sqrt(parmList$varR[asample,cC,acluster]), log=TRUE))
##                 }else{0} +
##                     ## integer covariates
##                     if(length(dC)>0){
##                         colSums(dbinom(t(X[,dC]), prob=parmList$probI[asample,dC,acluster], size=parmList$sizeI[asample,dC,acluster], log=TRUE))
##                     }else{0}
##             }, numeric(ndataz)))
##         )
##         ##
##         colSums(
##                 if(maincov %in% realCovs){
##                     parmList$meanR[asample,maincov,]
##                 }else{
##                     parmList$probI[asample,maincov,] * parmList$sizeI[asample,maincov,] 
##                 } * W)/colSums(W)
##     }/nsamples
## }
## ##
## ## Gives samples of mean of one covariate conditional on several feature values
## samplesmeanRgivenX <- function(maincov, X, parmList, nfsamples=NULL, inorder=FALSE){
##     realCovs <- dimnames(parmList$meanR)[[2]]
##     integerCovs <- dimnames(parmList$probI)[[2]]
##     cC <- setdiff(realCovs, maincov)
##     dC <- setdiff(integerCovs, maincov)
##     ndataz <- nrow(X)
##     q <- parmList$q
##     nclusters <- ncol(q)
##     if(is.numeric(nfsamples)){
##         fsubsamples <- seq(1, nrow(q), length.out=nfsamples)
##     }else{
##         nfsamples <- nrow(q)
##         fsubsamples <- seq_len(nfsamples)
##     }
##     ##
##     freqs <- foreach(asample=fsubsamples, .combine=c, .inorder=inorder)%dopar%{
##         ## W: rows=clusters, cols=datapoints
##         W <- exp(
##             log(q[asample,]) +
##             t(vapply(seq_len(nclusters), function(acluster){
##                 ## real covariates
##                 if(length(cC)>0){
##                     colSums(dnorm(t(X[,cC]), mean=parmList$meanR[asample,cC,acluster], sd=sqrt(parmList$varR[asample,cC,acluster]), log=TRUE))
##                 }else{0} +
##                     ## integer covariates
##                     if(length(dC)>0){
##                         colSums(dbinom(t(X[,dC]), prob=parmList$probI[asample,dC,acluster], size=parmList$sizeI[asample,dC,acluster], log=TRUE))
##                     }else{0}
##             }, numeric(ndataz)))
##         )
##         ## pR: rows=clusters, cols= P at grid points
##         ##
##                 ( (if(maincov %in% realCovs){
##                     parmList$meanR[asample,maincov,]
##                 }else{
##                     parmList$probI[asample,maincov,] * parmList$sizeI[asample,maincov,] 
##                 }) %*% W)/colSums(W)
##     }
##     dim(freqs) <- c(ndataz, nfsamples)
##     freqs
## }
## ##
## ## Gives samples of variance of one covariate conditional on several feature values
## samplesvarRgivenX <- function(maincov, X, parmList, nfsamples=NULL, inorder=FALSE){
##     realCovs <- dimnames(parmList$meanR)[[2]]
##     integerCovs <- dimnames(parmList$probI)[[2]]
##     cC <- setdiff(realCovs, maincov)
##     dC <- setdiff(integerCovs, maincov)
##     ndataz <- nrow(X)
##     q <- parmList$q
##     nclusters <- ncol(q)
##     if(is.numeric(nfsamples)){
##         fsubsamples <- seq(1, nrow(q), length.out=nfsamples)
##     }else{
##         nfsamples <- nrow(q)
##         fsubsamples <- seq_len(nfsamples)
##     }
##     ##
##     freqs <- foreach(asample=fsubsamples, .combine=c, .inorder=inorder)%dopar%{
##         ## W: rows=clusters, cols=datapoints
##         W <- exp(
##             log(q[asample,]) +
##             t(vapply(seq_len(nclusters), function(acluster){
##                 ## real covariates
##                 if(length(cC)>0){
##                     colSums(dnorm(t(X[,cC]), mean=parmList$meanR[asample,cC,acluster], sd=sqrt(parmList$varR[asample,cC,acluster]), log=TRUE))
##                 }else{0} +
##                     ## integer covariates
##                     if(length(dC)>0){
##                         colSums(dbinom(t(X[,dC]), prob=parmList$probI[asample,dC,acluster], size=parmList$sizeI[asample,dC,acluster], log=TRUE))
##                     }else{0}
##             }, numeric(ndataz)))
##         )
##         ## pR: rows=clusters, cols= P at grid points
##         ##
##                 ( (if(maincov %in% realCovs){
##                     parmList$varR[asample,maincov,]
##                 }else{
##                     (1-parmList$probI[asample,maincov,]) * parmList$probI[asample,maincov,] * parmList$sizeI[asample,maincov,]
##                 }) %*% W)/colSums(W)
##     }
##     dim(freqs) <- c(ndataz, nfsamples)
##     freqs
## }
## ##
## ## Gives samples of cumulative frequency distributions of log_RMSD conditional on several feature values
## samplescdfRgivenX <- function(maincov, X, parmList, covgrid, nfsamples=NULL, inorder=FALSE){
##     realCovs <- dimnames(parmList$meanR)[[2]]
##     integerCovs <- dimnames(parmList$probI)[[2]]
##     cC <- setdiff(realCovs, maincov)
##     dC <- setdiff(integerCovs, maincov)
##     ndataz <- nrow(X)
##     lcovgrid <- length(covgrid)
##     q <- parmList$q
##     nclusters <- ncol(q)
##     if(is.numeric(nfsamples)){
##         fsubsamples <- seq(1, nrow(q), length.out=nfsamples)
##     }else{
##         nfsamples <- nrow(q)
##         fsubsamples <- seq_len(nfsamples)
##     }
##     ##
##     freqs <- foreach(asample=fsubsamples, .combine=c, .inorder=inorder)%dopar%{
##         ## W: rows=clusters, cols=datapoints
##         W <- exp(
##             log(q[asample,]) +
##             t(vapply(seq_len(nclusters), function(acluster){
##                 ## real covariates
##                 if(length(cC)>0){
##                     colSums(dnorm(t(X[,cC]), mean=parmList$meanR[asample,cC,acluster], sd=sqrt(parmList$varR[asample,cC,acluster]), log=TRUE))
##                 }else{0} +
##                     ## integer covariates
##                     if(length(dC)>0){
##                         colSums(dbinom(t(X[,dC]), prob=parmList$probI[asample,dC,acluster], size=parmList$sizeI[asample,dC,acluster], log=TRUE))
##                     }else{0}
##             }, numeric(ndataz)))
##         )
##         W <- matrix(rep(W, lcovgrid), nrow=nclusters) # strings copies of W column-wise
##         ## pR: rows=clusters, cols= P at grid points
##         pR <- t(vapply(seq_len(nclusters), function(acluster){
##             rep(
##                 if(maincov %in% realCovs){
##                     pnorm(q=covgrid, mean=parmList$meanR[asample,maincov,acluster], sd=sqrt(parmList$varR[asample,maincov,acluster]))
##                 }else{
##                     pbinom(q=covgrid, prob=parmList$probI[asample,maincov,acluster], size=parmList$sizeI[asample,maincov,acluster])
##                 } , each=ndataz)
##         }, numeric(lcovgrid * ndataz)))
##         ##
##         colSums(pR * W)/colSums(W)
##     }
##     dim(freqs) <- c(ndataz, lcovgrid, nfsamples)
##     freqs
## }
## ##
## ## Gives samples of frequency distributions of log_RMSD conditional on several feature values
## samplesfRgivenX <- function(maincov, X, parmList, covgrid, nfsamples=NULL, inorder=FALSE){
##     realCovs <- dimnames(parmList$meanR)[[2]]
##     integerCovs <- dimnames(parmList$probI)[[2]]
##     cC <- setdiff(realCovs, maincov)
##     dC <- setdiff(integerCovs, maincov)
##     ndataz <- nrow(X)
##     lcovgrid <- length(covgrid)
##     q <- parmList$q
##     nclusters <- ncol(q)
##     if(is.numeric(nfsamples)){
##         fsubsamples <- seq(1, nrow(q), length.out=nfsamples)
##     }else{
##         nfsamples <- nrow(q)
##         fsubsamples <- seq_len(nfsamples)
##     }
##     ##
##     freqs <- foreach(asample=fsubsamples, .combine=c, .inorder=inorder)%dopar%{
##         ## W: rows=clusters, cols=datapoints
##         W <- exp(
##             log(q[asample,]) +
##             t(vapply(seq_len(nclusters), function(acluster){
##                 ## real covariates
##                 if(length(cC)>0){
##                     colSums(dnorm(t(X[,cC]), mean=parmList$meanR[asample,cC,acluster], sd=sqrt(parmList$varR[asample,cC,acluster]), log=TRUE))
##                 }else{0} +
##                     ## integer covariates
##                     if(length(dC)>0){
##                         colSums(dbinom(t(X[,dC]), prob=parmList$probI[asample,dC,acluster], size=parmList$sizeI[asample,dC,acluster], log=TRUE))
##                     }else{0}
##             }, numeric(ndataz)))
##         )
##         W <- matrix(rep(W, lcovgrid), nrow=nclusters) # strings copies of W column-wise
##         ## pR: rows=clusters, cols= P at grid points
##         pR <- t(vapply(seq_len(nclusters), function(acluster){
##             rep(
##                 if(maincov %in% realCovs){
##                     dnorm(x=covgrid, mean=parmList$meanR[asample,maincov,acluster], sd=sqrt(parmList$varR[asample,maincov,acluster]))
##                 }else{
##                     dbinom(x=covgrid, prob=parmList$probI[asample,maincov,acluster], size=parmList$sizeI[asample,maincov,acluster])
##                 } , each=ndataz)
##         }, numeric(lcovgrid * ndataz)))
##         ##
##         colSums(pR * W)/colSums(W)
##     }
##     dim(freqs) <- c(ndataz, lcovgrid, nfsamples)
##     freqs
## }
## ##
## ## Gives samples of frequency distributions of log_RMSD conditional on one feature value
## samplesfRgivenX1 <- function(X, parmList, RMSDgrid, nfsamples=NULL){
##     realCovs <- dimnames(parmList$meanR)[[2]]
##     integerCovs <- dimnames(parmList$probI)[[2]]
##     cC <- setdiff(realCovs, 'log_RMSD')
##     X <- as.matrix(X)[1,]
##     lRMSDgrid <- length(RMSDgrid)
##     q <- parmList$q
##     nclusters <- ncol(q)
##     if(is.numeric(nfsamples)){
##         fsubsamples <- round(seq(1, nrow(q), length.out=nfsamples))
##     }else{
##         nfsamples <- nrow(q)
##         fsubsamples <- seq_len(nfsamples)
##     }
    
##     ##
##         ## W: rows=samples, cols=clusters
##         W <- exp(
##             log(q[fsubsamples,]) +
##             vapply(seq_len(nclusters), function(acluster){
##                 ## real covariates
##                 colSums(dnorm(X[cC], mean=aperm(parmList$meanR[fsubsamples,cC,acluster,drop=FALSE], c(2,1,3)), sd=sqrt(aperm(parmList$varR[fsubsamples,cC,acluster,drop=FALSE], c(2,1,3))), log=TRUE)) +
##                     ## integer covariates
##                     colSums(dbinom(X[integerCovs], prob=aperm(parmList$probI[fsubsamples,,acluster,drop=FALSE], c(2,1,3)), size=aperm(parmList$sizeI[fsubsamples,,acluster,drop=FALSE], c(2,1,3)), log=TRUE))
##             }, numeric(nfsamples))
##         )
##         ## W: rows=clusters*gridpoints, cols=samples
##         W <- matrix(rep(W, lRMSDgrid), ncol=nfsamples, byrow=TRUE) # strings copies of transposed W row-wise
##         dim(W) <- c(nclusters, lRMSDgrid, nfsamples)
##         ## pR: rows=clusters, cols= P at grid points
##         pR <- vapply(RMSDgrid, function(aRvalue){
##            dnorm(x=aRvalue, mean=parmList$meanR[fsubsamples,'log_RMSD',], sd=sqrt(parmList$varR[fsubsamples,'log_RMSD',]))
##         }, numeric(nclusters * nfsamples))
##     dim(pR) <- c(nfsamples, nclusters, lRMSDgrid)
##     pR <- aperm(pR, c(2,3,1))
##         ##
##     colSums(pR * W)/colSums(W)
## }
## ##
## ## Calculates the predictive distribution on a grid of log_RMSD conditional on several feature values
## pRgivenX <- function(X, parmList, RMSDgrid){
##     realCovs <- dimnames(parmList$meanR)[[2]]
##     integerCovs <- dimnames(parmList$probI)[[2]]
##     cC <- setdiff(realCovs, 'log_RMSD')
##     ndataz <- nrow(X)
##     q <- parmList$q
##     nsamples <- nrow(q)
##     nclusters <- ncol(q)
##     lRMSDgrid <- length(RMSDgrid)
##     ##
##     freqs <- foreach(asample=seq_len(nsamples), .combine='+', .inorder=FALSE)%dopar%{
##         ## W: rows=clusters, cols=datapoints
##         W <- exp(
##             log(q[asample,]) +
##             t(vapply(seq_len(nclusters), function(acluster){
##                 ## real covariates
##                 colSums(dnorm(t(X[,cC]), mean=parmList$meanR[asample,cC,acluster], sd=sqrt(parmList$varR[asample,cC,acluster]), log=TRUE)) +
##                     ## integer covariates
##                     colSums(dbinom(t(X[,integerCovs]), prob=parmList$probI[asample,,acluster], size=parmList$sizeI[asample,,acluster], log=TRUE))
##             }, numeric(ndataz)))
##         )
##         W <- matrix(rep(W,lRMSDgrid),nrow=nrow(W))
##         ##
##         pR <- t(vapply(seq_len(nclusters), function(acluster){
##             ## real covariates
##             dnorm(x=rep(RMSDgrid, each=ndataz), mean=parmList$meanR[asample,'log_RMSD',acluster], sd=sqrt(parmList$varR[asample,'log_RMSD',acluster]))
##         }, numeric(ndataz*lRMSDgrid)))
##         ##
##         colSums(pR * W)/colSums(W)
##     }/nsamples
##     dim(freqs) <- c(ndataz, lRMSDgrid)
##     freqs
## }
## ##
## ## Gives samples of marginal frequency distributions of a covariate
## samplesfX <- function(acov, parmList, acovgrid, nfsamples=100){
##     realCovs <- dimnames(parmList$meanR)[[2]]
##     integerCovs <- dimnames(parmList$probI)[[2]]
##     lacovgrid <- length(acovgrid)
##     q <- parmList$q
##     nclusters <- ncol(q)
##     if(is.numeric(nfsamples)){
##         fsubsamples <- seq(1, nrow(q), length.out=nfsamples)
##     }else{
##         nfsamples <- nrow(q)
##         fsubsamples <- seq_len(nfsamples)
##     }
##     ##
##     if(acov %in% realCovs){ ## real covariates
##     freqs <- foreach(asample=fsubsamples, .combine=cbind, .inorder=FALSE)%dopar%{
##         ## W: rows=clusters, cols=datapoints
##         colSums( exp(
##             log(q[asample,]) +
##             t(vapply(seq_len(nclusters), function(acluster){
##                 dnorm(acovgrid, mean=parmList$meanR[asample,acov,acluster], sd=sqrt(parmList$varR[asample,acov,acluster]), log=TRUE)
##             }, numeric(lacovgrid)))
##             ) )
##     }
##     }else{ ## integer covariates
##     freqs <- foreach(asample=fsubsamples, .combine=cbind, .inorder=FALSE)%dopar%{
##         ## W: rows=clusters, cols=datapoints
##         colSums( exp(
##             log(q[asample,]) +
##             t(vapply(seq_len(nclusters), function(acluster){
##                 dbinom(acovgrid, prob=parmList$probI[asample,acov,acluster], size=parmList$sizeI[asample,acov,acluster], log=TRUE)
##                             }, numeric(lacovgrid)))
##             ) )
##     }
##     }
##     freqs
## }






## transfold <- function(x, varinfo, Iout='init', Sout='data', Tout='data', variates=NULL){ # 'in' 'data' 'aux'
##     x <- cbind(data.matrix(x))
##     if(!is.null(variates)){colnames(x) <- variates}
##     matrix(sapply(colnames(x), function(v){
##         datum <- x[,v,drop=F]
##         info <- varinfo[v,]
##         type <- info['type']
##         loc <- info['location']
##         sca <- info['scale']
##         ##
##         if(type %in% variatetypes[c('R','B')]){ # strictly positive or binary
##             datum <- (datum-loc)/sca
##             ##
##         } else if(type == variatetypes['L']){ # logarithmic
##             datum <- (log(datum)-loc)/sca
##             ##
##         } else if(type == variatetypes['I']){ # integer, discrete ordinal
##             n <- info['n'] 
##             if(Iout == 'init'){ # in sampling functions or init MCMC
##                 datum[is.na(datum)] <- (info['min']+info['max'])/2
##                 datum <- qnorm((round((datum-loc)/sca)+0.5)/n)
##             } else if(Iout == 'left'){ # as left for MCMC
##                 datum <- qnorm(pmax(0,round((datum-loc)/sca))/n)
##                 datum[is.na(datum)] <- -Inf
##             } else if(Iout == 'right'){ # as right for MCMC
##                 datum <- qnorm(pmin(n,round((datum-loc)/sca)+1)/n)
##                 datum[is.na(datum)] <- +Inf
##             } else if(Iout == 'index'){ # in output functions
##                 datum <- round((datum-loc)/sca)+1L
##             }
##             ##
##         } else if(type == variatetypes['S']){ # logarithmic censored
##             datum <- (log(datum)-loc)/sca
##             if(Sout == 'left'){ # in MCMC
##                 sel <- is.na(datum) | (x[,v] < info['max'])
##                 datum[sel] <- -Inf
##                 datum[!sel] <- (log(info['max'])-loc)/sca
##             } else if(Sout == 'data'){ # data in MCMC
##                 sel <- is.na(datum) | (x[,v] >= info['max'])
##                 datum[sel] <- NA
##             } else if(Sout == 'init'){ #init in MCMC
##                 datum[is.na(datum)] <- 0L
##                 datum[x[,v] >= info['max']] <- (log(info['max'])-loc)/sca+0.125
##             } else if(Sout == 'index'){ #in sampling functions
##                 datum[x[,v] >= info['max']] <- +Inf
##             }
##             ##
##         } else if(type == variatetypes['T']){ # continuous doubly-bounded
##             datum <- qnorm((datum-loc)/sca)
##             if(Tout == 'left'){ # in sampling functions
##                 sel <- is.na(datum) | (x[,v] < info['max'])
##                 datum[sel] <- -Inf
##                 datum[!sel] <- qnorm(1-info['n'])
##             } else if(Tout == 'right'){ # in sampling functions
##                 sel <- is.na(datum) | (x[,v] > info['min'])
##                 datum[sel] <- +Inf
##                 datum[!sel] <- qnorm(info['n'])
##             } else if(Tout == 'data'){ # as init for MCMC
##                 datum[is.na(datum) | (x[,v] >= info['max']) | (x[,v] <= info['min'])] <- NA
##             } else if(Tout == 'init'){ # as init for MCMC
##                 datum[is.na(datum)] <- 0L
##                 datum[(x[,v] >= info['max'])] <- qnorm(1-info['n'])+0.125
##                 datum[(x[,v] <= info['min'])] <- qnorm(info['n'])-0.125
##             } else if(Tout == 'index'){ # in sampling functions
##                 datum[(x[,v] >= info['max'])] <- +Inf
##                 datum[(x[,v] <= info['min'])] <- -Inf
##             }
##         }
##         datum
##     }),ncol=ncol(x),dimnames=list(NULL,colnames(x)))
## }

## }

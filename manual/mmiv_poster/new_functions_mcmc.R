## Extend range
extendrange <- function(vec){
    1.125 * range(vec, na.rm=TRUE)- 0.125 * mean(range(vec, na.rm=TRUE))
}
## Normalize vector
normalize <- function(freqs){freqs/sum(freqs)}
## Normalize rows of matrix
normalizerows <- function(freqs){freqs/rowSums(freqs)}
##
## Construct a list of parameter samples from the raw MCMC samples
mcsamples2parmlist <- function(mcsamples){
    parmNames <- c('q', 'meanR', 'tauR', 'probI', 'sizeI', 'probB')
    nclusters <- sum(grepl('^q\\[', colnames(mcsamples)))
    nrcovs <- sum(grepl('^meanR\\[[^,]*, 1]', colnames(mcsamples)))
    if(nrcovs != length(realCovs)){
        warning('**WARNING: some problems with real variates**')
        }
    nicovs <- sum(grepl('^probI\\[[^,]*, 1]', colnames(mcsamples)))
    if(nicovs != length(integerCovs)){
        warning('**WARNING: some problems with integer variates**')
        }
    nbcovs <- sum(grepl('^probB\\[[^,]*, 1]', colnames(mcsamples)))
    if(nbcovs != length(binaryCovs)){
        warning('**WARNING: some problems with binary variates**')
        }
    ##
    parmList <- foreach(var=parmNames)%do%{
        out <- mcsamples[,grepl(paste0('^',var,'\\['), colnames(mcsamples))]
        if((var=='meanR'||var=='tauR')){
            dim(out) <- c(nrow(mcsamples), nrcovs, nclusters)
            dimnames(out) <- list(NULL, realCovs, NULL)
        } else if((var=='probI'||var=='sizeI')){
            dim(out) <- c(nrow(mcsamples), nicovs, nclusters)
            dimnames(out) <- list(NULL, integerCovs, NULL)
        } else if((var=='probB')){
            dim(out) <- c(nrow(mcsamples), nbcovs, nclusters)
            dimnames(out) <- list(NULL, binaryCovs, NULL)
        } else if(var=='q'){
            dim(out) <- c(nrow(mcsamples), nclusters)
        } else {NULL}
            out
    }
    names(parmList) <- parmNames
    parmList
}
##
## Construct a list of parameter samples from the raw MCMC samples for the second monitored set
finalstate2list <- function(mcsamples){
    if(!is.vector(mcsamples)){print('ERROR!')}
    parmNames <- c('q', 'meanR', 'tauR', 'probI', 'sizeI', 'probB', 'C')
    nclusters <- sum(grepl('^q\\[', names(mcsamples)))
    nrcovs <- sum(grepl('^meanR\\[[^,]*, 1]', names(mcsamples)))
    nicovs <- sum(grepl('^probI\\[[^,]*, 1]', names(mcsamples)))
    nbcovs <- sum(grepl('^probB\\[[^,]*, 1]', names(mcsamples)))
    ##
    parmList <- foreach(var=parmNames)%dopar%{
        out <- mcsamples[grepl(paste0('^',var,'\\['), names(mcsamples))]
        if(var=='meanR'||var=='tauR'){
            dim(out) <- c(nrcovs, nclusters)
            dimnames(out) <- list(realCovs, NULL)
        } else if(var=='probI'||var=='sizeI'){
            dim(out) <- c(nicovs, nclusters)
            dimnames(out) <- list(integerCovs, NULL)
        } # 'q' and 'C' are vectors with no names
            out
    }
    names(parmList) <- parmNames
    parmList
}
##
## Calculates the probability of the data (likelihood of parameters) for several MCMC samples
llSamples <- function(dat, parmList){
    realCovs <- dimnames(parmList$meanR)[[2]]
    integerCovs <- dimnames(parmList$probI)[[2]]
    binaryCovs <- dimnames(parmList$probB)[[2]]
    covNames <- c(realCovs, integerCovs, binaryCovs)
    nrcovs <- length(realCovs)
    nicovs <- length(integerCovs)
    nbcovs <- length(binaryCovs)
    ndata <- nrow(dat[[which.max(c(nrcovs,nicovs,nbcovs))]])
    q <- parmList$q
    ##
    foreach(asample=seq_len(nrow(q)), .combine=c, .inorder=TRUE)%dopar%{
        sum( log( colSums(
            exp(
                log(q[asample,]) +
                t(vapply(seq_len(ncol(q)), function(acluster){
                    ## real covariates
                    if(nrcovs>0){
                        colSums(dnorm(t(dat$X), mean=parmList$meanR[asample,,acluster], sd=1/sqrt(parmList$tauR[asample,,acluster]), log=TRUE))
                        }else{0} +
                        ## integer covariates
                        if(nicovs>0){
                            colSums(dbinom(t(dat$Y), prob=parmList$probI[asample,,acluster], size=parmList$sizeI[asample,,acluster], log=TRUE))
                            }else{0} +
                    ## binary covariates
                    if(nbcovs>0){
                            colSums(log(
                                parmList$probB[asample,,acluster] * t(dat$Z) +
                                (1-parmList$probB[asample,,acluster]) * (1-t(dat$Z))
                                ))
                    }else{0}
    }, numeric(ndata)))
            )
        ) ) )
    }
}
## tottime <- Sys.time()
## for(i in 1:1e6){extraDistr::dbern(x=testx, prob=testp, log=TRUE)}
## Sys.time()-tottime
## ## Time difference of 56.0237 secs
## tottime <- Sys.time()
## for(i in 1:1e6){log(testp*testx + (1-testp)*(1L-testx))}
## Sys.time()-tottime
## ## Time difference of 31.73911 secs
## tottime <- Sys.time()
## for(i in 1:1e6){log(testx + (-1L)^testx*(1-testp))}
## Sys.time()-tottime
## ## Time difference of 51.8646 secs
##
## Calculates means and covariances of the sampled frequency distributions
moments12Samples <- function(parmList){
    realCovs <- dimnames(parmList$meanR)[[2]]
    integerCovs <- dimnames(parmList$probI)[[2]]
    binaryCovs <- dimnames(parmList$probB)[[2]]
    covNames <- c(realCovs, integerCovs, binaryCovs)
    nrcovs <- length(realCovs)
    nicovs <- length(integerCovs)
    nbcovs <- length(binaryCovs)
    ncovs <- length(covNames)
    q <- t(parmList$q)
    ##
if(length(realCovs)>0){
    meansr <- aperm(parmList$meanR, c(3, 1, 2))
    quadr <- aperm(1/parmList$tauR + parmList$meanR * parmList$meanR, c(3, 1, 2))
}else{
    meansr <- NULL
    quadr <- NULL
}
    if(length(integerCovs)>0){
    meansi <- aperm(parmList$probI * parmList$sizeI, c(3, 1, 2))
    quadi <- aperm(parmList$probI * parmList$sizeI *
                    (1 + parmList$probI * (parmList$sizeI - 1)), c(3, 1, 2))
    }else{
        meansi <- NULL
        quadi <- NULL
}
    if(length(binaryCovs)>0){
    meansb <- aperm(parmList$probB, c(3, 1, 2))
    quadb <- aperm(parmList$probB * parmList$probB, c(3, 1, 2))
    }else{
        meansb <- NULL
        quadb <- NULL
}
    clustermeans <- c(meansr, meansi, meansb)
    dim(clustermeans) <- c(dim((q)), ncovs)
    mixmeans <- colSums(c(q) * clustermeans)
    dimnames(mixmeans) <- list(NULL, paste0('MEAN_', covNames))
    ##
    mixvars <- c(quadr, quadi, quadb)
    dim(mixvars) <- c(dim((q)), ncovs)
    mixvars <- colSums(c(q) * mixvars) - mixmeans*mixmeans
    dimnames(mixvars) <- list(NULL, paste0('VAR_', covNames))
    ##
    mixcovars <- foreach(cov1=seq_len(ncovs-1), .combine=cbind)%:%foreach(cov2=(cov1+1):ncovs, .combine=cbind)%do%{
       colSums(c(q)*clustermeans[,,cov1,drop=FALSE]*clustermeans[,,cov2,drop=FALSE]) - mixmeans[,cov1,drop=FALSE]*mixmeans[,cov2,drop=FALSE]
    }
    colnames(mixcovars) <- foreach(cov1=seq_len(ncovs-1), .combine=c)%:%foreach(cov2=(cov1+1):ncovs, .combine=c)%do%{paste0('COV_',covNames[cov1],'|',covNames[cov2])}
    ##
    list(
        Dcovars=((matrix(colSums(c(q) * apply(
        (clustermeans -
         array(rep(c(mixmeans), each=nrow(q)), dim=dim(clustermeans)))/array(rep(sqrt(c(mixvars)), each=nrow(q)), dim=dim(clustermeans)),
        c(1,2), prod)),
        ncol=1, dimnames=list(NULL, 'Dcov')))),
        means=mixmeans,
        vars=(mixvars),
        covars=mixcovars
        )
}
##
## Gives samples of frequency distributions of any set of Y conditional on any set of X
samplesF <- function(Y, X=NULL, parmList, nfsamples=NULL, inorder=FALSE){
    rCovs <- dimnames(parmList$meanR)[[2]]
    iCovs <- dimnames(parmList$probI)[[2]]
    bCovs <- dimnames(parmList$probB)[[2]]
    covNames <- c(rCovs, iCovs, bCovs)
    nrcovs <- length(rCovs)
    nicovs <- length(iCovs)
    nbcovs <- length(bCovs)
    ##
    Y <- data.matrix(rbind(Y))
    rY <- colnames(Y)[colnames(Y) %in% rCovs]
    iY <- colnames(Y)[colnames(Y) %in% iCovs]
    bY <- colnames(Y)[colnames(Y) %in% bCovs]
    ##
    if(!is.null(X)){
        X <- data.matrix(rbind(X))
        rX <- colnames(X)[colnames(X) %in% rCovs]
        if(length(intersect(rX,rY))>0){
            warning('*WARNING: predictor and predictand have real variates in common. Removing from predictor*')
            rX <- setdiff(rX,rY)
        }
        iX <- colnames(X)[colnames(X) %in% iCovs]
        if(length(intersect(iX,iY))>0){
            warning('*WARNING: predictor and predictand have integer variates in common. Removing from predictor*')
            iX <- setdiff(iX,iY)
        }
        bX <- colnames(X)[colnames(X) %in% bCovs]
        if(length(intersect(bX,bY))>0){
            warning('*WARNING: predictor and predictand have binary variates in common. Removing from predictor*')
            bX <- setdiff(bX,bY)
        }
        if(length(c(rX,iX,bX))==0){X <- NULL}
    }
    ##
    if(!is.null(X)){
        if(nrow(X) < nrow(Y)){
        warning('*Note: X has fewer data than Y. Recycling*')
        X <- t(matrix(rep(t(X), ceiling(nrow(Y)/nrow(X))), nrow=ncol(X), dimnames=list(colnames(X),NULL)))[1:nrow(Y),,drop=FALSE]
    }
    if(nrow(X) > nrow(Y)){
        warning('*Note: X has more data than Y. Recycling*')
        Y <- t(matrix(rep(t(Y), ceiling(nrow(X)/nrow(Y))), nrow=ncol(Y), dimnames=list(colnames(Y),NULL)))[1:nrow(X),,drop=FALSE]
    }
    }
    ndata <- nrow(Y)
    ##
    q <- parmList$q
    nclusters <- ncol(q)
    if(is.numeric(nfsamples)){
        fsubsamples <- round(seq(1, nrow(q), length.out=nfsamples))
    }else{
        nfsamples <- nrow(q)
        fsubsamples <- seq_len(nfsamples)
    }
    ##
    if(!is.null(X)){
        freqs <- foreach(asample=fsubsamples, .combine=cbind, .inorder=inorder)%dopar%{
            ## pX: rows=clusters, cols=datapoints
            pX <- exp(
                log(q[asample,]) +
                t(rbind(vapply(seq_len(nclusters), function(acluster){
                    ## real covariates
                    if(length(rX)>0){
                        colSums(dnorm(x=t(X[,rX,drop=FALSE]), mean=parmList$meanR[asample,rX,acluster], sd=1/sqrt(parmList$tauR[asample,rX,acluster]), log=TRUE))
                    }else{0} +
                        ## integer covariates
                        if(length(iX)>0){
                            colSums(dbinom(x=t(X[,iX,drop=FALSE]), prob=parmList$probI[asample,iX,acluster], size=parmList$sizeI[asample,iX,acluster], log=TRUE))
                        }else{0} +
                        ## binary covariates
                        if(length(bX)>0){
                            colSums(log(
                                parmList$probB[asample,bX,acluster] * t(X[,bX,drop=FALSE]) +
                                (1-parmList$probB[asample,bX,acluster]) * (1-t(X[,bX,drop=FALSE]))
                            ))
                        }else{0}
                }, numeric(ndata))))
            )
            ## pY: rows=clusters, cols=datapoints
            pY <- exp(
                t(rbind(vapply(seq_len(nclusters), function(acluster){
                    ## real covariates
                    if(length(rY)>0){
                        colSums(dnorm(x=t(Y[,rY,drop=FALSE]), mean=parmList$meanR[asample,rY,acluster], sd=1/sqrt(parmList$tauR[asample,rY,acluster]), log=TRUE))
                    }else{0} +
                        ## integer covariates
                        if(length(iY)>0){
                            colSums(dbinom(x=t(Y[,iY,drop=FALSE]), prob=parmList$probI[asample,iY,acluster], size=parmList$sizeI[asample,iY,acluster], log=TRUE))
                        }else{0} +
                        ## binary covariates
                        if(length(bY)>0){
                            colSums(log(
                                parmList$probB[asample,bY,acluster] * t(Y[,bY,drop=FALSE]) +
                                (1-parmList$probB[asample,bY,acluster]) * (1-t(Y[,bY,drop=FALSE]))
                            ))
                        }else{0}
                }, numeric(ndata))))
            )
            ##
            colSums(pX * pY)/colSums(pX)
        }
    }else{
        freqs <- foreach(asample=fsubsamples, .combine=cbind, .inorder=inorder)%dopar%{
            ## pY: rows=clusters, cols=datapoints
            pY <- exp(
                log(q[asample,]) +
                t(rbind(vapply(seq_len(nclusters), function(acluster){
                    ## real covariates
                    if(length(rY)>0){
                        colSums(dnorm(x=t(Y[,rY,drop=FALSE]), mean=parmList$meanR[asample,rY,acluster], sd=1/sqrt(parmList$tauR[asample,rY,acluster]), log=TRUE))
                    }else{0} +
                        ## integer covariates
                        if(length(iY)>0){
                            colSums(dbinom(x=t(Y[,iY,drop=FALSE]), prob=parmList$probI[asample,iY,acluster], size=parmList$sizeI[asample,iY,acluster], log=TRUE))
                        }else{0} +
                        ## binary covariates
                        if(length(bY)>0){
                            colSums(log(
                                parmList$probB[asample,bY,acluster] * t(Y[,bY,drop=FALSE]) +
                                (1-parmList$probB[asample,bY,acluster]) * (1-t(Y[,bY,drop=FALSE]))
                            ))
                        }else{0}
                }, numeric(ndata))))
            )
            ##
            colSums(pY)
        }
    }
    freqs
}
##
## Gives samples of means of distributions of any set of Y conditional on any set of X
samplesMeans <- function(Y, X=NULL, parmList, nfsamples=NULL, inorder=FALSE){
    rCovs <- dimnames(parmList$meanR)[[2]]
    iCovs <- dimnames(parmList$probI)[[2]]
    bCovs <- dimnames(parmList$probB)[[2]]
    covNames <- c(rCovs, iCovs, bCovs)
    nrcovs <- length(rCovs)
    nicovs <- length(iCovs)
    nbcovs <- length(bCovs)
    ##
    if(length(Y)==0 | length(intersect(Y, covNames))==0){stop('invalid variates in Y')}
    Y <- intersect(Y, covNames)
    rY <- Y[Y %in% rCovs]
    iY <- Y[Y %in% iCovs]
    bY <- Y[Y %in% bCovs]
    ##
    if(!is.null(X)){
        X <- data.matrix(rbind(X))
        rX <- colnames(X)[colnames(X) %in% rCovs]
        if(length(intersect(rX,rY))>0){
            warning('*WARNING: predictor and predictand have real variates in common. Removing from predictor*')
            rX <- setdiff(rX,rY)
        }
        iX <- colnames(X)[colnames(X) %in% iCovs]
        if(length(intersect(iX,iY))>0){
            warning('*WARNING: predictor and predictand have integer variates in common. Removing from predictor*')
            iX <- setdiff(iX,iY)
        }
        bX <- colnames(X)[colnames(X) %in% bCovs]
        if(length(intersect(bX,bY))>0){
            warning('*WARNING: predictor and predictand have binary variates in common. Removing from predictor*')
            bX <- setdiff(bX,bY)
        }
        if(length(c(rX,iX,bX))==0){X <- NULL}
    }
    ##
    ## if(!is.null(X)){
    ##     if(nrow(X) < nrow(Y)){
    ##     warning('*Note: X has fewer data than Y. Recycling*')
    ##     X <- t(matrix(rep(t(X), ceiling(nrow(Y)/nrow(X))), nrow=ncol(X), dimnames=list(colnames(X),NULL)))[1:nrow(Y),,drop=FALSE]
    ## }
    ## if(nrow(X) > nrow(Y)){
    ##     warning('*Note: X has more data than Y. Recycling*')
    ##     Y <- t(matrix(rep(t(Y), ceiling(nrow(X)/nrow(Y))), nrow=ncol(Y), dimnames=list(Y,NULL)))[1:nrow(X),,drop=FALSE]
    ## }
    ## }
    ndata <- nrow(X)
    ##
    q <- parmList$q
    nclusters <- ncol(q)
    if(is.numeric(nfsamples)){
        fsubsamples <- round(seq(1, nrow(q), length.out=nfsamples))
    }else{
        nfsamples <- nrow(q)
        fsubsamples <- seq_len(nfsamples)
    }
    ##
    if(!is.null(X)){
        freqs <- foreach(asample=fsubsamples, .combine=cbind, .inorder=inorder)%dopar%{
            ## pX: rows=clusters, cols=datapoints
            pX <- exp(
                log(q[asample,]) +
                t(rbind(vapply(seq_len(nclusters), function(acluster){
                    ## real covariates
                    if(length(rX)>0){
                        colSums(dnorm(x=t(X[,rX,drop=FALSE]), mean=parmList$meanR[asample,rX,acluster], sd=1/sqrt(parmList$tauR[asample,rX,acluster]), log=TRUE))
                    }else{0} +
                        ## integer covariates
                        if(length(iX)>0){
                            colSums(dbinom(x=t(X[,iX,drop=FALSE]), prob=parmList$probI[asample,iX,acluster], size=parmList$sizeI[asample,iX,acluster], log=TRUE))
                        }else{0} +
                        ## binary covariates
                        if(length(bX)>0){
                            colSums(log(
                                parmList$probB[asample,bX,acluster] * t(X[,bX,drop=FALSE]) +
                                (1-parmList$probB[asample,bX,acluster]) * (1-t(X[,bX,drop=FALSE]))
                            ))
                        }else{0}
                }, numeric(ndata))))
            )
            ## pY: rows=means, cols=clusters
            pY <- rbind(
                    ## real covariates
                if(length(rY)>0){
                    parmList$meanR[asample,rY,]
                    },
                        ## integer covariates
                if(length(iY)>0){
                    parmList$probI[asample,iY,] * parmList$sizeI[asample,iY,]
                    },
                ## binary covariates
                if(length(bY)>0){
                    parmList$probB[asample,bY,]
                }
                )
            ##
            t(sapply(1:nrow(pY),function(amean){
                colSums(pY[amean,] * pX)/colSums(pX)
            }))
        }
    }else{
        freqs <- foreach(asample=fsubsamples, .combine=cbind, .inorder=inorder)%dopar%{
            ## pY: rows=means, cols=clusters
            pY <- rbind(
                    ## real covariates
                if(length(rY)>0){
                    parmList$meanR[asample,rY,]
                    },
                        ## integer covariates
                if(length(iY)>0){
                    parmList$probI[asample,iY,] * parmList$sizeI[asample,iY,]
                    },
                ## binary covariates
                if(length(bY)>0){
                    parmList$probB[asample,bY,]
                }
                )
            ##
            colSums(q[asample,]*t(pY))
        }
    }
    dim(freqs) <- c(length(Y), ndata, nfsamples)
    dimnames(freqs) <- list(Y,NULL,NULL)
    freqs
}
##
## Gives samples of variances of distributions of any set of Y conditional on any set of X
samplesVars <- function(Y, X=NULL, parmList, nfsamples=NULL, inorder=FALSE){
    rCovs <- dimnames(parmList$meanR)[[2]]
    iCovs <- dimnames(parmList$probI)[[2]]
    bCovs <- dimnames(parmList$probB)[[2]]
    covNames <- c(rCovs, iCovs, bCovs)
    nrcovs <- length(rCovs)
    nicovs <- length(iCovs)
    nbcovs <- length(bCovs)
    ##
    if(length(Y)==0 | length(intersect(Y, covNames))==0){stop('invalid variates in Y')}
    Y <- intersect(Y, covNames)
    rY <- Y[Y %in% rCovs]
    iY <- Y[Y %in% iCovs]
    bY <- Y[Y %in% bCovs]
    ##
    if(!is.null(X)){
        X <- data.matrix(rbind(X))
        rX <- colnames(X)[colnames(X) %in% rCovs]
        if(length(intersect(rX,rY))>0){
            warning('*WARNING: predictor and predictand have real variates in common. Removing from predictor*')
            rX <- setdiff(rX,rY)
        }
        iX <- colnames(X)[colnames(X) %in% iCovs]
        if(length(intersect(iX,iY))>0){
            warning('*WARNING: predictor and predictand have integer variates in common. Removing from predictor*')
            iX <- setdiff(iX,iY)
        }
        bX <- colnames(X)[colnames(X) %in% bCovs]
        if(length(intersect(bX,bY))>0){
            warning('*WARNING: predictor and predictand have binary variates in common. Removing from predictor*')
            bX <- setdiff(bX,bY)
        }
        if(length(c(rX,iX,bX))==0){X <- NULL}
    }
    ##
    ## if(!is.null(X)){
    ##     if(nrow(X) < nrow(Y)){
    ##     warning('*Note: X has fewer data than Y. Recycling*')
    ##     X <- t(matrix(rep(t(X), ceiling(nrow(Y)/nrow(X))), nrow=ncol(X), dimnames=list(colnames(X),NULL)))[1:nrow(Y),,drop=FALSE]
    ## }
    ## if(nrow(X) > nrow(Y)){
    ##     warning('*Note: X has more data than Y. Recycling*')
    ##     Y <- t(matrix(rep(t(Y), ceiling(nrow(X)/nrow(Y))), nrow=ncol(Y), dimnames=list(Y,NULL)))[1:nrow(X),,drop=FALSE]
    ## }
    ## }
    ndata <- nrow(X)
    ##
    q <- parmList$q
    nclusters <- ncol(q)
    if(is.numeric(nfsamples)){
        fsubsamples <- round(seq(1, nrow(q), length.out=nfsamples))
    }else{
        nfsamples <- nrow(q)
        fsubsamples <- seq_len(nfsamples)
    }
    ##
    if(!is.null(X)){
        freqs <- foreach(asample=fsubsamples, .combine=cbind, .inorder=inorder)%dopar%{
            ## pX: rows=clusters, cols=datapoints
            pX <- exp(
                log(q[asample,]) +
                t(rbind(vapply(seq_len(nclusters), function(acluster){
                    ## real covariates
                    if(length(rX)>0){
                        colSums(dnorm(x=t(X[,rX,drop=FALSE]), mean=parmList$meanR[asample,rX,acluster], sd=1/sqrt(parmList$tauR[asample,rX,acluster]), log=TRUE))
                    }else{0} +
                        ## integer covariates
                        if(length(iX)>0){
                            colSums(dbinom(x=t(X[,iX,drop=FALSE]), prob=parmList$probI[asample,iX,acluster], size=parmList$sizeI[asample,iX,acluster], log=TRUE))
                        }else{0} +
                        ## binary covariates
                        if(length(bX)>0){
                            colSums(log(
                                parmList$probB[asample,bX,acluster] * t(X[,bX,drop=FALSE]) +
                                (1-parmList$probB[asample,bX,acluster]) * (1-t(X[,bX,drop=FALSE]))
                            ))
                        }else{0}
                }, numeric(ndata))))
            )
            ## pY: rows=means, cols=clusters
            pY <- rbind(
                    ## real covariates
                if(length(rY)>0){
                    1/parmList$tauR[asample,rY,]
                    },
                        ## integer covariates
                if(length(iY)>0){
                    (1-parmList$probI[asample,iY,]) * parmList$probI[asample,iY,] * parmList$sizeI[asample,iY,]
                    },
                ## binary covariates
                if(length(bY)>0){
                    (1-parmList$probB[asample,bY,]) * parmList$probB[asample,bY,]
                }
                )
            ##
            t(sapply(1:nrow(pY),function(amean){
                colSums(pY[amean,] * pX)/colSums(pX)
            }))
        }
    }else{
        freqs <- foreach(asample=fsubsamples, .combine=cbind, .inorder=inorder)%dopar%{
            ## pY: rows=means, cols=clusters
            pY <- rbind(
                    ## real covariates
                if(length(rY)>0){
                    1/parmList$tauR[asample,rY,]
                    },
                        ## integer covariates
                if(length(iY)>0){
                    (1-parmList$probI[asample,iY,]) * parmList$probI[asample,iY,] * parmList$sizeI[asample,iY,]
                    },
                ## binary covariates
                if(length(bY)>0){
                    (1-parmList$probB[asample,bY,]) * parmList$probB[asample,bY,]
                }
                )
            ##
            colSums(q[asample,]*t(pY))
        }
    }
    dim(freqs) <- c(length(Y), ndata, nfsamples)
    dimnames(freqs) <- list(Y,NULL,NULL)
    freqs
}




####################################################
#### Functions below are not used at the moment ####
####################################################
##
## Calculates the median and IQR of each sample frequency
calcSampleMQ <- function(parmList, maxD=1000){
    realCovs <- dimnames(parmList$meanR)[[2]]
    integerCovs <- dimnames(parmList$probI)[[2]]
    covNames <- c(realCovs, integerCovs)
    ncovs <- length(covNames)
    q <- parmList$q
    nsamples <- nrow(q)
    ##
    quants <- foreach(asample=seq_len(nsamples), .combine=c)%:%foreach(acov=covNames, .combine=c)%dopar%{
        if(acov %in% realCovs){
                mixq <- function(x){sum(q[asample,] * pnorm(x, mean=parmList$meanR[asample,acov,], sd=1/sqrt(parmList$tauR[asample,acov,])))}
                fn <- function(par){
                    (mixq(par[1]) - 0.5)^2 +
                        (mixq(par[2]) - 0.25)^2 +
                        (mixq(par[3]) - 0.75)^2
                }
                out <- myoptim(par=rep(q[asample,]%*%parmList$meanR[asample,acov,],3), fn=fn)$par
            ##     out <- sapply(c(0.005, 0.995), function(border){
            ## optim(0, #rep(q[asample,] %*% parmList$meanR[asample,acov,], 2),
            ##              fn=fn,
            ##              gr = function(x) pracma::grad(fn, x), 
            ##              method = "L-BFGS-B",
            ##              lower = -Inf, upper = Inf,
            ##              control = list(factr = 1e-10, pgtol = 0, maxit = 100))
        }else{
            searchgrid <- 0:maxD
            dq <- colSums(c(q[asample,]) * pbinom(matrix(searchgrid, ncol=length(searchgrid), nrow=ncol(q), byrow=TRUE), prob=parmList$probI[asample,acov,], size=parmList$sizeI[asample,acov,]))
            out <- c(
                which.min(abs(dq - 0.5))-1,
                which.min(abs(dq - 0.25))-1,
                which.min(abs(dq - 0.75))-1
            )
        }
        out
    }
    ##
    dim(quants) <- c(3, ncovs, nsamples)
    quants <- aperm(quants)
    dimnames(quants) <- list(NULL, covNames, c('50%', '25%', '75%'))
    quants
}
##
## Calculates the probability of some datapoints for the MCMC samples
probValuesSamples <- function(X, parmList){
    realCovs <- dimnames(parmList$meanR)[[2]]
    integerCovs <- dimnames(parmList$probI)[[2]]
    covNames <- c(realCovs, integerCovs)
    ncovs <- length(covNames)
    q <- parmList$q
    ndataz <- nrow(X)
    ##
    (foreach(asample=seq_len(nrow(q)), .combine=cbind, .inorder=TRUE)%dopar%{
        colSums(
            exp(
                log(q[asample,]) +
                t(vapply(seq_len(ncol(q)), function(acluster){
                    ## real covariates
                    colSums(dnorm(t(X[,realCovs]), mean=parmList$meanR[asample,,acluster,drop=FALSE], sd=1/sqrt(parmList$tauR[asample,,acluster,drop=FALSE]), log=TRUE)) +
                        ## integer covariates
                    colSums(dbinom(t(X[,integerCovs]), prob=parmList$probI[asample,,acluster,drop=FALSE], size=parmList$sizeI[asample,,acluster,drop=FALSE], log=TRUE))
    }, numeric(ndataz)))
            )
        )
    })
}
##
## Improved optimization functions
myoptim <- function(par, fn){
    resu0 <- list(par=par)
    resu <- optim(par=resu0$par, fn=fn, control=list(factr = 1e-10, maxit=10000))
    while(any(resu$par!=resu0$par)){
        resu0 <- resu
        resu <- optim(par=resu0$par, fn=fn, control=list(factr = 1e-10, maxit=10000))
    }
    resu}
myoptimbounds <- function(par, fn, lower, upper, maxit=100){
    optim(par, fn=fn,
          gr = function(x) pracma::grad(fn, x), 
          method = "L-BFGS-B",
          lower = lower, upper = upper,
          control = list(factr = 1e-10, maxit = maxit))
}
##
## Calculates the 0.5% and 99.5% quantiles of each sample frequency
calcSampleQuantiles <- function(parmList){
    realCovs <- dimnames(parmList$meanR)[[2]]
    integerCovs <- dimnames(parmList$probI)[[2]]
    covNames <- c(realCovs, integerCovs)
    ncovs <- length(covNames)
    q <- parmList$q
    nsamples <- nrow(q)
    ##
    quants <- foreach(asample=seq_len(nsamples), .combine=c)%:%foreach(acov=covNames, .combine=c)%dopar%{
        if(acov %in% realCovs){
                mixq <- function(x){sum(q[asample,] * pnorm(x, mean=parmList$meanR[asample,acov,], sd=1/sqrt(parmList$tauR[asample,acov,])))}
                fn <- function(par){(mixq(par[1]) - 0.005)^2 +
                                        (mixq(par[2]) - 0.995)^2}
                out <- myoptim(par=rep(q[asample,]%*%parmList$meanR[asample,acov,],2), fn=fn)$par
            ##     out <- sapply(c(0.005, 0.995), function(border){
            ## optim(0, #rep(q[asample,] %*% parmList$meanR[asample,acov,], 2),
            ##              fn=fn,
            ##              gr = function(x) pracma::grad(fn, x), 
            ##              method = "L-BFGS-B",
            ##              lower = -Inf, upper = Inf,
            ##              control = list(factr = 1e-10, pgtol = 0, maxit = 100))
        }else{
            searchgrid <- 0:max(parmList$sizeI[asample,acov,])
            dq <- colSums(c(q[asample,]) * pbinom(matrix(searchgrid, ncol=length(searchgrid), nrow=ncol(q), byrow=TRUE), prob=parmList$probI[asample,acov,], size=parmList$sizeI[asample,acov,]))
            out <- c(which.min(abs(dq - 0.005))-1, which.min(abs(dq - 0.995))-1)
        }
        out
    }
    ##
    dim(quants) <- c(2, ncovs, nsamples)
    quants <- aperm(quants)
    dimnames(quants) <- list(NULL, covNames, c('0.5%', '99.5%'))
    quants
}
##
## Calculates the probability of the data (likelihood of parameters) for one MCMC samples
## if(!exists('calcLL')){
## calcLL <- nimbleFunction( run=function(
##     X=double(2), Y=double(2), Q=double(1),
##     MeanC=double(2), TauC=double(2),
##     ProbD=double(2), SizeD=double(2)
##     ){
##     returnType(double(0))
##     Nclusters <- length(Q)
##     Ndata <- dim(X)[1]
##     Nrcovs <- dim(X)[2]
##     Nicovs <- dim(Y)[2]
##     LL <- 0
##     for(adatum in 1:Ndata){
##         clustersum <- log(Q)
##         for(acov in 1:Nrcovs){
##             clustersum <- clustersum +
##                 dnorm(x=X[adatum,acov], mean=MeanC[acov,], sd=1/sqrt(TauC[acov,]), log=TRUE)
##         }
##         for(acov in 1:Nicovs){
##             clustersum <- clustersum +
##                 dbinom(x=Y[adatum,acov], prob=ProbD[acov,], size=SizeD[acov,], log=TRUE)
##         }
##         LL <- LL + log(sum(exp(clustersum)))
##     }
##     return(LL)
## } )
## CcalcLL <- compileNimble(calcLL)
## assign('CcalcLL', CcalcLL, envir = .GlobalEnv)
## assign('calcLL', calcLL, envir = .GlobalEnv)
## }
##
## Calculates the probability of several datapoints for several MCMC samples
probJointSamples <- function(dat, parmList, log=FALSE, inorder=FALSE){
    ndataz <- nrow(dat$X)
    q <- parmList$q
    ##
    freqs <- foreach(asample=seq_len(nrow(q)), .combine=cbind, .inorder=inorder)%dopar%{
        colSums(
            exp(
                log(q[asample,]) +
                t(vapply(seq_len(ncol(q)), function(acluster){
                    ## real covariates
                    colSums(dnorm(t(dat$X), mean=parmList$meanR[asample,,acluster], sd=1/sqrt(parmList$tauR[asample,,acluster]), log=TRUE)) +
                        ## integer covariates
                    colSums(dbinom(t(dat$Y), prob=parmList$probI[asample,,acluster], size=parmList$sizeI[asample,,acluster], log=TRUE))
    }, numeric(ndataz)))
            )
        )
    }
    if(!log){freqs} else {log(freqs)}
}
## Calculates the MCMC posterior probability of several datapoints
probJointMean <- function(X, parmList){
    realCovs <- dimnames(parmList$meanR)[[2]]
    integerCovs <- dimnames(parmList$probI)[[2]]
    if(is.list(X)){ X <- cbind(X$X, X$Y) }
    ndataz <- nrow(X)
    q <- parmList$q
    nsamples <- nrow(q)
    ##
    foreach(asample=seq_len(nsamples), .combine='+', .inorder=FALSE)%dopar%{
        colSums(
            exp(
                log(q[asample,]) +
                t(vapply(seq_len(ncol(q)), function(acluster){
                    ## real covariates
                    colSums(dnorm(t(X[,realCovs]), mean=parmList$meanR[asample,,acluster], sd=1/sqrt(parmList$tauR[asample,,acluster]), log=TRUE)) +
                        ## integer covariates
                    colSums(dbinom(t(X[,integerCovs]), prob=parmList$probI[asample,,acluster], size=parmList$sizeI[asample,,acluster], log=TRUE))
    }, numeric(ndataz)))
            )
        )
    }/nsamples
}
##
## Produces multidimensional samples from several MCMC distribution-samples
options(doFuture.rng.onMisuse = "ignore")
samplesFsamples <- function(parmList, nxsamples=1000, nfsamples=NULL, seed=149){
    realCovs <- dimnames(parmList$meanR)[[2]]
    integerCovs <- dimnames(parmList$probI)[[2]]
    ncC <- length(realCovs)
    ndC <- length(integerCovs)
    q <- parmList$q
        if(is.numeric(nfsamples)){
        fsubsamples <- seq(1, nrow(q), length.out=nfsamples)
    }else{
        nfsamples <- nrow(q)
        fsubsamples <- seq_len(nfsamples)
    }
    ##
    rng <- RNGseq( nfsamples * nxsamples, seed)
    allsamples <- foreach(afsample=fsubsamples, .combine=cbind, .inorder=FALSE)%:%foreach(axsample=seq_len(nxsamples), r=rng[(afsample-1)*nxsamples + 1:nxsamples], .combine=c, .inorder=FALSE)%dopar%{
        rngtools::setRNG(r)
        acluster <- rcat(n=1, prob=q[afsample,])
        c(
            ## real covariates
            rnorm(n=ncC, mean=parmList$meanR[afsample,realCovs,acluster], sd=1/sqrt(parmList$tauR[afsample,realCovs,acluster])),
            ## integer covariates
            rbinom(n=ndC, prob=parmList$probI[afsample,integerCovs,acluster], size=parmList$sizeI[afsample,integerCovs,acluster])
        )
    }
    dim(allsamples) <- c(ncC+ndC, nxsamples, nfsamples)
    dimnames(allsamples) <- list(c(realCovs,integerCovs), NULL, NULL)
    allsamples
}
##
## Calculates the predictive probability for several log_RMSD values conditional on several feature values
expeRgivenX <- function(maincov, X, parmList){
    realCovs <- dimnames(parmList$meanR)[[2]]
    integerCovs <- dimnames(parmList$probI)[[2]]
    cC <- setdiff(realCovs, maincov)
    dC <- setdiff(integerCovs, maincov)
    ndataz <- nrow(X)
    q <- parmList$q
    nsamples <- nrow(q)
    ##
    foreach(asample=seq_len(nsamples), .combine='+', .inorder=FALSE)%dopar%{
        ## W: rows=clusters, cols=datapoints
        W <- exp(
            log(q[asample,]) +
            t(vapply(seq_len(ncol(q)), function(acluster){
                ## real covariates
                if(length(cC)>0){
                    colSums(dnorm(t(X[,cC]), mean=parmList$meanR[asample,cC,acluster], sd=1/sqrt(parmList$tauR[asample,cC,acluster]), log=TRUE))
                }else{0} +
                    ## integer covariates
                    if(length(dC)>0){
                        colSums(dbinom(t(X[,dC]), prob=parmList$probI[asample,dC,acluster], size=parmList$sizeI[asample,dC,acluster], log=TRUE))
                    }else{0}
            }, numeric(ndataz)))
        )
        ##
        colSums(
                if(maincov %in% realCovs){
                    parmList$meanR[asample,maincov,]
                }else{
                    parmList$probI[asample,maincov,] * parmList$sizeI[asample,maincov,] 
                } * W)/colSums(W)
    }/nsamples
}
##
## Gives samples of mean of one covariate conditional on several feature values
samplesmeanRgivenX <- function(maincov, X, parmList, nfsamples=NULL, inorder=FALSE){
    realCovs <- dimnames(parmList$meanR)[[2]]
    integerCovs <- dimnames(parmList$probI)[[2]]
    cC <- setdiff(realCovs, maincov)
    dC <- setdiff(integerCovs, maincov)
    ndataz <- nrow(X)
    q <- parmList$q
    nclusters <- ncol(q)
    if(is.numeric(nfsamples)){
        fsubsamples <- seq(1, nrow(q), length.out=nfsamples)
    }else{
        nfsamples <- nrow(q)
        fsubsamples <- seq_len(nfsamples)
    }
    ##
    freqs <- foreach(asample=fsubsamples, .combine=c, .inorder=inorder)%dopar%{
        ## W: rows=clusters, cols=datapoints
        W <- exp(
            log(q[asample,]) +
            t(vapply(seq_len(nclusters), function(acluster){
                ## real covariates
                if(length(cC)>0){
                    colSums(dnorm(t(X[,cC]), mean=parmList$meanR[asample,cC,acluster], sd=1/sqrt(parmList$tauR[asample,cC,acluster]), log=TRUE))
                }else{0} +
                    ## integer covariates
                    if(length(dC)>0){
                        colSums(dbinom(t(X[,dC]), prob=parmList$probI[asample,dC,acluster], size=parmList$sizeI[asample,dC,acluster], log=TRUE))
                    }else{0}
            }, numeric(ndataz)))
        )
        ## pR: rows=clusters, cols= P at grid points
        ##
                ( (if(maincov %in% realCovs){
                    parmList$meanR[asample,maincov,]
                }else{
                    parmList$probI[asample,maincov,] * parmList$sizeI[asample,maincov,] 
                }) %*% W)/colSums(W)
    }
    dim(freqs) <- c(ndataz, nfsamples)
    freqs
}
##
## Gives samples of variance of one covariate conditional on several feature values
samplesvarRgivenX <- function(maincov, X, parmList, nfsamples=NULL, inorder=FALSE){
    realCovs <- dimnames(parmList$meanR)[[2]]
    integerCovs <- dimnames(parmList$probI)[[2]]
    cC <- setdiff(realCovs, maincov)
    dC <- setdiff(integerCovs, maincov)
    ndataz <- nrow(X)
    q <- parmList$q
    nclusters <- ncol(q)
    if(is.numeric(nfsamples)){
        fsubsamples <- seq(1, nrow(q), length.out=nfsamples)
    }else{
        nfsamples <- nrow(q)
        fsubsamples <- seq_len(nfsamples)
    }
    ##
    freqs <- foreach(asample=fsubsamples, .combine=c, .inorder=inorder)%dopar%{
        ## W: rows=clusters, cols=datapoints
        W <- exp(
            log(q[asample,]) +
            t(vapply(seq_len(nclusters), function(acluster){
                ## real covariates
                if(length(cC)>0){
                    colSums(dnorm(t(X[,cC]), mean=parmList$meanR[asample,cC,acluster], sd=1/sqrt(parmList$tauR[asample,cC,acluster]), log=TRUE))
                }else{0} +
                    ## integer covariates
                    if(length(dC)>0){
                        colSums(dbinom(t(X[,dC]), prob=parmList$probI[asample,dC,acluster], size=parmList$sizeI[asample,dC,acluster], log=TRUE))
                    }else{0}
            }, numeric(ndataz)))
        )
        ## pR: rows=clusters, cols= P at grid points
        ##
                ( (if(maincov %in% realCovs){
                    1/parmList$tauR[asample,maincov,]
                }else{
                    (1-parmList$probI[asample,maincov,]) * parmList$probI[asample,maincov,] * parmList$sizeI[asample,maincov,]
                }) %*% W)/colSums(W)
    }
    dim(freqs) <- c(ndataz, nfsamples)
    freqs
}
##
## Gives samples of cumulative frequency distributions of log_RMSD conditional on several feature values
samplescdfRgivenX <- function(maincov, X, parmList, covgrid, nfsamples=NULL, inorder=FALSE){
    realCovs <- dimnames(parmList$meanR)[[2]]
    integerCovs <- dimnames(parmList$probI)[[2]]
    cC <- setdiff(realCovs, maincov)
    dC <- setdiff(integerCovs, maincov)
    ndataz <- nrow(X)
    lcovgrid <- length(covgrid)
    q <- parmList$q
    nclusters <- ncol(q)
    if(is.numeric(nfsamples)){
        fsubsamples <- seq(1, nrow(q), length.out=nfsamples)
    }else{
        nfsamples <- nrow(q)
        fsubsamples <- seq_len(nfsamples)
    }
    ##
    freqs <- foreach(asample=fsubsamples, .combine=c, .inorder=inorder)%dopar%{
        ## W: rows=clusters, cols=datapoints
        W <- exp(
            log(q[asample,]) +
            t(vapply(seq_len(nclusters), function(acluster){
                ## real covariates
                if(length(cC)>0){
                    colSums(dnorm(t(X[,cC]), mean=parmList$meanR[asample,cC,acluster], sd=1/sqrt(parmList$tauR[asample,cC,acluster]), log=TRUE))
                }else{0} +
                    ## integer covariates
                    if(length(dC)>0){
                        colSums(dbinom(t(X[,dC]), prob=parmList$probI[asample,dC,acluster], size=parmList$sizeI[asample,dC,acluster], log=TRUE))
                    }else{0}
            }, numeric(ndataz)))
        )
        W <- matrix(rep(W, lcovgrid), nrow=nclusters) # strings copies of W column-wise
        ## pR: rows=clusters, cols= P at grid points
        pR <- t(vapply(seq_len(nclusters), function(acluster){
            rep(
                if(maincov %in% realCovs){
                    pnorm(q=covgrid, mean=parmList$meanR[asample,maincov,acluster], sd=1/sqrt(parmList$tauR[asample,maincov,acluster]))
                }else{
                    pbinom(q=covgrid, prob=parmList$probI[asample,maincov,acluster], size=parmList$sizeI[asample,maincov,acluster])
                } , each=ndataz)
        }, numeric(lcovgrid * ndataz)))
        ##
        colSums(pR * W)/colSums(W)
    }
    dim(freqs) <- c(ndataz, lcovgrid, nfsamples)
    freqs
}
##
## Gives samples of frequency distributions of log_RMSD conditional on several feature values
samplesfRgivenX <- function(maincov, X, parmList, covgrid, nfsamples=NULL, inorder=FALSE){
    realCovs <- dimnames(parmList$meanR)[[2]]
    integerCovs <- dimnames(parmList$probI)[[2]]
    cC <- setdiff(realCovs, maincov)
    dC <- setdiff(integerCovs, maincov)
    ndataz <- nrow(X)
    lcovgrid <- length(covgrid)
    q <- parmList$q
    nclusters <- ncol(q)
    if(is.numeric(nfsamples)){
        fsubsamples <- seq(1, nrow(q), length.out=nfsamples)
    }else{
        nfsamples <- nrow(q)
        fsubsamples <- seq_len(nfsamples)
    }
    ##
    freqs <- foreach(asample=fsubsamples, .combine=c, .inorder=inorder)%dopar%{
        ## W: rows=clusters, cols=datapoints
        W <- exp(
            log(q[asample,]) +
            t(vapply(seq_len(nclusters), function(acluster){
                ## real covariates
                if(length(cC)>0){
                    colSums(dnorm(t(X[,cC]), mean=parmList$meanR[asample,cC,acluster], sd=1/sqrt(parmList$tauR[asample,cC,acluster]), log=TRUE))
                }else{0} +
                    ## integer covariates
                    if(length(dC)>0){
                        colSums(dbinom(t(X[,dC]), prob=parmList$probI[asample,dC,acluster], size=parmList$sizeI[asample,dC,acluster], log=TRUE))
                    }else{0}
            }, numeric(ndataz)))
        )
        W <- matrix(rep(W, lcovgrid), nrow=nclusters) # strings copies of W column-wise
        ## pR: rows=clusters, cols= P at grid points
        pR <- t(vapply(seq_len(nclusters), function(acluster){
            rep(
                if(maincov %in% realCovs){
                    dnorm(x=covgrid, mean=parmList$meanR[asample,maincov,acluster], sd=1/sqrt(parmList$tauR[asample,maincov,acluster]))
                }else{
                    dbinom(x=covgrid, prob=parmList$probI[asample,maincov,acluster], size=parmList$sizeI[asample,maincov,acluster])
                } , each=ndataz)
        }, numeric(lcovgrid * ndataz)))
        ##
        colSums(pR * W)/colSums(W)
    }
    dim(freqs) <- c(ndataz, lcovgrid, nfsamples)
    freqs
}
##
## Gives samples of frequency distributions of log_RMSD conditional on one feature value
samplesfRgivenX1 <- function(X, parmList, RMSDgrid, nfsamples=NULL){
    realCovs <- dimnames(parmList$meanR)[[2]]
    integerCovs <- dimnames(parmList$probI)[[2]]
    cC <- setdiff(realCovs, 'log_RMSD')
    X <- as.matrix(X)[1,]
    lRMSDgrid <- length(RMSDgrid)
    q <- parmList$q
    nclusters <- ncol(q)
    if(is.numeric(nfsamples)){
        fsubsamples <- round(seq(1, nrow(q), length.out=nfsamples))
    }else{
        nfsamples <- nrow(q)
        fsubsamples <- seq_len(nfsamples)
    }
    
    ##
        ## W: rows=samples, cols=clusters
        W <- exp(
            log(q[fsubsamples,]) +
            vapply(seq_len(nclusters), function(acluster){
                ## real covariates
                colSums(dnorm(X[cC], mean=aperm(parmList$meanR[fsubsamples,cC,acluster,drop=FALSE], c(2,1,3)), sd=1/sqrt(aperm(parmList$tauR[fsubsamples,cC,acluster,drop=FALSE], c(2,1,3))), log=TRUE)) +
                    ## integer covariates
                    colSums(dbinom(X[integerCovs], prob=aperm(parmList$probI[fsubsamples,,acluster,drop=FALSE], c(2,1,3)), size=aperm(parmList$sizeI[fsubsamples,,acluster,drop=FALSE], c(2,1,3)), log=TRUE))
            }, numeric(nfsamples))
        )
        ## W: rows=clusters*gridpoints, cols=samples
        W <- matrix(rep(W, lRMSDgrid), ncol=nfsamples, byrow=TRUE) # strings copies of transposed W row-wise
        dim(W) <- c(nclusters, lRMSDgrid, nfsamples)
        ## pR: rows=clusters, cols= P at grid points
        pR <- vapply(RMSDgrid, function(aRvalue){
           dnorm(x=aRvalue, mean=parmList$meanR[fsubsamples,'log_RMSD',], sd=1/sqrt(parmList$tauR[fsubsamples,'log_RMSD',]))
        }, numeric(nclusters * nfsamples))
    dim(pR) <- c(nfsamples, nclusters, lRMSDgrid)
    pR <- aperm(pR, c(2,3,1))
        ##
    colSums(pR * W)/colSums(W)
}
##
## Calculates the predictive distribution on a grid of log_RMSD conditional on several feature values
pRgivenX <- function(X, parmList, RMSDgrid){
    realCovs <- dimnames(parmList$meanR)[[2]]
    integerCovs <- dimnames(parmList$probI)[[2]]
    cC <- setdiff(realCovs, 'log_RMSD')
    ndataz <- nrow(X)
    q <- parmList$q
    nsamples <- nrow(q)
    nclusters <- ncol(q)
    lRMSDgrid <- length(RMSDgrid)
    ##
    freqs <- foreach(asample=seq_len(nsamples), .combine='+', .inorder=FALSE)%dopar%{
        ## W: rows=clusters, cols=datapoints
        W <- exp(
            log(q[asample,]) +
            t(vapply(seq_len(nclusters), function(acluster){
                ## real covariates
                colSums(dnorm(t(X[,cC]), mean=parmList$meanR[asample,cC,acluster], sd=1/sqrt(parmList$tauR[asample,cC,acluster]), log=TRUE)) +
                    ## integer covariates
                    colSums(dbinom(t(X[,integerCovs]), prob=parmList$probI[asample,,acluster], size=parmList$sizeI[asample,,acluster], log=TRUE))
            }, numeric(ndataz)))
        )
        W <- matrix(rep(W,lRMSDgrid),nrow=nrow(W))
        ##
        pR <- t(vapply(seq_len(nclusters), function(acluster){
            ## real covariates
            dnorm(x=rep(RMSDgrid, each=ndataz), mean=parmList$meanR[asample,'log_RMSD',acluster], sd=1/sqrt(parmList$tauR[asample,'log_RMSD',acluster]))
        }, numeric(ndataz*lRMSDgrid)))
        ##
        colSums(pR * W)/colSums(W)
    }/nsamples
    dim(freqs) <- c(ndataz, lRMSDgrid)
    freqs
}
##
## Gives samples of marginal frequency distributions of a covariate
samplesfX <- function(acov, parmList, acovgrid, nfsamples=100){
    realCovs <- dimnames(parmList$meanR)[[2]]
    integerCovs <- dimnames(parmList$probI)[[2]]
    lacovgrid <- length(acovgrid)
    q <- parmList$q
    nclusters <- ncol(q)
    if(is.numeric(nfsamples)){
        fsubsamples <- seq(1, nrow(q), length.out=nfsamples)
    }else{
        nfsamples <- nrow(q)
        fsubsamples <- seq_len(nfsamples)
    }
    ##
    if(acov %in% realCovs){ ## real covariates
    freqs <- foreach(asample=fsubsamples, .combine=cbind, .inorder=FALSE)%dopar%{
        ## W: rows=clusters, cols=datapoints
        colSums( exp(
            log(q[asample,]) +
            t(vapply(seq_len(nclusters), function(acluster){
                dnorm(acovgrid, mean=parmList$meanR[asample,acov,acluster], sd=1/sqrt(parmList$tauR[asample,acov,acluster]), log=TRUE)
            }, numeric(lacovgrid)))
            ) )
    }
    }else{ ## integer covariates
    freqs <- foreach(asample=fsubsamples, .combine=cbind, .inorder=FALSE)%dopar%{
        ## W: rows=clusters, cols=datapoints
        colSums( exp(
            log(q[asample,]) +
            t(vapply(seq_len(nclusters), function(acluster){
                dbinom(acovgrid, prob=parmList$probI[asample,acov,acluster], size=parmList$sizeI[asample,acov,acluster], log=TRUE)
                            }, numeric(lacovgrid)))
            ) )
    }
    }
    freqs
}







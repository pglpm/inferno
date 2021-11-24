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
    parmNames <- c('q', 'meanC', 'tauC', 'probD', 'sizeD')
    nclusters <- sum(grepl('^q\\[', colnames(mcsamples)))
    nccovs <- sum(grepl('^meanC\\[[^,]*, 1]', colnames(mcsamples)))
    ndcovs <- sum(grepl('^probD\\[[^,]*, 1]', colnames(mcsamples)))
    ##
    parmList <- foreach(var=parmNames)%do%{
        out <- mcsamples[,grepl(paste0('^',var,'\\['), colnames(mcsamples))]
        if(var=='meanC'||var=='tauC'){
            dim(out) <- c(nrow(mcsamples), nccovs, nclusters)
            dimnames(out) <- list(NULL, continuousCovs, NULL)
        } else if(var=='probD'||var=='sizeD'){
            dim(out) <- c(nrow(mcsamples), ndcovs, nclusters)
            dimnames(out) <- list(NULL, discreteCovs, NULL)
        } else if(var=='q'){
            dim(out) <- c(nrow(mcsamples), nclusters)
        }
            out
    }
    names(parmList) <- parmNames
    parmList
}
##
## Construct a list of parameter samples from the raw MCMC samples for the second monitored set
finalstate2list <- function(mcsamples){
    if(!is.vector(mcsamples)){print('ERROR!')}
    parmNames <- c('q', 'meanC', 'tauC', 'probD', 'sizeD', 'C')
    nclusters <- sum(grepl('^q\\[', names(mcsamples)))
    nccovs <- sum(grepl('^meanC\\[[^,]*, 1]', names(mcsamples)))
    ndcovs <- sum(grepl('^probD\\[[^,]*, 1]', names(mcsamples)))
    ##
    parmList <- foreach(var=parmNames)%dopar%{
        out <- mcsamples[grepl(paste0('^',var,'\\['), names(mcsamples))]
        if(var=='meanC'||var=='tauC'){
            dim(out) <- c(nccovs, nclusters)
            dimnames(out) <- list(continuousCovs, NULL)
        } else if(var=='probD'||var=='sizeD'){
            dim(out) <- c(ndcovs, nclusters)
            dimnames(out) <- list(discreteCovs, NULL)
        } # 'q' and 'C' are vectors with no names
            out
    }
    names(parmList) <- parmNames
    parmList
}
## ##
## ## percentile functions for cont. and discr. covars
## qpdf <- function(parmprob){
##     myoptim(rep(0,length(prob)),
##             function(xe){
##                 sum((sapply(xe,function(i){
##                     sum(qs * pnorm(i, mean=Kmeans, sd=Ksds))
##                 })/prob-1)^2)
##                 })$par
## }
## ##
## qpmf <- function(qs, Kprobs, Ksizes, prob){
##     sapply(prob,function(aprob){
##         x <- 0
##         while((sum(qs * pbinom(x, prob=Kprobs, size=Ksizes)) - aprob) < 0){
##             x <- x + 1
##         }
##     x - 1
##     })
## }
##
## Calculates the median and IQR of each sample frequency
calcSampleMQ <- function(parmList, maxD=1000){
    continuousCovs <- dimnames(parmList$meanC)[[2]]
    discreteCovs <- dimnames(parmList$probD)[[2]]
    covNames <- c(continuousCovs, discreteCovs)
    ncovs <- length(covNames)
    q <- parmList$q
    nsamples <- nrow(q)
    ##
    quants <- foreach(asample=seq_len(nsamples), .combine=c)%:%foreach(acov=covNames, .combine=c)%dopar%{
        if(acov %in% continuousCovs){
                mixq <- function(x){sum(q[asample,] * pnorm(x, mean=parmList$meanC[asample,acov,], sd=1/sqrt(parmList$tauC[asample,acov,])))}
                fn <- function(par){
                    (mixq(par[1]) - 0.5)^2 +
                        (mixq(par[2]) - 0.25)^2 +
                        (mixq(par[3]) - 0.75)^2
                }
                out <- myoptim(par=rep(q[asample,]%*%parmList$meanC[asample,acov,],3), fn=fn)$par
            ##     out <- sapply(c(0.005, 0.995), function(border){
            ## optim(0, #rep(q[asample,] %*% parmList$meanC[asample,acov,], 2),
            ##              fn=fn,
            ##              gr = function(x) pracma::grad(fn, x), 
            ##              method = "L-BFGS-B",
            ##              lower = -Inf, upper = Inf,
            ##              control = list(factr = 1e-10, pgtol = 0, maxit = 100))
        }else{
            searchgrid <- 0:maxD
            dq <- colSums(c(q[asample,]) * pbinom(matrix(searchgrid, ncol=length(searchgrid), nrow=ncol(q), byrow=TRUE), prob=parmList$probD[asample,acov,], size=parmList$sizeD[asample,acov,]))
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
## Calculates means and covariances of the sampled frequency distributions
moments12Samples <- function(parmList){
    continuousCovs <- dimnames(parmList$meanC)[[2]]
    discreteCovs <- dimnames(parmList$probD)[[2]]
    covNames <- c(continuousCovs, discreteCovs)
    ncovs <- length(covNames)
    q <- t(parmList$q)
    ##
    meansc <- aperm(parmList$meanC, c(3, 1, 2))
    meansd <- aperm(parmList$probD * parmList$sizeD, c(3, 1, 2))
    clustermeans <- c(meansc, meansd)
    dim(clustermeans) <- c(dim(meansc)[-3], ncovs)
    mixmeans <- colSums(c(q) * clustermeans)
    dimnames(mixmeans) <- list(NULL, paste0('MEAN_', covNames))
    ##
    quadrc <- aperm(1/parmList$tauC + parmList$meanC * parmList$meanC, c(3, 1, 2))
    quadrd <- aperm(parmList$probD * parmList$sizeD *
                    (1 + parmList$probD * (parmList$sizeD - 1)), c(3, 1, 2))
    mixvars <- c(quadrc, quadrd)
    dim(mixvars) <- c(dim(quadrc)[-3], ncovs)
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
    continuousCovs <- dimnames(parmList$meanC)[[2]]
    discreteCovs <- dimnames(parmList$probD)[[2]]
    covNames <- c(continuousCovs, discreteCovs)
    ncovs <- length(covNames)
    q <- parmList$q
    nsamples <- nrow(q)
    ##
    quants <- foreach(asample=seq_len(nsamples), .combine=c)%:%foreach(acov=covNames, .combine=c)%dopar%{
        if(acov %in% continuousCovs){
                mixq <- function(x){sum(q[asample,] * pnorm(x, mean=parmList$meanC[asample,acov,], sd=1/sqrt(parmList$tauC[asample,acov,])))}
                fn <- function(par){(mixq(par[1]) - 0.005)^2 +
                                        (mixq(par[2]) - 0.995)^2}
                out <- myoptim(par=rep(q[asample,]%*%parmList$meanC[asample,acov,],2), fn=fn)$par
            ##     out <- sapply(c(0.005, 0.995), function(border){
            ## optim(0, #rep(q[asample,] %*% parmList$meanC[asample,acov,], 2),
            ##              fn=fn,
            ##              gr = function(x) pracma::grad(fn, x), 
            ##              method = "L-BFGS-B",
            ##              lower = -Inf, upper = Inf,
            ##              control = list(factr = 1e-10, pgtol = 0, maxit = 100))
        }else{
            searchgrid <- 0:max(parmList$sizeD[asample,acov,])
            dq <- colSums(c(q[asample,]) * pbinom(matrix(searchgrid, ncol=length(searchgrid), nrow=ncol(q), byrow=TRUE), prob=parmList$probD[asample,acov,], size=parmList$sizeD[asample,acov,]))
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
## Calculates the probability of some datapoints for the MCMC samples
probValuesSamples <- function(X, parmList){
    continuousCovs <- dimnames(parmList$meanC)[[2]]
    discreteCovs <- dimnames(parmList$probD)[[2]]
    covNames <- c(continuousCovs, discreteCovs)
    ncovs <- length(covNames)
    q <- parmList$q
    ndataz <- nrow(X)
    ##
    (foreach(asample=seq_len(nrow(q)), .combine=cbind, .inorder=TRUE)%dopar%{
        colSums(
            exp(
                log(q[asample,]) +
                t(vapply(seq_len(ncol(q)), function(acluster){
                    ## continuous covariates
                    colSums(dnorm(t(X[,continuousCovs]), mean=parmList$meanC[asample,,acluster,drop=FALSE], sd=1/sqrt(parmList$tauC[asample,,acluster,drop=FALSE]), log=TRUE)) +
                        ## discrete covariates
                    colSums(dbinom(t(X[,discreteCovs]), prob=parmList$probD[asample,,acluster,drop=FALSE], size=parmList$sizeD[asample,,acluster,drop=FALSE], log=TRUE))
    }, numeric(ndataz)))
            )
        )
    })
}
##
## Calculates the probability of the data (likelihood of parameters) for several MCMC samples
llSamples <- function(dat, parmList){
    ndataz <- nrow(dat$X)
    q <- parmList$q
    ##
    foreach(asample=seq_len(nrow(q)), .combine=c, .inorder=TRUE)%dopar%{
        sum( log( colSums(
            exp(
                log(q[asample,]) +
                t(vapply(seq_len(ncol(q)), function(acluster){
                    ## continuous covariates
                    colSums(dnorm(t(dat$X), mean=parmList$meanC[asample,,acluster], sd=1/sqrt(parmList$tauC[asample,,acluster]), log=TRUE)) +
                        ## discrete covariates
                    colSums(dbinom(t(dat$Y), prob=parmList$probD[asample,,acluster], size=parmList$sizeD[asample,,acluster], log=TRUE))
    }, numeric(ndataz)))
            )
        ) ) )
    }
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
##     Nccovs <- dim(X)[2]
##     Ndcovs <- dim(Y)[2]
##     LL <- 0
##     for(adatum in 1:Ndata){
##         clustersum <- log(Q)
##         for(acov in 1:Nccovs){
##             clustersum <- clustersum +
##                 dnorm(x=X[adatum,acov], mean=MeanC[acov,], sd=1/sqrt(TauC[acov,]), log=TRUE)
##         }
##         for(acov in 1:Ndcovs){
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
                    ## continuous covariates
                    colSums(dnorm(t(dat$X), mean=parmList$meanC[asample,,acluster], sd=1/sqrt(parmList$tauC[asample,,acluster]), log=TRUE)) +
                        ## discrete covariates
                    colSums(dbinom(t(dat$Y), prob=parmList$probD[asample,,acluster], size=parmList$sizeD[asample,,acluster], log=TRUE))
    }, numeric(ndataz)))
            )
        )
    }
    if(!log){freqs} else {log(freqs)}
}
## Calculates the MCMC posterior probability of several datapoints
probJointMean <- function(X, parmList){
    continuousCovs <- dimnames(parmList$meanC)[[2]]
    discreteCovs <- dimnames(parmList$probD)[[2]]
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
                    ## continuous covariates
                    colSums(dnorm(t(X[,continuousCovs]), mean=parmList$meanC[asample,,acluster], sd=1/sqrt(parmList$tauC[asample,,acluster]), log=TRUE)) +
                        ## discrete covariates
                    colSums(dbinom(t(X[,discreteCovs]), prob=parmList$probD[asample,,acluster], size=parmList$sizeD[asample,,acluster], log=TRUE))
    }, numeric(ndataz)))
            )
        )
    }/nsamples
}
##
## Produces multidimensional samples from several MCMC distribution-samples
options(doFuture.rng.onMisuse = "ignore")
samplesFsamples <- function(parmList, nxsamples=1000, nfsamples=NULL, seed=149){
    continuousCovs <- dimnames(parmList$meanC)[[2]]
    discreteCovs <- dimnames(parmList$probD)[[2]]
    ncC <- length(continuousCovs)
    ndC <- length(discreteCovs)
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
            ## continuous covariates
            rnorm(n=ncC, mean=parmList$meanC[afsample,continuousCovs,acluster], sd=1/sqrt(parmList$tauC[afsample,continuousCovs,acluster])),
            ## discrete covariates
            rbinom(n=ndC, prob=parmList$probD[afsample,discreteCovs,acluster], size=parmList$sizeD[afsample,discreteCovs,acluster])
        )
    }
    dim(allsamples) <- c(ncC+ndC, nxsamples, nfsamples)
    dimnames(allsamples) <- list(c(continuousCovs,discreteCovs), NULL, NULL)
    allsamples
}
##
## Calculates the predictive probability for several log_RMSD values conditional on several feature values
expeRgivenX <- function(maincov, X, parmList){
    continuousCovs <- dimnames(parmList$meanC)[[2]]
    discreteCovs <- dimnames(parmList$probD)[[2]]
    cC <- setdiff(continuousCovs, maincov)
    dC <- setdiff(discreteCovs, maincov)
    ndataz <- nrow(X)
    q <- parmList$q
    nsamples <- nrow(q)
    ##
    foreach(asample=seq_len(nsamples), .combine='+', .inorder=FALSE)%dopar%{
        ## W: rows=clusters, cols=datapoints
        W <- exp(
            log(q[asample,]) +
            t(vapply(seq_len(ncol(q)), function(acluster){
                ## continuous covariates
                if(length(cC)>0){
                    colSums(dnorm(t(X[,cC]), mean=parmList$meanC[asample,cC,acluster], sd=1/sqrt(parmList$tauC[asample,cC,acluster]), log=TRUE))
                }else{0} +
                    ## discrete covariates
                    if(length(dC)>0){
                        colSums(dbinom(t(X[,dC]), prob=parmList$probD[asample,dC,acluster], size=parmList$sizeD[asample,dC,acluster], log=TRUE))
                    }else{0}
            }, numeric(ndataz)))
        )
        ##
        colSums(
                if(maincov %in% continuousCovs){
                    parmList$meanC[asample,maincov,]
                }else{
                    parmList$probD[asample,maincov,] * parmList$sizeD[asample,maincov,] 
                } * W)/colSums(W)
    }/nsamples
}
##
## Gives samples of mean of one covariate conditional on several feature values
samplesmeanRgivenX <- function(maincov, X, parmList, nfsamples=NULL, inorder=FALSE){
    continuousCovs <- dimnames(parmList$meanC)[[2]]
    discreteCovs <- dimnames(parmList$probD)[[2]]
    cC <- setdiff(continuousCovs, maincov)
    dC <- setdiff(discreteCovs, maincov)
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
                ## continuous covariates
                if(length(cC)>0){
                    colSums(dnorm(t(X[,cC]), mean=parmList$meanC[asample,cC,acluster], sd=1/sqrt(parmList$tauC[asample,cC,acluster]), log=TRUE))
                }else{0} +
                    ## discrete covariates
                    if(length(dC)>0){
                        colSums(dbinom(t(X[,dC]), prob=parmList$probD[asample,dC,acluster], size=parmList$sizeD[asample,dC,acluster], log=TRUE))
                    }else{0}
            }, numeric(ndataz)))
        )
        ## pR: rows=clusters, cols= P at grid points
        ##
                ( (if(maincov %in% continuousCovs){
                    parmList$meanC[asample,maincov,]
                }else{
                    parmList$probD[asample,maincov,] * parmList$sizeD[asample,maincov,] 
                }) %*% W)/colSums(W)
    }
    dim(freqs) <- c(ndataz, nfsamples)
    freqs
}
##
## Gives samples of variance of one covariate conditional on several feature values
samplesvarRgivenX <- function(maincov, X, parmList, nfsamples=NULL, inorder=FALSE){
    continuousCovs <- dimnames(parmList$meanC)[[2]]
    discreteCovs <- dimnames(parmList$probD)[[2]]
    cC <- setdiff(continuousCovs, maincov)
    dC <- setdiff(discreteCovs, maincov)
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
                ## continuous covariates
                if(length(cC)>0){
                    colSums(dnorm(t(X[,cC]), mean=parmList$meanC[asample,cC,acluster], sd=1/sqrt(parmList$tauC[asample,cC,acluster]), log=TRUE))
                }else{0} +
                    ## discrete covariates
                    if(length(dC)>0){
                        colSums(dbinom(t(X[,dC]), prob=parmList$probD[asample,dC,acluster], size=parmList$sizeD[asample,dC,acluster], log=TRUE))
                    }else{0}
            }, numeric(ndataz)))
        )
        ## pR: rows=clusters, cols= P at grid points
        ##
                ( (if(maincov %in% continuousCovs){
                    1/parmList$tauC[asample,maincov,]
                }else{
                    (1-parmList$probD[asample,maincov,]) * parmList$probD[asample,maincov,] * parmList$sizeD[asample,maincov,]
                }) %*% W)/colSums(W)
    }
    dim(freqs) <- c(ndataz, nfsamples)
    freqs
}
##
## Gives samples of cumulative frequency distributions of log_RMSD conditional on several feature values
samplescdfRgivenX <- function(maincov, X, parmList, covgrid, nfsamples=NULL, inorder=FALSE){
    continuousCovs <- dimnames(parmList$meanC)[[2]]
    discreteCovs <- dimnames(parmList$probD)[[2]]
    cC <- setdiff(continuousCovs, maincov)
    dC <- setdiff(discreteCovs, maincov)
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
                ## continuous covariates
                if(length(cC)>0){
                    colSums(dnorm(t(X[,cC]), mean=parmList$meanC[asample,cC,acluster], sd=1/sqrt(parmList$tauC[asample,cC,acluster]), log=TRUE))
                }else{0} +
                    ## discrete covariates
                    if(length(dC)>0){
                        colSums(dbinom(t(X[,dC]), prob=parmList$probD[asample,dC,acluster], size=parmList$sizeD[asample,dC,acluster], log=TRUE))
                    }else{0}
            }, numeric(ndataz)))
        )
        W <- matrix(rep(W, lcovgrid), nrow=nclusters) # strings copies of W column-wise
        ## pR: rows=clusters, cols= P at grid points
        pR <- t(vapply(seq_len(nclusters), function(acluster){
            rep(
                if(maincov %in% continuousCovs){
                    pnorm(q=covgrid, mean=parmList$meanC[asample,maincov,acluster], sd=1/sqrt(parmList$tauC[asample,maincov,acluster]))
                }else{
                    pbinom(q=covgrid, prob=parmList$probD[asample,maincov,acluster], size=parmList$sizeD[asample,maincov,acluster])
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
    continuousCovs <- dimnames(parmList$meanC)[[2]]
    discreteCovs <- dimnames(parmList$probD)[[2]]
    cC <- setdiff(continuousCovs, maincov)
    dC <- setdiff(discreteCovs, maincov)
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
                ## continuous covariates
                if(length(cC)>0){
                    colSums(dnorm(t(X[,cC]), mean=parmList$meanC[asample,cC,acluster], sd=1/sqrt(parmList$tauC[asample,cC,acluster]), log=TRUE))
                }else{0} +
                    ## discrete covariates
                    if(length(dC)>0){
                        colSums(dbinom(t(X[,dC]), prob=parmList$probD[asample,dC,acluster], size=parmList$sizeD[asample,dC,acluster], log=TRUE))
                    }else{0}
            }, numeric(ndataz)))
        )
        W <- matrix(rep(W, lcovgrid), nrow=nclusters) # strings copies of W column-wise
        ## pR: rows=clusters, cols= P at grid points
        pR <- t(vapply(seq_len(nclusters), function(acluster){
            rep(
                if(maincov %in% continuousCovs){
                    dnorm(x=covgrid, mean=parmList$meanC[asample,maincov,acluster], sd=1/sqrt(parmList$tauC[asample,maincov,acluster]))
                }else{
                    dbinom(x=covgrid, prob=parmList$probD[asample,maincov,acluster], size=parmList$sizeD[asample,maincov,acluster])
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
    continuousCovs <- dimnames(parmList$meanC)[[2]]
    discreteCovs <- dimnames(parmList$probD)[[2]]
    cC <- setdiff(continuousCovs, 'log_RMSD')
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
                ## continuous covariates
                colSums(dnorm(X[cC], mean=aperm(parmList$meanC[fsubsamples,cC,acluster,drop=FALSE], c(2,1,3)), sd=1/sqrt(aperm(parmList$tauC[fsubsamples,cC,acluster,drop=FALSE], c(2,1,3))), log=TRUE)) +
                    ## discrete covariates
                    colSums(dbinom(X[discreteCovs], prob=aperm(parmList$probD[fsubsamples,,acluster,drop=FALSE], c(2,1,3)), size=aperm(parmList$sizeD[fsubsamples,,acluster,drop=FALSE], c(2,1,3)), log=TRUE))
            }, numeric(nfsamples))
        )
        ## W: rows=clusters*gridpoints, cols=samples
        W <- matrix(rep(W, lRMSDgrid), ncol=nfsamples, byrow=TRUE) # strings copies of transposed W row-wise
        dim(W) <- c(nclusters, lRMSDgrid, nfsamples)
        ## pR: rows=clusters, cols= P at grid points
        pR <- vapply(RMSDgrid, function(aRvalue){
           dnorm(x=aRvalue, mean=parmList$meanC[fsubsamples,'log_RMSD',], sd=1/sqrt(parmList$tauC[fsubsamples,'log_RMSD',]))
        }, numeric(nclusters * nfsamples))
    dim(pR) <- c(nfsamples, nclusters, lRMSDgrid)
    pR <- aperm(pR, c(2,3,1))
        ##
    colSums(pR * W)/colSums(W)
}
##
## Calculates the predictive distribution on a grid of log_RMSD conditional on several feature values
pRgivenX <- function(X, parmList, RMSDgrid){
    continuousCovs <- dimnames(parmList$meanC)[[2]]
    discreteCovs <- dimnames(parmList$probD)[[2]]
    cC <- setdiff(continuousCovs, 'log_RMSD')
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
                ## continuous covariates
                colSums(dnorm(t(X[,cC]), mean=parmList$meanC[asample,cC,acluster], sd=1/sqrt(parmList$tauC[asample,cC,acluster]), log=TRUE)) +
                    ## discrete covariates
                    colSums(dbinom(t(X[,discreteCovs]), prob=parmList$probD[asample,,acluster], size=parmList$sizeD[asample,,acluster], log=TRUE))
            }, numeric(ndataz)))
        )
        W <- matrix(rep(W,lRMSDgrid),nrow=nrow(W))
        ##
        pR <- t(vapply(seq_len(nclusters), function(acluster){
            ## continuous covariates
            dnorm(x=rep(RMSDgrid, each=ndataz), mean=parmList$meanC[asample,'log_RMSD',acluster], sd=1/sqrt(parmList$tauC[asample,'log_RMSD',acluster]))
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
    continuousCovs <- dimnames(parmList$meanC)[[2]]
    discreteCovs <- dimnames(parmList$probD)[[2]]
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
    if(acov %in% continuousCovs){ ## continuous covariates
    freqs <- foreach(asample=fsubsamples, .combine=cbind, .inorder=FALSE)%dopar%{
        ## W: rows=clusters, cols=datapoints
        colSums( exp(
            log(q[asample,]) +
            t(vapply(seq_len(nclusters), function(acluster){
                dnorm(acovgrid, mean=parmList$meanC[asample,acov,acluster], sd=1/sqrt(parmList$tauC[asample,acov,acluster]), log=TRUE)
            }, numeric(lacovgrid)))
            ) )
    }
    }else{ ## discrete covariates
    freqs <- foreach(asample=fsubsamples, .combine=cbind, .inorder=FALSE)%dopar%{
        ## W: rows=clusters, cols=datapoints
        colSums( exp(
            log(q[asample,]) +
            t(vapply(seq_len(nclusters), function(acluster){
                dbinom(acovgrid, prob=parmList$probD[asample,acov,acluster], size=parmList$sizeD[asample,acov,acluster], log=TRUE)
                            }, numeric(lacovgrid)))
            ) )
    }
    }
    freqs
}

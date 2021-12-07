## Normalize vector
normalize <- function(freqs){freqs/sum(freqs)}
## Normalize rows of matrix
normalizerows <- function(freqs){freqs/rowSums(freqs)}
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
                    colSums(dnbinom(t(dat$Y), prob=parmList$probD[asample,,acluster], size=parmList$sizeD[asample,,acluster], log=TRUE))
    }, numeric(ndataz)))
            )
        ) ) )
    }
}
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
                    colSums(dnbinom(t(dat$Y), prob=parmList$probD[asample,,acluster], size=parmList$sizeD[asample,,acluster], log=TRUE))
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
                    colSums(dnbinom(t(X[,discreteCovs]), prob=parmList$probD[asample,,acluster], size=parmList$sizeD[asample,,acluster], log=TRUE))
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
            rnbinom(n=ndC, prob=parmList$probD[afsample,discreteCovs,acluster], size=parmList$sizeD[afsample,discreteCovs,acluster])
        )
    }
    dim(allsamples) <- c(ncC+ndC, nxsamples, nfsamples)
    dimnames(allsamples) <- list(c(continuousCovs,discreteCovs), NULL, NULL)
    allsamples
}
##
## Calculates the predictive probability for several log_RMSD values conditional on several feature values
expeRgivenX <- function(X, parmList){
    continuousCovs <- dimnames(parmList$meanC)[[2]]
    discreteCovs <- dimnames(parmList$probD)[[2]]
    cC <- setdiff(continuousCovs, 'log_RMSD')
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
                colSums(dnorm(t(X[,cC]), mean=parmList$meanC[asample,cC,acluster], sd=1/sqrt(parmList$tauC[asample,cC,acluster]), log=TRUE)) +
                    ## discrete covariates
                    colSums(dnbinom(t(X[,discreteCovs]), prob=parmList$probD[asample,,acluster], size=parmList$sizeD[asample,,acluster], log=TRUE))
            }, numeric(ndataz)))
        )
        ##
        colSums(parmList$meanC[asample,'log_RMSD',] * W)/colSums(W)
    }/nsamples
}
##
## Gives samples of frequency distributions of log_RMSD conditional on several feature values
samplesfRgivenX <- function(X, parmList, RMSDgrid, nfsamples=NULL){
    continuousCovs <- dimnames(parmList$meanC)[[2]]
    discreteCovs <- dimnames(parmList$probD)[[2]]
    cC <- setdiff(continuousCovs, 'log_RMSD')
    ndataz <- nrow(X)
    lRMSDgrid <- length(RMSDgrid)
    q <- parmList$q
    nclusters <- ncol(q)
    if(is.numeric(nfsamples)){
        fsubsamples <- seq(1, nrow(q), length.out=nfsamples)
    }else{
        nfsamples <- nrow(q)
        fsubsamples <- seq_len(nfsamples)
    }
    ##
    freqs <- foreach(asample=fsubsamples, .combine=c, .inorder=FALSE)%dopar%{
        ## W: rows=clusters, cols=datapoints
        W <- exp(
            log(q[asample,]) +
            t(vapply(seq_len(nclusters), function(acluster){
                ## continuous covariates
                colSums(dnorm(t(X[,cC]), mean=parmList$meanC[asample,cC,acluster], sd=1/sqrt(parmList$tauC[asample,cC,acluster]), log=TRUE)) +
                    ## discrete covariates
                    colSums(dnbinom(t(X[,discreteCovs]), prob=parmList$probD[asample,,acluster], size=parmList$sizeD[asample,,acluster], log=TRUE))
            }, numeric(ndataz)))
        )
        W <- matrix(rep(W, lRMSDgrid), nrow=nclusters) # strings copies of W column-wise
        ## pR: rows=clusters, cols= P at grid points
        pR <- t(vapply(seq_len(nclusters), function(acluster){
            rep( dnorm(x=RMSDgrid, mean=parmList$meanC[asample,'log_RMSD',acluster], sd=1/sqrt(parmList$tauC[asample,'log_RMSD',acluster])) , each=ndataz)
        }, numeric(lRMSDgrid * ndataz)))
        ##
        colSums(pR * W)/colSums(W)
    }
    dim(freqs) <- c(ndataz, lRMSDgrid, nfsamples)
    freqs
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
                    colSums(dnbinom(t(X[,discreteCovs]), prob=parmList$probD[asample,,acluster], size=parmList$sizeD[asample,,acluster], log=TRUE))
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

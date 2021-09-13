## Author: PGL  Porta Mana
## Created: 2021-03-20T10:07:17+0100
## Last-Updated: 2021-09-13T14:55:34+0200
################
## Script for direct regression, continuous RMSD
################

#### Custom setup ####
## Colour-blind friendly palettes, from https://personal.sron.nl/~pault/
library('khroma')
palette(colour('bright')())
## palette(colour('muted')())
library('data.table')
library('ggplot2')
library('ggthemes')
theme_set(theme_bw(base_size=18))
scale_colour_discrete <- scale_colour_bright
#library('cowplot')
library('png')
library('foreach')
library('doFuture')
library('doRNG')
registerDoFuture()
#library('LaplacesDemon') # used for Dirichlet generator
library('ash')
library('extraDistr')
library('PReMiuM')
library('mvtnorm')
options(bitmapType='cairo')
pdff <- function(filename){pdf(file=paste0(filename,'.pdf'),paper='a4r',height=11.7,width=16.5)} # to output in pdf format
pngf <- function(filename,res=300){png(file=paste0(filename,'.png'),height=11.7*1.2,width=16.5,units='in',res=res,pointsize=36)} # to output in pdf format
#### End custom setup ####

#######################################
#### FUNCTION TO CALCULATE MUTUAL INFO FROM FREQUENCY PAIRS
## freqs[,S] = response freqs for stimulus S: one column per stimulus
## assumes all stimuli equally probable
mutualinfo <- function(freqs,base=2){##in bits by default
    freqs1 <- rowSums(freqs)
    freqs2 <- colSums(freqs)
    sum(freqs *
        log2(freqs/outer(freqs1,freqs2)), na.rm=TRUE)/log2(base)
}
condentropy21 <- function(freqs,base=2){##in bits by default
    freqs1 <- rowSums(freqs)
    freqs2 <- colSums(freqs)
    -sum(freqs *
        log2(freqs/outer(freqs1,rep(1,length(freqs2)))), na.rm=TRUE)/log2(base)
}
entropy <- function(freqs,base=2){##in bits by default
    -sum(freqs * log2(freqs), na.rm=TRUE)/log2(base)
}
## function to normalize absolute frequencies
normalize <- function(freqs){freqs/sum(freqs)}
## for rows of frequency distributions
normalizem <- function(freqs){freqs/rowSums(freqs)}
##

## Read and reorganize data
rm(alldata)
alldata <- fread('../data_processed_transformed_rescaled_shuffled.csv', sep=' ')
nameFeatures <- names(alldata)
nSamples <- nrow(alldata)
nFeatures <- ncol(alldata)
##
##
#set.seed(222)
covNames <-  c('log_RMSD',
               'log_mcs_unbonded_polar_sasa',
               'logit_ec_tanimoto_similarity',
               'mcs_NumHeteroAtoms',
               ##'scale_fc_tanimoto_similarity'
               'docked_HeavyAtomCount',
               'mcs_RingCount',
               'docked_NumRotatableBonds'
               )
discreteCovs <- covNames[sapply(covNames, function(x){is.integer(alldata[[x]])})]
continuousCovs <- covNames[sapply(covNames, function(x){is.double(alldata[[x]])})]
covNames <- c(continuousCovs, discreteCovs)
##
##

library('nimble')
##

rm(constants, dat, inits, bayesnet, model, Cmodel, confmodel, mcmcsampler, Cmcmcsampler)
gc()
##
nclusters <- 100
ndata <- 6000 # nSamples = 37969
ncvars <- length(continuousCovs)
ndvars <- length(discreteCovs)
##
pSamples <- nimbleFunction(
    run = function(X=double(2),
                   Y=double(2),
                   q=double(2),
                   meanC=double(3),
                   tauC=double(3),
                   probD=double(3),
                   sizeD=double(3),
                   log=integer(0, default=1)
                   ){
        ##
        returnType(double(2))
        ndataz <- dim(X)[1]
        nsamplesz <- dim(q)[1]
        nclustersz <- dim(q)[2]
        ncvarx <- dim(X)[2]
        ndvarx <- dim(Y)[2]
        pout <- nimMatrix(nrow=ndataz, ncol=nsamplesz, init=FALSE)
        for(idat in 1:ndataz){
            for(isam in 1:nsamplesz){
                sumclusters <- 0
                for(iclu in 1:nclustersz){
                sumclusters <- sumclusters +
                    exp(
                        log(q[isam, iclu]) +
                        sum( dnorm(x=X[idat, 1:ncvarx], mean=meanC[isam, 1:ncvarx, iclu], sd=1/sqrt(tauC[isam, 1:ncvarx, iclu]), log=TRUE)) + 
                        sum( dnbinom(x=Y[idat, 1:ndvarx], prob=probD[isam, 1:ndvarx, iclu], size=sizeD[isam, 1:ndvarx, iclu], log=TRUE))
                    )
                }
                pout[idat, isam] <- sumclusters
            }
        }
        ##
        if(log) return( log(pout))
        else return(pout)
})
CpSamples <- compileNimble(pSamples)
##
pData <- nimbleFunction(
    run = function(X=double(2),
                   Y=double(2),
                   q=double(2),
                   meanC=double(3),
                   tauC=double(3),
                   probD=double(3),
                   sizeD=double(3),
                   log=integer(0, default=1)
                   ){
        ##
        returnType(double(1))
        ndataz <- dim(X)[1]
        nsamplesz <- dim(q)[1]
        nclustersz <- dim(q)[2]
        ncvarx <- dim(X)[2]
        ndvarx <- dim(Y)[2]
        pout <- numeric(length=nsamplesz, init=FALSE)
        for(isam in 1:nsamplesz){
            pouttemp <- 0
            for(idat in 1:ndataz){
                sumclusters <- 0
                for(iclu in 1:nclustersz){
                sumclusters <- sumclusters +
                    exp(
                        log(q[isam, iclu]) +
                        sum( dnorm(x=X[idat, 1:ncvarx], mean=meanC[isam, 1:ncvarx, iclu], sd=1/sqrt(tauC[isam, 1:ncvarx, iclu]), log=TRUE)) + 
                        sum( dnbinom(x=Y[idat, 1:ndvarx], prob=probD[isam, 1:ndvarx, iclu], size=sizeD[isam, 1:ndvarx, iclu], log=TRUE))
                    )
                }
                pouttemp <- pouttemp + log(sumclusters)
            }
            pout[isam] <- pouttemp
        }
        ##
        if(log) return( log(pout))
        else return(pout)
})
CpData <- compileNimble(pData)
## testtime <- Sys.time()
## testpout2 <- CpData(X=dat$X, Y=dat$Y, q=parmList$q, meanC=parmList$meanC, tauC=parmList$tauC, probD=parmList$probD, sizeD=parmList$sizeD, log=T)
## Sys.time() - testtime
##
pData2 <- function(X, Y, q, meanC, tauC, probD, sizeD){
    ndataz <- dim(X)[1]
    nsamplesz <- dim(q)[1]
    nclustersz <- dim(q)[2]
    ncvarx <- dim(X)[2]
    ndvarx <- dim(Y)[2]
    foreach(idat=1:ndataz, .combine=rbind)%:%foreach(isam=1:nsamplesz, .combine=cbind)%dopar%{
        log(sum(exp( sapply(1:nclustersz, function(iclu){
            log(q[isam, iclu]) +
                sum( dnorm(x=X[idat, 1:ncvarx], mean=meanC[isam, 1:ncvarx, iclu], sd=1/sqrt(tauC[isam, 1:ncvarx, iclu]), log=TRUE)) + 
                sum( dnbinom(x=Y[idat, 1:ndvarx], prob=probD[isam, 1:ndvarx, iclu], size=sizeD[isam, 1:ndvarx, iclu], log=TRUE))
        }) )))
    }
}
## testtime <- Sys.time()
## plan(sequential)
## plan(multisession, workers = 6L)
## testpout <- pData2(X=dat$X, Y=dat$Y, q=parmList$q, meanC=parmList$meanC, tauC=parmList$tauC, probD=parmList$probD, sizeD=parmList$sizeD)
## plan(sequential)
## Sys.time() - testtime
##
constants <- list(
    nClusters=nclusters,
    nData=ndata,
    nCvars=ncvars,
    nDvars=ndvars,
    alpha0=rep(1,nclusters)*4/nclusters,
    meanC0=0,
    tauC0=1/3^2,
    shapeC0=1,
    rateC0=1,
    shape1D0=1,
    shape2D0=1,
    shapeD0=1,
    rateD0=1
)
##
dat <- list(
    X=as.matrix(alldata[1:ndata, ..continuousCovs]),
    Y=as.matrix(alldata[1:ndata, ..discreteCovs])
)
##
initsFunction <- function(){
    list( q=rdirch(n=1, alpha=rep(1,nclusters)/nclusters),
    meanC=matrix(rnorm(n=ncvars*nclusters, mean=0, sd=10), nrow=ncvars, ncol=nclusters),
    tauC=matrix(rgamma(n=ncvars*nclusters, shape=1, rate=1), nrow=ncvars, ncol=nclusters),
    probD=matrix(rbeta(n=ndvars*nclusters, shape1=1, shape2=2), nrow=ndvars, ncol=nclusters),
    sizeD=matrix(rgamma(n=ndvars*nclusters, shape=1, rate=1), nrow=ndvars, ncol=nclusters),
    C=rcat(n=ndata, prob=rep(1/nclusters,nclusters))
    )
}
##
inits <- list(
    q=rep(1,nclusters)/nclusters,
    meanC=matrix(0, nrow=length(continuousCovs), ncol=nclusters),
    tauC=matrix(1, nrow=length(continuousCovs), ncol=nclusters),
    probD=matrix(0.5, nrow=length(discreteCovs), ncol=nclusters),
    sizeD=matrix(50, nrow=length(discreteCovs), ncol=nclusters),
    C=rcat(n=ndata, prob=rep(1,nclusters)/nclusters)
)
##
bayesnet <- nimbleCode({
    q[1:nClusters] ~ ddirch(alpha=alpha0[1:nClusters])
    for(i in 1:nClusters){
        for(j in 1:nCvars){
            meanC[j,i] ~ dnorm(mean=meanC0, tau=tauC0)
            tauC[j,i] ~ dgamma(shape=shapeC0, rate=rateC0)
        }
        for(j in 1:nDvars){
            probD[j,i] ~ dbeta(shape1=shape1D0, shape2=shape2D0)
            sizeD[j,i] ~ dgamma(shape=shapeD0, rate=rateD0)
        }
    }
    ##
    for(i in 1:nData){
        C[i] ~ dcat(prob=q[1:nClusters])
    }
    for(i in 1:nData){
        for(j in 1:nCvars){
            X[i,j] ~ dnorm(mean=meanC[j,C[i]], tau=tauC[j,C[i]])
        }
        for(j in 1:nDvars){
            Y[i,j] ~ dnbinom(prob=probD[j,C[i]], size=sizeD[j,C[i]])
        }
    }
})

model <- nimbleModel(code=bayesnet, name='model1', constants=constants, inits=inits, data=dat)
Cmodel <- compileNimble(model, showCompilerOutput=TRUE)

confmodel <- configureMCMC(Cmodel)
## confmodel$removeSamplers(paste0('sizeD'))
## for(i in 1:nclusters){ for(j in 1:length(discreteCovs)){
##                            confmodel$addSampler(target=paste0('sizeD[',j,', ',i,']'), type='slice', control=list(adaptInterval=100))
##                        } }
## print(confmodel)

mcmcsampler <- buildMCMC(confmodel)
Cmcmcsampler <- compileNimble(mcmcsampler, resetFunctions = TRUE)

totaltime <- Sys.time()
## NB: putting all data in one cluster at start leads to slow convergence
mcsamples <- runMCMC(Cmcmcsampler, nburnin=0, niter=2000, thin=1, inits=initsFunction, setSeed=123)
totaltime <- Sys.time() - totaltime
print(totaltime)
## 7 vars, 6000 data, 100 cl, 200 iter: 12.52 mins
saveRDS(mcsamples,file=paste0('_mcsamples',length(covNames),'-d',ndata,'-c',nclusters,'-i',nrow(mcsamples),'.rds'))
##
parmNames <- c('q', 'meanC', 'tauC', 'probD', 'sizeD')
##
parmList <- foreach(var=parmNames)%dopar%{
    out <- mcsamples[,grepl(paste0(var,'\\['), colnames(mcsamples))]
    if(grepl('C', var)){
        dim(out) <- c(nrow(mcsamples), ncvars, nclusters)
    } else if(grepl('D', var)){
        dim(out) <- c(nrow(mcsamples), ndvars, nclusters)
    } else {dim(out) <- c(nrow(mcsamples), nclusters) }
    out
}
names(parmList) <- parmNames
##

lpdat <- CpData(X=dat$X, Y=dat$Y, q=parmList$q, meanC=parmList$meanC, tauC=parmList$tauC, probD=parmList$probD, sizeD=parmList$sizeD, log=T)
##


pdff(paste0('mcsummary',length(covNames),'-d',ndata,'-c',nclusters,'-i',nrow(mcsamples),'.rds'))
matplot(lpdat[1000:2000], type='l', lty=1, col=palette()[2], main='logprobData')
matplot(log(t(apply(parmList$q,1,range))),type='l',lty=1, main='range p-clusters')
##
for(j in c(1,nclusters)){
    vcol <- paste0('q[',j,']')
    matplot(log(mcsamples[,vcol]),type='l',lty=1, main=vcol)
}
for(j in c(1,nclusters)){
    for(i in c(1,ncvars)){
        vcol <- paste0('meanC[',i,', ',j,']')
        matplot(mcsamples[,vcol], type='l', lty=1, main=vcol)
        vcol <- paste0('tauC[',i,', ',j,']')
        matplot(log(mcsamples[,vcol]), type='l', lty=1, main=vcol)
    }
    for(i in c(1,ndvars)){
        vcol <- paste0('probD[',i,', ',j,']')
        matplot(log(mcsamples[,vcol]), type='l', lty=1, main=vcol)
        vcol <- paste0('sizeD[',i,', ',j,']')
        matplot(log(mcsamples[,vcol]), type='l', lty=1, main=vcol)
    }
}
dev.off()
















indq <- grepl('logProb_q\\[', colnames(mcsamples))
matplot(identity(mcsamples[,indq]),type='l',lty=1)

indq <- grepl('meanC\\[1, 1]', colnames(mcsamples)) || grepl('meanC\\[1, 1]', colnames(mcsamples))

totaltime <- Sys.time()
mcsamplesb <- runMCMC(Cmcmcsampler, nburnin=0, niter=2000, thin=1, nchains=2, inits=initsFunction, setSeed=123)
totaltime <- Sys.time() - totaltime
totaltime
saveRDS(mcsamplesb,file=paste0('_mcsampleswburnin2_v',length(covNames),'-d',ndata,'-c',nclusters,'.rds'))

indq <- grepl('q\\[', colnames(mcsamplesb))
matplot(identity(mcsamplesb[,indq][,1]),type='l',lty=1)


#####################################################################
#####################################################################
#####################################################################

testnf <- nimbleFunction(
    run = function(x=double(1),
                   ##meanC=double(2), sdC=double(2),
                   log=integer(0, default=0)){
        returnType(double(1))
        prob <- numeric(length=length(x), init=FALSE)
        ##        prob <- dnorm(x, mean=meanC, sd=sdC)
        prob[1:1] <- x[1:1] * 2
        prob[2:length(x)] <- x[2:length(x)] * 3
        if(log) return(log(prob))
        else return(prob)
        })
Ctestnf <- compileNimble(testnf)

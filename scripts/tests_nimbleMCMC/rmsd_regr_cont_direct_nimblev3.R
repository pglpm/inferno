## Author: PGL  Porta Mana
## Created: 2021-03-20T10:07:17+0100
## Last-Updated: 2021-09-17T09:06:00+0200
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
##normalizem <- function(freqs){freqs/rowSums(freqs)}
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
##
constants <- list(
    nClusters=nclusters,
    nData=ndata,
    nCvars=ncvars,
    nDvars=ndvars
    ## alpha0=rep(1,nclusters)*4/nclusters,
    ## meanC0=0,
    ## tauC0=1/3^2,
    ## shapeC0=1,
    ## rateC0=1,
    ## shape1D0=1,
    ## shape2D0=1,
    ## shapeD0=1,
    ## rateD0=1
)
##
dat <- list(
    X=as.matrix(alldata[1:ndata, ..continuousCovs]),
    Y=as.matrix(alldata[1:ndata, ..discreteCovs])
)
##
inits <- list( alpha0=rep(10/nclusters, nclusters),
         meanC0=0,
         tauC0=1/(0.5^2),
         shapeC0=0.9, #3, #7, #0.5, #0.6,
         rateC0=0.1, #2, #6, #0.03, #0.1,
         shape1D0=1,
         shape2D0=1,
         shapeD0=2,
         rateD0=1/2,
         ##
                  q=rdirch(n=1, alpha=rep(10/nclusters, nclusters)),
    meanC=matrix(rnorm(n=ncvars*nclusters, mean=0, sd=0.5), nrow=ncvars, ncol=nclusters),
    tauC=matrix(rgamma(n=ncvars*nclusters, shape=0.9, rate=0.1), nrow=ncvars, ncol=nclusters),
    probD=matrix(rbeta(n=ndvars*nclusters, shape1=1, shape2=1), nrow=ndvars, ncol=nclusters),
    sizeD=matrix(rgamma(n=ndvars*nclusters, shape=2, rate=1/2), nrow=ndvars, ncol=nclusters),
    C=rcat(n=ndata, prob=rep(1/nclusters,nclusters))
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
##
mcmcsampler <- buildMCMC(confmodel)
Cmcmcsampler <- compileNimble(mcmcsampler, resetFunctions = TRUE)


version <- 3
initsFunction <- function(){
    list( alpha0=rep(10/nclusters, nclusters),
         meanC0=0,
         tauC0=1/(0.5^2),
         shapeC0=0.9, #3, #7, #0.5, #0.6,
         rateC0=0.1, #2, #6, #0.03, #0.1,
         shape1D0=1,
         shape2D0=1,
         shapeD0=2,
         rateD0=1/2,
         ##
         q=rdirch(n=1, alpha=rep(10/nclusters, nclusters)),
    meanC=matrix(rnorm(n=ncvars*nclusters, mean=0, sd=0.5), nrow=ncvars, ncol=nclusters),
    tauC=matrix(rgamma(n=ncvars*nclusters, shape=0.9, rate=0.1), nrow=ncvars, ncol=nclusters),
    probD=matrix(rbeta(n=ndvars*nclusters, shape1=1, shape2=1), nrow=ndvars, ncol=nclusters),
    sizeD=matrix(rgamma(n=ndvars*nclusters, shape=2, rate=1/2), nrow=ndvars, ncol=nclusters),
    C=rcat(n=ndata, prob=rep(1/nclusters,nclusters))
         )
}
##
totaltime <- Sys.time()
## NB: putting all data in one cluster at start leads to slow convergence
mcsamples <- runMCMC(Cmcmcsampler, nburnin=2000, niter=10000, thin=10, inits=initsFunction, setSeed=149)
#Cmcmcsampler$run(niter=10000, thin=10, reset=FALSE, resetMV=TRUE)
totaltime <- Sys.time() - totaltime
print(totaltime)
mcsamples <- as.matrix(Cmcmcsampler$mvSamples)
## 7 vars, 6000 data, 100 cl, 200 iter: 12.52 mins
## 7 vars, 6000 data, 100 cl, 2000 iter: 1.88 hours
saveRDS(mcsamples,file=paste0('_mcsamples-run',version,'-v',length(covNames),'-d',ndata,'-c',nclusters,'-i',nrow(mcsamples),'.rds'))
save(model,Cmodel,confmodel,mcmcsampler,Cmcmcsampler, file=paste0('_model-',version,'-v',length(covNames),'-d',ndata,'-c',nclusters,'-i',nrow(mcsamples),'.RData'))
##
parmNames <- c('q', 'meanC', 'tauC', 'probD', 'sizeD')
parmList <- foreach(var=parmNames)%dopar%{
    out <- mcsamples[,grepl(paste0(var,'\\['), colnames(mcsamples))]
    if(grepl('C', var)){
        dim(out) <- c(nrow(mcsamples), ncvars, nclusters)
        dimnames(out) <- list(NULL, continuousCovs, NULL)
    } else if(grepl('D', var)){
        dim(out) <- c(nrow(mcsamples), ndvars, nclusters)
        dimnames(out) <- list(NULL, discreteCovs, NULL)
    } else {dim(out) <- c(nrow(mcsamples), nclusters) }
    out
}
names(parmList) <- parmNames
##
timecount <- Sys.time()
plan(sequential)
plan(multisession, workers = 6L)
loglikelihood <- llSamples(dat=dat, parmList=parmList)
plan(sequential)
print(Sys.time()-timecount)
##
pdff(paste0('mcsummary-run',version,'-v',length(covNames),'-d',ndata,'-c',nclusters,'-i',nrow(mcsamples),'.rds'))
matplot(loglikelihood, type='l', lty=1, col=palette()[2], main='logprobData')
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
##
save.image(file=paste0('_nimbleoutput-run',version,'.RData'))





subsamplep <- round(seq(1, nrow(parmList$q), length.out=100))
redparmList <- list(
    q=parmList$q[subsamplep,],
    meanC=parmList$meanC[subsamplep,,],
    tauC=parmList$tauC[subsamplep,,],
    probD=parmList$probD[subsamplep,,],
    sizeD=parmList$sizeD[subsamplep,,]
)

nxsamples <- 1000
##
timecount <- Sys.time()
plan(sequential)
plan(multisession, workers = 6L)
xsamples <- samplesFsamples(varNames=covNames, parmList=redparmList, nxsamples=nxsamples)
plan(sequential)
print(Sys.time()-timecount)

##
plotvarRanges <- plotvarQs <- xlim <- ylim <- list()
for(var in covNames){
    plotvarQs[[var]] <- quantile(alldata[[var]], prob=c(0.05,0.95))
    plotvarRanges[[var]] <- thisrange <- range(alldata[[var]])
    xlim[[var]] <- c( min(thisrange, quantile(c(xsamples[var,,]),prob=0.05)),
                     max(thisrange, quantile(c(xsamples[var,,]),prob=0.95)) )
}
##
subsamplex <- 1:nxsamples #round(seq(1, dim(xsamples)[2], length.out=1000))
pdff(paste0('post_samplesvars2D'))#'.pdf'), height=11.7, width=16.5)
par(mfrow = c(2, 3))
for(addvar in setdiff(covNames, 'log_RMSD')){
    matplot(x=alldata[['log_RMSD']][seq(1,nrow(alldata),length.out=1000)],
            y=alldata[[addvar]][seq(1,nrow(alldata),length.out=1000)],
            xlim=xlim[['log_RMSD']],
            ylim=xlim[[addvar]],
            xlab='log_RMSD',
            ylab=addvar,
            type='p', pch=1, cex=0.2, lwd=1, col=palette()[2])
        matlines(x=c(rep(plotvarRanges[['log_RMSD']], each=2), plotvarRanges[['log_RMSD']][1]),
             y=c(plotvarRanges[[addvar]], rev(plotvarRanges[[addvar]]), plotvarRanges[[addvar]][1]),
                 lwd=2, col=paste0(palette()[2],'88'))
        matlines(x=c(rep(plotvarQs[['log_RMSD']], each=2), plotvarQs[['log_RMSD']][1]),
             y=c(plotvarQs[[addvar]], rev(plotvarQs[[addvar]]), plotvarQs[[addvar]][1]),
                 lwd=2, col=paste0(palette()[4],'88'))
}
for(asample in 1:length(subsamplep)){
par(mfrow = c(2, 3))
for(addvar in setdiff(covNames, 'log_RMSD')){
    matplot(x=xsamples['log_RMSD',, asample],
            y=xsamples[addvar,, asample],
            xlim=xlim[['log_RMSD']],
            ylim=xlim[[addvar]],
            xlab='log_RMSD',
            ylab=addvar,
            type='p', pch=1, cex=0.2, lwd=1, col=palette()[1])
    matlines(x=c(rep(plotvarRanges[['log_RMSD']], each=2), plotvarRanges[['log_RMSD']][1]),
             y=c(plotvarRanges[[addvar]], rev(plotvarRanges[[addvar]]), plotvarRanges[[addvar]][1]),
             lwd=2, col=paste0(palette()[2],'88'))
            matlines(x=c(rep(plotvarQs[['log_RMSD']], each=2), plotvarQs[['log_RMSD']][1]),
             y=c(plotvarQs[[addvar]], rev(plotvarQs[[addvar]]), plotvarQs[[addvar]][1]),
                 lwd=2, col=paste0(palette()[4],'88'))
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


###################################################
###### not used ###################################
###################################################

probJointSamples2 <- nimbleFunction(
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
CprobJointSamples2 <- compileNimble(probJointSamples2)
##


#####################################################################
#####################################################################
#####################################################################

testnf <- nimbleFunction(
    run = function(x=double(2), y=double(1)
                   ##meanC=double(2), sdC=double(2),
                   ){
        returnType(double(2))
        nr <- dim(x)[1]
        nc <- dim(x)[2]
        out <- nimMatrix(0, nrow=nr, ncol=nc)
        out <- out + x
        return(out[1:nr,1:nc] + y[1:nr])
        })
Ctestnf <- compileNimble(testnf)

#!/usr/bin/env Rscript

## Author: PGL  Porta Mana
## Created: 2021-03-20T10:07:17+0100
## Last-Updated: 2021-10-01T14:51:15+0200
################
## Script for direct regression, continuous RMSD
################

.libPaths(c("/cluster/home/pglpm/R",.libPaths()))
#### Custom setup ####
## Colour-blind friendly palettes, from https://personal.sron.nl/~pault/
library('khroma')
palette(colour('bright')())
## palette(colour('muted')())
library('data.table')
## library('ggplot2')
## library('ggthemes')
## theme_set(theme_bw(base_size=18))
scale_colour_discrete <- scale_colour_bright
#library('cowplot')
library('png')
library('foreach')
## library('doFuture')
## library('doRNG')
## registerDoFuture()
## library('LaplacesDemon') # used for Dirichlet generator
## library('ash')
## library('extraDistr')
## library('PReMiuM')
## library('mvtnorm')
options(bitmapType='cairo')
pdff <- function(filename){pdf(file=paste0(filename,'.pdf'),paper='a4r',height=11.7,width=16.5)} # to output in pdf format
pngf <- function(filename,res=300){png(file=paste0(filename,'.png'),height=11.7*1.2,width=16.5,units='in',res=res,pointsize=36)} # to output in pdf format
library('nimble')
library('coda')
#### End custom setup ####

#######################################
## Read and reorganize data
alldata <- fread('data_id_processed_transformed_rescaled_shuffled.csv', sep=',')
nameFeatures <- names(which(sapply(alldata, is.numeric)==TRUE))
nSamples <- nrow(alldata)
nFeatures <- length(nameFeatures)
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

gc()
##
nclusters <- 100
ndata <- 6000 # nSamples = 37969
nccovs <- length(continuousCovs)
ndcovs <- length(discreteCovs)
meansccovs <- apply(alldata[1:ndata,..continuousCovs],2,mean)
varsccovs <- apply(alldata[1:ndata,..continuousCovs],2,function(x)var(x, na.rm=T))
meansdcovs <- apply(alldata[1:ndata,..discreteCovs],2,mean)
varsdcovs <- apply(alldata[1:ndata,..discreteCovs],2,function(x)var(x, na.rm=T))
maxdcovs <- apply(alldata[1:ndata,..discreteCovs],2,max)
tauQccovs <- sapply(continuousCovs, function(acov){
    resu <- list(par=c(0.1,0.1))
    for(i in 1:1000){
        resu <- optim(par=resu$par,
                      fn=function(parms){
                          (pinvgamma(varsccovs[acov]/100, shape=parms[1], scale=parms[2]) - 0.05)^2 +
                              (pinvgamma(varsccovs[acov]*100, shape=parms[1], scale=parms[2]) - 0.95)^2
                      }
                    , control=list(maxit=1000000)
                      )
    }
    resu$par})
##
constants <- list(
    nClusters=nclusters,
    nData=ndata,
    nCcovs=nccovs,
    nDcovs=ndcovs
)
##
dat <- list(
    X=as.matrix(alldata[1:ndata, ..continuousCovs]),
    Y=as.matrix(alldata[1:ndata, ..discreteCovs])
)
##
inits <- list(
    alphaK=rep(5/nclusters, nclusters),
    meanCmean=meansccovs,
    meanCtau=0.5/varsccovs,
    tauCshape=tauQccovs[1,],
    tauCrate=tauQccovs[2,],
    sizeDpar1=1/(1+maxdcovs),
    sizeDpar2=1+0*maxdcovs,
    ##
    q=rep(1/nclusters, nclusters),
    meanC=matrix(meansccovs, nrow=nccovs, ncol=nclusters),
    tauC=matrix(1/varsccovs, nrow=nccovs, ncol=nclusters),
    probD=matrix(meansdcovs/maxdcovs, nrow=ndcovs, ncol=nclusters),
    sizeD=matrix(maxdcovs, nrow=ndcovs, ncol=nclusters),
    ## meanC=matrix(rnorm(n=nccovs*nclusters, mean=meansccovs, sd=sqrt(varsccovs)), nrow=nccovs, ncol=nclusters),
    ## tauC=matrix(rgamma(n=nccovs*nclusters, shape=tauQccovs[1,], rate=tauQccovs[2,]), nrow=nccovs, ncol=nclusters),
    ## probD=matrix(rbeta(n=ndcovs*nclusters, shape1=1, shape2=1), nrow=ndcovs, ncol=nclusters),
    ## sizeD=matrix(rnbinom(n=ndcovs*nclusters, prob=1/(1+maxdcovs), size=maxdcovs), nrow=ndcovs, ncol=nclusters),
    C=rep(1,ndata) # rcat(n=ndata, prob=rep(1/nclusters,nclusters))
         )
##
bayesnet <- nimbleCode({
    q[1:nClusters] ~ ddirch(alpha=alphaK[1:nClusters])
    for(acluster in 1:nClusters){
        for(acov in 1:nCcovs){
            meanC[acov,acluster] ~ dnorm(mean=meanCmean[acov], tau=meanCtau[acov])
            tauC[acov,acluster] ~ dgamma(shape=tauCshape[acov], rate=tauCrate[acov])
        }
        for(acov in 1:nDcovs){
            probD[acov,acluster] ~ dbeta(shape1=1, shape2=1)
            sizeD[acov,acluster] ~ dnbinom(prob=sizeDpar1[acov], size=sizeDpar2[acov])
        }
    }
    if(posterior){
        for(adatum in 1:nData){
            C[adatum] ~ dcat(prob=q[1:nClusters])
        }            ##
        for(adatum in 1:nData){
            for(acov in 1:nCcovs){
                X[adatum,acov] ~ dnorm(mean=meanC[acov,C[adatum]], tau=tauC[acov,C[adatum]])
            }
            for(acov in 1:nDcovs){
                Y[adatum,acov] ~ dbinom(prob=probD[acov,C[adatum]], size=sizeD[acov,C[adatum]])
            }
        }
    }
})
##

posterior <- TRUE
if(posterior){
    model <- nimbleModel(code=bayesnet, name='model1', constants=constants, inits=inits, data=dat, check=FALSE)
}else{
    model <- nimbleModel(code=bayesnet, name='model1', constants=constants, inits=inits, data=list(), check=FALSE)
    }
Cmodel <- compileNimble(model)
gc()
##
confmodel <- configureMCMC(Cmodel, monitors=c('q','meanC', 'tauC', 'probD', 'sizeD')) #, control=list(adaptive=FALSE))
## confmodel$removeSamplers(paste0('sizeD'))
## for(i in 1:nclusters){
##     confmodel$addSampler(target=paste0('sizeD[',1:ndcovs,', ',i,']'), type='AF_slice', control=list(sliceAdaptFactorInterval=100))
## }
## print(confmodel)
##
## samplerConfList <- confmodel$getSamplers()

mcmcsampler <- buildMCMC(confmodel)
Cmcmcsampler <- compileNimble(mcmcsampler, resetFunctions = TRUE)
gc()

##
source('functions_rmsdregr_nimble_binom.R')
initsFunction <- function(){
inits <- list(
    alphaK=rep(5/nclusters, nclusters),
    meanCmean=meansccovs,
    meanCtau=0.5/varsccovs,
    tauCshape=tauQccovs[1,],
    tauCrate=tauQccovs[2,],
    sizeDpar1=1/(1+maxdcovs),
    sizeDpar2=1+0*maxdcovs,
    ##
    q=rep(1/nclusters, nclusters),
    meanC=matrix(meansccovs, nrow=nccovs, ncol=nclusters),
    tauC=matrix(1/varsccovs, nrow=nccovs, ncol=nclusters),
    probD=matrix(meansdcovs/maxdcovs, nrow=ndcovs, ncol=nclusters),
    sizeD=matrix(maxdcovs, nrow=ndcovs, ncol=nclusters),
    ## meanC=matrix(rnorm(n=nccovs*nclusters, mean=meansccovs, sd=sqrt(varsccovs)), nrow=nccovs, ncol=nclusters),
    ## tauC=matrix(rgamma(n=nccovs*nclusters, shape=tauQccovs[1,], rate=tauQccovs[2,]), nrow=nccovs, ncol=nclusters),
    ## probD=matrix(rbeta(n=ndcovs*nclusters, shape1=1, shape2=1), nrow=ndcovs, ncol=nclusters),
    ## sizeD=matrix(rnbinom(n=ndcovs*nclusters, prob=1/(1+maxdcovs), size=maxdcovs), nrow=ndcovs, ncol=nclusters),
    C=rep(1,ndata) # rcat(n=ndata, prob=rep(1/nclusters,nclusters))
         )
}

##
version <- 'post1slurm'
gc()
totalruntime <- Sys.time()
mcsamples <- runMCMC(Cmcmcsampler, nburnin=0, niter=10, thin=1, inits=initsFunction, setSeed=149)
## Cmcmcsampler$run(niter=1000, thin=1, reset=FALSE, resetMV=FALSE)
## mcsamples <- as.matrix(Cmcmcsampler$mvSamples)
totalruntime <- Sys.time() - totalruntime
print(totalruntime)
## 7 vars, 6000 data, 100 cl, 1000 iter, slice: 1.15 h
## 7 vars, 6000 data, 100 cl, 4000 iter, slice: 4.66 h
## 7 vars, 6000 data, 100 cl, 4000 iter, slice: 4.58 h
## 7 vars, 6000 data, 100 cl, 6000 iter, slice: 7.58 h
##
saveRDS(mcsamples,file=paste0('_mcsamples-R',version,'-V',length(covNames),'-D',ndata,'-K',nclusters,'-I',nrow(mcsamples),'.rds'))
## save(model,Cmodel,confmodel,mcmcsampler,Cmcmcsampler, file=paste0('_model-',version,'-v',length(covNames),'-d',ndata,'-c',nclusters,'-i',nrow(mcsamples),'.RData'))
##
parmNames <- c('q', 'meanC', 'tauC', 'probD', 'sizeD')
parmList <- foreach(var=parmNames)%do%{
    out <- mcsamples[,grepl(paste0(var,'\\['), colnames(mcsamples))]
    if(grepl('C', var)){
        dim(out) <- c(nrow(mcsamples), nccovs, nclusters)
        dimnames(out) <- list(NULL, continuousCovs, NULL)
    } else if(grepl('D', var)){
        dim(out) <- c(nrow(mcsamples), ndcovs, nclusters)
        dimnames(out) <- list(NULL, discreteCovs, NULL)
    } else {dim(out) <- c(nrow(mcsamples), nclusters) }
    out
}
names(parmList) <- parmNames
##
ess <- effectiveSize(as.mcmc(mcsamples))
print(summary(ess))
## ##
## if(posterior){
##     timecount <- Sys.time()
##     plan(sequential)
##     plan(multisession, workers = 6L)
##     loglikelihood <- llSamples(dat=dat, parmList=parmList)
##     plan(sequential)
##     print(Sys.time()-timecount)
##     lless <- effectiveSize(as.mcmc(loglikelihood))
## ##    print(lless)
##     ##
##     pdff(paste0('mcsummary-R',version,'-V',length(covNames),'-D',ndata,'-K',nclusters,'-I',nrow(mcsamples)))
##     matplot(1:2,1:2, type='l', lty=1, col='white')
##     legend('topleft',legend=c(paste0('time: ',signif(totalruntime,3)), paste0('LL ESS: ', signif(lless,3)), paste(c('ESS ',names(summary(ess))),collapse=' '), paste(c('ESS ',signif(summary(ess),3)),collapse=' ')),bty='n',cex=1)
##     matplot(loglikelihood, type='l', lty=1, col=palette()[2], main='logprobData', ylab='logprobData')
##     matplot(apply((parmList$q),1,function(x){mean(log(x[x>0]))}),type='l',lty=1, main='mean log-q', ylab='mean log-q')
##     matplot(log(apply((parmList$q),1,function(x){sd(log(x[x>0]))})), type='l',lty=1, main='log-SD log-q', ylab='log-SD log-q')
##     ##
##     for(i in continuousCovs){
##         matplot(rowMeans((parmList$meanC[,i,])),type='l',lty=1, main=paste0('mean means ', i), ylab=paste0('mean means ', i))
##         matplot(log(apply(parmList$meanC[,i,],1,sd)), type='l',lty=1, main=paste0('log-SD means ', i), ylab=paste0('log-SD means ', i))
##         matplot(rowMeans(log(parmList$tauC[,i,])),type='l',lty=1, main=paste0('mean taus ', i), ylab=paste0('mean log-taus ', i))
##         matplot(log(apply(log(parmList$tauC[,i,]),1,sd)),type='l',lty=1, main=paste0('log-SD taus ', i), ylab=paste0('log-SD log-taus ', i))
##     }
##     for(i in discreteCovs){
##         matplot(rowMeans(log(parmList$probD[,i,])),type='l',lty=1, main=paste0('mean probs ', i), ylab=paste0('mean log-probs ', i))
##         matplot(log(apply(log(parmList$probD[,i,]),1,sd)),type='l',lty=1, main=paste0('log-SD probs ', i), ylab=paste0('log-SD log-probs ', i))
##         matplot(rowMeans((parmList$sizeD[,i,])),type='l',lty=1, main=paste0('mean sizes ', i), ylab=paste0('mean sizes ', i))
##         matplot(log(apply(parmList$sizeD[,i,],1,sd)),type='l',lty=1, main=paste0('log-SD sizes ', i), ylab=paste0('log-SD sizes ', i))
##     }
##     dev.off()
## }
##
##save.image(file=paste0('_nimbleoutput-run',version,'.RData'))


## ##
## nxsamples <- 1000
## ##
## timecount <- Sys.time()
## plan(sequential)
## plan(multisession, workers = 6L)
## xsamples <- samplesFsamples(parmList=parmList, nxsamples=nxsamples, nfsamples=NULL)
## plan(sequential)
## print(Sys.time()-timecount)
## ##
## plotvarRanges <- plotvarQs <- xlim <- ylim <- list()
## for(var in covNames){
##     plotvarQs[[var]] <- quantile(alldata[[var]], prob=c(0.05,0.95))
##     plotvarRanges[[var]] <- thisrange <- range(alldata[[var]])
##     xlim[[var]] <- c( min(thisrange, quantile(xsamples[var,,],prob=0.05)),
##                      max(thisrange, quantile(xsamples[var,,],prob=0.95)) )
## }
## ##
## subsamplep <- round(seq(1, dim(xsamples)[3], length.out=100))
## subsamplex <- round(seq(1, dim(xsamples)[2], length.out=1000))
## ##
## pdff(paste0('hypersamplesvars2D-R',version,'-V',length(covNames),'-D',ndata,'-K',nclusters,'-I',nrow(mcsamples)))
## par(mfrow = c(2, 3))
## for(addvar in setdiff(covNames, 'log_RMSD')){
##     matplot(x=alldata[['log_RMSD']][subsamplex],
##             y=alldata[[addvar]][subsamplex],
##             xlim=xlim[['log_RMSD']],
##             ylim=xlim[[addvar]],
##             xlab='log_RMSD',
##             ylab=addvar,
##             type='p', pch=1, cex=0.2, lwd=1, col=palette()[2])
##         matlines(x=c(rep(plotvarRanges[['log_RMSD']], each=2), plotvarRanges[['log_RMSD']][1]),
##              y=c(plotvarRanges[[addvar]], rev(plotvarRanges[[addvar]]), plotvarRanges[[addvar]][1]),
##                  lwd=2, col=paste0(palette()[2],'88'))
##         matlines(x=c(rep(plotvarQs[['log_RMSD']], each=2), plotvarQs[['log_RMSD']][1]),
##              y=c(plotvarQs[[addvar]], rev(plotvarQs[[addvar]]), plotvarQs[[addvar]][1]),
##                  lwd=2, col=paste0(palette()[4],'88'))
## }
## for(asample in subsamplep){
## par(mfrow = c(2, 3))
## for(addvar in setdiff(covNames, 'log_RMSD')){
##     matplot(x=xsamples['log_RMSD', subsamplex, asample][subsamplex],
##             y=xsamples[addvar, subsamplex, asample][subsamplex],
##             xlim=xlim[['log_RMSD']],
##             ylim=xlim[[addvar]],
##             xlab='log_RMSD',
##             ylab=addvar,
##             type='p', pch=1, cex=0.2, lwd=1, col=palette()[1])
##     matlines(x=c(rep(plotvarRanges[['log_RMSD']], each=2), plotvarRanges[['log_RMSD']][1]),
##              y=c(plotvarRanges[[addvar]], rev(plotvarRanges[[addvar]]), plotvarRanges[[addvar]][1]),
##              lwd=2, col=paste0(palette()[2],'88'))
##             matlines(x=c(rep(plotvarQs[['log_RMSD']], each=2), plotvarQs[['log_RMSD']][1]),
##              y=c(plotvarQs[[addvar]], rev(plotvarQs[[addvar]]), plotvarQs[[addvar]][1]),
##                  lwd=2, col=paste0(palette()[4],'88'))
## }
## }
## dev.off()




## #######################################################################


## subsamplep <- round(seq(1, nrow(parmList$q), length.out=100))
## redparmList <- list(
##     q=parmList$q[subsamplep,],
##     meanC=parmList$meanC[subsamplep,,],
##     tauC=parmList$tauC[subsamplep,,],
##     probD=parmList$probD[subsamplep,,],
##     sizeD=parmList$sizeD[subsamplep,,]
## )

## nxsamples <- 1000
## ##
## timecount <- Sys.time()
## plan(sequential)
## plan(multisession, workers = 6L)
## xsamples <- samplesFsamples(varNames=covNames, parmList=redparmList, nxsamples=nxsamples)
## plan(sequential)
## print(Sys.time()-timecount)

## ##
## plotvarRanges <- plotvarQs <- xlim <- ylim <- list()
## for(var in covNames){
##     plotvarQs[[var]] <- quantile(alldata[[var]], prob=c(0.05,0.95))
##     plotvarRanges[[var]] <- thisrange <- range(alldata[[var]])
##     xlim[[var]] <- c( min(thisrange, quantile(c(xsamples[var,,]),prob=0.05)),
##                      max(thisrange, quantile(c(xsamples[var,,]),prob=0.95)) )
## }
## ##
## subsamplex <- 1:nxsamples #round(seq(1, dim(xsamples)[2], length.out=1000))
## pdff(paste0('post_samplesvars2D'))#'.pdf'), height=11.7, width=16.5)
## par(mfrow = c(2, 3))
## for(addvar in setdiff(covNames, 'log_RMSD')){
##     matplot(x=alldata[['log_RMSD']][seq(1,nrow(alldata),length.out=1000)],
##             y=alldata[[addvar]][seq(1,nrow(alldata),length.out=1000)],
##             xlim=xlim[['log_RMSD']],
##             ylim=xlim[[addvar]],
##             xlab='log_RMSD',
##             ylab=addvar,
##             type='p', pch=1, cex=0.2, lwd=1, col=palette()[2])
##         matlines(x=c(rep(plotvarRanges[['log_RMSD']], each=2), plotvarRanges[['log_RMSD']][1]),
##              y=c(plotvarRanges[[addvar]], rev(plotvarRanges[[addvar]]), plotvarRanges[[addvar]][1]),
##                  lwd=2, col=paste0(palette()[2],'88'))
##         matlines(x=c(rep(plotvarQs[['log_RMSD']], each=2), plotvarQs[['log_RMSD']][1]),
##              y=c(plotvarQs[[addvar]], rev(plotvarQs[[addvar]]), plotvarQs[[addvar]][1]),
##                  lwd=2, col=paste0(palette()[4],'88'))
## }
## for(asample in 1:length(subsamplep)){
## par(mfrow = c(2, 3))
## for(addvar in setdiff(covNames, 'log_RMSD')){
##     matplot(x=xsamples['log_RMSD',, asample],
##             y=xsamples[addvar,, asample],
##             xlim=xlim[['log_RMSD']],
##             ylim=xlim[[addvar]],
##             xlab='log_RMSD',
##             ylab=addvar,
##             type='p', pch=1, cex=0.2, lwd=1, col=palette()[1])
##     matlines(x=c(rep(plotvarRanges[['log_RMSD']], each=2), plotvarRanges[['log_RMSD']][1]),
##              y=c(plotvarRanges[[addvar]], rev(plotvarRanges[[addvar]]), plotvarRanges[[addvar]][1]),
##              lwd=2, col=paste0(palette()[2],'88'))
##             matlines(x=c(rep(plotvarQs[['log_RMSD']], each=2), plotvarQs[['log_RMSD']][1]),
##              y=c(plotvarQs[[addvar]], rev(plotvarQs[[addvar]]), plotvarQs[[addvar]][1]),
##                  lwd=2, col=paste0(palette()[4],'88'))
## }
## }
## dev.off()

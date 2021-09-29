## Author: PGL  Porta Mana
## Created: 2021-03-20T10:07:17+0100
## Last-Updated: 2021-09-28T12:21:20+0200
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
library('nimble')
library('coda')
#### End custom setup ####

#######################################
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

rm(constants, dat, inits, bayesnet, model, Cmodel, confmodel, mcmcsampler, Cmcmcsampler)
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
    q=rdirch(n=1, alpha=rep(5/nclusters, nclusters)),
    meanC=matrix(rnorm(n=nccovs*nclusters, mean=meansccovs, sd=sqrt(varsccovs)), nrow=nccovs, ncol=nclusters),
    tauC=matrix(rgamma(n=nccovs*nclusters, shape=tauQccovs[1,], rate=tauQccovs[2,]), nrow=nccovs, ncol=nclusters),
    probD=matrix(rbeta(n=ndcovs*nclusters, shape1=1, shape2=1), nrow=ndcovs, ncol=nclusters),
    sizeD=matrix(rnbinom(n=ndcovs*nclusters, prob=1/(1+maxdcovs), size=maxdcovs), nrow=ndcovs, ncol=nclusters),
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
    model <- nimbleModel(code=bayesnet, name='model1', constants=constants, inits=inits, data=dat)
}else{
    model <- nimbleModel(code=bayesnet, name='model1', constants=constants, inits=inits, data=list())
    }
Cmodel <- compileNimble(model, showCompilerOutput=TRUE)
gc()

##
confmodel <- configureMCMC(Cmodel, monitors=c('q','meanC', 'tauC', 'probD', 'sizeD')) #, control=list(adaptive=FALSE))
confmodel$removeSamplers(paste0('sizeD'))
for(i in 1:nclusters){
    confmodel$addSampler(target=paste0('sizeD[',1:ndcovs,', ',i,']'), type='AF_slice', control=list(maxContractionsWarning=FALSE))
}
print(confmodel)
##
## samplerConfList <- confmodel$getSamplers()
##
mcmcsampler <- buildMCMC(confmodel)
Cmcmcsampler <- compileNimble(mcmcsampler, resetFunctions = TRUE)
gc()
##
source('functions_rmsdregr_nimble_binom.R')
initsFunction <- function(){
inits <- list(
    alphaK=rep(5/nclusters, nclusters),
    meanCmean=meansccovs,
    meanCtau=1/(4*varsccovs),
    tauCshape=tauQccovs[1,],
    tauCrate=tauQccovs[2,],
    sizeDpar1=1/(1+maxdcovs),
    sizeDpar2=1+0*maxdcovs,
    ##
    q=rdirch(n=1, alpha=rep(5/nclusters, nclusters)),
    meanC=matrix(rnorm(n=nccovs*nclusters, mean=meansccovs, sd=sqrt(4*varsccovs)), nrow=nccovs, ncol=nclusters),
    tauC=matrix(rgamma(n=nccovs*nclusters, shape=tauQccovs[1,], rate=tauQccovs[2,]), nrow=nccovs, ncol=nclusters),
    probD=matrix(rbeta(n=ndcovs*nclusters, shape1=1, shape2=1), nrow=ndcovs, ncol=nclusters),
    sizeD=matrix(rnbinom(n=ndcovs*nclusters, prob=1/(1+maxdcovs), size=maxdcovs), nrow=ndcovs, ncol=nclusters),
    C=rep(1,ndata) # rcat(n=ndata, prob=rep(1/nclusters,nclusters))
         )
}
##
version <- 'post6AS'
totaltime <- Sys.time()
## runMCMC(Cmcmcsampler, nburnin=1, niter=3, thin=1, inits=initsFunction, setSeed=149)
Cmcmcsampler$run(niter=5000, thin=1, reset=TRUE)
## Cmcmcsampler$run(niter=1000, thin=1, reset=FALSE, resetMV=TRUE)
totaltime <- Sys.time() - totaltime
print(totaltime)
##
mcsamples <- as.matrix(Cmcmcsampler$mvSamples)
## 7 vars, 6000 data, 100 cl, 5000 iter, AFslice: 2.19 days
saveRDS(mcsamples,file=paste0('_mcsamples-R',version,'-V',length(covNames),'-D',ndata,'-K',nclusters,'-I',nrow(mcsamples),'.rds'))
## save(model,Cmodel,confmodel,mcmcsampler,Cmcmcsampler, file=paste0('_model-',version,'-v',length(covNames),'-d',ndata,'-c',nclusters,'-i',nrow(mcsamples),'.RData'))
##
parmNames <- c('q', 'meanC', 'tauC', 'probD', 'sizeD')
parmList <- foreach(var=parmNames)%dopar%{
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

##
if(posterior){
    timecount <- Sys.time()
    plan(sequential)
    plan(multisession, workers = 6L)
    loglikelihood <- llSamples(dat=dat, parmList=parmList)
    plan(sequential)
    print(Sys.time()-timecount)
    ##
    pdff(paste0('mcsummary-R',version,'-V',length(covNames),'-D',ndata,'-K',nclusters,'-I',nrow(mcsamples)))
    matplot(loglikelihood, type='l', lty=1, col=palette()[2], main='logprobData', ylab='logprobData')
    matplot(apply((parmList$q),1,function(x){mean(log(x[x>0]))}),type='l',lty=1, main='mean log-q', ylab='mean log-q')
    matplot(log(apply((parmList$q),1,function(x){sd(log(x[x>0]))})), type='l',lty=1, main='log-SD log-q', ylab='log-SD log-q')
    ##
    for(i in 1:nccovs){
        matplot(rowMeans((parmList$meanC[,i,])),type='l',lty=1, main=paste0('mean means ', i), ylab=paste0('mean means ', i))
        matplot(log(apply(parmList$meanC[,i,],1,sd)), type='l',lty=1, main=paste0('log-SD means ', i), ylab=paste0('log-SD means ', i))
        matplot(rowMeans(log(parmList$tauC[,i,])),type='l',lty=1, main=paste0('mean taus ', i), ylab=paste0('mean taus ', i))
        matplot(log(apply(log(parmList$tauC[,i,]),1,sd)),type='l',lty=1, main=paste0('log-SD taus ', i), ylab=paste0('log-SD taus ', i))
    }
    for(i in 1:ndcovs){
        matplot(rowMeans(log(parmList$probD[,i,])),type='l',lty=1, main=paste0('mean probs ', i), ylab=paste0('mean probs ', i))
        matplot(log(apply(log(parmList$probD[,i,]),1,sd)),type='l',lty=1, main=paste0('log-SD probs ', i), ylab=paste0('log-SD probs ', i))
        matplot(rowMeans((parmList$sizeD[,i,])),type='l',lty=1, main=paste0('mean sizes ', i), ylab=paste0('mean sizes ', i))
        matplot(log(apply(parmList$sizeD[,i,],1,sd)),type='l',lty=1, main=paste0('log-SD sizes ', i), ylab=paste0('log-SD sizes ', i))
    }
    dev.off()
}
##
##save.image(file=paste0('_nimbleoutput-run',version,'.RData'))


##
nxsamples <- 1000
##
timecount <- Sys.time()
plan(sequential)
plan(multisession, workers = 6L)
xsamples <- samplesFsamples(parmList=parmList, nxsamples=nxsamples, nfsamples=NULL)
plan(sequential)
print(Sys.time()-timecount)
##
plotvarRanges <- plotvarQs <- xlim <- ylim <- list()
for(var in covNames){
    plotvarQs[[var]] <- quantile(alldata[[var]], prob=c(0.05,0.95))
    plotvarRanges[[var]] <- thisrange <- range(alldata[[var]])
    xlim[[var]] <- c( min(thisrange, quantile(xsamples[var,,],prob=0.05)),
                     max(thisrange, quantile(xsamples[var,,],prob=0.95)) )
}
##
subsamplep <- round(seq(1, dim(xsamples)[3], length.out=100))
subsamplex <- round(seq(1, dim(xsamples)[2], length.out=1000))
##
pdff(paste0('hypersamplesvars2D-R',version,'-V',length(covNames),'-D',ndata,'-K',nclusters,'-I',nrow(mcsamples)))
par(mfrow = c(2, 3))
for(addvar in setdiff(covNames, 'log_RMSD')){
    matplot(x=alldata[['log_RMSD']][subsamplex],
            y=alldata[[addvar]][subsamplex],
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
for(asample in subsamplep){
par(mfrow = c(2, 3))
for(addvar in setdiff(covNames, 'log_RMSD')){
    matplot(x=xsamples['log_RMSD', subsamplex, asample][subsamplex],
            y=xsamples[addvar, subsamplex, asample][subsamplex],
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




#######################################################################


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

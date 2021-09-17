## Author: PGL  Porta Mana
## Created: 2021-03-20T10:07:17+0100
## Last-Updated: 2021-09-17T11:51:54+0200
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
library('nimble')
##

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
normalizerows <- function(freqs){freqs/rowSums(freqs)}
## for columns of frequency distributions
normalizecols <- function(freqs){t(t(freqs)/colSums(freqs))}
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


rm(constants, dat, inits, bayesnet, model, Cmodel, confmodel, mcmcsampler, Cmcmcsampler)
gc()
##
nclusters <- 100
ndata <- 6000 # nSamples = 37969
ncvars <- length(continuousCovs)
ndvars <- length(discreteCovs)
##
constants <- list(
    nClusters=nclusters,
    nCvars=ncvars,
    nDvars=ndvars
)
##
dat <- list(
    X=as.matrix(alldata[1:ndata, ..continuousCovs]),
    Y=as.matrix(alldata[1:ndata, ..discreteCovs])
)
##
inits <- list(
    qalpha <- rinvgamma(n=1, shape=0.5, scale=0.5),
    meanCmean <- rnorm(n=ncvars, mean=0, sd=1),
    meanCtau <- rgamma(n=ncvars, shape=0.5, rate=0.5),
    tauCshape <- rinvgamma(n=ncvars, shape=0.5, scale=0.5),
    tauCrate <- rgamma(n=ncvars, shape=0.5, rate=0.5),
    probDshape1 <- rinvgamma(n=ndvars, shape=0.5, scale=0.5),
    probDshape2 <- rinvgamma(n=ndvars, shape=0.5, scale=0.5),
    sizeDshape <- rinvgamma(n=ndvars, shape=0.5, scale=0.5),
    sizeDrate <- rgamma(n=ndvars, shape=0.5, rate=0.5),
    ##
    q=rdirch(n=1, alpha=rep(1, nclusters)),
    meanC=matrix(rnorm(n=ncvars*nclusters, mean=0, sd=1), nrow=ncvars, ncol=nclusters),
    tauC=matrix(rgamma(n=ncvars*nclusters, shape=0.5, rate=0.5), nrow=ncvars, ncol=nclusters),
    probD=matrix(rbeta(n=ndvars*nclusters, shape1=1, shape2=1), nrow=ndvars, ncol=nclusters),
    sizeD=matrix(rgamma(n=ndvars*nclusters, shape=0.5, rate=0.5), nrow=ndvars, ncol=nclusters),
    C=rcat(n=ndata, prob=rep(1/nclusters,nclusters))
)
##
priorbayesnet <- nimbleCode({
    qalpha ~ dinvgamma(shape=0.5, scale=0.5)
    alphad[1:nClusters] <- qalpha
    q[1:nClusters] ~ ddirch(alpha=alphad[1:nClusters])
    for(i in 1:nClusters){
        for(j in 1:nCvars){
            meanC[j,i] ~ dnorm(mean=meanCmean[j], tau=meanCtau[j])
            tauC[j,i] ~ dgamma(shape=tauCshape[j], rate=tauCrate[j])
        }
        for(j in 1:nDvars){
            probD[j,i] ~ dbeta(shape1=probDshape1[j], shape2=probDshape2[j])
            sizeD[j,i] ~ dgamma(shape=sizeDshape[j], rate=sizeDrate[j])
        }
    }
    for(j in 1:nCvars){
        meanCmean[j] ~ dnorm(mean=0, sd=1)
        meanCtau[j] ~ dgamma(shape=0.5, rate=0.5)
        tauCshape[j] ~ dinvgamma(shape=0.5, scale=0.5)
        tauCrate[j] ~ dgamma(shape=0.5, rate=0.5)
    }
    for(j in 1:nDvars){
        probDshape1[j] ~ dinvgamma(shape=0.5, scale=0.5)
        probDshape2[j] ~ dinvgamma(shape=0.5, scale=0.5)
        sizeDshape[j] ~ dinvgamma(shape=0.5, scale=0.5)
        sizeDrate[j] ~ dgamma(shape=0.5, rate=0.5)
    }
    ##
    ## for(i in 1:nData){
    ##     C[i] ~ dcat(prob=q[1:nClusters])
    ## }
    ## for(i in 1:nData){
    ##     for(j in 1:nCvars){
    ##         X[i,j] ~ dnorm(mean=meanC[j,C[i]], tau=tauC[j,C[i]])
    ##     }
    ##     for(j in 1:nDvars){
    ##         Y[i,j] ~ dnbinom(prob=probD[j,C[i]], size=sizeD[j,C[i]])
    ##     }
    ## }
})
##
source('functions_rmsdregr_nimble.R')


priormodel <- nimbleModel(code=priorbayesnet, name='priormodel1', constants=constants, inits=list(), data=list())
Cpriormodel <- compileNimble(priormodel, showCompilerOutput=TRUE)
##
confpriormodel <- configureMCMC(Cpriormodel, monitors=c('q','meanC', 'tauC', 'probD', 'sizeD'))
## confpriormodel$removeSamplers(paste0('sizeD'))
## for(i in 1:nclusters){ for(j in 1:length(discreteCovs)){
##                            confpriormodel$addSampler(target=paste0('sizeD[',j,', ',i,']'), type='slice', control=list(adaptInterval=100))
##                        } }
## print(confpriormodel)
##
priormcmcsampler <- buildMCMC(confpriormodel)
Cpriormcmcsampler <- compileNimble(priormcmcsampler, resetFunctions = TRUE)


initsFunction <- function(){
inits <- list(
    qalpha=rinvgamma(n=1, shape=0.5, scale=0.5),
    meanCmean=rnorm(n=ncvars, mean=0, sd=1),
    meanCtau=rgamma(n=ncvars, shape=0.5, rate=0.5),
    tauCshape=rinvgamma(n=ncvars, shape=0.5, scale=0.5),
    tauCrate=rgamma(n=ncvars, shape=0.5, rate=0.5),
    probDshape1=rinvgamma(n=ndvars, shape=0.5, scale=0.5),
    probDshape2=rinvgamma(n=ndvars, shape=0.5, scale=0.5),
    sizeDshape=rinvgamma(n=ndvars, shape=0.5, scale=0.5),
    sizeDrate=rgamma(n=ndvars, shape=0.5, rate=0.5),
    ##
    q=rdirch(n=1, alpha=rep(1, nclusters)),
    meanC=matrix(rnorm(n=ncvars*nclusters, mean=0, sd=1), nrow=ncvars, ncol=nclusters),
    tauC=matrix(rgamma(n=ncvars*nclusters, shape=0.5, rate=0.5), nrow=ncvars, ncol=nclusters),
    probD=matrix(rbeta(n=ndvars*nclusters, shape1=1, shape2=1), nrow=ndvars, ncol=nclusters),
    sizeD=matrix(rgamma(n=ndvars*nclusters, shape=0.5, rate=0.5), nrow=ndvars, ncol=nclusters)
)
}
##
totaltime <- Sys.time()
## NB: putting all data in one cluster at start leads to slow convergence
priormcsamples <- runMCMC(Cpriormcmcsampler, nburnin=0, niter=101, thin=1, inits=initsFunction, setSeed=123)
## Cpriormcmcsampler$run(niter=2000, thin=1, reset=TRUE, setSeed=123, init=initsFunction)
## priormcsamples <- as.matrix(Cpriormcmcsampler$mvSamples)
totaltime <- Sys.time() - totaltime
print(totaltime)
##
## Cpriormcmcsampler$run(niter=2000, thin=1, reset=FALSE, resetMV=TRUE)
## priormcsamples <- as.matrix(Cmcmcsampler$mvSamples)
## 7 vars, 6000 data, 100 cl, 200 iter: 12.52 mins
## 7 vars, 6000 data, 100 cl, 2000 iter: 1.88 hours
## saveRDS(priormcsamples,file=paste0('_priormcsamples',length(covNames),'-c',nclusters,'-i',nrow(priormcsamples),'.rds'))
## save(priormodel,Cpriormodel,confpriormodel,priormcmcsampler,Cpriormcmcsampler, file=paste0('_priormodel',length(covNames),'-c',nclusters,'-i',nrow(priormcsamples),'.RData'))
##
parmNames <- c('q', 'meanC', 'tauC', 'probD', 'sizeD')
parmList <- foreach(var=parmNames)%dopar%{
    out <- priormcsamples[,grepl(paste0(var,'\\['), colnames(priormcsamples))]
    if(grepl('C', var)){
        dim(out) <- c(nrow(priormcsamples), ncvars, nclusters)
        dimnames(out) <- list(NULL, continuousCovs, NULL)
    } else if(grepl('D', var)){
        dim(out) <- c(nrow(priormcsamples), ndvars, nclusters)
        dimnames(out) <- list(NULL, discreteCovs, NULL)
    } else {dim(out) <- c(nrow(priormcsamples), nclusters) }
    out
}
names(parmList) <- parmNames

##
nxsamples <- 1000
##
timecount <- Sys.time()
plan(sequential)
plan(multisession, workers = 6L)
xsamples <- samplesFsamples(parmList=parmList, nxsamples=nxsamples)
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
pdff(paste0('priorsamplesvars2D_run2'))#'.pdf'), height=11.7, width=16.5)
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

matplot(x=1:2,y=1:2)










subsamplep <- round(seq(1, nrow(parmList$q), length.out=100))
redparmList <- list(
    q=parmList$q[subsamplep,],
    meanC=parmList$meanC[subsamplep,,],
    tauC=parmList$tauC[subsamplep,,],
    probD=parmList$probD[subsamplep,,],
    sizeD=parmList$sizeD[subsamplep,,]
)
##
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
pdff(paste0('prior_run2_samplesvars2D'))#'.pdf'), height=11.7, width=16.5)
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














    ## pdff('samplesvars2D')
    for(asample in subsamplep){
        matplot(x=xsamples[1,subsamplex,asample],y=xsamples[2,subsamplex,asample], type='p', pch=1, lwd=1, xlab=rownames(xsamples)[1], ylab=rownames(xsamples)[2], xlim=xlim, ylim=ylim)
    }
    dev.off()
}






## rm(alldata)
## alldata <- fread('../processed_data_scaled.csv', sep=' ')
plotVars <- c(
    "log_RMSD"
    ,"log_mcs_unbonded_polar_sasa"
    ## ,"logit_ec_tanimoto_similarity",
    ## ,"mcs_NumHeteroAtoms",
    ## ,"docked_HeavyAtomCount",
    ## ,"mcs_RingCount",
    ## ,"docked_NumRotatableBonds"
  )
##
ngridpoints <- 100
plotgrids <- lapply(plotVars,function(var){
    varrange <- summary(alldata[[var]])
    varrange <- c(diff(varrange[c(4,1)]), diff(varrange[c(4,6)]))*1+varrange[4]
    seq(varrange[1], varrange[2], length.out=ngridpoints)
    ##
    ##seq(-10, 10, length.out=ngridpoints)
})
names(plotgrids) <- plotVars
##
##gridpoints2 <- mesh(plotgrids[[1]], plotgrids[[2]])
gridpoints <- cbind(rep(plotgrids[[1]],each=length(plotgrids[[2]])),
                    rep(plotgrids[[2]],length(plotgrids[[1]])))
colnames(gridpoints) <- plotVars

##
plan(sequential)
plan(multisession, workers = 6L)
zvalues <- Fsamples(X=gridpoints, parmList=parmList)
plan(sequential)
origz <- zvalues
dim(zvalues) <- c(ngridpoints, ngridpoints, ncol(zvalues))

##
indir <- ''
pdff(paste0(indir, 'samplesvars'))
subsample <- round(seq(1, ncol(zvalues), length.out=20))
for(asample in subsample){
persp(plotgrids[[2]],plotgrids[[1]],zvalues[,,asample],zlim=c(0,max(zvalues[,,subsample])),ticktype='detailed',theta = 45, phi = 15,xlab=plotVars[2],ylab=plotVars[1],zlab='freq. density')
}
dev.off()




q1 <- parmList$q[1,]
meanC1 <- parmList$meanC[1,,]
tauC1 <- parmList$tauC[1,,]
probD1 <- parmList$probD[1,,]
sizeD1 <- parmList$sizeD[1,,]

xx <- gridpoints[1,]

testd <- foreach(asample=1:2, .combine=cbind)%:%foreach(i=1:nrow(gridpoints), .combine=c)%do%{
    xx <- gridpoints[i,]
    sum(parmList$q[asample,] * exp(rowSums(t(dnorm(xx, mean=parmList$meanC[asample,plotVars,], sd=1/sqrt(parmList$tauC[asample,plotVars,]), log=TRUE)))))
    }

testd2 <- exp(rowSums(t(dnorm(xx, mean=meanC1[plotVars,], sd=1/sqrt(tauC1[plotVars,]), log=TRUE))))

testd <- dnorm(xx, mean=meanC1[plotVars,], sd=1/sqrt(tauC1[plotVars,]))

prod(dnorm(xx, mean=meanC1[plotVars,3], sd=1/sqrt(tauC1[plotVars,3])))

foreach(var=plotVars, .combine=c)%do%{dnorm(xx[var], mean=meanC1[var,3], sd=1/sqrt(tauC1[var,3]))}


timecount <- Sys.time()
plan(sequential)
plan(multisession, workers = 6L)
loglikelihood <- llSamples(dat=dat, parmList=parmList)
plan(sequential)
print(Sys.time()-timecount)
##
pdff(paste0('mcsummary2',length(covNames),'-d',ndata,'-c',nclusters,'-i',nrow(mcsamples),'.rds'))
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

## sampleX <- nimbleFunction(
##     run = function( niCparms=integer(0),
##                    niDparms=integer(0),
##                    q=double(2),
##                    meanC=double(3),
##                    tauC=double(3),
##                    probD=double(3),
##                    sizeD=double(3)
##                    ){
##         ##
##         returnType(double(2))
##         nsamplesz <- dim(q)[1]
##         nclustersz <- dim(q)[2]
##         niCparmsd <- niCparms+niDparms
##         Xout <- nimMatrix(nrow=niCparmsd, ncol=nsamplesz, init=FALSE)
##         for(asample in 1:nsamplesz){
##             acluster <- rcat(n=1, prob=q[asample,1:nclustersz])
##             if(niCparms > 0){
##                 Xout[1:niCparms, asample] <- rnorm(n=niCparms, mean=meanC[asample, 1:niCparms, acluster], sd=1/sqrt(tauC[asample, 1:niCparms, acluster]))
##                 }
##             if(niDparms > 0){
##                 niCparms1 <- niCparms+1
##                 Xout[niCparms1:niCparmsd, asample] <- rnbinom(n=niDparms, prob=probD[asample, 1:niDparms, acluster], size=sizeD[asample, 1:niDparms, acluster])
##             }
##         }
##         return(Xout)
## })
## CsampleX <- compileNimble(sampleX, showCompilerOutput = TRUE)
## assign('CsampleX', CsampleX, envir = .GlobalEnv)
## assign('sampleX', sampleX, envir = .GlobalEnv)
## ##

## tests <- CsampleX(niCparms=1, niDparms=0, q=parmList$q, meanC=parmList$meanC[,1,,drop=FALSE], tauC=parmList$tauC[,1,,drop=FALSE], probD=parmList$probD[,1,,drop=FALSE], sizeD=parmList$sizeD[,1,,drop=FALSE])

## samplesFsamples2 <- function(varNames, parmList, nxsamples){
##     cC <- varNames[varNames %in% continuousCovs]
##     ncC <- length(cC)
##     dC <- varNames[varNames %in% discreteCovs]
##     ndC <- length(dC)
##     q <- parmList$q
##     if(ncC > 0){
##     meanC <- parmList$meanC[,cC,,drop=FALSE]
##     tauC <- parmList$tauC[,cC,,drop=FALSE]
##     } else { meanC <- tauC <- array(0, dim=rep(2,3)) }
##     if(ndC > 0){
##     probD <- parmList$probD[,dC,,drop=FALSE]
##     sizeD <- parmList$sizeD[,dC,,drop=FALSE]
##     } else { probD <- sizeD <- array(0, dim=rep(2,3)) }
##     ##
##     ## rng <- RNGseq( npsamples * nxsamples * (ncC+ndC), seed)
##     allsamples <- foreach(asample=seq_len(nxsamples), .combine=cbind, .inorder=FALSE)%dorng%{ 
## ##        CsampleX(niCparms=1, niDparms=0, q=parmList$q, meanC=parmList$meanC[,1,,drop=FALSE], tauC=parmList$tauC[,1,,drop=FALSE], probD=parmList$probD[,1,,drop=FALSE], sizeD=parmList$sizeD[,1,,drop=FALSE])
##          CsampleX(niCparms=ncC, niDparms=ndC, q=q, meanC=meanC, tauC=tauC, probD=probD, sizeD=sizeD)
##     }
##     dim(allsamples) <- c(ncC+ndC, nxsamples, nrow(q))
##     dimnames(allsamples) <- list(c(cC,dC), NULL, NULL)
##     allsamples
## }
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

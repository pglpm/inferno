## Author: PGL  Porta Mana
## Created: 2021-03-20T10:07:17+0100
## Last-Updated: 2021-09-17T14:04:25+0200
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
    sizeDscale <- rgamma(n=ndvars, shape=0.5, rate=0.5),
    ##
    q=rdirch(n=1, alpha=rep(1, nclusters)),
    meanC=matrix(rnorm(n=ncvars*nclusters, mean=0, sd=1), nrow=ncvars, ncol=nclusters),
    tauC=matrix(rgamma(n=ncvars*nclusters, shape=0.5, rate=0.5), nrow=ncvars, ncol=nclusters),
    probD=matrix(rbeta(n=ndvars*nclusters, shape1=1, shape2=1), nrow=ndvars, ncol=nclusters),
    sizeD=matrix(rinvgamma(n=ndvars*nclusters, shape=0.5, scale=0.5), nrow=ndvars, ncol=nclusters)
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
            sizeD[j,i] ~ dinvgamma(shape=sizeDshape[j], scale=sizeDscale[j])
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
        sizeDscale[j] ~ dgamma(shape=0.5, rate=0.5)
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
##
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
##
initsFunction <- function(){
    list(
    qalpha=rinvgamma(n=1, shape=0.5, scale=0.5),
    meanCmean=rnorm(n=ncvars, mean=0, sd=1),
    meanCtau=rgamma(n=ncvars, shape=0.5, rate=0.5),
    tauCshape=rinvgamma(n=ncvars, shape=0.5, scale=0.5),
    tauCrate=rgamma(n=ncvars, shape=0.5, rate=0.5),
    probDshape1=rinvgamma(n=ndvars, shape=0.5, scale=0.5),
    probDshape2=rinvgamma(n=ndvars, shape=0.5, scale=0.5),
    sizeDshape=rinvgamma(n=ndvars, shape=0.5, scale=0.5),
    sizeDscale=rgamma(n=ndvars, shape=0.5, rate=1),
    ##
    q=rdirch(n=1, alpha=rep(1, nclusters)),
    meanC=matrix(rnorm(n=ncvars*nclusters, mean=0, sd=1), nrow=ncvars, ncol=nclusters),
    tauC=matrix(rgamma(n=ncvars*nclusters, shape=0.5, rate=0.5), nrow=ncvars, ncol=nclusters),
    probD=matrix(rbeta(n=ndvars*nclusters, shape1=1, shape2=1), nrow=ndvars, ncol=nclusters),
    sizeD=matrix(rinvgamma(n=ndvars*nclusters, shape=0.5, scale=0.5), nrow=ndvars, ncol=nclusters)
)
}
##
totaltime <- Sys.time()
## NB: putting all data in one cluster at start leads to slow convergence
priormcsamples <- runMCMC(Cpriormcmcsampler, nburnin=1, niter=101, thin=1, inits=initsFunction, setSeed=123)
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
pdff(paste0('hyerpriorsamplesvars2D_run2'))#'.pdf'), height=11.7, width=16.5)
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

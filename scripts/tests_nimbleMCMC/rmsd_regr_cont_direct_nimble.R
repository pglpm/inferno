## Author: PGL  Porta Mana
## Created: 2021-03-20T10:07:17+0100
## Last-Updated: 2021-09-12T17:29:26+0200
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
rm(constants, dat, inits, bayesnet, model, Cmodel, confmodel, mcmcsampler, Cmcmcsampler)
gc()

library('nimble')
##
nclusters <- 100
ndata <- 6000 # nSamples = 37969
##
constants <- list(
    nClusters=nclusters,
    nData=ndata,
    nCvars=length(continuousCovs),
    nDvars=length(discreteCovs),
    alpha0=rep(1,nclusters)*4/nclusters,
    meanC0=0,
    tauC0=1/10^2,
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
            Y[i,j] ~ dnegbin(prob=probD[j,C[i]], size=sizeD[j,C[i]])
        }
    }
})

model <- nimbleModel(code=bayesnet, name='model1', constants=constants, inits=inits, data=dat)

Cmodel <- compileNimble(model, showCompilerOutput=TRUE)

confmodel <- configureMCMC(Cmodel)

confmodel$removeSamplers(paste0('sizeD'))
for(i in 1:nclusters){ for(j in 1:length(discreteCovs)){
        confmodel$addSampler(target=paste0('sizeD[',j,', ',i,']'), type='slice', control=list(adaptInterval=100))
                       } }
confmodel$addMonitors('logProb_X')

mcmcsampler <- buildMCMC(confmodel)
Cmcmcsampler <- compileNimble(mcmcsampler, resetFunctions = TRUE)

totaltime <- Sys.time()
mcsamples <- runMCMC(Cmcmcsampler, nburnin=1000, niter=2000, thin=1, setSeed=123)
totaltime <- Sys.time() - totaltime
totaltime
## 7 vars, 1000 data, 100 cl: 38 min
## 7 vars, 2000 data, 100 cl: 38.37 min
## 7 vars, 4000 data, 100 cl: 1.26 hours
saveRDS(mcsamples,file=paste0('_mcsamples_v',length(covNames),'-d',ndata,'-c',nclusters,'.rds'))

indq <- grepl('meanC\\[', colnames(mcsamples))
matplot(identity(mcsamples[,indq][,1]),type='l',lty=1)




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

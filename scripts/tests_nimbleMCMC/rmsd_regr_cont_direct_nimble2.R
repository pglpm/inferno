## Author: PGL  Porta Mana
## Created: 2021-03-20T10:07:17+0100
## Last-Updated: 2021-09-13T14:52:55+0200
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

rm(constants2, dat2, inits2, bayesnet2, model2, Cmodel2, confmodel2, mcmcsampler2, Cmcmcsampler2)
gc()
nclusters <- 100
ndata <- 20 # nSamples = 37969
## RMSD variable
##indc <- which(grepl('log_RMSD', covNames))
indc <- which(covNames %in% continuousCovs)
ncvars <- length(indc)
icvars <- min(indc)
fcvars <- max(indc)
## ## 0-inf variables
## indc0 <- which(grepl('sasa', covNames))
## nc0vars <- length(indc0)
## ic0vars <- min(indc0)
## fc0vars <- max(indc0)
## ## 0-1 variables
## indc1 <- which(grepl('tanimoto', covNames))
## nc1vars <- length(indc1)
## ic1vars <- min(indc1)
## fc1vars <- max(indc1)
## discrete variables
indd <- which(covNames %in% discreteCovs)
ndvars <- length(indd)
idvars <- min(indd)
fdvars <- max(indd)
##
icparms1 <- icvars
fcparms1 <- fcvars
icparms2 <- ncvars + ndvars + icvars
fcparms2 <- ncvars + ndvars + fcvars
##
idparms1 <- idvars
fdparms1 <- fdvars
idparms2 <- ncvars + ndvars + idvars
fdparms2 <- ncvars + ndvars + fdvars

##
## testnf <- nimbleFunction(
##     run = function(x=double(1),
##                    meanC=double(2), sdC=double(2),
##                    log=integer(0, default=0)){
##         returnType(double(0))
##         prob <- x[1]
## ##        prob <- dnorm(x, mean=meanC, sd=sdC)
##         if(log) return(log(prob))
##         else return(prob)
##         })
## Ctestnf <- compileNimble(testnf)
dMix <- nimbleFunction(
    run = function(x=double(1),
                   logq=double(1), parms=double(2),
                   ##
                   nClusters=integer(0, default=nclusters),
                   iCvars=integer(0, default=icvars),
                   fCvars=integer(0, default=fcvars),
                   iCparms1=integer(0, default=icparms1),
                   fCparms1=integer(0, default=fcparms1),
                   iCparms2=integer(0, default=icparms2),
                   fCparms2=integer(0, default=fcparms2),
                   iDvars=integer(0, default=idvars),
                   fDvars=integer(0, default=fdvars),
                   iDparms1=integer(0, default=idparms1),
                   fDparms1=integer(0, default=fdparms1),
                   iDparms2=integer(0, default=idparms2),
                   fDparms2=integer(0, default=fdparms2),
                   ##
                   log=integer(0, default=0)){
        returnType(double(0))
        sumclusters <- 0
        for(i in 1:nClusters){
            sumclusters <- sumclusters + exp(
                                             logq[i] +
                                             sum(dnorm(x=x[iCvars:fCvars], mean=parms[iCparms1:fCparms1,i], sd=exp(parms[iCparms2:fCparms2,i]), log=TRUE)) +
                                             ## sum(dgamma(x=x[iC0vars:fC0vars], shape=exp(logshapeC[1:nC0vars,i]), scale=exp(logscaleC[1:nC0vars,i]), log=TRUE)) +
                                             ## sum(dbeta(x=x[iC1vars:fC1vars], shape1=exp(logshape1C[1:nC1vars,i]), shape2=exp(logshape2C[1:nC1vars,i]), log=TRUE)) +
                                             sum(dnbinom(x=x[iDvars:fDvars], prob=1/(1+exp(-parms[iDparms1:fDparms1,i])), size=exp(parms[iDparms2:fDparms2,i]), log=TRUE))
                                         )
        }
        if(log) return( log(sumclusters) - log(sum(exp(logq))) )
            else return( sumclusters/sum(exp(logq)) )
        })
##CdMix <- compileNimble(dMix)
##
rMix <- nimbleFunction(
    run = function(n=integer(0, default=1),
                   logq=double(1), parms=double(2),
                   ##
                   nClusters=integer(0, default=nclusters),
                   iCvars=integer(0, default=icvars),
                   fCvars=integer(0, default=fcvars),
                   iCparms1=integer(0, default=icparms1),
                   fCparms1=integer(0, default=fcparms1),
                   iCparms2=integer(0, default=icparms2),
                   fCparms2=integer(0, default=fcparms2),
                   iDvars=integer(0, default=idvars),
                   fDvars=integer(0, default=fdvars),
                   iDparms1=integer(0, default=idparms1),
                   fDparms1=integer(0, default=fdparms1),
                   iDparms2=integer(0, default=idparms2),
                   fDparms2=integer(0, default=fdparms2)
                   ){
        ##
        returnType(double(1))
        nCv <- fCvars - iCvars + 1
        nDv <- fDvars - iDvars + 1
        xout <- numeric(length=nCv+nDv, init=FALSE)
        cluster <- rcat(n=1, prob=exp(logq[1:nClusters])/sum(exp(logq[1:nClusters])))
        ##
        xout[iCvars:fCvars] <- rnorm(n=nCv, mean=parms[iCparms1:fCparms1,cluster], sd=exp(parms[iCparms2:fCparms2,cluster]))
        ## xout[iC0vars:fC0vars] <- rgamma(n=nC0vars, shape=exp(logshapeC[1:nC0vars,cluster]), scale=exp(logscaleC[1:nC0vars,cluster]))
        ## xout[iC1vars:fC1vars] <- rbeta(n=nC1vars, shape1=exp(logshape1C[1:nC1vars,cluster]), shape2=exp(logshape2C[1:nC1vars,cluster]))
        xout[iDvars:fDvars] <- rnbinom(n=nDv, prob=1/(1+exp(-parms[iDparms1:fDparms1,cluster])), size=exp(parms[iDparms2:fDparms2,cluster]))
        ##
        return(xout)
})
##CrMix <- compileNimble(rMix)
##
dparmsPrior <- nimbleFunction(
    run = function(x=double(2),
                   meanCmean=double(0),
                   meanCsd=double(0),
                   sdCshape=double(0),
                   sdCrate=double(0),
                   probDshape1=double(0),
                   probDshape2=double(0),
                   sizeDshape=double(0),
                   sizeDrate=double(0),
                   nClusters=integer(0),
                   ##
                   iCparms1=integer(0, default=icparms1),
                   fCparms1=integer(0, default=fcparms1),
                   iCparms2=integer(0, default=icparms2),
                   fCparms2=integer(0, default=fcparms2),
                   iDparms1=integer(0, default=idparms1),
                   fDparms1=integer(0, default=fdparms1),
                   iDparms2=integer(0, default=idparms2),
                   fDparms2=integer(0, default=fdparms2),
                   ##
                   log=integer(0, default=0)){
        returnType(double(0))
        ##
        lprob <- sum( -(x[iCparms1:fCparms1,1:nClusters] - meanCmean)^2/(2*meanCsd^2) ) +
            sum( -sdCshape * 2 * x[iCparms2:fCparms2,1:nClusters] - sdCrate * exp(-2*x[iCparms2:fCparms2,1:nClusters]) ) +
            sum( probDshape1 * x[iDparms1:fDparms1,1:nClusters] -
                 (probDshape1+probDshape2) * log1p(exp(x[iDparms1:fDparms1,1:nClusters])) ) +
            ## sum( -(probDshape1-1) * log1p(exp(-x[iDparms1:fDparms1])) -
            ##      (probDshape2+1) * log1p(exp(x[iDparms1:fDparms1])) +
            ##      x[iDparms1:fDparms1] ) +
            sum( sizeDshape * x[iDparms2:fDparms2,1:nClusters] - sizeDrate * exp(x[iDparms2:fDparms2,1:nClusters]) )
        ## lprob <- sum(dnorm(x=x[iCparms1:fCparms1], mean=meanCmean, sd=meanCsd, log=TRUE)) +
        ##     sum(dgamma(x=exp(-2*x[iCparms2:fCparms2]), shape=sdCshape, rate=sdCrate, log=TRUE) - 2*x[iCparms2:fCparms2]) +
        ##     sum(dbeta(x=1/(x[iDparms1:fDparms1])^2, shape=probDshape1, rate=probDshape2, log=TRUE) +  )
        ## sum(dgamma(x=exp(x[iDparms2:fDparms2]), shape=sizeDshape, rate=sizeDrate, log=TRUE) + x[iDparms2:fDparms2])
        ##
        if(log) return( lprob )
        else return( exp(lprob) )
    })
##CdparmsPrior <- compileNimble(dparmsPrior)
##
rparmsPrior <- nimbleFunction(
    run = function(n=integer(0, default=1),
                   meanCmean=double(0),
                   meanCsd=double(0),
                   sdCshape=double(0),
                   sdCrate=double(0),
                   probDshape1=double(0),
                   probDshape2=double(0),
                   sizeDshape=double(0),
                   sizeDrate=double(0),
                   nClusters=integer(0),
                   ##
                   iCparms1=integer(0, default=icparms1),
                   fCparms1=integer(0, default=fcparms1),
                   iCparms2=integer(0, default=icparms2),
                   fCparms2=integer(0, default=fcparms2),
                   iDparms1=integer(0, default=idparms1),
                   fDparms1=integer(0, default=fdparms1),
                   iDparms2=integer(0, default=idparms2),
                   fDparms2=integer(0, default=fdparms2)
                   ){
        returnType(double(2))
        ##
        nCv <- fCparms1 - iCparms1 + 1
        nDv <- fDparms1 - iDparms1 + 1
        parmsout <- matrix(nrow=2*(nCv+nDv), ncol=nClusters, init=FALSE)
        for(i in 1:nClusters){
            parmsout[iCparms1:fCparms1,i] <- rnorm(n=nCv, mean=meanCmean, sd=meanCsd)
            parmsout[iCparms2:fCparms2,i] <- rgamma(n=nCv, shape=sdCshape, rate=sdCrate)
            parmsout[iDparms1:fDparms1,i] <- rbeta(n=nDv, shape1=probDshape1, shape2=probDshape2)
            parmsout[iDparms2:fDparms2,i] <- rgamma(n=nDv, shape=sizeDshape, rate=sizeDrate)
        }
        return( parmsout )
    })
##CrparmsPrior <- compileNimble(rparmsPrior)
##
dlogqPrior <- nimbleFunction(
    run = function(x=double(1),
                   alpha=double(1),
                   ##
                   log=integer(0, default=0)){
        returnType(double(0))
        ##
        lprob <- sum( alpha * x ) - 1000 * (sum(exp(x))-1)^2
        ##
        if(log) return( lprob )
        else return( exp(lprob) )
    })
##CdlogqPrior <- compileNimble(dlogqPrior)
##
rlogqPrior <- nimbleFunction(
    run = function(n=integer(0, default=1),
                   alpha=double(1)
                   ){
        returnType(double(1))
        ##
        logqout <- rdirch(n=1, alpha)
        return( logqout )
    })
##CrlogqPrior <- compileNimble(rlogqPrior)
##
constants2 <- list(
    nData=ndata,
    nClusters=nclusters,
    nCvars=ncvars,
    iCvars=icvars,
    fCvars=fcvars,
    nDvars=ndvars,
    iDvars=idvars,
    fDvars=fdvars,
    iCparms1=icparms1,
    fCparms1=fcparms1,
    iCparms2=icparms2,
    fCparms2=fcparms2,
    iDparms1=idparms1,
    fDparms1=fdparms1,
    iDparms2=idparms2,
    fDparms2=fdparms2,
    nVars=ncvars+ndvars,
    nParms=2*(ncvars+ndvars),
    ##
    alpha=rep(1,nclusters)*4/nclusters,
    meanCmean=0,
    meanCsd=2,
    sdCshape=1,
    sdCrate=1,
    probDshape1=1,
    probDshape2=1,
    sizeDshape=1,
    sizeDrate=1
)
##
dat2 <- list(
    X=as.matrix(alldata[1:ndata, ..covNames])
)
##
inits2 <- list(
    logq=log(rep(1/nclusters, nclusters)),
    Parms=rbind( matrix(0, nrow=ncvars, ncol=nclusters),
                qlogis(matrix(0.5, nrow=ndvars, ncol=nclusters)),
                log(matrix(1, nrow=ncvars, ncol=nclusters)),
                log(matrix(10, nrow=ndvars, ncol=nclusters)) )
    ## meanC=matrix(0, nrow=ncvars, ncol=nclusters),
    ## logsdC=log(matrix(1, nrow=ncvars, ncol=nclusters)),
    ## ## logshapeC=log(matrix(1, nrow=nc0vars, ncol=nclusters)),
    ## ## logscaleC=log(matrix(1, nrow=nc0vars, ncol=nclusters)),
    ## ## logshape1C=log(matrix(1, nrow=nc1vars, ncol=nclusters)),
    ## ## logshape2C=log(matrix(1, nrow=nc1vars, ncol=nclusters)),
    ## logitprobD=plogis(matrix(0.5, nrow=ndvars, ncol=nclusters)),
    ## logsizeD=log(matrix(10, nrow=ndvars, ncol=nclusters))
)
##
bayesnet2 <- nimbleCode({
    logq[1:nClusters] ~ dlogqPrior(alpha=alpha[1:nClusters])
    Parms[1:nParms, 1:nClusters] ~ dparmsPrior( meanCmean=meanCmean,
                                         meanCsd=meanCsd,
                                         sdCshape=sdCshape,
                                         sdCrate=sdCrate,
                                         probDshape1=probDshape1,
                                         probDshape2=probDshape2,
                                         sizeDshape=sizeDshape,
                                         sizeDrate=sizeDrate,
                                         nClusters=nClusters,
                                         iCparms1=iCparms1, fCparms1=fCparms1, 
                                         iCparms2=iCparms2, fCparms2=fCparms2,
                                         iDparms1=iDparms1, fDparms1=fDparms1, 
                                         iDparms2=iDparms2, fDparms2=fDparms2 )
    ##
    for(i in 1:nData){
        X[i, 1:nVars] ~ dMix( logq=logq[1:nClusters],
                           parms=Parms[1:nParms, 1:nClusters],
                           nClusters=nClusters,
                           iCvars=iCvars, fCvars=fCvars,
                           iCparms1=iCparms1, fCparms1=fCparms1, 
                           iCparms2=iCparms2, fCparms2=fCparms2,
                           iDvars=iDvars, fDvars=fDvars,
                           iDparms1=iDparms1, fDparms1=fDparms1, 
                           iDparms2=iDparms2, fDparms2=fDparms2 )
    }
})

model2 <- nimbleModel(code=bayesnet2, name='model2', constants=constants2, inits=inits2, data=dat2)

Cmodel2 <- compileNimble(model2, showCompilerOutput=TRUE)

confmodel2 <- configureMCMC(Cmodel2, nodes=NULL)
    confmodel2$addSampler(target=c('logq','Parms'), type='AF_slice', control=list(sliceAdaptFactorMaxIter=1000, sliceAdaptFactorInterval=100, sliceAdaptWidthMaxIter=100, sliceMaxSteps=100, maxContractions=100))
## confmodel2$addSampler(target=paste0('Parms'), type='AF_slice', control=list(sliceAdaptFactorMaxIter=1000, sliceAdaptFactorInterval=100, sliceAdaptWidthMaxIter=100, sliceMaxSteps=100, maxContractions=100))
## for(i in 1:nclusters){
##     confmodel2$addSampler(target=paste0('Parms[1:',2*(ncvars+ndvars),', ',i,']'), type='AF_slice', control=list(sliceAdaptFactorMaxIter=1000, sliceAdaptFactorInterval=100, sliceAdaptWidthMaxIter=100, sliceMaxSteps=100, maxContractions=100))
## }
confmodel2

mcmcsampler2 <- buildMCMC(confmodel2)
Cmcmcsampler2 <- compileNimble(mcmcsampler2, resetFunctions = TRUE)

totaltime <- Sys.time()
mcsamples2 <- runMCMC(Cmcmcsampler2, nburnin=2, niter=4, thin=1, setSeed=123)
totaltime <- Sys.time() - totaltime
totaltime

indq <- grepl('logq\\[', colnames(mcsamples2))
testqs <- normalizem(exp(mcsamples2[,indq]))
matplot(identity(testqs),type='l',lty=1)
sum(testqs[2,])

###########################################

library('nimble')

bayesnet3 <- nimbleCode({
    m ~ dnorm(mean=0, sd=1)
    x ~ dnorm(mean=m, sd=1)
})
##
constants <- list()
dat <- list(x=3)
inits <- list(m=-5)

model3 <- nimbleModel(code=bayesnet3, name='model3', constants=constants, inits=inits, data=dat)
Cmodel3 <- compileNimble(model3, showCompilerOutput=TRUE)

confmodel3 <- configureMCMC(Cmodel3)

mcmcsampler3 <- buildMCMC(confmodel3)
Cmcmcsampler3 <- compileNimble(mcmcsampler3)

Cmcmcsampler3$run(niter=4, nburnin=2)

Cmcmcsampler3$run(niter=4, nburnin=0, reset=FALSE)

Cmcmcsampler3$run(niter=4, nburnin=1, reset=FALSE)

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
CpData <- compileNimble(pData)

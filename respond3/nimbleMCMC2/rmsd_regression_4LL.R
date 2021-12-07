## Author: PGL  Porta Mana
## Created: 2021-03-20T10:07:17+0100
## Last-Updated: 2021-10-05T13:34:09+0200
################
## Script for direct regression, continuous RMSD
################

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
library('doFuture')
library('doRNG')
registerDoFuture()
#library('LaplacesDemon') # used for Dirichlet generator
## library('ash')
## library('extraDistr')
## library('PReMiuM')
## library('mvtnorm')
options(bitmapType='cairo')
pdff <- function(filename){pdf(file=paste0(filename,'.pdf'),paper='a4r',height=11.7,width=16.5)} # to output in pdf format
pngf <- function(filename,res=300){png(file=paste0(filename,'.png'),height=11.7*1.2,width=16.5,units='in',res=res,pointsize=36)} # to output in pdf format
library('nimble')
## library('coda')
#### End custom setup ####

#######################################
## Read and reorganize data
if(exists('alldata')){rm(alldata)}
alldata <- fread('../data_id_processed_transformed_rescaled_shuffled.csv', sep=',')
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

for(obj in c('constants', 'dat', 'inits', 'bayesnet', 'model', 'Cmodel', 'confmodel', 'mcmcsampler', 'Cmcmcsampler')){if(exists(obj)){do.call(rm,list(obj))}}
gc()
##
nclusters <- 100
ndata <- 6000 # nSamples = 37969
nccovs <- length(continuousCovs)
ndcovs <- length(discreteCovs)
##
source('functions_rmsdregr_nimble_binom.R')
meansccovs <- apply(alldata[1:ndata,..continuousCovs],2,mean)
varsccovs <- apply(alldata[1:ndata,..continuousCovs],2,function(x)var(x, na.rm=T))
## shape & scale parameters for the gamma distribution for tau
tauQccovs <- sapply(continuousCovs, function(acov){
    fn <- function(parms){
        (pinvgamma(varsccovs[acov]/sqrt(10), shape=parms[1], scale=parms[2]) - 0.005)^2 +
            (pinvgamma(varsccovs[acov]*sqrt(10), shape=parms[1], scale=parms[2]) - 0.995)^2
    }
    myoptimbounds(c(1, 1), fn=fn, lower=0, upper=Inf)$par
})
## print(sapply(continuousCovs, function(acov){
##     fn <- function(parms){c(
##          qinvgamma(0.005, shape=parms[1], scale=parms[2]),
##          varsccovs[acov]/sqrt(10),
##          qinvgamma(0.995, shape=parms[1], scale=parms[2]),
##          varsccovs[acov]*sqrt(10)
##          ) }
##     fn(tauQccovs[,acov])
## }))
##
meansdcovs <- apply(alldata[1:ndata,..discreteCovs],2,mean)
varsdcovs <- apply(alldata[1:ndata,..discreteCovs],2,function(x)var(x, na.rm=T))
maxdcovs <- apply(alldata[1:ndata,..discreteCovs],2,max)
## prob and size parameters for the neg-binomial distribution for sizeD
sizeQdcovs <- sapply(discreteCovs, function(acov){
    fn <- function(parms){
        (pnbinom(round(maxdcovs[acov]/sqrt(10)), prob=parms[1], size=parms[2]) - 0.005)^2 +
            (pnbinom(round(maxdcovs[acov]*sqrt(10)), prob=parms[1], size=parms[2]) - 0.995)^2
    }
    myoptimbounds(c(0.5, 1), fn=fn, lower = c(0.001,0), upper = c(1,Inf))$par
})
## print(sapply(discreteCovs, function(acov){
##     fn <- function(parms){c(
##          qnbinom(0.005, prob=parms[1], size=parms[2]),
##          round(maxdcovs[acov]/sqrt(10)),
##          qnbinom(0.995, prob=parms[1], size=parms[2]),
##          round(maxdcovs[acov]*sqrt(10))
##          ) }
##     fn(sizeQdcovs[,acov])
## }))
##
## shape1 and shape2 parameters for the beta distribution for probD
shapesratio <- (maxdcovs-meansdcovs)/meansdcovs
alphadcovs <- sapply(discreteCovs, function(acov){
    fn <- function(parms){
        parms2 <- parms*shapesratio[acov]
        -(lbeta(parms,parms2) - (parms-1)*digamma(parms) - (parms2-1)*digamma(parms2) + (parms+parms2-2)*digamma(parms+parms2))
    }
    myoptimbounds(1, fn=fn, lower = 0, upper = Inf)$par
})
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
    probDshape1=alphadcovs,
    probDshape2=alphadcovs*shapesratio,
    sizeDpar1=sizeQdcovs[1,], # 1/(maxdcovs), # 1/(1+2*meansdcovs),
    sizeDpar2=sizeQdcovs[2,], # maxdcovs/(maxdcovs-1), # 1+0*maxdcovs,
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
    C=rep(1,ndata), # rcat(n=ndata, prob=rep(1/nclusters,nclusters))
    ML=0
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
            probD[acov,acluster] ~ dbeta(shape1=probDshape1[acov], shape2=probDshape2[acov])
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
        ML <- calcLL(X=X[1:nData,1:nCcovs], Y=Y[1:nData,1:nDcovs],
                      Q=q[1:nClusters],
                      MeanC=meanC[1:nCcovs,1:nClusters], TauC=tauC[1:nCcovs,1:nClusters],
                      ProbD=probD[1:nDcovs,1:nClusters], SizeD=sizeD[1:nDcovs,1:nClusters])
    }
})
##

posterior <- TRUE
if(posterior){
    model <- nimbleModel(code=bayesnet, name='model1', constants=constants, inits=inits, data=dat)
}else{
    model <- nimbleModel(code=bayesnet, name='model1', constants=constants, inits=inits, data=list())
    }
Cmodel <- compileNimble(model, showCompilerOutput=FALSE)
gc()
##
if(posterior){
    confmodel <- configureMCMC(Cmodel, monitors=c('q','meanC', 'tauC', 'probD', 'sizeD', 'ML')) #, control=list(adaptive=FALSE))
}else{
    confmodel <- configureMCMC(Cmodel, monitors=c('q','meanC', 'tauC', 'probD', 'sizeD')) #, control=list(adaptive=FALSE))
}
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
list(
    alphaK=rep(5/nclusters, nclusters),
    meanCmean=meansccovs,
    meanCtau=0.5/varsccovs,
    tauCshape=tauQccovs[1,],
    tauCrate=tauQccovs[2,],
    probDshape1=alphadcovs,
    probDshape2=alphadcovs*shapesratio,
    sizeDpar1=sizeQdcovs[1,], # 1/(maxdcovs), # 1/(1+2*meansdcovs),
    sizeDpar2=sizeQdcovs[2,], # maxdcovs/(maxdcovs-1), # 1+0*maxdcovs,
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
version <- 'testLLpostHM'
gc()
totalruntime <- Sys.time()
mcsamples <- runMCMC(Cmcmcsampler, nburnin=1, niter=101, thin=1, inits=initsFunction, setSeed=149)
## Cmcmcsampler$run(niter=2000, thin=1, reset=FALSE, resetMV=TRUE)
## mcsamples <- as.matrix(Cmcmcsampler$mvSamples)
totalruntime <- Sys.time() - totalruntime
print(totalruntime)
## 7 vars, 6000 data, 100 cl, 2000 iter, slice: 2.48 h
## 7 vars, 6000 data, 100 cl, 5001 iter, slice: 6.84 h
## 7 vars, 6000 data, 100 cl: rougly 8.2 min/(100 iterations)




##
saveRDS(mcsamples,file=paste0('_mcsamples-R',version,'-V',length(covNames),'-D',ndata,'-K',nclusters,'-I',nrow(mcsamples),'.rds'))
##
## mcsamples <- readRDS(file='_mcsamples-RpriorHM-V7-D6000-K100-I2000.rds')
parmList <- mcsamples2parmlist(mcsamples, parmNames=c('q', 'meanC', 'tauC', 'probD', 'sizeD'))
momentstraces <- moments12Samples(parmList)
allmomentstraces <- cbind(Dcov=plogis(momentstraces$Dcovars, scale=1/10), momentstraces$means, log(momentstraces$vars), momentstraces$covars)
##
diagnESS <- LaplacesDemon::ESS(allmomentstraces)
diagnBMK <- LaplacesDemon::BMK.Diagnostic(allmomentstraces, batches=2)[,1]
diagnMCSE <- 100*LaplacesDemon::MCSE(allmomentstraces)/apply(allmomentstraces, 2, sd)
diagnStat <- apply(allmomentstraces, 2, function(x){LaplacesDemon::is.stationary(as.matrix(x,ncol=1))})
diagnBurn <- apply(allmomentstraces, 2, function(x){LaplacesDemon::burnin(x)})
##
##
timecount <- Sys.time()
plan(sequential)
plan(multisession, workers = 6L)
samplesQuantiles <- calcSampleQuantiles(parmList)
plan(sequential)
print(Sys.time()-timecount)
## 7 covs, 2000 samples, serial: 1.722 min 
##
alldataRanges <- dataQuantiles <- xlimits <- list()
for(acov in covNames){
    dataQuantiles[[acov]] <- quantile(alldata[[acov]], prob=c(0.005,0.995))
    alldataRanges[[acov]] <- range(alldata[[acov]])
    xlimits[[acov]] <- range(c(alldataRanges[[acov]], samplesQuantiles[,acov,]))
}
##
##
pdff(paste0('mcsummary-R',version,'-V',length(covNames),'-D',ndata,'-K',nclusters,'-I',nrow(mcsamples)))
matplot(1:2, type='l', col='white', main='Stats')
legend(x='topleft', bty='n', cex=2, legend=c( paste0('min ESS = ', min(diagnESS)),
                                        paste0('max BMK = ', max(diagnBMK)),
                                        paste0('max MCSE = ', max(diagnMCSE)),
                                        paste0('all stationary: ', all(diagnStat)),
                                        paste0('burn: ', max(diagnBurn))))
##
par(mfrow = c(2, 3))
for(addvar in setdiff(covNames, 'log_RMSD')){
    matplot(x=c(rep(alldataRanges[['log_RMSD']], each=2),
                alldataRanges[['log_RMSD']][1]),
            y=c(alldataRanges[[addvar]], rev(alldataRanges[[addvar]]),
                alldataRanges[[addvar]][1]),
            type='l', lwd=2, col=paste0(palette()[2], '88'),
             xlim=xlimits[['log_RMSD']],
            ylim=xlimits[[addvar]],
            xlab='log_RMSD',
            ylab=addvar
            )
    matlines(x=c(rep(dataQuantiles[['log_RMSD']], each=2),
                 dataQuantiles[['log_RMSD']][1]),
             y=c(dataQuantiles[[addvar]], rev(dataQuantiles[[addvar]]),
                 dataQuantiles[[addvar]][1]),
             lwd=2, col=paste0(palette()[4], '88'))
}
##
par(mfrow=c(1,1))
for(acov in colnames(allmomentstraces)){
    matplot(allmomentstraces[, acov], type='l', lty=1,
            col=palette()[if(grepl('^MEAN_', acov)){1}else if(grepl('^VAR_', acov)){3}else if(acov=='Dcov'){2}else{4}],
            main=paste0(acov,
                        '\nESS = ', signif(diagnESS[acov], 3),
                        ' | BMK = ', signif(diagnBMK[acov], 3),
                        ' | MCSE(6.27) = ', signif(diagnMCSE[acov], 3),
                        ' | stat: ', diagnStat[acov],
                        ' | burn: ', diagnBurn[acov]
                        ),
            ylab=acov)
}
dev.off()
##
##
nxsamples <- 1000
nfsamples <- 100
##
timecount <- Sys.time()
plan(sequential)
plan(multisession, workers = 6L)
xsamples <- samplesFsamples(parmList=parmList, nxsamples=nxsamples, nfsamples=nfsamples)
plan(sequential)
print(Sys.time()-timecount)
##
subsamplep <- round(seq(1, dim(xsamples)[3], length.out=100))
subsamplex <- round(seq(1, dim(xsamples)[2], length.out=1000))
##
pdff(paste0('hypersamplesvars2D-R',version,'-V',length(covNames),'-D',ndata,'-K',nclusters,'-I',nrow(mcsamples)))
par(mfrow = c(2, 3))
for(addvar in setdiff(covNames, 'log_RMSD')){
    matplot(x=alldata[['log_RMSD']][subsamplex],
            y=alldata[[addvar]][subsamplex],
            xlim=xlimits[['log_RMSD']],
            ylim=xlimits[[addvar]],
            xlab='log_RMSD',
            ylab=addvar,
            type='p', pch=1, cex=0.2, lwd=1, col=palette()[3])
    matlines(x=c(rep(alldataRanges[['log_RMSD']], each=2), alldataRanges[['log_RMSD']][1]),
             y=c(alldataRanges[[addvar]], rev(alldataRanges[[addvar]]), alldataRanges[[addvar]][1]),
             lwd=2, col=paste0(palette()[2],'88'))
    matlines(x=c(rep(dataQuantiles[['log_RMSD']], each=2), dataQuantiles[['log_RMSD']][1]),
             y=c(dataQuantiles[[addvar]], rev(dataQuantiles[[addvar]]), dataQuantiles[[addvar]][1]),
             lwd=2, col=paste0(palette()[4],'88'))
}
for(asample in subsamplep){
    par(mfrow = c(2, 3))
    for(addvar in setdiff(covNames, 'log_RMSD')){
        matplot(x=xsamples['log_RMSD', subsamplex, asample][subsamplex],
                y=xsamples[addvar, subsamplex, asample][subsamplex],
                xlim=xlimits[['log_RMSD']],
                ylim=xlimits[[addvar]],
                xlab='log_RMSD',
                ylab=addvar,
                type='p', pch=1, cex=0.2, lwd=1, col=palette()[1])
        matlines(x=c(rep(alldataRanges[['log_RMSD']], each=2), alldataRanges[['log_RMSD']][1]),
                 y=c(alldataRanges[[addvar]], rev(alldataRanges[[addvar]]), alldataRanges[[addvar]][1]),
                 lwd=2, col=paste0(palette()[2],'88'))
        matlines(x=c(rep(dataQuantiles[['log_RMSD']], each=2), dataQuantiles[['log_RMSD']][1]),
                 y=c(dataQuantiles[[addvar]], rev(dataQuantiles[[addvar]]), dataQuantiles[[addvar]][1]),
                 lwd=2, col=paste0(palette()[4],'88'))
    }
}
dev.off()

##
if(TRUE){
timecount <- Sys.time()
plan(sequential)
plan(multisession, workers = 6L)
loglikelihood <- loglikelihoodorig <- matrix(llSamples(dat=dat, parmList=parmList), ncol=1)
plan(sequential)
print(Sys.time()-timecount)
loglikelihood[is.infinite(loglikelihood)] <- NA
allmomentstraces <- cbind(matrix(loglikelihood, ncol=1, dimnames=list(NULL, 'Log-likelihood')), Dcov=plogis(momentstraces$Dcovars, scale=1/10), momentstraces$means, log(momentstraces$vars), momentstraces$covars)
##
diagnESS <- LaplacesDemon::ESS(allmomentstraces)
diagnBMK <- LaplacesDemon::BMK.Diagnostic(allmomentstraces, batches=2)[,1]
diagnMCSE <- 100*apply(allmomentstraces, 2, function(x){LaplacesDemon::MCSE(x[!is.na(x)])/sd(x, na.rm=TRUE)})
diagnStat <- apply(allmomentstraces, 2, function(x){LaplacesDemon::is.stationary(as.matrix(x,ncol=1))})
diagnBurn <- apply(allmomentstraces, 2, function(x){LaplacesDemon::burnin(x)})
##
pdff(paste0('mcsummary-R',version,'-V',length(covNames),'-D',ndata,'-K',nclusters,'-I',nrow(mcsamples)))
matplot(1:2, type='l', col='white', main='Stats')
legend(x='topleft', bty='n', cex=2, legend=c( paste0('min ESS = ', min(diagnESS)),
                                        paste0('max BMK = ', max(diagnBMK)),
                                        paste0('max MCSE = ', max(diagnMCSE)),
                                        paste0('all stationary: ', all(diagnStat)),
                                        paste0('burn: ', max(diagnBurn))))
for(acov in colnames(allmomentstraces)){
    matplot(allmomentstraces[, acov], type='l', lty=1,
            col=palette()[if(grepl('^MEAN_', acov)){1}else if(grepl('^VAR_', acov)){3}else if(acov=='Dcov'){2}else{4}],
            main=paste0(acov,
                        '\nESS = ', signif(diagnESS[acov], 3),
                        ' | BMK = ', signif(diagnBMK[acov], 3),
                        ' | MCSE(6.27) = ', signif(diagnMCSE[acov], 3),
                        ' | stat: ', diagnStat[acov],
                        ' | burn: ', diagnBurn[acov]
                        ),
            ylab=acov)
}
dev.off()
}


#######################################################################
#######################################################################
#######################################################################
#######################################################################

mcsamples <- readRDS(file='_mcsamples-RpriorHM-V7-D6000-K100-I2000.rds')[1:1001,]
parmList <- mcsamples2parmlist(mcsamples, parmNames=c('q', 'meanC', 'tauC', 'probD', 'sizeD'))

gc()
timecount <- Sys.time()
plan(sequential)
plan(multisession, workers = 6L)
testll2 <- llSamples(dat=dat, parmList=parmList)
plan(sequential)
print(Sys.time()-timecount)
## Time difference of 2.336992 mins



calcLL <- nimbleFunction( run=function(
    X=double(2), Y=double(2), Q=double(1),
    MeanC=double(2), TauC=double(2),
    ProbD=double(2), SizeD=double(2)
    ){
    returnType(double(0))
    Nclusters <- length(Q)
    Ndata <- dim(X)[1]
    Nccovs <- dim(X)[2]
    Ndcovs <- dim(Y)[2]
    LL <- 0
    for(adatum in 1:Ndata){
        clustersum <- log(Q)
        for(acov in 1:Nccovs){
            clustersum <- clustersum +
                dnorm(x=X[adatum,acov], mean=MeanC[acov,], sd=1/sqrt(TauC[acov,]), log=TRUE)
        }
        for(acov in 1:Ndcovs){
            clustersum <- clustersum +
                dbinom(x=Y[adatum,acov], prob=ProbD[acov,], size=SizeD[acov,], log=TRUE)
        }
        LL <- LL + log(sum(exp(clustersum)))
    }
    return(LL)
} )
CcalcLL <- compileNimble(calcLL)
assign('CcalcLL', CcalcLL, envir = .GlobalEnv)
assign('calcLL', calcLL, envir = .GlobalEnv)

gc()
timecount <- Sys.time()
plan(sequential)
plan(multisession, workers = 6L)
testll <- foreach(asample=1:10, .combine=c, .export='CcalcLL', .packages='nimble', .inorder=TRUE)%do%{
CcalcLL(X=dat$X, Y=dat$Y, Q=parmList$q[asample,],
       MeanC=parmList$meanC[asample,,], TauC=parmList$tauC[asample,,],
       ProbD=parmList$probD[asample,,], SizeD=parmList$sizeD[asample,,])
}
plan(sequential)
print(Sys.time()-timecount)

calcLLall <- nimbleFunction( run=function(
    X=double(2), Y=double(2), Q=double(2),
    MeanC=double(3), TauC=double(3),
    ProbD=double(3), SizeD=double(3)
    ){
    returnType(double(1))
    Nsamples <- dim(Q)[1]
    Nclusters <- dim(Q)[2]
    Ndata <- dim(X)[1]
    Nccovs <- dim(X)[2]
    Ndcovs <- dim(Y)[2]
    LL <- numeric(value=0, length=Nclusters)
    for(asample in 1:Nsamples){
        pLL <- 0
    for(adatum in 1:Ndata){
        clustersum <- log(Q[asample,])
        for(acov in 1:Nccovs){
            clustersum <- clustersum +
                dnorm(x=X[adatum,acov], mean=MeanC[asample,acov,], sd=1/sqrt(TauC[asample,acov,]), log=TRUE)
        }
        for(acov in 1:Ndcovs){
            clustersum <- clustersum +
                dbinom(x=Y[adatum,acov], prob=ProbD[asample,acov,], size=SizeD[asample,acov,], log=TRUE)
        }
        pLL <- pLL + log(sum(exp(clustersum)))
    }
        LL[asample] <- pLL
    }
    return(LL)
} )
CcalcLLall <- compileNimble(calcLLall)
assign('CcalcLLall', CcalcLLall, envir = .GlobalEnv)
assign('calcLLall', calcLLall, envir = .GlobalEnv)


gc()
timecount <- Sys.time()
testll3 <- CcalcLLall(X=dat$X, Y=dat$Y, Q=parmList$q,
       MeanC=parmList$meanC, TauC=parmList$tauC,
       ProbD=parmList$probD, SizeD=parmList$sizeD)
print(Sys.time()-timecount)




CcalcLL(X=dat$X, Y=dat$Y, Q=parmList$q[1,],
       MeanC=parmList$meanC[1,,], TauC=parmList$tauC[1,,],
       ProbD=parmList$probD[1,,], SizeD=parmList$sizeD[1,,])


nsamples <- nrow(parmList$q)
chunksize <- 2
nchunks <- nsamples/chunksize
gc()
timecount <- Sys.time()
plan(sequential)
##plan(multisession, workers = 6L)
testll <- foreach(achunk=0:(nchunks-1), .combine=c)%dopar%{
    subsamples <- chunksize*achunk+(1:chunksize)
    probs <- rowSums( sapply(continuousCovs, function(acov){
        dnorm(rep(dat$X[,acov], each=chunksize*nclusters), mean=c(t(parmList$meanC[subsamples,acov,])), sd=1/sqrt(c(t(parmList$tauC[subsamples,acov,]))), log=TRUE)
    }) ) +
        rowSums( sapply(discreteCovs, function(acov){
            dbinom(rep(dat$Y[,acov], each=chunksize*nclusters), prob=c(t(parmList$probD[subsamples,acov,])), size=c(t(parmList$sizeD[subsamples,acov,])), log=TRUE)
        }) ) +
        log(c(t(parmList$q[subsamples,])))
    dim(probs) <- c(nclusters, chunksize, ndata)
    rowSums(log(colSums(exp(probs))))
}
plan(sequential)
print(Sys.time()-timecount)


nsamples <- nrow(parmList$q)
nchunks <- nsamples/chunksize
gc()
timecount <- Sys.time()
plan(sequential)
plan(multisession, workers = 6L)
testll3 <- foreach(asample=1:nsamples, .combine=c, .inorder=TRUE)%:%foreach(adatum=1:ndata, .combine='+', .inorder=FALSE)%dopar%{
    log(sum(exp(colSums(
        dnorm(dat$X[adatum,], mean=parmList$meanC[asample,,], sd=1/sqrt(parmList$tauC[asample,,]), log=TRUE) ) +
        colSums( 
            dbinom(dat$Y[adatum,], prob=parmList$probD[asample,,], size=parmList$sizeD[asample,,], log=TRUE) ) +
        log(parmList$q[asample,]))))
}
plan(sequential)
print(Sys.time()-timecount)


asample <- 4
adatum <- 1
log(sum(foreach(acluster=1:nclusters, .combine=c)%do%{
    parmList$q[asample,acluster] *
        prod(dnorm(dat$X[adatum,],
                   mean=parmList$meanC[asample,,acluster],
                   sd=1/sqrt(parmList$tauC[asample,,acluster]))) *
        prod(dbinom(dat$Y[adatum,], prob=parmList$probD[asample,,acluster],
                    size=parmList$sizeD[asample,,acluster]))
}))


testme <- meansdcovs[1]
testma <- maxdcovs[1]
xgrid <- 0:ceiling(1.5*testma)
matplot(xgrid, dnbinom(xgrid, prob=))



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
alldataRanges <- dataQuantiles <- xlim <- ylim <- list()
for(var in covNames){
    dataQuantiles[[var]] <- quantile(alldata[[var]], prob=c(0.05,0.95))
    alldataRanges[[var]] <- thisrange <- range(alldata[[var]])
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
        matlines(x=c(rep(alldataRanges[['log_RMSD']], each=2), alldataRanges[['log_RMSD']][1]),
             y=c(alldataRanges[[addvar]], rev(alldataRanges[[addvar]]), alldataRanges[[addvar]][1]),
                 lwd=2, col=paste0(palette()[2],'88'))
        matlines(x=c(rep(dataQuantiles[['log_RMSD']], each=2), dataQuantiles[['log_RMSD']][1]),
             y=c(dataQuantiles[[addvar]], rev(dataQuantiles[[addvar]]), dataQuantiles[[addvar]][1]),
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
    matlines(x=c(rep(alldataRanges[['log_RMSD']], each=2), alldataRanges[['log_RMSD']][1]),
             y=c(alldataRanges[[addvar]], rev(alldataRanges[[addvar]]), alldataRanges[[addvar]][1]),
             lwd=2, col=paste0(palette()[2],'88'))
            matlines(x=c(rep(dataQuantiles[['log_RMSD']], each=2), dataQuantiles[['log_RMSD']][1]),
             y=c(dataQuantiles[[addvar]], rev(dataQuantiles[[addvar]]), dataQuantiles[[addvar]][1]),
                 lwd=2, col=paste0(palette()[4],'88'))
}
}
dev.off()











indq <- grepl('logProb_q\\[', colnames(mcsamples))
matplot(identity(mcsamples[,indq]),type='l',lty=1)

indq <- grepl('meanC\\[1, 1]', colnames(mcsamples)) || grepl('meanC\\[1, 1]', colnames(mcsamples))

totalruntime <- Sys.time()
mcsamplesb <- runMCMC(Cmcmcsampler, nburnin=0, niter=2000, thin=1, nchains=2, inits=initsFunction, setSeed=123)
totalruntime <- Sys.time() - totalruntime
totalruntime
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

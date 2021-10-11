## Author: PGL  Porta Mana
## Created: 2021-03-20T10:07:17+0100
## Last-Updated: 2021-10-11T11:04:55+0200
################
## Batch script for direct regression, continuous RMSD
################
if(file.exists("/cluster/home/pglpm/R")){
    .libPaths(c("/cluster/home/pglpm/R",.libPaths()))
}

seed <- 149
baseversion <- 'checksHMUp_'
nclusters <- 2L^6
ndata <- 2L^12 # nSamples = 37969
niter <- 2L^10
nstages <- 31L
ncheckpoints <- 8L
covNames <-  c('log_RMSD'
               ,'log_mcs_unbonded_polar_sasa'
               ,'logit_ec_tanimoto_similarity'
               ,'mcs_NumHeteroAtoms'
               ## ,'scale_fc_tanimoto_similarity'
               ,'docked_HeavyAtomCount'
               ,'mcs_RingCount'
               ,'docked_NumRotatableBonds'
               ,'mcs_NOCount'
               )
## pdff('check_mutualinfo')
## for(i in 1:(length(covNames)-1)){
##     for(j in (i+1):length(covNames)){
##         matplot(alldata[[covNames[i]]], alldata[[covNames[j]]], type='p', pch='.', col=paste0('#000000','88'), xlab=covNames[i], ylab=covNames[j])
##     }
## }
## dev.off()


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
## library('LaplacesDemon')
## library('extraDistr')
## library('mvtnorm')
options(bitmapType='cairo')
pdff <- function(filename){pdf(file=paste0(filename,'.pdf'),paper='a4r',height=11.7,width=16.5)} # to output in pdf format
pngf <- function(filename,res=300){png(file=paste0(filename,'.png'),height=11.7*1.2,width=16.5,units='in',res=res,pointsize=36)} # to output in pdf format
library('nimble')
## library('coda')
#### End custom setup ####

##
if(exists('alldata')){rm(alldata)}
datafile <- 'data_id_processed_transformed_rescaled_shuffled.csv'
if(!file.exists(datafile)){
    datafile <- paste0('../', datafile)
}
alldata <- fread(datafile, sep=',')
discreteCovs <- covNames[sapply(covNames, function(x){is.integer(alldata[[x]])})]
continuousCovs <- covNames[sapply(covNames, function(x){is.double(alldata[[x]])})]
covNames <- c(continuousCovs, discreteCovs)
nccovs <- length(continuousCovs)
ndcovs <- length(discreteCovs)
rm(alldata)

for(obj in c('constants', 'dat', 'inits', 'bayesnet', 'model', 'Cmodel', 'confmodel', 'mcmcsampler', 'Cmcmcsampler')){if(exists(obj)){do.call(rm,list(obj))}}
gc()
##
#################
## Create test dataset
#################
## Load last state from a past semi-converged sampling and calculate params
source('functions_rmsdregr_nimble_binom.R')
alldataRanges <- readRDS(file='alldataRanges.rds')
dataQuantiles <- readRDS(file='dataQuantiles.rds')
meansccovs <- readRDS(file='meansccovs.rds')
varsccovs <- readRDS(file='varsccovs.rds')
tauQccovs <- readRDS(file='tauQccovs.rds')
meansdcovs <- readRDS(file='meansdcovs.rds')
varsdcovs <- readRDS(file='varsdcovs.rds')
maxdcovs <- readRDS(file='maxdcovs.rds')
sizeQdcovs <- readRDS(file='sizeQdcovs.rds')
shapesratiodcovs <- readRDS(file='shapesratiodcovs.rds')
alphadcovs <- readRDS(file='alphadcovs.rds')
set.seed(222)
parmListTest <- list()
##
parmListTest$q <- matrix(rdirch(n=1, alpha=rep(1/nclusters, nclusters)), nrow=1)
##
parmListTest$meanC <- array(rnorm(n=nccovs*nclusters, mean=meansccovs[continuousCovs], sd=sqrt(sqrt(10)*varsccovs[continuousCovs])), dim=c(1, nccovs, nclusters))
parmListTest$tauC <- array(rgamma(n=nccovs*nclusters, shape=tauQccovs[1,continuousCovs], rate=tauQccovs[2,continuousCovs]), dim=c(1, nccovs, nclusters))
dimnames(parmListTest$meanC) <- dimnames(parmListTest$tauC) <- list(NULL, continuousCovs, NULL)
##
parmListTest$probD <- array(rbeta(n=ndcovs*nclusters, shape1=alphadcovs[discreteCovs], shape2=alphadcovs[discreteCovs]*shapesratiodcovs[discreteCovs]), dim=c(1, ndcovs, nclusters))
parmListTest$sizeD <- array(rnbinom(n=ndcovs*nclusters, prob=sizeQdcovs[1,discreteCovs], size=sizeQdcovs[2,discreteCovs]), dim=c(1, ndcovs, nclusters))
dimnames(parmListTest$probD) <- dimnames(parmListTest$sizeD) <- list(NULL, discreteCovs, NULL)
##
saveRDS(parmListTest, file=paste0('_parmListTest-R',baseversion,'-V',length(covNames),'-K',nclusters,'.rds'))

##
Cs <- rcat(n=ndata, prob=parmListTest$q[1,])
alldata <- cbind(
    data.table(matrix(rnorm(n=nccovs*ndata, mean=t(parmListTest$meanC[1,continuousCovs,Cs]), sd=1/sqrt(t(parmListTest$tauC[1,continuousCovs,Cs]))), nrow=ndata, ncol=nccovs)),
    data.table(matrix(rbinom(n=ndcovs*ndata, prob=t(parmListTest$probD[1,discreteCovs,Cs]), size=t(parmListTest$sizeD[1,discreteCovs,Cs])), nrow=ndata, ncol=ndcovs))
    )
colnames(alldata) <- covNames
saveRDS(alldata,file=paste0('_alldataTest-R',baseversion,'-V',length(covNames),'-D',ndata,'-K',nclusters,'.rds'))
##
##
dat <- list(
    X=alldata[, ..continuousCovs],
    Y=alldata[, ..discreteCovs]
)
##

##
source('functions_rmsdregr_nimble_binom.R')
alldataRanges <- dataQuantiles <- list()
for(acov in covNames){
        dataQuantiles[[acov]] <- quantile(alldata[[acov]], prob=c(0.005,0.995))
        alldataRanges[[acov]] <- range(alldata[[acov]])
}
meansccovs <- apply(alldata[1:ndata,..continuousCovs],2,mean)
varsccovs <- apply(alldata[1:ndata,..continuousCovs],2,function(x)var(x, na.rm=T))
## shape & scale parameters for the gamma distribution for tau
tauQccovs <- sapply(continuousCovs, function(acov){
    fn <- function(parms){
        parms <- exp(parms)
        (qinvgamma(0.005, shape=parms[1], scale=parms[2])/(varsccovs[acov]/128)-1)^2 +
            (qinvgamma(0.995, shape=parms[1], scale=parms[2])/(varsccovs[acov]*128)-1)^2
    }
    exp(myoptim(c(0, 0), fn=fn)$par)
})
print('Diagnostics tauQccovs')
fn <- function(parms){cbind(
                          qinvgamma(0.005, shape=parms[1,], scale=parms[2,]),
                          varsccovs/128,
                          qinvgamma(0.995, shape=parms[1,], scale=parms[2,]),
                          varsccovs*128
                      ) }
fn(tauQccovs)
##
meansdcovs <- apply(alldata[1:ndata,..discreteCovs],2,mean)
varsdcovs <- apply(alldata[1:ndata,..discreteCovs],2,function(x)var(x, na.rm=T))
maxdcovs <- apply(alldata[1:ndata,..discreteCovs],2,max)
## prob and size parameters for the neg-binomial distribution for sizeD
sizeQdcovs <- sapply(discreteCovs, function(acov){
    fn <- function(parms){
        parms <- c(plogis(parms[1]), exp(parms[2]))
        (qnbinom(0.005, prob=parms[1], size=parms[2])/ceiling(maxdcovs[acov]/128)-1)^2 +
            (qnbinom(0.995, prob=parms[1], size=parms[2])/round(maxdcovs[acov]*128)-1)^2
    }
    out <- myoptim(c(0, 1), fn=fn)$par
    c(plogis(out[1]), exp(out[2]))
})
print('Diagnostics sizeQdcovs')
fn <- function(parms){cbind(
                          qnbinom(0.005, prob=parms[1,], size=parms[2,]),
                          ceiling(maxdcovs/128),
                          qnbinom(0.995, prob=parms[1,], size=parms[2,]),
                          round(maxdcovs*128)
                      ) }
fn(sizeQdcovs)
##
## shape1 and shape2 parameters for the beta distribution for probD
shapesratiodcovs <- (maxdcovs-meansdcovs)/meansdcovs
alphadcovs <- sapply(discreteCovs, function(acov){
    fn <- function(parms){
        parms2 <- parms*shapesratiodcovs[acov]
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
## dat <- list(
##     X=as.matrix(alldata[1:ndata, ..continuousCovs]),
##     Y=as.matrix(alldata[1:ndata, ..discreteCovs])
## )
##
inits <- list(
    alphaK=rep(1/nclusters, nclusters),
    meanCmean=meansccovs,
    meanCtau=1/(128*varsccovs),
    tauCshape=tauQccovs[1,],
    tauCrate=tauQccovs[2,],
    probDshape1=alphadcovs,
    probDshape2=alphadcovs*shapesratiodcovs,
    sizeDpar1=sizeQdcovs[1,], # 1/(maxdcovs), # 1/(1+2*meansdcovs),
    sizeDpar2=sizeQdcovs[2,], # maxdcovs/(maxdcovs-1), # 1+0*maxdcovs,
    ##
    q=rep(1/nclusters, nclusters),    
    ## meanC=matrix(meansccovs, nrow=nccovs, ncol=nclusters),
    ## tauC=matrix(1/varsccovs, nrow=nccovs, ncol=nclusters),
    ## probD=matrix(meansdcovs/maxdcovs, nrow=ndcovs, ncol=nclusters),
    ## sizeD=matrix(maxdcovs, nrow=ndcovs, ncol=nclusters),
    meanC=matrix(rnorm(n=nccovs*nclusters, mean=meansccovs, sd=sqrt(128*varsccovs)), nrow=nccovs, ncol=nclusters),
    tauC=matrix(rgamma(n=nccovs*nclusters, shape=tauQccovs[1,], rate=tauQccovs[2,]), nrow=nccovs, ncol=nclusters),
    probD=matrix(rbeta(n=ndcovs*nclusters, shape1=alphadcovs, shape2=alphadcovs*shapesratiodcovs), nrow=ndcovs, ncol=nclusters),
    sizeD=apply(matrix(rnbinom(n=ndcovs*nclusters, prob=sizeQdcovs[1,], size=sizeQdcovs[2,]), nrow=ndcovs, ncol=nclusters), 2, function(x){maxdcovs*(x<maxdcovs)+x*(x>=maxdcovs)}),
    ## C=rep(1,ndata)
    C=rcat(n=ndata, prob=rep(1/nclusters,nclusters))
         )
##
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
    }
})
##

timecount <- Sys.time()

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
confmodel <- configureMCMC(Cmodel,
                           monitors=c('q','meanC', 'tauC', 'probD', 'sizeD'),
                           monitors2=c('C')) #, control=list(adaptive=FALSE))
}else{
confmodel <- configureMCMC(Cmodel,
                           monitors=c('q','meanC', 'tauC', 'probD', 'sizeD'))
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
## source('functions_rmsdregr_nimble_binom.R')
initsFunction <- function(){
list(
    alphaK=rep(1/nclusters, nclusters),
    meanCmean=meansccovs,
    meanCtau=1/(128*varsccovs),
    tauCshape=tauQccovs[1,],
    tauCrate=tauQccovs[2,],
    probDshape1=alphadcovs,
    probDshape2=alphadcovs*shapesratiodcovs,
    sizeDpar1=sizeQdcovs[1,], # 1/(maxdcovs), # 1/(1+2*meansdcovs),
    sizeDpar2=sizeQdcovs[2,], # maxdcovs/(maxdcovs-1), # 1+0*maxdcovs,
    ##
    q=rep(1/nclusters, nclusters),    
    ## meanC=matrix(meansccovs, nrow=nccovs, ncol=nclusters),
    ## tauC=matrix(1/varsccovs, nrow=nccovs, ncol=nclusters),
    ## probD=matrix(meansdcovs/maxdcovs, nrow=ndcovs, ncol=nclusters),
    ## sizeD=matrix(maxdcovs, nrow=ndcovs, ncol=nclusters),
    meanC=matrix(rnorm(n=nccovs*nclusters, mean=meansccovs, sd=sqrt(128*varsccovs)), nrow=nccovs, ncol=nclusters),
    tauC=matrix(rgamma(n=nccovs*nclusters, shape=tauQccovs[1,], rate=tauQccovs[2,]), nrow=nccovs, ncol=nclusters),
    probD=matrix(rbeta(n=ndcovs*nclusters, shape1=alphadcovs, shape2=alphadcovs*shapesratiodcovs), nrow=ndcovs, ncol=nclusters),
    sizeD=apply(matrix(rnbinom(n=ndcovs*nclusters, prob=sizeQdcovs[1,], size=sizeQdcovs[2,]), nrow=ndcovs, ncol=nclusters), 2, function(x){maxdcovs*(x<maxdcovs)+x*(x>=maxdcovs)}),
    ## C=rep(1,ndata)
    C=rcat(n=ndata, prob=rep(1/nclusters,nclusters))
         )
}
set.seed(941)
checkpoints <- rbind(
    c(meansccovs, round(meansdcovs)),
    c(meansccovs+sqrt(varsccovs), round(meansdcovs+sqrt(varsdcovs))),
    c(meansccovs-sqrt(varsccovs), sapply(round(meansdcovs-sqrt(varsdcovs)), function(x){max(0,x)})),
    as.matrix(alldata[sample(1:ndata, size=ncheckpoints), ..covNames])
)
rownames(checkpoints) <- c('Pdatamean', 'PdatacornerHi', 'PdatacornerLo', paste0('Pdatum',1:ncheckpoints))
saveRDS(checkpoints,file=paste0('_checkpoints-R',baseversion,'-V',length(covNames),'-D',ndata,'-K',nclusters,'.rds'))

print('Setup time:')
print(Sys.time() - timecount)

##
for(stage in 0:nstages){
    totalruntime <- Sys.time()

    print(paste0('==== STAGE ', stage, ' ===='))
    version <- paste0(baseversion, stage)
    gc()
    if(stage==0){
        mcsamples <- runMCMC(Cmcmcsampler, nburnin=1, niter=niter+1, thin=1, thin2=niter, inits=initsFunction, setSeed=seed)
    }else{
        Cmcmcsampler$run(niter=niter, thin=1, thin2=niter, reset=FALSE, resetMV=TRUE)
    }
    mcsamples <- as.matrix(Cmcmcsampler$mvSamples)
    print('MCMC time:')
    print(Sys.time() - totalruntime)
    ## 7 vars, 6000 data, 100 cl, 2000 iter, slice: 2.48 h
    ## 7 vars, 6000 data, 100 cl, 5001 iter, slice: 6.84 h
    ## 7 vars, 6000 data, 100 cl: rougly 8.2 min/(100 iterations)
    ##
    saveRDS(mcsamples,file=paste0('_mcsamples-R',version,'-V',length(covNames),'-D',ndata,'-K',nclusters,'-I',nrow(mcsamples),'.rds'))
    ## save final state of MCMC chain
    finalstate <- as.matrix(Cmcmcsampler$mvSamples2)
    usedclusters <- length(unique(c(finalstate)))
    print(paste0('OCCUPIED CLUSTERS: ', usedclusters, ' OF ', nclusters))
    finalstate <- c(mcsamples[nrow(mcsamples),], finalstate[nrow(finalstate),])
    saveRDS(finalstate2list(finalstate),file=paste0('_finalstate-R',version,'-V',length(covNames),'-D',ndata,'-K',nclusters,'-I',nrow(mcsamples),'.rds'))
    ##
    parmList <- mcsamples2parmlist(mcsamples)
    ## Traces to follow for diagnostics
    ll <- llSamples(dat, parmList)
    llTest <- llSamples(dat, parmListTest)
    momentstraces <- moments12Samples(parmList)
    probCheckpoints <- t(probValuesSamples(checkpoints, parmList))
    traces <- cbind(probCheckpoints, do.call(cbind, momentstraces))
    saveRDS(traces,file=paste0('_traces-R',version,'-V',length(covNames),'-D',ndata,'-K',nclusters,'-I',nrow(mcsamples),'.rds'))
    ##
    if(nrow(traces)>=1000){
        funMCSE <- function(x){LaplacesDemon::MCSE(x, method='batch.means')$se}
    }else{
        funMCSE <- function(x){LaplacesDemon::MCSE(x)}
    }
    diagnESS <- LaplacesDemon::ESS(traces)
    diagnBMK <- LaplacesDemon::BMK.Diagnostic(traces, batches=2)[,1]
    diagnMCSE <- 100*apply(traces, 2, function(x){funMCSE(x)/sd(x)})
    diagnStat <- apply(traces, 2, function(x){LaplacesDemon::is.stationary(as.matrix(x,ncol=1))})
    diagnBurn <- apply(traces, 2, function(x){LaplacesDemon::burnin(matrix(x[1:(10*trunc(length(x)/10))], ncol=1))})
    ##
    momentstracesTest <- moments12Samples(parmListTest)
    probCheckpointsTest <- t(probValuesSamples(checkpoints, parmListTest))
    tracesTest <- cbind(probCheckpointsTest, do.call(cbind, momentstracesTest))
    diagnSE <- vapply(colnames(traces), function(x){100*ecdf(traces[,x])(tracesTest[1,x])}, FUN.VALUE=numeric(1))
    ##
    tracegroups <- list('maxD'=1:(ncheckpoints+4),
                        '1D'=(ncheckpoints+4)+(1:(2*(nccovs+ndcovs))),
                        '2D'=((ncheckpoints+4)+2*(nccovs+ndcovs)+1):ncol(traces) )
    grouplegends <- foreach(agroup=1:length(tracegroups))%do%{
        c( paste0('-- STATS ', names(tracegroups)[agroup], ' --'),
          paste0('min ESS = ', min(diagnESS[tracegroups[[agroup]]])),
          paste0('max BMK = ', max(diagnBMK[tracegroups[[agroup]]])),
          paste0('max MCSE = ', max(diagnMCSE[tracegroups[[agroup]]])),
          paste0('all stationary: ', all(diagnStat[tracegroups[[agroup]]])),
          paste0('burn: ', max(diagnBurn[tracegroups[[agroup]]])),
          paste0('min percent.: ', min(c(diagnSE[tracegroups[[agroup]]],100-diagnSE[tracegroups[[agroup]]]))) )
    }
    ##
    ## plan(sequential)
    ## plan(multisession, workers = 6L)
    samplesQuantiles <- calcSampleQuantiles(parmList)
    ## plan(sequential)
    ## 7 covs, 2000 samples, serial: 1.722 min 
    ##
    xlimits <- list()
    for(acov in covNames){
        xlimits[[acov]] <- range(c(alldataRanges[[acov]], samplesQuantiles[,acov,]))
    }
    ##

    ##
    colpalette <- sapply(colnames(traces),function(acov){
        if(acov=='Pdatamean'){1}
        else if(grepl('^Pcorner', acov)){3}
        else if(grepl('^Pdatum', acov)){4}
        else if(grepl('^Dcov', acov)){2}
        else if(grepl('^MEAN_', acov)){5}
        else if(grepl('^VAR_', acov)){3}
        else{4}
    })
    names(colpalette) <- colnames(traces)
    ##
    pdff(paste0('mcsummary-R',version,'-V',length(covNames),'-D',ndata,'-K',nclusters,'-I',nrow(mcsamples)))
    matplot(1:2, type='l', col='white', main=paste0('Stats stage ',stage), axes=FALSE, ann=FALSE)
    legendpositions <- c('topleft','bottomleft','topright')
    for(alegend in 1:length(grouplegends)){
        legend(x=legendpositions[alegend], bty='n', cex=1.5,
               legend=grouplegends[[alegend]] )
    }
    legend(x='bottomright', bty='n', cex=1.5,
           legend=c(
               paste0('Occupied clusters: ', usedclusters, ' of ', nclusters),
               paste0('LL: ', signif(mean(ll),3), ' +- ', signif(sd(ll),3))
           ))
    ##
    par(mfrow=c(1,1))
    for(acov in continuousCovs){
        Xgrid <- seq(min(alldata[[acov]]), max(alldata[[acov]]), length.out=2^8)
        plotsamples <- samplesfX(acov, parmList, Xgrid, nfsamples=64)
        plotsamplesTest <- samplesfX(acov, parmListTest, Xgrid, nfsamples=1)
        matplot(Xgrid, plotsamples, type='l', col=paste0(palette()[7], '44'), lty=1, lwd=2, xlab=acov, ylab='probability density', cex.axis=1.5, cex.lab=1.5)
        matlines(Xgrid,plotsamplesTest, col=palette()[2], lty=1, lwd=5)
    }
    for(acov in discreteCovs){
        Xgrid <- seq(min(alldata[[acov]]), max(alldata[[acov]]), by=1)
        plotsamples <- samplesfX(acov, parmList, Xgrid, nfsamples=64)
        plotsamplesTest <- samplesfX(acov, parmListTest, Xgrid, nfsamples=1)
        matplot(Xgrid, plotsamples, type='l', col=paste0(palette()[7], '44'), lty=1, lwd=2, xlab=acov, ylab='probability density', cex.axis=1.5, cex.lab=1.5)
        matlines(Xgrid,plotsamplesTest, col=palette()[2], lty=1, lwd=5)
    }
    ##
    par(mfrow = c(2, 4))
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
    matplot(ll, type='l', col=palette()[3], lty=1, main='LL', ylab='LL', ylim=range(c(ll,llTest)))
    matlines(x=c(1,length(ll)), rep(llTest,2), lwd=5, lty=1, col=paste0('#000000','88'))
    for(acov in colnames(traces)){
        if(grepl('^[PDV]', acov)){transf <- function(x){log(abs(x))}
        }else{transf <- identity}
        matplot(transf(traces[,acov]), type='l', lty=1, col=colpalette[acov],
                main=paste0(acov,
                            '\nESS = ', signif(diagnESS[acov], 3),
                            ' | BMK = ', signif(diagnBMK[acov], 3),
                            ' | MCSE(6.27) = ', signif(diagnMCSE[acov], 3),
                            ' | stat: ', diagnStat[acov],
                            ' | burn: ', diagnBurn[acov],
                            ' | percentile: ', signif(diagnSE[acov], 3)
                            ),
                ylab=acov,
                ylim=range(c(transf(traces[,acov]), transf(tracesTest[,acov][abs(tracesTest[,acov])<Inf]))))
        matlines(x=c(1,nrow(traces)), rep(transf(tracesTest[,acov]),2), lwd=5, lty=1, col=paste0('#000000','88'))
    }
    dev.off()

    print('Total runtime:')
    print(Sys.time() - totalruntime)

}


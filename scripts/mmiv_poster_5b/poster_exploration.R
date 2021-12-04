## Author: PGL  Porta Mana
## Created: 2021-03-20T10:07:17+0100
## Last-Updated: 2021-12-02T09:22:01+0100
################
## Exploration for MMIV poster
################
if(file.exists("/cluster/home/pglpm/R")){
    .libPaths(c("/cluster/home/pglpm/R",.libPaths()))
}
#### Custom setup ####
## Colour-blind friendly palettes, from https://personal.sron.nl/~pault/
## library('khroma')
## palette(colour('bright')())
## scale_colour_discrete <- scale_colour_bright
## palette(colour('muted')())
library('data.table')
## library('ggplot2')
## library('ggthemes')
## theme_set(theme_bw(base_size=18))
#library('cowplot')
library('png')
library('foreach')
library('doFuture')
library('doRNG')
registerDoFuture()
print('availableCores:')
print(availableCores())
print('availableCores-multicore:')
print(availableCores('multicore'))
if(file.exists("/cluster/home/pglpm/R")){
    plan(multicore, workers=availableCores()-1)
}else{
    plan(multisession, workers=6)
}
## library('LaplacesDemon')
## library('extraDistr')
## library('mvtnorm')
## options(bitmapType='cairo')
## pdff <- function(filename){pdf(file=paste0(filename,'.pdf'),paper='a4r',height=11.7,width=16.5)} # to output in pdf format
## pngf <- function(filename,res=300){png(file=paste0(filename,'.png'),height=11.7*1.2,width=16.5,units='in',res=res,pointsize=36)} # to output in pdf format
library('nimble')
## library('coda')
#### End custom setup ####

oalldata <- data.table(
group=as.integer(c(0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0)),
y=c(64, 60, 59, 65, 64, 62, 54, 68, 67, 79, 45, 48, 59, 65, 87, 49, 46, 46, 97, 36, 67, 56, 62, 65, 65, 81, 83, 71),
x=c(74, 70, 65, 67, 62, 67, 51, 93, 56, 78, 58, 52, 60, 76, 74, 36, 42, 50, 79, 35, 60, 58, 57, 68, 60, 89, 58, 70)
)

## tplot(oalldata[group==1,x], oalldata[group==1,y], xlim=range(oalldata$x), ylim=range(oalldata$y), type='p')
## matplot(oalldata[group==2,x], oalldata[group==2,y], type='p', pch=0, col=2, add=T)

## hi <- foreach(g=1:2)%do%{
##     thist(oalldata[group==g,x], n=tticks(oalldata$x, n=8))
## }
## tplot(cbind(hi[[1]]$breaks, hi[[2]]$breaks), cbind(hi[[1]]$density, hi[[2]]$density))
## scatteraxis(1,oalldata[group==1,x], col=1)
## scatteraxis(3,oalldata[group==2,x], col=2)

## do this for paired analysis
## oalldata <- rbind(data.table(x=oalldata$x, group=0L),
##                   data.table(x=oalldata$y, group=1L))


pdff('new_test_dataplot', height=4, paper='special')
## par(mai=c(0.8,2.8,0.8,0))
##    par(mai=c(0.8,0.8,0.8,0))
## matplot(x=NA, y=NA, xlim=range(oalldata$x)+c(-10,0), ylim=c(0,2), axes=F, xlab='X = % tumour-volume increase after 6 months', ylab=NA, cex.lab=2)
## axis(1, tticks(range(oalldata$x)), #labels=sprintf('%2g', tticks(xlim)),
##          lwd=0, cex.axis=2, mgp=c(3,0.5,0), col.axis='#000000', gap.axis=0.25)
## axis(2, c(1.5,0.5), labels=c('group A', 'group B'), lwd=0, cex.axis=2, las=1, mgp=c(3,0.5,0), col.axis='black', gap.axis=0.25)
## text(x=-10 + min(oalldata$x), y=1.5, labels='group A', col=1, cex=2.5, xpd=NA)
## text(x=-10 + min(oalldata$x), y=0.5, labels='group B', col=2, cex=2.5, xpd=NA)
## for(i in pretty(range(oalldata$x))){abline(v=i, lty=1, lwd=3, col='#BBBBBB60')}
#add1 <- rnorm(length(oalldata[group==0,x]),0,0.5)
#add2 <- rnorm(length(oalldata[group==1,x]),0,0.5)
## matlines(x=rbind(oalldata[group==0,x]+add1,oalldata[group==0,x]+add1), y=cbind(c(1,1.9)), lwd=5, lty=1, col=paste0(palette()[1],'88'))
tplot(x=oalldata[group==0,x],
##      y=runif(n=length(oalldata[group==0,x]), 1.25,1.75),
      y=sample(seq(1.25, 1.75, length.out=length(oalldata[group==0,x]))),
      type='p', pch=16, col=1, cex=2.25, alpha=0.25,
      xlim=range(oalldata$x)+c(0,0),
      ylim=c(0,2),
      xlab='X = % tumour-volume increase after 6 months',
      yticks=-30, 
      mar=c(3.25, 12, 1, 1)+c(1,1,1,1)
      )
##
tplot(x=oalldata[group==1,x],
#      y=runif(n=length(oalldata[group==1,x]), 0.25, 0.75),
      y=sample(seq(0.25, 0.75, length.out=length(oalldata[group==1,x]))),
      type='p', pch=15, col=2, cex=2.25, alpha=0.25,
      add=T)
## matlines(x=rbind(oalldata[group==1,x]+add2,oalldata[group==1,x]+add2), y=cbind(c(0.1,1)), lwd=5, lty=1, col=paste0(palette()[2],'88'))
 text(x=-10 + min(oalldata$x), y=1.5, labels='group A', col=1, cex=2.5, xpd=NA)
 text(x=-10 + min(oalldata$x), y=0.5, labels='group B', col=2, cex=2.5, xpd=NA)
dev.off()


#################################
## Setup for Monte Carlo sampling
#################################

seed <- 149
baseversion <- 'mmivexample_X_'
nclusters <- as.integer(16) #as.integer(2^6)
niter <- as.integer(2^11) #as.integer(2^11)
niter0 <- as.integer(2^10) #as.integer(2^10)
thin <- as.integer(10)
nstages <- as.integer(1)
ncheckpoints <- as.integer(4)
covNames <-  c('x', 'group')
posterior <- TRUE

discreteCovs <- covNames[sapply(covNames, function(x){is.integer(oalldata[[x]])})]
continuousCovs <- covNames[sapply(covNames, function(x){is.double(oalldata[[x]])})]
nccovs <- length(continuousCovs)
ndcovs <- length(discreteCovs)
##
alldata <- oalldata[,..covNames]
ndata <- as.integer(nrow(alldata))

for(obj in c('constants', 'dat', 'inits', 'bayesnet', 'model', 'Cmodel', 'confmodel', 'mcmcsampler', 'Cmcmcsampler')){if(exists(obj)){do.call(rm,list(obj))}}
gc()
##
##
dat <- list(
    X=alldata[, ..continuousCovs],
    Y=alldata[, ..discreteCovs]
)
##
##
source('functions_rmsdregr_nimble_binom.R')
irq2sd <- 1/(2*qnorm(3/4))
alldataRanges <- dataQuantiles <- list()
for(acov in covNames){
        dataQuantiles[[acov]] <- quant(alldata[[acov]], prob=c(0.005,0.995))
        alldataRanges[[acov]] <- range(alldata[[acov]])
}

medianccovs <- apply(alldata[1:ndata,..continuousCovs],2,median)
widthccovs <- irq2sd * apply(alldata[1:ndata,..continuousCovs],2,IQR)
##
mediandcovs <- apply(alldata[1:ndata,..discreteCovs],2,function(x){max(median(x),1)})
widthdcovs <- ceiling(apply(alldata[1:ndata,..discreteCovs],2,IQR))
maxdcovs <- apply(alldata[1:ndata,..discreteCovs],2,max)
##
print('Creating and saving checkpoints')
checkpointsFile <- paste0('_checkpoints-',ncheckpoints,'-R',baseversion,'-V',length(covNames),'-D',ndata,'-K',nclusters,'.rds')
checkpoints <- rbind(
    c(medianccovs, round(mediandcovs)),
    c(medianccovs+widthccovs, round(mediandcovs+widthdcovs)),
    c(medianccovs-widthccovs, sapply(round(mediandcovs-widthdcovs), function(x){max(0,x)})),
    as.matrix(alldata[sample(1:ndata, size=ncheckpoints), ..covNames])
)
rownames(checkpoints) <- c('Pdatamedians', 'PdatacornerHi', 'PdatacornerLo', paste0('Pdatum',1:ncheckpoints))
saveRDS(checkpoints,file=checkpointsFile)
##
##
constants <- list(
    nClusters=nclusters,
    nData=ndata,
    nCcovs=nccovs,
    nDcovs=ndcovs
)

##
initsFunction <- function(){
list(
    alphaK=rep(1/nclusters, nclusters),
    meanCmean=medianccovs,
    meanCshape1=rep(1/2, nccovs),
    meanCrate2=1/(widthccovs/2)^2, # dims = inv. variance
    ##
    tauCshape1=rep(1, nccovs),
    tauCrate2=1/(widthccovs/2)^2, # dims = inv. variance
    ##
    probDa1=rep(100, ndcovs),
    probDb1=rep(100, ndcovs),
    sizeDprob1=rep(1/2,2),
    sizeDsize1=2,
    ##
    ##
    meanCtau1=1/(widthccovs/2)^2, # dims = inv. variance
    meanCrate1=(widthccovs/2)^2, # dims = variance
    tauCrate1=(widthccovs/2)^2, # dims = variance
    ##
    ##
    q=rep(1/nclusters, nclusters),
    meanC=matrix(rnorm(n=nccovs*(nclusters), mean=medianccovs, sd=1/sqrt(rgamma(n=nccovs, shape=1/2, rate=rgamma(n=nccovs, shape=1/2, rate=1/(widthccovs/2)^2)))), nrow=nccovs, ncol=nclusters),
    tauC=matrix(rgamma(n=nccovs*(nclusters), shape=1/2, rate=rgamma(n=nccovs, shape=1/2, rate=1/(widthccovs/2)^2)), nrow=nccovs, ncol=nclusters),
    probD=matrix(rbeta(n=ndcovs*(nclusters), shape1=1, shape2=1), nrow=ndcovs, ncol=nclusters),
    sizeD=matrix(rcat(ndcovs*nclusters, prob=rep(1/2,2)), nrow=ndcovs, ncol=nclusters),
    ## sizeD=matrix(apply(matrix(rnbinom(n=ndcovs*(nclusters), prob=rbeta(n=ndcovs, shape1=16, shape2=16), size=maxdcovs), nrow=ndcovs, ncol=nclusters), 2, function(x){maxdcovs*(x<maxdcovs)+x*(x>=maxdcovs)}), nrow=ndcovs, ncol=nclusters),
    ##
    C=rep(1,ndata)
)
}

##
##
bayesnet <- nimbleCode({
    q[1:nClusters] ~ ddirch(alpha=alphaK[1:nClusters])
    for(acluster in 1:nClusters){
        for(acov in 1:nCcovs){
            meanC[acov,acluster] ~ dnorm(mean=meanCmean[acov], tau=meanCtau1[acov])
            tauC[acov,acluster] ~ dgamma(shape=tauCshape1[acov], rate=tauCrate1[acov])
        }
        for(acov in 1:nDcovs){
            probD[acov,acluster] ~ dbeta(shape1=probDa1[acov], shape2=probDb1[acov])
            sizeD[acov,acluster] ~ dcat(prob=sizeDprob1[1:2])
        }
    }
    ##
    for(acov in 1:nCcovs){
        meanCtau1[acov] ~ dgamma(shape=meanCshape1[acov], rate=meanCrate1[acov])
        meanCrate1[acov] ~ dgamma(shape=meanCshape1[acov], rate=meanCrate2[acov])
        tauCrate1[acov] ~ dgamma(shape=tauCshape1[acov], rate=tauCrate2[acov])
    }
    ##
    if(posterior){
        for(adatum in 1:nData){
            C[adatum] ~ dcat(prob=q[1:nClusters])
        }            ##
        for(adatum in 1:nData){
            for(acov in 1:nCcovs){
                X[adatum,acov] ~ dnorm(mean=meanC[acov,C[adatum]], tau=tauC[acov,C[adatum]])
            }
            for(acov in 1:nDcovs){
                Y[adatum,acov] ~ dbinom(prob=probD[acov,C[adatum]], size=sizeD[acov,C[adatum]]-1)
            }
        }
    }
})
##

timecount <- Sys.time()

if(posterior){
    model <- nimbleModel(code=bayesnet, name='model1', constants=constants, inits=initsFunction(), data=dat)
}else{
    model <- nimbleModel(code=bayesnet, name='model1', constants=constants, inits=initsFunction(), data=list())
    }
Cmodel <- compileNimble(model, showCompilerOutput=FALSE)
gc()
##

    confmodel <- configureMCMC(Cmodel,
                               monitors=c('q','meanC', 'tauC', 'probD', 'sizeD'),
                               monitors2=c('C', 'meanCtau1', 'meanCrate1', 'tauCrate1'))


if(posterior){
    confmodel <- configureMCMC(Cmodel, nodes=NULL,
                               monitors=c('q','meanC', 'tauC', 'probD', 'sizeD'),
                               monitors2=c('C', 'meanCtau1', 'meanCrate1', 'tauCrate1'))
    for(adatum in 1:ndata){
        confmodel$addSampler(target=paste0('C[', adatum, ']'), type='categorical')
    }
    confmodel$addSampler(target=paste0('q[1:', nclusters, ']'), type='conjugate')
    for(acluster in 1:nclusters){
        for(acov in 1:nccovs){
            confmodel$addSampler(target=paste0('meanC[', acov, ', ', acluster, ']'), type='conjugate')
            confmodel$addSampler(target=paste0('tauC[', acov, ', ', acluster, ']'), type='conjugate')
        }
        for(acov in 1:ndcovs){
            confmodel$addSampler(target=paste0('probD[', acov, ', ', acluster, ']'), type='conjugate')
            confmodel$addSampler(target=paste0('sizeD[', acov, ', ', acluster, ']'), type='categorical')
        }
    }
    for(acov in 1:nccovs){
        confmodel$addSampler(target=paste0('meanCtau1[', acov, ']'), type='conjugate')
        confmodel$addSampler(target=paste0('meanCrate1[', acov, ']'), type='conjugate')
        confmodel$addSampler(target=paste0('tauCrate1[', acov, ']'), type='conjugate')
    }
}else{
    confmodel <- configureMCMC(Cmodel,
                               monitors=c('q','meanC', 'tauC', 'probD', 'sizeD'))
}
## confmodel$printSamplers(executionOrder=TRUE)
print(confmodel)

mcmcsampler <- buildMCMC(confmodel)
Cmcmcsampler <- compileNimble(mcmcsampler, resetFunctions = TRUE)
gc()

print('Setup time:')
print(Sys.time() - timecount)

##################################################
## Monte Carlo sampler and plots of MC diagnostics
##################################################
nstages <- 3
for(stage in 0:nstages){
    totalruntime <- Sys.time()

    print(paste0('==== STAGE ', stage, ' ===='))
    version <- paste0(baseversion, stage)
    gc()
    if(stage==0){
        mcsamples <- runMCMC(Cmcmcsampler, nburnin=1, niter=niter0+1, thin=1, thin2=niter0, inits=initsFunction, setSeed=seed)
    }else{
        Cmcmcsampler$run(niter=niter*thin, thin=thin, thin2=niter*thin, reset=FALSE, resetMV=TRUE)
    }
    mcsamples <- as.matrix(Cmcmcsampler$mvSamples)
    print('MCMC time:')
    print(Sys.time() - totalruntime)
    ## 7 vars, 6000 data, 100 cl, 2000 iter, slice: 2.48 h
    ## 7 vars, 6000 data, 100 cl, 5001 iter, slice: 6.84 h
    ## 7 vars, 6000 data, 100 cl: rougly 8.2 min/(100 iterations)
    ##
    if(any(is.na(mcsamples))){print('WARNING: SOME NA OUTPUTS')}
    if(any(!is.finite(mcsamples))){print('WARNING: SOME INFINITE OUTPUTS')}
    saveRDS(mcsamples,file=paste0('_mcsamples-R',version,'-V',length(covNames),'-D',ndata,'-K',nclusters,'-I',nrow(mcsamples),'.rds'))
    ## save final state of MCMC chain
    finalstate <- as.matrix(Cmcmcsampler$mvSamples2)
    finalstate <- c(mcsamples[nrow(mcsamples),], finalstate[nrow(finalstate),])
    occupations <- finalstate[grepl('^C\\[', names(finalstate))]
    usedclusters <- length(unique(occupations))
    print(paste0('OCCUPIED CLUSTERS: ', usedclusters, ' OF ', nclusters))
    saveRDS(finalstate2list(finalstate),file=paste0('_finalstate-R',version,'-V',length(covNames),'-D',ndata,'-K',nclusters,'-I',nrow(mcsamples),'.rds'))
    ##
    parmList <- mcsamples2parmlist(mcsamples)
    parmList$sizeD <- parmList$sizeD-1
    ## Traces to follow for diagnostics
    ll <- llSamples(dat, parmList)
    momentstraces <- moments12Samples(parmList)
    probCheckpoints <- t(probValuesSamples(checkpoints, parmList))
    ## miqrtraces <- calcSampleMQ(parmList)
    ## medians <- miqrtraces[,,1]
    ## colnames(medians) <- paste0('MEDIAN_', colnames(miqrtraces))
    ## Q1s <- miqrtraces[,,2]
    ## colnames(Q1s) <- paste0('Q1_', colnames(miqrtraces))
    ## Q3s <- miqrtraces[,,3]
    ## colnames(Q3s) <- paste0('Q3_', colnames(miqrtraces))
    ## iqrs <- Q3s - Q1s
    ## colnames(iqrs) <- paste0('IQR_', colnames(miqrtraces))
    traces <- cbind(LL=ll, probCheckpoints, #medians, iqrs, Q1s, Q3s,
                    do.call(cbind, momentstraces))
    saveRDS(traces,file=paste0('_traces-R',version,'-V',length(covNames),'-D',ndata,'-K',nclusters,'-I',nrow(mcsamples),'.rds'))
    ##
    if(nrow(traces)>=1000){
        funMCSE <- function(x){LaplacesDemon::MCSE(x, method='batch.means')$se}
    }else{
        funMCSE <- function(x){LaplacesDemon::MCSE(x)}
    }
    diagnESS <- LaplacesDemon::ESS(traces * (abs(traces) < Inf))
    diagnBMK <- LaplacesDemon::BMK.Diagnostic(traces, batches=2)[,1]
    diagnMCSE <- 100*apply(traces, 2, function(x){funMCSE(x)/sd(x)})
    diagnStat <- apply(traces, 2, function(x){LaplacesDemon::is.stationary(as.matrix(x,ncol=1))})
    diagnBurn <- apply(traces, 2, function(x){LaplacesDemon::burnin(matrix(x[1:(10*trunc(length(x)/10))], ncol=1))})
    ##
    ##
    tracenames <- colnames(traces)
    tracegroups <- list(
        'maxD'=tracenames[grepl('^(Pdat|Dcov|LL)', tracenames)],
        '1D'=tracenames[grepl('^(MEDIAN|Q1|Q3|IQR|MEAN|VAR)_', tracenames)],
        '2D'=tracenames[grepl('^COV_', tracenames)]
    )
    grouplegends <- foreach(agroup=1:length(tracegroups))%do%{
        c( paste0('-- STATS ', names(tracegroups)[agroup], ' --'),
          paste0('min ESS = ', min(diagnESS[tracegroups[[agroup]]])),
          paste0('max BMK = ', max(diagnBMK[tracegroups[[agroup]]])),
          paste0('max MCSE = ', max(diagnMCSE[tracegroups[[agroup]]])),
          paste0('all stationary: ', all(diagnStat[tracegroups[[agroup]]])),
          paste0('burn: ', max(diagnBurn[tracegroups[[agroup]]]))
          )
    }
    colpalette <- sapply(tracenames, function(atrace){
        c(2, 3, 1) %*%
        sapply(tracegroups, function(agroup){atrace %in% agroup})
        })
    ##
    samplesQuantiles <- calcSampleQuantiles(parmList)
    ##
    xlimits <- list()
    for(acov in covNames){
        xlimits[[acov]] <- range(c(alldataRanges[[acov]], samplesQuantiles[,acov,]))
    }
    ##

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
        Xgrid <- seq(extendrange(alldata[[acov]])[1], extendrange(alldata[[acov]])[2], length.out=2^8)
        df <- 1/min(diff(Xgrid))/2
        plotsamples <- samplesfX(acov, parmList, Xgrid, nfsamples=64)
        tplot(Xgrid, plotsamples, type='l', col=paste0(palette()[7], '44'), lty=1, lwd=2, xlab=acov, ylab='probability density', cex.axis=1.5, cex.lab=1.5, ylim=c(0, max(plotsamples[plotsamples<df])))
    }
    for(acov in discreteCovs){
        Xgrid <- seq(max(0,min(alldata[[acov]])-1), max(alldata[[acov]])+1, by=1)
        plotsamples <- samplesfX(acov, parmList, Xgrid, nfsamples=64)
        tplot(Xgrid, plotsamples, type='l', col=paste0(palette()[7], '44'), lty=1, lwd=2, xlab=acov, ylab='probability density', cex.axis=1.5, cex.lab=1.5, ylim=c(0, max(plotsamples)))
    }
    ##
    par(mfrow = c(2, 2))
    for(maincov in covNames){
    for(addvar in setdiff(covNames, maincov)){
        tplot(x=c(rep(alldataRanges[[maincov]], each=2),
                    alldataRanges[[maincov]][1]),
                y=c(alldataRanges[[addvar]], rev(alldataRanges[[addvar]]),
                    alldataRanges[[addvar]][1]),
                type='l', lwd=2, col=paste0(palette()[2], '88'),
                xlim=xlimits[[maincov]],
                ylim=xlimits[[addvar]],
                xlab=maincov,
                ylab=addvar
                )
        matlines(x=c(rep(dataQuantiles[[maincov]], each=2),
                     dataQuantiles[[maincov]][1]),
                 y=c(dataQuantiles[[addvar]], rev(dataQuantiles[[addvar]]),
                     dataQuantiles[[addvar]][1]),
                 lwd=2, col=paste0(palette()[4], '88'))
    }}
    ##
    par(mfrow=c(1,1))
#    matplot(ll, type='l', col=palette()[3], lty=1, main='LL', ylab='LL', ylim=range(ll[abs(ll)<Inf]))
    for(acov in colnames(traces)){
        if(grepl('^[PDV]', acov)){transf <- function(x){log(abs(x)+1e-12)}
        }else{transf <- identity}
        tplot(y=transf(traces[,acov]), type='l', lty=1, col=colpalette[acov],
                main=paste0(acov,
                            '\nESS = ', signif(diagnESS[acov], 3),
                            ' | BMK = ', signif(diagnBMK[acov], 3),
                            ' | MCSE(6.27) = ', signif(diagnMCSE[acov], 3),
                            ' | stat: ', diagnStat[acov],
                            ' | burn: ', diagnBurn[acov]
                            ),
                ylab=acov
              #, ylim=range(c(transf(traces[,acov][abs(transf(traces[,acov]))<Inf])))
              )
    }
    dev.off()

    print('Total runtime:')
    print(Sys.time() - totalruntime)

}
############################################################
## End MCMC
############################################################


source('functions_rmsdregr_nimble_binom.R')
if(!exists('mcsamples')){
mcsamples <- readRDS('_mcsamples-Rmmivexample_M12S1_3-V2-D28-K16-I2048.rds')
parmList <- mcsamples2parmlist(mcsamples)
}


############################################################
## Example plots of various kinds of conditional predictions
############################################################

xrange <- extendrange(alldata$x)
xrangemin <- range(alldata$x)
xgrid <- seq(xrange[1], xrange[2], length.out=256)
xtest <- round(quant(alldata$x, (1:3)/4)*2)/2

fdists <- samplesfRgivenX(maincov='x', X=cbind(group=0:1), parmList=parmList, covgrid=xgrid, inorder=T)
fdists <- aperm(fdists)
preddists <- colMeans(fdists)
margdists1 <- apply(fdists,c(2,3),function(x){quant(x, c(1,3)/4)})
margdists2 <- apply(fdists,c(2,3),function(x){quant(x, c(1,7)/8)})

edists <- samplesmeanRgivenX(maincov='x', X=cbind(group=0:1), parmList=parmList, inorder=T)
edists <- aperm(edists)
meandiffs <- apply(edists,1,function(x){-diff(x)})
mdh1 <- thist(edists[,1])
mdh2 <- thist(edists[,2])
mdh <- thist(meandiffs, n=16)

vdists <- samplesvarRgivenX(maincov='x', X=cbind(group=0:1), parmList=parmList, inorder=T)
vdists <- aperm(vdists)
vardiffs <- apply(vdists,1,function(x){-diff(x)})
vdh1 <- thist(vdists[,1])
vdh2 <- thist(vdists[,2])
vdh <- thist(vardiffs)

## https://www.socscistatistics.com/confidenceinterval/default4.aspx
## μ1 - μ2 = (M1 - M2) = 15.11111, 95% CI [5.5625164, 24.6597036].
## You can be 95% confident that the difference between your two population means (μ1 - μ2) lies between 5.5625164 and 24.6597036.
## μ1 - μ2 = (M1 - M2) = 15.11111, 90% CI [7.1879676, 23.0342524].
## You can be 90% confident that the difference between your two population means (μ1 - μ2) lies between 7.1879676 and 23.0342524.

## https://stats.libretexts.org/Learning_Objects/02%3A_Interactive_Statistics/32%3A_Two_Independent_Samples_With_Statistics_Calculator
## ≠ p 0.002225456123561953
## < p 0.998887271938219
## > p 0.0011127280617809765
## t 3.4604299890224537
## LB90% 7.612632347927092
## UB90% 22.6095876520729
## LB95% 6.0548586063425205
## UB95% 24.167361393657472


subsample <- round(seq(1,nrow(fdists),length.out=4*8))
##
pdff('example_plot_2ndmeas_popsamples')
par(mfrow=c(4,8))
xticks <- tticks(xgrid)
yticks <- tticks(c(0,max(fdists[subsample,,])))
for(subs in subsample){
    par(mai=c(3,3,0,0)/12)
    matplot(x=xgrid, y=fdists[subs,,1:2], type='l',
            xlab='', ylab='', lwd=c(3,3), lty=c('solid','91'), xlim=range(xgrid), ylim=c(0,max(fdists[subsample,,])),axes=F)
        for(i in xticks){abline(v=i, lty=1, lwd=1, col='#BBBBBB60')}
    for(i in yticks){abline(h=i, lty=1, lwd=1, col='#BBBBBB60')}
    if(subs==subsample[3*8+1]){
    ##     title(xlab='x')
     axis(1, xticks, lwd=0, cex.axis=1, mgp=c(0,0,0), col.axis='#555555', gap.axis=0.25)
     axis(2, yticks, lwd=0, cex.axis=1, mgp=c(0,0,0), col.axis='#555555', gap.axis=0.25)
        text(x=xticks[length(xticks)/2],y=yticks[1]-1.75*diff(yticks)[1], 'x', cex=1.5, xpd=NA)
        text(y=yticks[length(yticks)/2],x=xticks[1]-1.3*diff(xticks)[1], 'distribution', cex=1.5, xpd=NA, srt=90)
     }
## matlines(x=xgrid, y=t(fdists[subsample,,1]), col=paste0(palette()[1],'40'),lty=1,lwd=1)
## matlines(x=xgrid, y=t(fdists[subsample,,2]), col=paste0(palette()[2],'40'),lty=1,lwd=1)
## scatteraxis(1,alldata[group==0,x],col=1)
## scatteraxis(1,alldata[group==1,x],col=2)
}
dev.off()
##

subsample2 <- round(seq(1,nrow(fdists),length.out=128*2))
pdff('example_plot_2ndmeas_popsamplestogether')
ytoplot <- rbind(t(fdists[subsample2,,1]), t(fdists[subsample2,,2]))
dim(ytoplot) <- c(length(xgrid), length(subsample2)*2)
tplot(x=xgrid, y=ytoplot, cex.lab=1.5, xlab='', ylab='', lwd=1, lty=1, col=1:2, alpha='40', xlim=range(xgrid), ylim=c(0,max(margdists2[2,,])),xaxis=F,yaxis=F)#, main='predicted distribution of full population')
#tplot(x=xgrid, y=,  lwd=1, lty=1, col=2, alpha='40', add=T)#, main='predicted distribution of full population')
dev.off()


##
pdff('example_plot_2ndmeas_A')
tplot(x=xgrid, y=preddists[,1], xlab='X', ylab='full-population distribution', lwd=4, ylim=c(0,max(margdists2)))#, main='predicted distribution of full population')
polygon(x=c(xgrid,rev(xgrid)), y=c(margdists1[1,,1],rev(margdists1[2,,1])), col=paste0(palette()[1],'40'), border=NA)
polygon(x=c(xgrid,rev(xgrid)), y=c(margdists2[1,,1],rev(margdists2[2,,1])), col=paste0(palette()[1],'40'), border=NA)
legend('topright', legend=c(
                       'estimate distribution A',
                       '50% uncertainty',
                       '75% uncertainty'
                   ),
       lty=c(1,1,1), lwd=c(3,10,10),
       col=paste0(palette()[1],c('ff','80','40')),
       bty='n', cex=1.25)
dev.off()
##
pdff('example_plot_2ndmeas_B')
tplot(x=xgrid, y=preddists[,2], xlab='X', ylab='full-population distribution', lwd=5, lty=2, col=2, ylim=c(0,max(margdists2)))#, main='predicted distribution of full population')
polygon(x=c(xgrid,rev(xgrid)), y=c(margdists1[1,,2],rev(margdists1[2,,2])), col=paste0(palette()[2],'40'), border=NA)
polygon(x=c(xgrid,rev(xgrid)), y=c(margdists2[1,,2],rev(margdists2[2,,2])), col=paste0(palette()[2],'40'), border=NA)
## scatteraxis(1,alldata[group==0,x],col=1)
## scatteraxis(1,alldata[group==1,x],col=2)
legend('topright', legend=c(
                       'estimate distribution B',
                       '50% uncertainty',
                       '75% uncertainty'
                   ),
       lty=c(2,1,1), lwd=c(3,10,10),
       col=paste0(palette()[2],c('ff','80','40')),
       bty='n', cex=1.25)
dev.off()
##

pdff('example_plot_2ndmeas_AB')
tplot(x=xgrid, y=preddists, xlab='X', ylab='distribution', lwd=4:5, lty=c('solid','91'), ylim=c(0,max(margdists2)), xlim=range(xgrid))#, main='predicted distribution of full population')
polygon(x=c(xgrid,rev(xgrid)), y=c(margdists1[1,,1],rev(margdists1[2,,1])), col=paste0(palette()[1],'40'), border=NA)
polygon(x=c(xgrid,rev(xgrid)), y=c(margdists2[1,,1],rev(margdists2[2,,1])), col=paste0(palette()[1],'40'), border=NA)
polygon(x=c(xgrid,rev(xgrid)), y=c(margdists1[1,,2],rev(margdists1[2,,2])), col=paste0(palette()[2],'40'), border=NA)
polygon(x=c(xgrid,rev(xgrid)), y=c(margdists2[1,,2],rev(margdists2[2,,2])), col=paste0(palette()[2],'40'), border=NA)
## scatteraxis(1,alldata[group==0,x],col=1)
## scatteraxis(1,alldata[group==1,x],col=2)
legend('topright', legend=c(
                       'estimate distribution A',
                       '50% uncertainty',
                       '75% uncertainty',
                       'estimate distribution B',
                       '50% uncertainty',
                       '75% uncertainty'
                   ),
       lty=c('solid','solid','solid','91','solid','solid'), lwd=c(3,10,10),
       col=paste0(palette()[c(1,1,1,2,2,2)],c('ff','80','40')),
       bty='n', cex=1.25)
##
dev.off()

pdff('new_example_plot_2ndmeas_means')
tplot(mdh1$mid, mdh1$density, col=5, alpha=0.25, xlab='mean', ylab='probability density', xlim=range(c(mdh1$mid,mdh2$mid)), lwd=6, border=NA, ylim=c(0,max(c(mdh1$density,mdh2$density))))#, main='predicted means of groups A and B')
tplot(mdh2$mid, mdh2$density, col=6, alpha=0.25, xlab='mean group B', ylab='probability', main='predicted mean of group A', xlim=range(c(mdh1$mid,mdh2$mid)), ylim=c(0,max(c(mdh1$density,mdh2$density))), lty='82', lwd=4, border=NA, add=T)
legend('topright', legend=c(
                       'uncertainty of mean for group A',
                       'uncertainty of mean for group B'
                   ),
       lty=c('solid','31'), lwd=c(5,4),
       col=c(5,6),
       bty='n', cex=1.25)
dev.off()

##
pdff('example_plot_2ndmeas_meansdiff')
tplot(mdh$mid, mdh$density, col=3, lwd=5, xlab='(mean group A) - (mean group B)', ylab='probability')#, main='predicted difference between means of groups A and B')
legend('topright', legend=paste0(
                       'estimated difference between means: ',signif(mean(meandiffs),2),
                       ## signif(quant(meandiffs,c(1)/8)*2,2),' < difference < ',signif(quant(meandiffs,c(7)/8)*2,2),'\nwith 75% probability','\n\n',
                                 '\n\n',signif(quant(meandiffs,c(2.5)/100),2),' < difference < ',signif(quant(meandiffs,c(97.5)/100),2),'\nwith 95% probability\n\n',
                       'difference > 0\nwith ',signif(sum(meandiffs>0)/length(meandiffs)*100,2),'% probability'
                                 ), bty='n', cex=1.5)
dev.off()
##
pdff('example_plot_2ndmeas_vars')
tplot(vdh1$mid, vdh1$density, col=1, xlab='variance', ylab='probability', main='predicted variances of groups A and B', #xlim=quant(vdists, c(2.5,97.5)/100),
      ylim=c(0,max(c(vdh1$density,vdh2$density))))
tplot(vdh2$mid, vdh2$density, col=2, xlab='mean group B', ylab='probability', main='predicted mean of group A', xlim=range(c(vdh1$mid,vdh2$mid)), ylim=c(0,max(c(vdh1$density,vdh2$density))), add=T)
dev.off()
##

vdh <- thist(vardiffs,n=256*2)
pdff('example_plot_2ndmeas_varsdiff')
tplot(vdh$mid, vdh$density, col=4, lwd=5, xlab='(variance group A) - (variance group B)', ylab='probability', xlim=quant(vardiffs,c(1,99)/100))#, main='predicted difference between variances of groups A and B') 
legend(x=-280,y=0.013, legend=paste0(
                       'estimated difference between variances: ',signif(mean(vardiffs),2),'\u00b1',signif(sd(vardiffs),2),
                       ## signif(quant(meandiffs,c(1)/8)*2,2),' < difference < ',signif(quant(meandiffs,c(7)/8)*2,2),'\nwith 75% probability','\n\n',
                                 '\n\n',signif(quant(vardiffs,c(2.5)/100),2),' < difference < ',signif(quant(vardiffs,c(97.5)/100),2),'\nwith 95% probability\n\n',
                       'difference > 0\nwith ',signif(sum(vardiffs>0)/length(vardiffs)*100,2),'% probability'
                                 ), bty='n', cex=1.5)
dev.off()


for(g in 1:2){
print(paste0(signif(quant(edists[,g],c(2.5)/100),2),' < ',mean(edists[,g]),' < ',signif(quant(edists[,g],c(97.5)/100),2)))
}


## hi <- foreach(g=1:2)%do%{
##     thist(edists[,g], n=tticks(c(edists),n=16))
## }
## tplot(sapply(hi,function(x){x$breaks}), sapply(hi, function(x){x$density}), xlim=range(c(edists)))

## cdists <- samplescdfRgivenX(maincov='x', X=cbind(group=0:1), parmList=parmList, covgrid=xgrid)
## cdists <- aperm(cdists)
## qinter <- foreach(g=1:2)%do%{ approxfun(colMeans(cdists[,,g]), xgrid) }




















xrange <- extendrange(alldata$x)
xrangemin <- range(alldata$x)
xgrid <- seq(xrange[1], xrange[2], length.out=256)
xtest <- round(quant(alldata$x, (1:3)/4)*2)/2

egx1 <- cbind(group=1)
fdists <- samplesfRgivenX(maincov='x', X=cbind(group=1:2), parmList=parmList, covgrid=xgrid, inorder=T)
fdists <- aperm(fdists)
preddists <- colMeans(fdists)
margdists1 <- apply(fdists,c(2,3),function(x){quant(x, c(1,3)/4)})
margdists2 <- apply(fdists,c(2,3),function(x){quant(x, c(1,7)/8)})

edists <- samplesmeanRgivenX(maincov='x', X=cbind(group=1:2), parmList=parmList, inorder=T)
edists <- aperm(edists)
meandiffs <- apply(edists,1,function(x){-diff(x)})
mdh <- thist(meandiffs)

vdists <- samplesvarRgivenX(maincov='x', X=cbind(group=1:2), parmList=parmList, inorder=T)
vdists <- aperm(vdists)
vardiffs <- apply(vdists,1,function(x){-diff(x)})
vdh <- thist(vardiffs)


hi <- foreach(g=1:2)%do%{
    thist(edists[,g], n=tticks(c(edists),n=16))
}
tplot(sapply(hi,function(x){x$mid}), sapply(hi, function(x){x$density}), xlim=range(c(edists)))

cdists <- samplescdfRgivenX(maincov='x', X=cbind(group=1:2), parmList=parmList, covgrid=xgrid)
cdists <- aperm(cdists)
qinter <- foreach(g=1:2)%do%{ approxfun(colMeans(cdists[,,g]), xgrid) }


subsample <- seq(1,nrow(fdists),length.out=64)

pdff('example_plot')
tplot(x=xgrid, y=preddists, xlab='x', ylab='full-population distribution', lwd=2:3, ylim=c(0,max(margdists)), main='predicted distribution of full population')
polygon(x=c(xgrid,rev(xgrid)), y=c(margdists1[1,,1],rev(margdists1[2,,1])), col=paste0(palette()[1],'40'), border=NA)
polygon(x=c(xgrid,rev(xgrid)), y=c(margdists2[1,,1],rev(margdists2[2,,1])), col=paste0(palette()[1],'40'), border=NA)
polygon(x=c(xgrid,rev(xgrid)), y=c(margdists1[1,,2],rev(margdists1[2,,2])), col=paste0(palette()[2],'40'), border=NA)
polygon(x=c(xgrid,rev(xgrid)), y=c(margdists2[1,,2],rev(margdists2[2,,2])), col=paste0(palette()[2],'40'), border=NA)
scatteraxis(1,alldata[group==1,x],col=1)
scatteraxis(1,alldata[group==2,x],col=2)
legend('topright', legend=c(
                   ))
##
tplot(mdh$breaks, mdh$density, col=3, xlab='(mean group A) - (mean group B)', ylab='probability', main='predicted difference between means of groups A and B')
legend('topright', legend=paste0(
                       'expected difference: ',signif(mean(meandiffs),2),
                       ## signif(quant(meandiffs,c(1)/8)*2,2),' < difference < ',signif(quant(meandiffs,c(7)/8)*2,2),'\nwith 75% probability','\n\n',
                                 '\n\n',signif(quant(meandiffs,c(5)/100),2),' < difference < ',signif(quant(meandiffs,c(95)/100),2),'\nwith 90% probability\n\n',
                       'difference > 0\nwith ',signif(sum(meandiffs>0)/length(meandiffs)*100,2),'% probability'
                                 ), bty='n', cex=1.5)
##
tplot(vdh$breaks, vdh$density, col=4, xlab='(variance group A) - (variance group B)', ylab='probability', main='predicted difference between variances of groups A and B')
legend('topright', legend=paste0(
                       'expected difference: ',signif(mean(vardiffs),2),'\u00b1',signif(sd(vardiffs),2),
                       ## signif(quant(meandiffs,c(1)/8)*2,2),' < difference < ',signif(quant(meandiffs,c(7)/8)*2,2),'\nwith 75% probability','\n\n',
                                 '\n\n',signif(quant(vardiffs,c(5)/100),2),' < difference < ',signif(quant(vardiffs,c(95)/100),2),'\nwith 90% probability\n\n',
                       'difference > 0\nwith ',signif(sum(vardiffs>0)/length(vardiffs)*100,2),'% probability'
                                 ), bty='n', cex=1.5)
##
tplot(x=xgrid, y=preddists, xlab='x', ylab='full-population distribution', lwd=2:3, ylim=c(0,max(fdists)/2))
matlines(x=xgrid, y=t(fdists[subsample,,1]), col=paste0(palette()[1],'40'),lty=1,lwd=1)
matlines(x=xgrid, y=t(fdists[subsample,,2]), col=paste0(palette()[2],'40'),lty=1,lwd=1)
scatteraxis(1,alldata[group==1,x],col=1)
scatteraxis(1,alldata[group==2,x],col=2)
dev.off()


pdff('example_statement')
mean1 <- xgrid %*% normalize(meandist1)
mean2 <- xgrid %*% normalize(meandist2)
    ##
    ymax=max(egygivenx1[,subsamples], egygivenx2[,subsamples])*1.2
myplot(ygrid*2+20, egygivenx1[,subsamples], type='l', col=paste0(palette()[5],'44'), lty=1, lwd=2, xlab='Y', ylab='predicted frequency in full population', ylim=c(0,ymax), cex.axis=1.5, cex.lab=1.5, bty="n")
    matplot(ygrid*2+20, egygivenx2[,subsamples], type='l', col=paste0(palette()[2],'44'), lty=1, lwd=2, xlab='Y', ylab='x', add=TRUE)
    matlines(ygrid*2+20, meandist1, type='l', col=1, lty=1, lwd=4)
matlines(ygrid*2+20, meandist2, type='l', col=6, lty=2, lwd=4)
seqpol <- seq(qt1[1],qt1[3],length.out=64)
polygon(x=c(qt1[1], seqpol, qt1[3])*2+20,
        y=c(0, pinter1(seqpol), 0), col=paste0(palette()[5],'33'), border=NA)
seqpol <- seq(qt2[1],qt2[3],length.out=64)
polygon(x=c(qt2[1], seqpol, qt2[3])*2+20,
        y=c(0, pinter2(seqpol), 0), col=paste0(palette()[2],'33'), border=NA)
##grid(lwd=0.5,lty=1, col=paste0(palette()[7],'88'))
    legend('topright', legend=paste0('x = ',xtest[3]*2+20), cex=2, bty='n')
    legend('topleft', legend=c(paste0('75% of full population A has ',signif(qt1[1]*2+20,3),' < Y < ',signif(qt1[3]*2+20,3)),
                               paste0('75% of full population B has ',signif(qt2[1]*2+20,3),' < Y < ',signif(qt2[3]*2+20,3))
                               ), lty=c(1,2), lwd=3, col=c(1,2), bty='n', cex=1.5)
##
hist((medians1-medians2)*2, breaks=10, xlab='(median Y group A) - (median Y group B)', ylab='probability', freq=F, col=3, cex.axis=1.5, cex.lab=1.5, border=NA, main='Forecast of difference between medians of full populations A and B', xlim=c(-4,4))#max(medians1-medians2)*2.1))
axis(1, at=seq(1,3,by=1), cex.axis=1.5)
grid(lwd=0.5,lty=1, col=paste0(palette()[7],'88'))
abline(v=0, lwd=2, col=paste0(palette()[7],'88'))
legend('topleft', legend=paste0(signif(quant(medians1-medians2,c(1)/8)*2,2),' < difference < ',signif(quant(medians1-medians2,c(7)/8)*2,2),'\nwith 75% probability','\n\n',
                signif(quant(medians1-medians2,c(5)/100)*2,2),' < difference < ',signif(quant(medians1-medians2,c(95)/100)*2,2),'\nwith 90% probability'                ), bty='n', cex=1.5)
dev.off()




## Calculations of quartiles and octiles
egx1 <- cbind(x=xtest[3],group=1)
egygivenx1 <- samplesfRgivenX(maincov='Y', X=egx1, parmList=parmList, covgrid=ygrid)[1,,]
pinter1 <- approxfun(ygrid, rowMeans(egygivenx1))
cdfygivenx1 <- samplescdfRgivenX(maincov='Y', X=egx1, parmList=parmList, covgrid=ygrid)[1,,]
qinter1 <- approxfun(rowMeans(cdfygivenx1), ygrid)
##
egx2 <- cbind(x=xtest[3],group=2)
egygivenx2 <- samplesfRgivenX(maincov='Y', X=egx2, parmList=parmList, covgrid=ygrid)[1,,]
pinter2 <- approxfun(ygrid, rowMeans(egygivenx2))
cdfygivenx2 <- samplescdfRgivenX(maincov='Y', X=egx2, parmList=parmList, covgrid=ygrid)[1,,]
qinter2 <- approxfun(rowMeans(cdfygivenx2), ygrid)

##############################################
## Given x, how distinguishable are the 75% ranges and the medians of Y for the two groups?
meandist1 <- rowMeans(egygivenx1)
qt1 <- qinter1(c(1,4,7)/8)
meandist2 <- rowMeans(egygivenx2)
qt2 <- qinter2(c(1,4,7)/8)
medians1 <- foreach(i=1:niter, .combine=c)%dopar%{
    approxfun(cdfygivenx1[,i], ygrid)(0.5)
}
medians2 <- foreach(i=1:niter, .combine=c)%dopar%{
    approxfun(cdfygivenx2[,i], ygrid)(0.5)
}
##
subsamples <- round(seq(1, niter, length.out=32))
pdff('few_example_statement')
mean1 <- ygrid %*% normalize(meandist1)
mean2 <- ygrid %*% normalize(meandist2)
    ##
    ymax=max(egygivenx1[,subsamples], egygivenx2[,subsamples])*1.2
myplot(ygrid*2+20, egygivenx1[,subsamples], type='l', col=paste0(palette()[5],'44'), lty=1, lwd=2, xlab='Y', ylab='predicted frequency in full population', ylim=c(0,ymax), cex.axis=1.5, cex.lab=1.5, bty="n")
    matplot(ygrid*2+20, egygivenx2[,subsamples], type='l', col=paste0(palette()[2],'44'), lty=1, lwd=2, xlab='Y', ylab='x', add=TRUE)
    matlines(ygrid*2+20, meandist1, type='l', col=1, lty=1, lwd=4)
matlines(ygrid*2+20, meandist2, type='l', col=6, lty=2, lwd=4)
seqpol <- seq(qt1[1],qt1[3],length.out=64)
polygon(x=c(qt1[1], seqpol, qt1[3])*2+20,
        y=c(0, pinter1(seqpol), 0), col=paste0(palette()[5],'33'), border=NA)
seqpol <- seq(qt2[1],qt2[3],length.out=64)
polygon(x=c(qt2[1], seqpol, qt2[3])*2+20,
        y=c(0, pinter2(seqpol), 0), col=paste0(palette()[2],'33'), border=NA)
##grid(lwd=0.5,lty=1, col=paste0(palette()[7],'88'))
    legend('topright', legend=paste0('x = ',xtest[3]*2+20), cex=2, bty='n')
    legend('topleft', legend=c(paste0('75% of full population A has ',signif(qt1[1]*2+20,3),' < Y < ',signif(qt1[3]*2+20,3)),
                               paste0('75% of full population B has ',signif(qt2[1]*2+20,3),' < Y < ',signif(qt2[3]*2+20,3))
                               ), lty=c(1,2), lwd=3, col=c(1,2), bty='n', cex=1.5)
##
hist((medians1-medians2)*2, breaks=10, xlab='(median Y group A) - (median Y group B)', ylab='probability', freq=F, col=3, cex.axis=1.5, cex.lab=1.5, border=NA, main='Forecast of difference between medians of full populations A and B', xlim=c(-4,4))#max(medians1-medians2)*2.1))
axis(1, at=seq(1,3,by=1), cex.axis=1.5)
grid(lwd=0.5,lty=1, col=paste0(palette()[7],'88'))
abline(v=0, lwd=2, col=paste0(palette()[7],'88'))
legend('topleft', legend=paste0(signif(quant(medians1-medians2,c(1)/8)*2,2),' < difference < ',signif(quant(medians1-medians2,c(7)/8)*2,2),'\nwith 75% probability','\n\n',
                signif(quant(medians1-medians2,c(5)/100)*2,2),' < difference < ',signif(quant(medians1-medians2,c(95)/100)*2,2),'\nwith 90% probability'                ), bty='n', cex=1.5)
dev.off()




##############################################
## Given x, how well does Y predict the group?
testx1 <- cbind(x=xtest,group=1)
samplesygivenx1 <- samplesfRgivenX(maincov='Y', X=testx1, parmList=parmList, covgrid=ygrid)
##
testx2 <- cbind(x=xtest,group=2)
samplesygivenx2 <- samplesfRgivenX(maincov='Y', X=testx2, parmList=parmList, covgrid=ygrid)
##
subsamples <- round(seq(1, niter, length.out=64))

pdff('few_testplots_ygivenx')
for(xvalue in 1:length(xtest)){
    meandist1 <- rowMeans(samplesygivenx1[xvalue,,])
    meandist2 <- rowMeans(samplesygivenx2[xvalue,,])
    mean1 <- ygrid %*% normalize(meandist1)
    mean2 <- ygrid %*% normalize(meandist2)
    ##
    ymax=max(samplesygivenx1[,,subsamples],samplesygivenx2[,,subsamples])
    myplot(ygrid*2+20, samplesygivenx1[xvalue,,subsamples], type='l', col=paste0(palette()[5],'44'), lty=1, lwd=2, xlab='Y', ylab='predicted frequency in full population', ylim=c(0,ymax), main='predicting Y given x')
    tplot(ygrid*2+20, samplesygivenx2[xvalue,,subsamples], type='l', col=paste0(palette()[2],'44'), lty=1, lwd=2, xlab='Y', ylab='x', add=TRUE)
    matlines(ygrid*2+20, meandist1, type='l', col=1, lty=1, lwd=4)
    matlines(ygrid*2+20, meandist2, type='l', col=6, lty=2, lwd=4)
#    grid(lwd=0.5,lty=1, col=paste0(palette()[7],'88'))
    legend('topright', legend=paste0('x = ', signif(testx1[xvalue]*2+20,3)), cex=2, bty='n')
    legend('topleft', legend=c(paste0('group A'),#, mean(y) = ', signif(mean1,2)),
                               paste0('group B')#', mean(y) = ', signif(mean2,2))
                               ), lty=c(1,2), lwd=3, col=c(1,2), bty='n', cex=1.5)
}
dev.off()

##
##
##############################################
## Given Y, how well does x predict the group?
testy1 <- cbind(Y=ytest,group=1)
samplesxgiveny1 <- samplesfRgivenX(maincov='x', X=testy1, parmList=parmList, covgrid=xgrid)
##
testy2 <- cbind(Y=ytest,group=2)
samplesxgiveny2 <- samplesfRgivenX(maincov='x', X=testy2, parmList=parmList, covgrid=xgrid)
##
subsamples <- round(seq(1, niter, length.out=64))
pdff('testplots_xgiveny')
for(yvalue in 1:length(ytest)){
    meandist1 <- rowMeans(samplesxgiveny1[yvalue,,])
    meandist2 <- rowMeans(samplesxgiveny2[yvalue,,])
    mean1 <- xgrid %*% normalize(meandist1)
    mean2 <- xgrid %*% normalize(meandist2)
    ##
    xmax=max(samplesxgiveny1[,,subsamples],samplesxgiveny2[,,subsamples])
    myplot(xgrid*2+20, samplesxgiveny1[yvalue,,subsamples], type='l', col=paste0(palette()[5],'44'), lty=1, lwd=2, xlab='x', ylab='predicted frequency in full population', ylim=c(0,xmax), main='predicting x given Y')
    matplot(xgrid*2+20, samplesxgiveny2[yvalue,,subsamples], type='l', col=paste0(palette()[2],'44'), lty=1, lwd=2, xlab='x', ylab='Y', add=TRUE)
    matlines(xgrid*2+20, meandist1, type='l', col=1, lty=1, lwd=4)
    matlines(xgrid*2+20, meandist2, type='l', col=6, lty=2, lwd=4)
##        grid(lwd=0.5,lty=1, col=paste0(palette()[7],'88'))
    legend('topright', legend=paste0('Y = ', signif(testy1[yvalue]*2+20,3)), cex=2, bty='n')
    legend('topleft', legend=c(paste0('group A'),#, mean(x) = ', signif(mean1,2)),
                               paste0('group B')#', mean(x) = ', signif(mean2,2))
                               ), lty=c(1,2), lwd=3, col=c(1,2), bty='n', cex=1.5)
}
dev.off()


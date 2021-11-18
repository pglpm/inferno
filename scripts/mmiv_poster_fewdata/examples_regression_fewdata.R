## Author: PGL  Porta Mana
## Created: 2021-03-20T10:07:17+0100
## Last-Updated: 2021-11-18T13:35:33+0100
################
## Exploration for MMIV poster
################
if(file.exists("/cluster/home/pglpm/R")){
    .libPaths(c("/cluster/home/pglpm/R",.libPaths()))
}
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
options(bitmapType='cairo')
pdff <- function(filename){pdf(file=paste0(filename,'.pdf'),paper='a4r',height=11.7,width=16.5)} # to output in pdf format
pngf <- function(filename,res=300){png(file=paste0(filename,'.png'),height=11.7*1.2,width=16.5,units='in',res=res,pointsize=36)} # to output in pdf format
library('nimble')
## library('coda')
#### End custom setup ####


## Create example data
set.seed(123)
noise <- 1
dxfromy <- function(x,y){dnorm((x-abs(y)^(1.5)*sign(y))/(noise+1/((y^2-1)^2+1)))}
rxfromy <- function(n,y){
    rnorm(n=n)*(noise+1/((y^2-1)^2+1)) + abs(y)^(1.5)*sign(y)
}
nsamples1 <- 100
ysamples1 <- rnorm(nsamples1, mean=0.5, sd=0.5)
xsamples1 <- rxfromy(nsamples1, ysamples1)
nsamples2 <- 75
ysamples2 <- rnorm(nsamples2, mean=-0.5, sd=0.5)
xsamples2 <- rxfromy(nsamples2, ysamples2)
## matplot(xsamples1, ysamples1, type='p', pch=1, cex=0.5, col=1, xlim=range(c(xsamples1,xsamples2)), ylim=range(c(ysamples1,ysamples2)))#, xlim=c(-2,3))
## matplot(xsamples2, ysamples2, type='p', pch=1, cex=0.5, col=2, add=T)
##
pdff('few_alldata')
matplot(c(xsamples1,xsamples2)*2+20, c(ysamples1,ysamples2)*2+20, type='p', pch=8, cex=0.5, cex.axis=1.5, cex.lab=1.5, col='black', xlim=range(c(xsamples1,xsamples2)*2+20), ylim=range(c(ysamples1,ysamples2)*2+20), xlab='T', ylab='Y', main='all data')
    grid(lwd=0.5,lty=1, col=paste0(palette()[7],'88'))
##
matplot(xsamples1*2+20, ysamples1*2+20, type='p', pch=1, cex=1, cex.axis=1.5, cex.lab=1.5, col=1, xlim=range(c(xsamples1,xsamples2)*2+20), ylim=range(c(ysamples1,ysamples2)*2+20), xlab='T', ylab='Y', main='data according to group')
    legend('topleft', legend=c('group A', 'group B'), pch=c(1,3), col=c(1,2), bty='n', cex=2)
matpoints(xsamples2*2+20, ysamples2*2+20, type='p', pch=3, cex=1, cex.axis=1.5, cex.lab=1.5, col=2, xlim=range(c(xsamples1,xsamples2)*2+20), ylim=range(c(ysamples1,ysamples2)*2+20), xlab='T', ylab='Y', main='all data')
    grid(lwd=0.5,lty=1, col=paste0(palette()[7],'88'))
matplot(xsamples1*2+20, ysamples1*2+20, type='p', pch=1, cex=1, cex.axis=1.5, cex.lab=1.5, col=1, xlim=range(c(xsamples1,xsamples2)*2+20), ylim=range(c(ysamples1,ysamples2)*2+20), xlab='T', ylab='Y')
    legend('topleft', legend=c('group A'), pch=c(1), col=c(1), bty='n', cex=2)
    grid(lwd=0.5,lty=1, col=paste0(palette()[7],'88'))
matplot(xsamples2*2+20, ysamples2*2+20, type='p', pch=3, cex=1, cex.axis=1.5, cex.lab=1.5, col=2, xlim=range(c(xsamples1,xsamples2)*2+20), ylim=range(c(ysamples1,ysamples2)*2+20), xlab='T', ylab='Y')
    legend('topleft', legend=c('group B'), pch=c(3), col=c(2), bty='n', cex=2)
    grid(lwd=0.5,lty=1, col=paste0(palette()[7],'88'))
dev.off()


## Put data together into a data table
oalldata <- rbind(
    data.table(T=xsamples1, Y=ysamples1, group=1L),
    data.table(T=xsamples2, Y=ysamples2, group=2L)
)


seed <- 149
baseversion <- 'regr-mmivposter-C21D11_'
nclusters <- as.integer(32) #as.integer(2^6)
ndata <- as.integer(nrow(alldata))
niter <- as.integer(2^11) #as.integer(2^11)
niter0 <- as.integer(2^10) #as.integer(2^10)
nstages <- as.integer(1)
ncheckpoints <- as.integer(8)
covNames <-  c('T'
              ,'Y'
              ,'group'
               )
posterior <- TRUE

## initial.options <- commandArgs(trailingOnly = FALSE)
## thisscriptname <- sub('--file=', "", initial.options[grep('--file=', initial.options)])
## file.copy(from=thisscriptname, to=paste0('script-R',baseversion,'-V',length(covNames),'-D',ndata,'-K',nclusters,'.Rscript'))

## pdff('check_mutualinfo')
## for(i in 1:(length(covNames)-1)){
##     for(j in (i+1):length(covNames)){
##         matplot(alldata[[covNames[i]]], alldata[[covNames[j]]], type='p', pch='.', col=paste0('#000000','88'), xlab=covNames[i], ylab=covNames[j])
##     }
## }
## dev.off()


##
## alldata <- fread(datafile, sep=',')
discreteCovs <- covNames[sapply(covNames, function(x){is.integer(oalldata[[x]])})]
continuousCovs <- covNames[sapply(covNames, function(x){is.double(oalldata[[x]])})]
covNames <- c(continuousCovs, discreteCovs)
nccovs <- length(continuousCovs)
ndcovs <- length(discreteCovs)
##
alldata <- oalldata[1:ndata,..covNames]


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
        dataQuantiles[[acov]] <- quantile(alldata[[acov]], prob=c(0.005,0.995))
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
    tauCshape1=rep(1/2, nccovs),
    tauCrate2=1/(widthccovs/2)^2, # dims = inv. variance
    ##
    probDa1=rep(1, ndcovs),
    probDb1=rep(1, ndcovs),
    sizeDprob1=rep(1/2, ndcovs),
    sizeDsize1=maxdcovs,
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
    sizeD=matrix(apply(matrix(rnbinom(n=ndcovs*(nclusters), prob=rbeta(n=ndcovs, shape1=16, shape2=16), size=maxdcovs), nrow=ndcovs, ncol=nclusters), 2, function(x){maxdcovs*(x<maxdcovs)+x*(x>=maxdcovs)}), nrow=ndcovs, ncol=nclusters),
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
            sizeD[acov,acluster] ~ dnbinom(prob=sizeDprob1[acov], size=sizeDsize1[acov])
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
                Y[adatum,acov] ~ dbinom(prob=probD[acov,C[adatum]], size=sizeD[acov,C[adatum]])
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
            confmodel$addSampler(target=paste0('sizeD[', acov, ', ', acluster, ']'), type='slice')
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

##
## source('functions_rmsdregr_nimble_binom.R')

print('Setup time:')
print(Sys.time() - timecount)

##
nstages <- 3
for(stage in 0:nstages){
    totalruntime <- Sys.time()

    print(paste0('==== STAGE ', stage, ' ===='))
    version <- paste0(baseversion, stage)
    gc()
    if(stage==0){
        mcsamples <- runMCMC(Cmcmcsampler, nburnin=1, niter=niter0+1, thin=1, thin2=niter0, inits=initsFunction, setSeed=seed)
    }else{
        Cmcmcsampler$run(niter=niter*2, thin=2, thin2=niter*2, reset=FALSE, resetMV=TRUE)
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
    finalstate <- c(mcsamples[nrow(mcsamples),], finalstate[nrow(finalstate),])
    occupations <- finalstate[grepl('^C\\[', names(finalstate))]
    usedclusters <- length(unique(occupations))
    print(paste0('OCCUPIED CLUSTERS: ', usedclusters, ' OF ', nclusters))
    saveRDS(finalstate2list(finalstate),file=paste0('_finalstate-R',version,'-V',length(covNames),'-D',ndata,'-K',nclusters,'-I',nrow(mcsamples),'.rds'))
    ##
    parmList <- mcsamples2parmlist(mcsamples)
    ## Traces to follow for diagnostics
    ll <- llSamples(dat, parmList)
    momentstraces <- moments12Samples(parmList)
    miqrtraces <- calcSampleMQ(parmList)
    probCheckpoints <- t(probValuesSamples(checkpoints, parmList))
    medians <- miqrtraces[,,1]
    colnames(medians) <- paste0('MEDIAN_', colnames(miqrtraces))
    Q1s <- miqrtraces[,,2]
    colnames(Q1s) <- paste0('Q1_', colnames(miqrtraces))
    Q3s <- miqrtraces[,,3]
    colnames(Q3s) <- paste0('Q3_', colnames(miqrtraces))
    iqrs <- Q3s - Q1s
    colnames(iqrs) <- paste0('IQR_', colnames(miqrtraces))
    traces <- cbind(LL=ll, probCheckpoints, medians, iqrs, Q1s, Q3s, do.call(cbind, momentstraces))
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
    ##
    ## tracegroups <- list(
    ##     'maxD'=(1:ncol(traces))[grepl('^(Pdat|Dcov|LL)', tracenames)],
    ##     '1D'=(1:ncol(traces))[grepl('^(MEDIAN|Q1|Q3|IQR|MEAN|VAR)_', tracenames)],
    ##     '2D'=(1:ncol(traces))[grepl('^COV_', tracenames)]
    ## )
    ## grouplegends <- foreach(agroup=1:length(tracegroups))%do%{
    ##     c( paste0('-- STATS ', names(tracegroups)[agroup], ' --'),
    ##       paste0('min ESS = ', min(diagnESS[tracegroups[[agroup]]])),
    ##       paste0('max BMK = ', max(diagnBMK[tracegroups[[agroup]]])),
    ##       paste0('max MCSE = ', max(diagnMCSE[tracegroups[[agroup]]])),
    ##       paste0('all stationary: ', all(diagnStat[tracegroups[[agroup]]])),
    ##       paste0('burn: ', max(diagnBurn[tracegroups[[agroup]]]))
    ##       )
    ## }
    ## ##
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
        matplot(Xgrid, plotsamples, type='l', col=paste0(palette()[7], '44'), lty=1, lwd=2, xlab=acov, ylab='probability density', cex.axis=1.5, cex.lab=1.5, ylim=c(0, max(plotsamples[plotsamples<df])))
    }
    for(acov in discreteCovs){
        Xgrid <- seq(min(alldata[[acov]]), max(alldata[[acov]]), by=1)
        plotsamples <- samplesfX(acov, parmList, Xgrid, nfsamples=64)
        matplot(Xgrid, plotsamples, type='l', col=paste0(palette()[7], '44'), lty=1, lwd=2, xlab=acov, ylab='probability density', cex.axis=1.5, cex.lab=1.5, ylim=c(0, max(plotsamples)))
    }
    ##
    par(mfrow = c(2, 2))
    for(maincov in c('T', 'Y')){
    for(addvar in setdiff(covNames, maincov)){
        matplot(x=c(rep(alldataRanges[[maincov]], each=2),
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
        if(grepl('^[PDV]', acov)){transf <- function(x){log(abs(x))}
        }else{transf <- identity}
        matplot(transf(traces[,acov]), type='l', lty=1, col=colpalette[acov],
                main=paste0(acov,
                            '\nESS = ', signif(diagnESS[acov], 3),
                            ' | BMK = ', signif(diagnBMK[acov], 3),
                            ' | MCSE(6.27) = ', signif(diagnMCSE[acov], 3),
                            ' | stat: ', diagnStat[acov],
                            ' | burn: ', diagnBurn[acov]
                            ),
                ylab=acov,
                ylim=range(c(transf(traces[,acov][abs(transf(traces[,acov]))<Inf]))))
    }
    dev.off()

    print('Total runtime:')
    print(Sys.time() - totalruntime)

}


xrange <- range(alldata$T)
xrangemin <- range(alldata$T)
xgrid <- seq(xrange[1], xrange[2], length.out=256)
xtest <- round(quantile(alldata$T, (1:3)/4)*2)/2
##
yrange <- range(alldata$Y)
yrangemin <- range(alldata$Y)
ygrid <- seq(yrange[1], yrange[2], length.out=256+1)
ytest <- round(quantile(alldata$Y, (1:3)/4)*2)/2
##

##
testx1 <- cbind(T=xtest,group=1)
samplesygivenx1 <- samplesfRgivenX(maincov='Y', X=testx1, parmList=parmList, covgrid=ygrid)
##
testx2 <- cbind(T=xtest,group=2)
samplesygivenx2 <- samplesfRgivenX(maincov='Y', X=testx2, parmList=parmList, covgrid=ygrid)
##
subsamples <- round(seq(1, niter, length.out=64))
pdff('testplots_ygivenx')
for(xvalue in 1:length(xtest)){
    meandist1 <- rowMeans(samplesygivenx1[xvalue,,])
    meandist2 <- rowMeans(samplesygivenx2[xvalue,,])
    mean1 <- ygrid %*% normalize(meandist1)
    mean2 <- ygrid %*% normalize(meandist2)
    ##
    ymax=max(samplesygivenx1[,,subsamples],samplesygivenx2[,,subsamples])
    matplot(ygrid*2+20, samplesygivenx1[xvalue,,subsamples], type='l', col=paste0(palette()[5],'44'), lty=1, lwd=2, xlab='Y', ylab='predicted frequency in full population', ylim=c(0,ymax), cex.axis=1.5, cex.lab=1.5, , bty="n", main='predicting Y given T')
    matplot(ygrid*2+20, samplesygivenx2[xvalue,,subsamples], type='l', col=paste0(palette()[2],'44'), lty=1, lwd=2, xlab='Y', ylab='T', add=TRUE)
    matlines(ygrid*2+20, meandist1, type='l', col=1, lty=1, lwd=4)
    matlines(ygrid*2+20, meandist2, type='l', col=6, lty=2, lwd=4)
    grid(lwd=0.5,lty=1, col=paste0(palette()[7],'88'))
    legend('topright', legend=paste0('T = ', signif(testx1[xvalue]*2+20,3)), cex=2, bty='n')
    legend('topleft', legend=c(paste0('group A'),#, mean(y) = ', signif(mean1,2)),
                               paste0('group B')#', mean(y) = ', signif(mean2,2))
                               ), lty=c(1,2), lwd=3, col=c(1,2), bty='n', cex=1.5)
}
dev.off()
##
##
testy1 <- cbind(Y=ytest,group=1)
samplesxgiveny1 <- samplesfRgivenX(maincov='T', X=testy1, parmList=parmList, covgrid=xgrid)
##
testy2 <- cbind(Y=ytest,group=2)
samplesxgiveny2 <- samplesfRgivenX(maincov='T', X=testy2, parmList=parmList, covgrid=xgrid)
##
subsamples <- round(seq(1, niter, length.out=64))
pdff('testplots_xgiveny')
for(yvalue in 1:length(ytest)){
    meandist1 <- rowMeans(samplesxgiveny1[yvalue,,])
    meandist2 <- rowMeans(samplesxgiveny2[yvalue,,])
    mean1 <- xgrid %*% normalize(meandist1)
    sd1 <- 
    mean2 <- xgrid %*% normalize(meandist2)
    ##
    xmax=max(samplesxgiveny1[,,subsamples],samplesxgiveny2[,,subsamples])
    matplot(xgrid*2+20, samplesxgiveny1[yvalue,,subsamples], type='l', col=paste0(palette()[5],'44'), lty=1, lwd=2, xlab='T', ylab='predicted frequency in full population', ylim=c(0,xmax), cex.axis=1.5, cex.lab=1.5, , bty="n", main='predicting T given Y')
    matplot(xgrid*2+20, samplesxgiveny2[yvalue,,subsamples], type='l', col=paste0(palette()[2],'44'), lty=1, lwd=2, xlab='T', ylab='Y', add=TRUE)
    matlines(xgrid*2+20, meandist1, type='l', col=1, lty=1, lwd=4)
    matlines(xgrid*2+20, meandist2, type='l', col=6, lty=2, lwd=4)
        grid(lwd=0.5,lty=1, col=paste0(palette()[7],'88'))
    legend('topright', legend=paste0('Y = ', signif(testy1[yvalue]*2+20,3)), cex=2, bty='n')
    legend('topleft', legend=c(paste0('group A'),#, mean(x) = ', signif(mean1,2)),
                               paste0('group B')#', mean(x) = ', signif(mean2,2))
                               ), lty=c(1,2), lwd=3, col=c(1,2), bty='n', cex=1.5)
}
dev.off()

egx1 <- cbind(T=xtest[3],group=1)
egygivenx1 <- samplesfRgivenX(maincov='Y', X=egx1, parmList=parmList, covgrid=ygrid)[1,,]
pinter1 <- approxfun(ygrid, rowMeans(egygivenx1))
cdfygivenx1 <- samplescdfRgivenX(maincov='Y', X=egx1, parmList=parmList, covgrid=ygrid)[1,,]
qinter1 <- approxfun(rowMeans(cdfygivenx1), ygrid)
##
egx2 <- cbind(T=xtest[3],group=2)
egygivenx2 <- samplesfRgivenX(maincov='Y', X=egx2, parmList=parmList, covgrid=ygrid)[1,,]
pinter2 <- approxfun(ygrid, rowMeans(egygivenx2))
cdfygivenx2 <- samplescdfRgivenX(maincov='Y', X=egx2, parmList=parmList, covgrid=ygrid)[1,,]
qinter2 <- approxfun(rowMeans(cdfygivenx2), ygrid)

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
matplot(ygrid*2+20, egygivenx1[,subsamples], type='l', col=paste0(palette()[5],'44'), lty=1, lwd=2, xlab='Y', ylab='predicted frequency in full population', ylim=c(0,ymax), cex.axis=1.5, cex.lab=1.5, bty="n")
    matplot(ygrid*2+20, egygivenx2[,subsamples], type='l', col=paste0(palette()[2],'44'), lty=1, lwd=2, xlab='Y', ylab='T', add=TRUE)
    matlines(ygrid*2+20, meandist1, type='l', col=1, lty=1, lwd=4)
matlines(ygrid*2+20, meandist2, type='l', col=6, lty=2, lwd=4)
seqpol <- seq(qt1[1],qt1[3],length.out=64)
polygon(x=c(qt1[1], seqpol, qt1[3])*2+20,
        y=c(0, pinter1(seqpol), 0), col=paste0(palette()[5],'33'), border=NA)
seqpol <- seq(qt2[1],qt2[3],length.out=64)
polygon(x=c(qt2[1], seqpol, qt2[3])*2+20,
        y=c(0, pinter2(seqpol), 0), col=paste0(palette()[2],'33'), border=NA)
grid(lwd=0.5,lty=1, col=paste0(palette()[7],'88'))
    legend('topright', legend=paste0('T = ',xtest[3]*2+20), cex=2, bty='n')
    legend('topleft', legend=c(paste0('75% of full population A has ',signif(qt1[1]*2+20,3),' < Y < ',signif(qt1[3]*2+20,3)),
                               paste0('75% of full population B has ',signif(qt2[1]*2+20,3),' < Y < ',signif(qt2[3]*2+20,3))
                               ), lty=c(1,2), lwd=3, col=c(1,2), bty='n', cex=1.5)
##
hist((medians1-medians2)*2, breaks=10, xlab='(median Y group A) - (median Y group B)', ylab='probability density', freq=F, col=3, cex.axis=1.5, cex.lab=1.5, border=NA, main='Forecast of difference between medians of full populations A and B', xlim=c(-4,4))#max(medians1-medians2)*2.1))
axis(1, at=seq(1,3,by=1), cex.axis=1.5)
grid(lwd=0.5,lty=1, col=paste0(palette()[7],'88'))
abline(v=0, lwd=2, col=paste0(palette()[7],'88'))
legend('topleft', legend=paste0(signif(quantile(medians1-medians2,c(1)/8)*2,2),' < difference < ',signif(quantile(medians1-medians2,c(7)/8)*2,2),'\nwith 75% probability','\n\n',
                signif(quantile(medians1-medians2,c(5)/100)*2,2),' < difference < ',signif(quantile(medians1-medians2,c(95)/100)*2,2),'\nwith 90% probability'                ), bty='n', cex=1.5)
dev.off()


## initsFunction <- function(){
## list(
##     alphaK=rep(4/nclusters, nclusters),
##     meanCmean=medianccovs,
##     meanCshape1=rep(1/2, nccovs),
##     meanCrate2=1/(widthccovs/2)^2, # dims = inv. variance
##     ##
##     tauCshape1=rep(1/2, nccovs),
##     tauCrate2=1/(widthccovs/2)^2, # dims = inv. variance
##     ##
##     probDa1=rep(1, ndcovs),
##     probDb1=rep(1, ndcovs),
##     sizeDsize1=maxdcovs,
##     sizeDa2=rep(32, ndcovs),
##     sizeDb2=rep(32, ndcovs),
##     ##
##     ##
##     meanCtau1=1/(widthccovs/2)^2, # dims = inv. variance
##     meanCrate1=(widthccovs/2)^2, # dims = variance
##     tauCrate1=(widthccovs/2)^2, # dims = variance
##     ##
##     sizeDprob1=rep(1/2, ndcovs),
##     ##
##     ##
##     q=rep(1/nclusters, nclusters),
##     ##
##     meanC=matrix(rnorm(n=nccovs*nclusters, mean=medianccovs, sd=1/sqrt(rgamma(n=nccovs,shape=1/2,rate=rgamma(n=nccovs,shape=1/2,rate=1/(widthccovs/2)^2)))), nrow=nccovs, ncol=nclusters),
##     tauC=matrix(rgamma(n=nccovs*nclusters, shape=1/2, rate=rgamma(n=nccovs,shape=1/2,rate=1/(widthccovs/2)^2)), nrow=nccovs, ncol=nclusters),
##     probD=matrix(rbeta(n=ndcovs*nclusters, shape1=1, shape2=1), nrow=ndcovs, ncol=nclusters),
##     sizeD=apply(matrix(rnbinom(n=ndcovs*nclusters, prob=rbeta(n=ndcovs,shape1=32,shape2=32), size=maxdcovs), nrow=ndcovs, ncol=nclusters), 2, function(x){maxdcovs*(x<maxdcovs)+x*(x>=maxdcovs)}),
##     ## C=rep(1,ndata)
##     C=rcat(n=ndata, prob=rep(1/nclusters,nclusters))
## )
## }
## temppdf <- function(x, acov){
##     inits <- initsFunction()
##     sum(inits$q * dnorm(x=x, mean=inits$meanC[acov,], sd=1/sqrt(inits$tauC[acov,])))
## }
## xgrid <- seq(-6,1,length.out=2^10)
## ygrid <- sapply(xgrid, function(x){temppdf(x,3)})
## matplot(xgrid,ygrid,type='l')
## pdff('hists')
## for(acov in covNames){
##     hist(alldata[[acov]], breaks=100, main=acov)
##     legend(x='top',legend=paste0('median: ', signif(c(medianccovs,mediandcovs)[acov],3), '  width: ', signif(c(widthccovs,widthdcovs)[acov],3)))
## }
## dev.off()

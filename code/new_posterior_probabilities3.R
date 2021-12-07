## Author: PGL  Porta Mana
## Created: 2021-11-25T14:52:14+0100
## Last-Updated: 2021-12-01T21:16:28+0100
################
## Prediction of population frequencies for Alzheimer study
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


seed <- 149
baseversion <- 'test_'
nclusters <- as.integer(64) #as.integer(2^6)
niter <- as.integer(256) #as.integer(2^11)
niter0 <- as.integer(256) #as.integer(2^10)
thin <- as.integer(1)
nstages <- as.integer(1)
ncheckpoints <- as.integer(4)
maincov <- 'Subgroup_num_'
family <- 'Palatino'
posterior <- TRUE
##

saveinfofile <- 'variates_info.csv'
variateinfo <- fread(saveinfofile, sep=',')
covNames <- variateinfo$variate
covTypes <- variateinfo$type
covMins <- variateinfo$min
covMaxs <- variateinfo$max
names(covTypes) <- names(covMins) <- names(covMaxs) <- covNames

datafile <- 'data_transformed_shuffled.csv'
alldata <- fread(datafile, sep=',')
alldata <- alldata[Usage_ == 'train']

#################################
## Setup for Monte Carlo sampling
#################################

realCovs <- covNames[sapply(covNames, function(acov){is.double(alldata[[acov]])})]
integerCovs <- covNames[sapply(covNames, function(acov){is.integer(alldata[[acov]]) && covMaxs[acov]>1})]
binaryCovs <- covNames[sapply(covNames, function(acov){is.integer(alldata[[acov]]) && covMaxs[acov]==1})]
covNames <- c(realCovs, integerCovs, binaryCovs)
nrcovs <- length(realCovs)
nicovs <- length(integerCovs)
nbcovs <- length(binaryCovs)
ncovs <- length(covNames)
ndata <- nrow(alldata)
##
dirname <- paste0(baseversion,'-V',length(covNames),'-D',ndata,'-K',nclusters,'-I',niter)
dir.create(dirname)
if(file.exists("/cluster/home/pglpm/R")){
initial.options <- commandArgs(trailingOnly = FALSE)
thisscriptname <- sub('--file=', "", initial.options[grep('--file=', initial.options)])
file.copy(from=thisscriptname, to=paste0(dirname,'/script-R',baseversion,'-V',length(covNames),'-D',ndata,'-K',nclusters,'.Rscript'))
}

for(obj in c('constants', 'dat', 'inits', 'bayesnet', 'model', 'Cmodel', 'confmodel', 'mcmcsampler', 'Cmcmcsampler')){if(exists(obj)){do.call(rm,list(obj))}}
gc()
##
##
dat <- list(
    X=as.matrix(alldata[1:ndata, ..realCovs]),
    Y=as.matrix(alldata[1:ndata, ..integerCovs]),
    Z=as.matrix(alldata[1:ndata, ..binaryCovs])
)
##
##
source('new_functions_mcmc.R')
alldataRanges <- dataQuantiles <- list()
for(acov in covNames){
        dataQuantiles[[acov]] <- quant(alldata[1:ndata][[acov]], prob=c(0.005,0.995))
        alldataRanges[[acov]] <- range(alldata[1:ndata][[acov]])
}
medianrcovs <- apply(alldata[1:ndata,..realCovs],2,median)
widthrcovs <- apply(alldata[1:ndata,..realCovs],2,IQR)*3/4
##
medianicovs <- apply(alldata[1:ndata,..integerCovs],2,median)
widthicovs <- ceiling(apply(alldata[1:ndata,..integerCovs],2,IQR))
maxicovs <- apply(alldata[1:ndata,..integerCovs],2,max)
thmaxicovs <- covMaxs[integerCovs]
matrixprobicovs <- matrix(0, nrow=nicovs, ncol=max(thmaxicovs), dimnames=list(integerCovs))
for(acov in integerCovs){
    matrixprobicovs[acov,1:thmaxicovs[acov]] <- (1:thmaxicovs[acov])/sum(1:thmaxicovs[acov])
}
##
medianbcovs <- apply(alldata[1:ndata,..binaryCovs],2,median)
##
print('Creating and saving checkpoints')
checkpointsFile <- paste0(dirname,'/_checkpoints-',ncheckpoints,'-R',baseversion,'-V',length(covNames),'-D',ndata,'-K',nclusters,'.rds')
checkpoints <- rbind(
    c(medianrcovs, medianicovs, medianbcovs),
    c(medianrcovs+widthrcovs, sapply(round(medianicovs+widthicovs), function(x){min(maxicovs,x)}), 1+0*medianbcovs),
    c(medianrcovs-widthrcovs, sapply(round(medianicovs-widthicovs), function(x){max(0,x)}), 0*medianbcovs),
    as.matrix(alldata[sample(1:ndata, size=ncheckpoints), ..covNames])
)
rownames(checkpoints) <- c('Pdatamedians', 'PdatacornerHi', 'PdatacornerLo', paste0('Pdatum',1:ncheckpoints))
saveRDS(checkpoints,file=checkpointsFile)
##
##
constants <- c(
    list(nClusters=nclusters),
    if(nrcovs>0){list(nRcovs=nrcovs)},
    if(nicovs>0){list(nIcovs=nicovs, maxIcovs=ncol(matrixprobicovs))},
    if(nbcovs>0){list(nBcovs=nbcovs)},
    if(posterior){list(nData=ndata)}
)

##
initsFunction <- function(){c(list(
       qalpha=rep(1/nclusters, nclusters),
       q=rep(1/nclusters, nclusters)
       ),
       if(nrcovs>0){list(
                        meanRmean=medianrcovs,
                        meanRshape1=rep(1/2, nrcovs),
                        meanRrate2=1/(widthrcovs/2)^2, # dims = inv. variance
                        ##
                        tauRshape1=rep(1, nrcovs),
                        tauRrate2=1/(widthrcovs/2)^2, # dims = inv. variance
                        ##
                        meanRtau1=1/(widthrcovs/2)^2, # dims = inv. variance
                        meanRrate1=(widthrcovs/2)^2, # dims = variance
                        tauRrate1=(widthrcovs/2)^2, # dims = variance
                        ##
                        meanR=matrix(rnorm(n=nrcovs*(nclusters), mean=medianrcovs, sd=1/sqrt(rgamma(n=nrcovs, shape=1/2, rate=rgamma(n=nrcovs, shape=1/2, rate=1/(widthrcovs/2)^2)))), nrow=nrcovs, ncol=nclusters),
                        tauR=matrix(rgamma(n=nrcovs*(nclusters), shape=1/2, rate=rgamma(n=nrcovs, shape=1/2, rate=1/(widthrcovs/2)^2)), nrow=nrcovs, ncol=nclusters)
                        )},
       if(nicovs>0){list(
                        probIa1=rep(2, nicovs),
                        probIb1=rep(1, nicovs),
                        sizeIprob1=matrixprobicovs,
                        ##
                        probI=matrix(rbeta(n=nicovs*(nclusters), shape1=1, shape2=1), nrow=nicovs, ncol=nclusters),
                        sizeI=matrix(maxicovs, nrow=nicovs, ncol=nclusters)
                        )},
       if(nbcovs>0){list(
                        probBa1=rep(1,nbcovs),
                        probBa2=rep(1,nbcovs),
                        ##
                        probB=matrix(0.5, nrow=nbcovs, ncol=nclusters)
                        )},
       if(posterior){list(C=rep(1,ndata))}
)}

##
##
bayesnet <- nimbleCode({
    q[1:nClusters] ~ ddirch(alpha=qalpha[1:nClusters])
    for(acluster in 1:nClusters){
        if(nrcovs>0){
            for(acov in 1:nRcovs){
                meanR[acov,acluster] ~ dnorm(mean=meanRmean[acov], tau=meanRtau1[acov])
                tauR[acov,acluster] ~ dgamma(shape=tauRshape1[acov], rate=tauRrate1[acov])
            }
        }
        if(nicovs>0){
            for(acov in 1:nIcovs){
                probI[acov,acluster] ~ dbeta(shape1=probIa1[acov], shape2=probIb1[acov])
                sizeI[acov,acluster] ~ dcat(prob=sizeIprob1[acov,1:maxIcovs])
            }
        }
        if(nbcovs>0){
            for(acov in 1:nBcovs){
                probB[acov,acluster] ~ dbeta(shape1=probBa1[acov], shape2=probBa2[acov])
            }
        }
    }
    ##
    if(nrcovs>0){
        for(acov in 1:nRcovs){
            meanRtau1[acov] ~ dgamma(shape=meanRshape1[acov], rate=meanRrate1[acov])
            meanRrate1[acov] ~ dgamma(shape=meanRshape1[acov], rate=meanRrate2[acov])
            tauRrate1[acov] ~ dgamma(shape=tauRshape1[acov], rate=tauRrate2[acov])
        }
    }
    ##
    if(posterior){
        for(adatum in 1:nData){
            C[adatum] ~ dcat(prob=q[1:nClusters])
        }            ##
        for(adatum in 1:nData){
            if(nrcovs>0){
                for(acov in 1:nRcovs){
                    X[adatum,acov] ~ dnorm(mean=meanR[acov,C[adatum]], tau=tauR[acov,C[adatum]])
                }
            }
            if(nicovs>0){
                for(acov in 1:nIcovs){
                    Y[adatum,acov] ~ dbinom(prob=probI[acov,C[adatum]], size=sizeI[acov,C[adatum]])
                }
            }
            if(nbcovs>0){
                for(acov in 1:nBcovs){
                    Z[adatum,acov] ~ dbern(prob=probB[acov,C[adatum]])
                }
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
                               monitors=c('q',
                                          if(nrcovs>0){c('meanR', 'tauR')},
                                          if(nicovs>0){c('probI', 'sizeI')},
                                          if(nbcovs>0){c('probB')}
                                          ),
                               monitors2=c('C',
                                           if(nrcovs>0){c('meanRtau1', 'meanRrate1', 'tauRrate1')}
                                           ))
    ##
    for(adatum in 1:ndata){
        confmodel$addSampler(target=paste0('C[', adatum, ']'), type='categorical')
    }
    confmodel$addSampler(target=paste0('q[1:', nclusters, ']'), type='conjugate')
    for(acluster in 1:nclusters){
        if(nrcovs>0){
            for(acov in 1:nrcovs){
                confmodel$addSampler(target=paste0('meanR[', acov, ', ', acluster, ']'), type='conjugate')
                confmodel$addSampler(target=paste0('tauR[', acov, ', ', acluster, ']'), type='conjugate')
            }
        }
        if(nicovs>0){
            for(acov in 1:nicovs){
                confmodel$addSampler(target=paste0('probI[', acov, ', ', acluster, ']'), type='conjugate')
                confmodel$addSampler(target=paste0('sizeI[', acov, ', ', acluster, ']'), type='categorical')
            }
        }
        if(nbcovs>0){
            for(acov in 1:nbcovs){
                confmodel$addSampler(target=paste0('probB[', acov, ', ', acluster, ']'), type='conjugate')
            }
        }
    }
    if(nrcovs>0){
        for(acov in 1:nrcovs){
            confmodel$addSampler(target=paste0('meanRtau1[', acov, ']'), type='conjugate')
            confmodel$addSampler(target=paste0('meanRrate1[', acov, ']'), type='conjugate')
            confmodel$addSampler(target=paste0('tauRrate1[', acov, ']'), type='conjugate')
        }
    }
}else{
    confmodel <- configureMCMC(Cmodel, 
                               monitors=c('q',
                                          if(nrcovs>0){c('meanR', 'tauR')},
                                          if(nicovs>0){c('probI', 'sizeI')},
                                          if(nbcovs>0){c('probB')}
                                          ))
}
print(confmodel)
## confmodel$printSamplers(executionOrder=TRUE)

mcmcsampler <- buildMCMC(confmodel)
Cmcmcsampler <- compileNimble(mcmcsampler, resetFunctions = TRUE)
gc()

print('Setup time:')
print(Sys.time() - timecount)

##################################################
## Monte Carlo sampler and plots of MC diagnostics
##################################################
## if(!posterior){
##     nstages <- 1
##     niter0 <- niter <- 64
## }
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
    saveRDS(mcsamples,file=paste0(dirname,'/_mcsamples-R',version,'-V',length(covNames),'-D',ndata,'-K',nclusters,'-I',nrow(mcsamples),'.rds'))
    ## save final state of MCMC chain
    finalstate <- as.matrix(Cmcmcsampler$mvSamples2)
    finalstate <- c(mcsamples[nrow(mcsamples),], finalstate[nrow(finalstate),])
    occupations <- finalstate[grepl('^C\\[', names(finalstate))]
    usedclusters <- length(unique(occupations))
    if(usedclusters > nclusters-5){print('WARNING: TOO MANY CLUSTERS OCCUPIED')}
    print(paste0('OCCUPIED CLUSTERS: ', usedclusters, ' OF ', nclusters))
    saveRDS(finalstate2list(finalstate),file=paste0(dirname,'/_finalstate-R',version,'-V',length(covNames),'-D',ndata,'-K',nclusters,'-I',nrow(mcsamples),'.rds'))
    ##
    parmList <- mcsamples2parmlist(mcsamples)
    saveRDS(parmList,file=paste0(dirname,'/_frequencies-R',version,'-V',length(covNames),'-D',ndata,'-K',nclusters,'-I',nrow(mcsamples),'.rds'))
    ## Traces to follow for diagnostics
    ll <- llSamples(dat, parmList)
    flagll <- FALSE
    if(!posterior && !any(is.finite(ll))){
        flagll <- TRUE
        ll <- rep(0, length(ll))}
    momentstraces <- moments12Samples(parmList)
    probCheckpoints <- samplesF(Y=checkpoints, parmList=parmList, inorder=TRUE)
    ## miqrtraces <- calcSampleMQ(parmList)
    ## medians <- miqrtraces[,,1]
    ## colnames(medians) <- paste0('MEDIAN_', colnames(miqrtraces))
    ## Q1s <- miqrtraces[,,2]
    ## colnames(Q1s) <- paste0('Q1_', colnames(miqrtraces))
    ## Q3s <- miqrtraces[,,3]
    ## colnames(Q3s) <- paste0('Q3_', colnames(miqrtraces))
    ## iqrs <- Q3s - Q1s
    ## colnames(iqrs) <- paste0('IQR_', colnames(miqrtraces))
    traces <- cbind(LL=ll, t(probCheckpoints), #medians, iqrs, Q1s, Q3s,
                    do.call(cbind, momentstraces))
    badcols <- foreach(i=1:ncol(traces), .combine=c)%do%{if(all(is.na(traces[,i]))){i}else{NULL}}
    if(!is.null(badcols)){traces <- traces[,-badcols]}
    saveRDS(traces,file=paste0(dirname,'/_traces-R',version,'-V',length(covNames),'-D',ndata,'-K',nclusters,'-I',nrow(mcsamples),'.rds'))
    
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
    ## samplesQuantiles <- calcSampleQuantiles(parmList)
    ##
    ## xlimits <- list()
    ## for(acov in covNames){
    ##     xlimits[[acov]] <- range(c(alldataRanges[[acov]], samplesQuantiles[,acov,]))
    ## }
    ##

    ##
    pdff(paste0(dirname,'/mcsummary-R',version,'-V',length(covNames),'-D',ndata,'-K',nclusters,'-I',nrow(mcsamples)))
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
    legend(x='center', bty='n', cex=1,
           legend=c('WARNINGS:',
                    if(any(is.na(mcsamples))){'some NA MC outputs'},
                    if(any(!is.finite(mcsamples))){'some infinite MC outputs'},
                    if(usedclusters > nclusters-5){'too many clusters occupied'},
                    if(flagll){'infinite values in likelihood'}
           ))
    ##
    par(mfrow=c(1,1))
    for(acov in covNames){
        datum <- alldata[1:ndata][[acov]]
        if(acov %in% realCovs){
            rg <- range(datum)+c(-1,1)*IQR(datum)
            Xgrid <- seq(rg[1], rg[2], length.out=256)
            tpar <- unlist(variateinfo[variate==acov,c('transfM','transfW')])
            if(!any(is.na(tpar))){
                Ogrid <- pretty(exp(tpar['transfW']*Xgrid + tpar['transfM']),n=10)
            }
        }else{
            rg <- range(datum)
            rg <- round(c((covMins[acov]+7*rg[1])/8, (covMaxs[acov]+7*rg[2])/8))
            Xgrid <- rg[1]:rg[2]
            tpar <- NA
        }
        Xgrid <- cbind(Xgrid)
        colnames(Xgrid) <- acov
        plotsamples <- samplesF(Y=Xgrid, parmList=parmList, nfsamples=min(64,nrow(mcsamples)), inorder=FALSE)
        ymax <- quant(apply(plotsamples,2,function(x){quant(x,99/100)}),99/100)
        ## ymax <- quant(apply(plotsamples,2,max),99/100)
        tplot(x=Xgrid, y=plotsamples, type='l', col=paste0(palette()[7], '44'), lty=1, lwd=2, xlab=acov, ylab='probability density', ylim=c(0, ymax), family=family)#max(plotsamples[plotsamples<df])))
        if(!any(is.na(tpar))){
            axis(3,at=(log(Ogrid)-tpar['transfM'])/tpar['transfW'],labels=Ogrid,lwd=0,lwd.ticks=1,col.ticks='#bbbbbb80')
        }
        if(acov %in% binaryCovs){
            histo <- thist(plotsamples[2,])
            tplot(histo$breaks, histo$density, col=7, xlab=paste0('P(',acov,' = 1)'), ylab='probability density', ylim=c(0, max(histo$density)), xlim=c(0,1), family=family)
        }
    }
    ## ##
    ## par(mfrow = rep(ceiling(sqrt(nicovs+nrcovs)), 2))
    ## for(addvar in setdiff(covNames, maincov)){
    ##     tplot(x=c(rep(alldataRanges[[maincov]], each=2),
    ##                 alldataRanges[[maincov]][1]),
    ##             y=c(alldataRanges[[addvar]], rev(alldataRanges[[addvar]]),
    ##                 alldataRanges[[addvar]][1]),
    ##             type='l', lwd=2, col=paste0(palette()[2], '88'),
    ##             xlim=xlimits[[maincov]],
    ##             ylim=xlimits[[addvar]],
    ##             xlab=maincov,
    ##             ylab=addvar
    ##             )
    ##     matlines(x=c(rep(dataQuantiles[[maincov]], each=2),
    ##                  dataQuantiles[[maincov]][1]),
    ##              y=c(dataQuantiles[[addvar]], rev(dataQuantiles[[addvar]]),
    ##                  dataQuantiles[[addvar]][1]),
    ##              lwd=2, col=paste0(palette()[4], '88'))
    ## }
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
                ylab=acov, family=family
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


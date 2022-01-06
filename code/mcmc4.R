## Author: PGL  Porta Mana
## Created: 2021-11-25T14:52:14+0100
## Last-Updated: 2022-01-06T09:24:14+0100
################
## Template code for model-free probabilistic analysis and prediction of data
## Works with continuous and discrete (categorical, binary, integer) variables
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


seed <- 707
baseversion <- 'testmcmc4_'
nclusters <- 64L
niter <- 128L
niter0 <- 128L
thin <- 1L
nstages <- 0L
ncheckprobs1 <- 16L
ncheckprobs2 <- 8L
maincov <- 'Vb'
family <- 'Palatino'
##ndata <- 128L
posterior <- FALSE
##

saveinfofile <- 'testdata_variate_info.csv'
variateinfo <- fread(saveinfofile, sep=',')
covNames <- variateinfo$variate
covTypes <- variateinfo$type
covMins <- variateinfo$min
covMaxs <- variateinfo$max
names(covTypes) <- names(covMins) <- names(covMaxs) <- covNames

datafile <- 'testdata.csv'
odata <- fread(datafile, sep=',')

#################################
## Setup for Monte Carlo sampling
#################################

realCovs <- covNames[covTypes=='double']
integerCovs <- covNames[covTypes=='integer']
binaryCovs <- covNames[covTypes=='binary']
covNames <- c(realCovs, integerCovs, binaryCovs)
nrcovs <- length(realCovs)
nicovs <- length(integerCovs)
nbcovs <- length(binaryCovs)
ncovs <- length(covNames)
if(!exists('ndata') || is.null(ndata) || is.na(ndata)){ndata <- nrow(odata)}
alldata <- odata[1:ndata, ..covNames]

##
for(obj in c('constants', 'dat', 'inits', 'bayesnet', 'model', 'Cmodel', 'confmodel', 'mcmcsampler', 'Cmcmcsampler')){if(exists(obj)){do.call(rm,list(obj))}}
gc()

dat <- list()
if(nrcovs>0){ dat$X=data.matrix(alldata[, ..realCovs])}
if(nicovs>0){ dat$Y=data.matrix(alldata[, ..integerCovs])}
if(nbcovs>0){ dat$Z=data.matrix(alldata[, ..binaryCovs])}
##
##
dirname <- paste0(baseversion,'-V',length(covNames),'-D',ndata,'-K',nclusters,'-I',niter)
dir.create(dirname)
if(file.exists("/cluster/home/pglpm/R")){
initial.options <- commandArgs(trailingOnly = FALSE)
thisscriptname <- sub('--file=', "", initial.options[grep('--file=', initial.options)])
file.copy(from=thisscriptname, to=paste0(dirname,'/script-R',baseversion,'-V',length(covNames),'-D',ndata,'-K',nclusters,'.Rscript'))
}

##
##
source('functions_mcmc.R')
## alldataRanges <- dataQuantiles <- list()
## for(avar in covNames){
##         dataQuantiles[[avar]] <- quant(alldata[1:ndata][[avar]], prob=c(0.005,0.995))
##         alldataRanges[[avar]] <- range(alldata[1:ndata][[avar]])
## }
##
if(nrcovs>0){
    medianrcovs <- apply(alldata[1:ndata,..realCovs],2,function(x)median(x, na.rm=TRUE))
    widthrcovs <- apply(alldata[1:ndata,..realCovs],2,function(x)IQR(x, na.rm=TRUE, type=8))
}
##
if(length(integerCovs)>0){
    medianicovs <- round(apply(alldata[1:ndata,..integerCovs],2,function(x)median(x, na.rm=TRUE)))
    widthicovs <- ceiling(apply(alldata[1:ndata,..integerCovs],2,function(x)IQR(x, na.rm=TRUE, type=8)))
    maxicovs <- apply(alldata[1:ndata,..integerCovs],2,function(x)max(x, na.rm=T))
    thmaxicovs <- covMaxs[integerCovs]
    matrixprobicovs <- matrix(0, nrow=nicovs, ncol=max(thmaxicovs), dimnames=list(integerCovs))
    for(avar in integerCovs){
        matrixprobicovs[avar,1:thmaxicovs[avar]] <- (1:thmaxicovs[avar])/sum(1:thmaxicovs[avar])
    }
}
##
if(nbcovs>0){
    medianbcovs <- apply(alldata[1:ndata,..binaryCovs],2,function(x)median(x, na.rm=TRUE))
}

##
##
if(posterior){
print('Creating and saving checkpoints')
checkprobsFile <- paste0(dirname,'/_checkprobs-R',baseversion,'-V',length(covNames),'-D',ndata,'-K',nclusters,'.rds')
if(exists('continue') && is.character(continue)){
    checkprobs <- readRDS(checkprobsFile)
}else{
    valuelist <- data.matrix(na.omit(alldata))
    nvalues <- nrow(valuelist)
    ncheckprobs1 <- min(ncheckprobs1, nvalues)
    ncheckprobs2 <- min(ncheckprobs2, nvalues)
    psamples0 <- sample(1:nvalues, ncheckprobs1)
    psamples1 <- sample(1:nvalues, ncheckprobs2)
    csamples1 <- sample(setdiff(covNames,maincov), ncheckprobs2, replace=(ncheckprobs2>ncovs-1))
    psamples2 <- sample(1:nvalues, ncheckprobs2)
    csamples2 <- sample(setdiff(covNames,maincov), ncheckprobs2, replace=(ncheckprobs2>ncovs-1))
    checkprobs <- c( list(
        list(y=valuelist[psamples0,,drop=F], x=NULL,
             names=c(
    paste0('joint_',psamples0),
    paste0(maincov,'|rest_',psamples0),
    paste0('rest|',maincov,'_',psamples0),
    paste0(csamples1,'|rest_',psamples1),
    paste0('rest|',csamples2,'_',psamples2) )
    ),
        list(y=valuelist[psamples0, maincov, drop=F], x=valuelist[psamples0, setdiff(covNames,maincov), drop=F]),
        list(y=valuelist[psamples0, setdiff(covNames,maincov), drop=F], x=valuelist[psamples0, maincov, drop=F])
    ),
    lapply(1:ncheckprobs2,
           function(apoint){list(y=valuelist[psamples1[apoint], csamples1[apoint], drop=F],
                                 x=valuelist[psamples1[apoint], setdiff(covNames,csamples1[apoint]), drop=F] )}),
    lapply(1:ncheckprobs2,
           function(apoint){list(y=valuelist[psamples2[apoint], setdiff(covNames,csamples2[apoint]), drop=F],
                                 x=valuelist[psamples2[apoint], csamples2[apoint], drop=F] )})
    )
    saveRDS(checkprobs,file=checkprobsFile)
}
}
##
##
constants <- list(nClusters=nclusters)
if(nrcovs>0){constants$nRcovs <- nrcovs}
if(nicovs>0){constants$nIcovs <- nicovs
    constants$maxIcovs <- ncol(matrixprobicovs)}
if(nbcovs>0){constants$nBcovs <- nbcovs}
if(posterior){constants$nData <- ndata}

##
initsFunction <- function(){
    c(list(# cluster probabilities
        qalpha0=rep(1/nclusters, nclusters) # hyperparameters
#        q=rep(1/nclusters, nclusters) # variables
    ),
    if(nrcovs>0){# real variates
        list(# hyperparameters
            meanRmean0=medianrcovs,
            meanRshape0=rep(1/2, nrcovs),
            meanRrate0=(widthrcovs^2)/2, # dims = variance
            tauRshape0=rep(1/2, nrcovs),
            tauRrate0=1/(widthrcovs/(2*qnorm(3/4)))^2 # dims = inv. variance
            ## integrated parameters
#            meanRtau=1/(widthrcovs/(2*qnorm(3/4)))^2, # dims = inv. variance
#            tauRrate=(widthrcovs/(2*qnorm(3/4)))^2, # dims = variance
            ## variables
#            meanR=matrix(rnorm(n=nrcovs*(nclusters), mean=medianrcovs, sd=1/sqrt(rgamma(n=nrcovs, shape=1/2, rate=(widthrcovs^2)/2))), nrow=nrcovs, ncol=nclusters),
#            tauR=matrix(rgamma(n=nrcovs*(nclusters), shape=1/2, rate=rgamma(n=nrcovs, shape=1/2, rate=1/(widthrcovs/(2*qnorm(3/4)))^2)), nrow=nrcovs, ncol=nclusters)
        )},
    if(nicovs>0){# integer variates
        list(# hyperparameters
            probIa0=rep(2, nicovs),
            probIb0=rep(1, nicovs),
            sizeIprob0=matrixprobicovs,
            ## variables
#            probI=matrix(rbeta(n=nicovs*(nclusters), shape1=2, shape2=1), nrow=nicovs, ncol=nclusters),
            sizeI=matrix(maxicovs, nrow=nicovs, ncol=nclusters)
        )},
    if(nbcovs>0){# binary variates
        list(# hyperparameters
            probBa0=rep(1,nbcovs),
            probBb0=rep(1,nbcovs)
            ## variables
#            probB=matrix(0.5, nrow=nbcovs, ncol=nclusters)
        )},
    if(posterior){list(C=rep(1, ndata))} # cluster occupations
)}

##
##
bayesnet <- nimbleCode({
    q[1:nClusters] ~ ddirch(alpha=qalpha0[1:nClusters])
    for(acluster in 1:nClusters){
        if(nrcovs>0){# real variates
            for(avar in 1:nRcovs){
                meanR[avar,acluster] ~ dnorm(mean=meanRmean0[avar], tau=meanRtau[avar])
                tauR[avar,acluster] ~ dgamma(shape=tauRshape0[avar], rate=tauRrate[avar])
            }
        }
        if(nicovs>0){# integer variates
            for(avar in 1:nIcovs){
                probI[avar,acluster] ~ dbeta(shape1=probIa0[avar], shape2=probIb0[avar])
                sizeI[avar,acluster] ~ dcat(prob=sizeIprob0[avar,1:maxIcovs])
            }
        }
        if(nbcovs>0){# binary variates
            for(avar in 1:nBcovs){
                probB[avar,acluster] ~ dbeta(shape1=probBa0[avar], shape2=probBb0[avar])
            }
        }
    }
    ##
    if(nrcovs>0){# integrated real variables
        for(avar in 1:nRcovs){
            meanRtau[avar] ~ dgamma(shape=meanRshape0[avar], rate=meanRrate0[avar])
            tauRrate[avar] ~ dgamma(shape=tauRshape0[avar], rate=tauRrate0[avar])
        }
    }
    ##
    if(posterior){# cluster occupations
        for(adatum in 1:nData){
            C[adatum] ~ dcat(prob=q[1:nClusters])
        }            ##
        for(adatum in 1:nData){
            if(nrcovs>0){
                for(avar in 1:nRcovs){
                    X[adatum,avar] ~ dnorm(mean=meanR[avar,C[adatum]], tau=tauR[avar,C[adatum]])
                }
            }
            if(nicovs>0){
                for(avar in 1:nIcovs){
                    Y[adatum,avar] ~ dbinom(prob=probI[avar,C[adatum]], size=sizeI[avar,C[adatum]])
                }
            }
            if(nbcovs>0){
                for(avar in 1:nBcovs){
                    Z[adatum,avar] ~ dbern(prob=probB[avar,C[adatum]])
                }
            }
        }
    }
})

##
##
timecount <- Sys.time()

if(posterior){
    model <- nimbleModel(code=bayesnet, name='model1', constants=constants, inits=initsFunction(), data=dat, dimensions=list(q=nclusters, meanRtau=nrcovs, tauRrate=nrcovs, meanR=c(nrcovs,nclusters), tauR=c(nrcovs,nclusters), probI=c(nicovs,nclusters), probB=c(nbcovs,nclusters), C=ndata) )
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
                                           if(nrcovs>0){c('meanRtau', 'tauRrate')}
                                           ))
    ##
    for(adatum in 1:ndata){
        confmodel$addSampler(target=paste0('C[', adatum, ']'), type='categorical')
    }
    for(acluster in 1:nclusters){
        if(nrcovs>0){
            for(avar in 1:nrcovs){
                confmodel$addSampler(target=paste0('meanR[', avar, ', ', acluster, ']'), type='conjugate')
                confmodel$addSampler(target=paste0('tauR[', avar, ', ', acluster, ']'), type='conjugate')
            }
        }
        if(nicovs>0){
            for(avar in 1:nicovs){
                confmodel$addSampler(target=paste0('probI[', avar, ', ', acluster, ']'), type='conjugate')
                confmodel$addSampler(target=paste0('sizeI[', avar, ', ', acluster, ']'), type='categorical')
            }
        }
        if(nbcovs>0){
            for(avar in 1:nbcovs){
                confmodel$addSampler(target=paste0('probB[', avar, ', ', acluster, ']'), type='conjugate')
            }
        }
    }
    if(nrcovs>0){
        for(avar in 1:nrcovs){
            confmodel$addSampler(target=paste0('meanRtau[', avar, ']'), type='conjugate')
            confmodel$addSampler(target=paste0('tauRrate[', avar, ']'), type='conjugate')
        }
    }
    confmodel$addSampler(target=paste0('q[1:', nclusters, ']'), type='conjugate')
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
for(stage in 0:nstages){
    totalruntime <- Sys.time()

    print(paste0('==== STAGE ', stage, ' ===='))
    version <- paste0(baseversion, stage)
    gc()
    if(stage==0){
        if(exists('continue') && is.character(continue)){
            initsc <- readRDS(continue)
            inits0 <- initsFunction()
            for(aname in names(initsc)){inits0[[aname]] <- initsc[[aname]]}
        }else{inits0 <- initsFunction}
        mcsamples <- runMCMC(Cmcmcsampler, nburnin=1, niter=niter0+1, thin=1, thin2=niter0, inits=inits0, setSeed=seed)
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
    saveRDS(parmList,file=paste0(dirname,'/_frequencies-R',version,'-V',length(covNames),'-D',ndata,'-K',nclusters,'-I',nrow(parmList$q),'.rds'))
    ## Traces to follow for diagnostics
    ll <- llSamples(dat, parmList)
    flagll <- FALSE
    if(!posterior && !any(is.finite(ll))){
        flagll <- TRUE
        ll <- rep(0, length(ll))}
    ##momentstraces <- moments12Samples(parmList)
if(posterior){
    probCheckprobs <- foreach(apoint=checkprobs, .combine=rbind)%do%{
        samplesF(Y=apoint$y, X=apoint$x, parmList=parmList, inorder=TRUE)
    }
    rownames(probCheckprobs) <- checkprobs[[1]]$names
    }
    ## miqrtraces <- calcSampleMQ(parmList)
    ## medians <- miqrtraces[,,1]
    ## colnames(medians) <- paste0('MEDIAN_', colnames(miqrtraces))
    ## Q1s <- miqrtraces[,,2]
    ## colnames(Q1s) <- paste0('Q1_', colnames(miqrtraces))
    ## Q3s <- miqrtraces[,,3]
    ## colnames(Q3s) <- paste0('Q3_', colnames(miqrtraces))
    ## iqrs <- Q3s - Q1s
    ## colnames(iqrs) <- paste0('IQR_', colnames(miqrtraces))
    traces <- cbind(loglikelihood=ll, t(log(probCheckprobs))) #medians, iqrs, Q1s, Q3s,
                    ##do.call(cbind, momentstraces))
    badcols <- foreach(i=1:ncol(traces), .combine=c)%do%{if(all(is.na(traces[,i]))){i}else{NULL}}
    if(!is.null(badcols)){traces <- traces[,-badcols]}
    saveRDS(traces,file=paste0(dirname,'/_probtraces-R',version,'-V',length(covNames),'-D',ndata,'-K',nclusters,'-I',nrow(parmList$q),'.rds'))
    
    ##
    if(nrow(traces)>=1000){
        funMCSE <- function(x){LaplacesDemon::MCSE(x, method='batch.means')$se}
    }else{
        funMCSE <- function(x){LaplacesDemon::MCSE(x)}
    }
    diagnESS <- LaplacesDemon::ESS(traces * (abs(traces) < Inf))
    diagnIAT <- apply(traces, 2, function(x){LaplacesDemon::IAT(x[is.finite(x)])})
    diagnBMK <- LaplacesDemon::BMK.Diagnostic(traces, batches=2)[,1]
    diagnMCSE <- 100*apply(traces, 2, function(x){funMCSE(x)/sd(x)})
    diagnStat <- apply(traces, 2, function(x){LaplacesDemon::is.stationary(as.matrix(x,ncol=1))})
    diagnBurn <- apply(traces, 2, function(x){LaplacesDemon::burnin(matrix(x[1:(10*trunc(length(x)/10))], ncol=1))})
    ##
    ##
    ## tracenames <- colnames(traces)
    tracegroups <- list(loglikelihood=1,
                        'joint probs'=c(1:ncheckprobs1)+1,
                        'single given rest'=c(ncheckprobs1+(1:ncheckprobs1),
                                                ncheckprobs1*3+(1:ncheckprobs2))+1,
                        'rest given single'=c(ncheckprobs1*2+(1:ncheckprobs1),
                                              ncheckprobs1*3+ncheckprobs2+(1:ncheckprobs2))+1
                        )
    ##     'maxD'=tracenames[grepl('^(Pdat|Dcov|LL)', tracenames)],
    ##     '1D'=tracenames[grepl('^(MEDIAN|Q1|Q3|IQR|MEAN|VAR)_', tracenames)],
    ##     '2D'=tracenames[grepl('^COV_', tracenames)]
    ## )
    grouplegends <- foreach(agroup=1:length(tracegroups))%do%{
        c( paste0('-- STATS ', names(tracegroups)[agroup], ' --'),
          paste0('min ESS = ', signif(min(diagnESS[tracegroups[[agroup]]]),6)),
          paste0('max IAT = ', signif(max(diagnIAT[tracegroups[[agroup]]]),6)),
          paste0('max BMK = ', signif(max(diagnBMK[tracegroups[[agroup]]]),6)),
          paste0('max MCSE = ', signif(max(diagnMCSE[tracegroups[[agroup]]]),6)),
          paste0('stationary: ', sum(diagnStat[tracegroups[[agroup]]]),'/',length(diagnStat[tracegroups[[agroup]]])),
          paste0('burn: ', signif(max(diagnBurn[tracegroups[[agroup]]]),6))
          )
    }
    colpalette <- c(7,
                    rep(2,ncheckprobs1), rep(1,ncheckprobs1), rep(3,ncheckprobs1),
                    rep(1,ncheckprobs2), rep(3,ncheckprobs2))
    names(colpalette) <- colnames(traces)
    ##
    ## samplesQuantiles <- calcSampleQuantiles(parmList)
    ##
    ## xlimits <- list()
    ## for(avar in covNames){
    ##     xlimits[[avar]] <- range(c(alldataRanges[[avar]], samplesQuantiles[,avar,]))
    ## }
    ##

    ##
    pdff(paste0(dirname,'/mcsummary2-R',version,'-V',length(covNames),'-D',ndata,'-K',nclusters,'-I',nrow(parmList$q)))
    matplot(1:2, type='l', col='white', main=paste0('Stats stage ',stage), axes=FALSE, ann=FALSE)
    legendpositions <- c('topleft','topright','bottomleft','bottomright')
    for(alegend in 1:length(grouplegends)){
        legend(x=legendpositions[alegend], bty='n', cex=1.5,
               legend=grouplegends[[alegend]] )
    }
    legend(x='center', bty='n', cex=1,
           legend=c(
               paste0('Occupied clusters: ', usedclusters, ' of ', nclusters),
               paste0('LL: ', signif(mean(ll),3), ' +- ', signif(sd(ll),3)),
               'WARNINGS:',
                    if(any(is.na(mcsamples))){'some NA MC outputs'},
                    if(any(!is.finite(mcsamples))){'some infinite MC outputs'},
                    if(usedclusters > nclusters-5){'too many clusters occupied'},
                    if(flagll){'infinite values in likelihood'}
           ))
    ##
    par(mfrow=c(1,1))
    for(avar in covNames){
        datum <- alldata[1:ndata][[avar]]
        if(avar %in% realCovs){
            rg <- range(datum, na.rm=T)+c(-1,1)*IQR(datum, type=8, na.rm=T)
            Xgrid <- seq(rg[1], rg[2], length.out=256)
            tpar <- unlist(variateinfo[variate==avar,c('transfM','transfW')])
            if(!any(is.na(tpar))){
                Ogrid <- pretty(exp(tpar['transfW']*Xgrid + tpar['transfM']),n=10)
            }
        }else{
            rg <- range(datum, na.rm=T)
            rg <- round(c((covMins[avar]+7*rg[1])/8, (covMaxs[avar]+7*rg[2])/8))
            Xgrid <- rg[1]:rg[2]
            tpar <- NA
        }
        Xgrid <- cbind(Xgrid)
        colnames(Xgrid) <- avar
        plotsamples <- samplesF(Y=Xgrid, parmList=parmList, nfsamples=min(64,nrow(mcsamples)), inorder=FALSE)
        ymax <- quant(apply(plotsamples,2,function(x){quant(x,99/100)}),99/100, na.rm=T)
        ## ymax <- quant(apply(plotsamples,2,max),99/100)
        tplot(x=Xgrid, y=plotsamples, type='l', col=paste0(palette()[7], '44'), lty=1, lwd=2, xlab=avar, ylab='probability density', ylim=c(0, ymax), family=family)#max(plotsamples[plotsamples<df])))
        if(!any(is.na(tpar))){
            axis(3,at=(log(Ogrid)-tpar['transfM'])/tpar['transfW'],labels=Ogrid,lwd=0,lwd.ticks=1,col.ticks='#bbbbbb80')
        }
        if(avar %in% binaryCovs){
            histo <- thist(plotsamples[2,])
            tplot(histo$breaks, histo$density, col=7, xlab=paste0('P(',avar,' = 1)'), ylab='probability density', ylim=c(0, max(histo$density)), xlim=c(0,1), family=family)
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
        transf <- identity
    for(avar in colnames(traces)){
        ## if(grepl('^[PDV]', avar)){transf <- function(x){log(abs(x)+1e-12)}
        ## if(grepl('^[PDV]', avar)){transf <- function(x){log(abs(x)+1e-12)}
        ## }else{transf <- identity}
        tplot(y=transf(traces[,avar]), type='l', lty=1, col=colpalette[avar],
                main=paste0(avar,
                            '\nESS = ', signif(diagnESS[avar], 3),
                            ' | IAT = ', signif(diagnIAT[avar], 3),
                            ' | BMK = ', signif(diagnBMK[avar], 3),
                            ' | MCSE(6.27) = ', signif(diagnMCSE[avar], 3),
                            ' | stat: ', diagnStat[avar],
                            ' | burn: ', diagnBurn[avar]
                            ),
                ylab=avar, family=family
              #, ylim=range(c(transf(traces[,avar][abs(transf(traces[,avar]))<Inf])))
              )
    }
    dev.off()

    print('Total runtime:')
    print(Sys.time() - totalruntime)

}
############################################################
## End MCMC
############################################################

xrg <- 10
xgrid <- seq(-xrg,xrg,length.out=512)
plotgrid <- c(8,8)
vrsamples <- samplesF(Y=cbind(Vr=xgrid), parmList=parmList, nfsamples=prod(plotgrid), inorder=T)
vi1samples <- samplesF(Y=cbind(Vi1=0:thmaxicovs['Vi1']), parmList=parmList, nfsamples=prod(plotgrid), inorder=T)
vi2samples <- samplesF(Y=cbind(Vi2=0:thmaxicovs['Vi2']), parmList=parmList, nfsamples=prod(plotgrid), inorder=T)
pdff('plots_vr')
par(mfrow=plotgrid,mar=c(0,0,0,0))
for(i in 1:prod(plotgrid)){
    tplot(x=xgrid, y=vrsamples[,i], lty=1, lwd=2, ylim=c(0,NA),
          xlab=NA, ylab=NA, xlabels=F, ylabels=F, mar=c(0,0,0,0)+1)
    }
for(i in 1:prod(plotgrid)){
    tplot(x=0:thmaxicovs['Vi1'], y=vi1samples[,i], lty=1, lwd=2, ylim=c(0,NA),
          xlab=NA, ylab=NA, xlabels=F, ylabels=F, mar=c(0,0,0,0)+1)
    }
for(i in 1:prod(plotgrid)){
    tplot(x=0:thmaxicovs['Vi2'], y=vi2samples[,i], lty=1, lwd=2, ylim=c(0,NA),
          xlab=NA, ylab=NA, xlabels=F, ylabels=F, mar=c(0,0,0,0)+1)
    }
dev.off()


nn <- 32L;
sz <- 50L;
nk <- 64L;
shape1 <- 2;
shape2 <- 1;
prob <- (1:sz);
##alph <- 1/nk;
bsamples <- dbinom(array(rep(0:sz,nk*nn),dim=c(sz+1,nn,nk)), prob=rbeta(nk*nn, shape1=shape1, shape2=shape2), size=sample(x=1:sz,size=nk*nn,replace=TRUE,prob=prob))
qsamples <- LaplacesDemon::rdirichlet(n=nn, alph=rep(1/nk,nk));
isamples <- t(apply(bsamples,1,function(x){rowSums(qsamples*x)}))
##
par(mfrow=c(4,8),mar=c(0,0,0,0)); for(i in 1:nn){tplot(x=0:sz,isamples[,i],lty=1,lwd=2,ylim=c(0,NA),xlim=c(-1,sz+1),xlab=NA,ylab=NA,xlabels=F,ylabels=F,mar=c(0,0,0,0)+1)}



nn <- 32L;
nk <- 64L;
shapem <- 1/2
shapev <- 1/2
xgrid <- seq(-10,10,length.out=256)
##alph <- 1/nk;
bsamples <- dnorm(array(rep(xgrid,nk*nn),dim=c(256,nn,nk)), mean=rnorm(nk*nn, mean=0, sd=1/sqrt(rgamma(nk*nn,shape=shapem, rate=1))), sd=1/sqrt(rgamma(nk*nn,shape=shapev,rate=rgamma(nk*nn,shape=shapev,rate=1))))
qsamples <- LaplacesDemon::rdirichlet(n=nn, alph=rep(1/nk,nk));
isamples <- t(apply(bsamples,1,function(x){rowSums(qsamples*x)}))
##
par(mfrow=c(4,8),mar=c(0,0,0,0));
for(i in 1:nn){tplot(x=xgrid,y=isamples[,i],lty=1,lwd=2,ylim=c(0,NA),xlab=NA,ylab=NA,xlabels=F,ylabels=F,mar=c(0,0,0,0)+1)}
dev.off()
tplot(x=xgrid,y=isamples[,1],lty=1,lwd=2,ylim=c(0,NA),xlab=NA,ylab=NA)

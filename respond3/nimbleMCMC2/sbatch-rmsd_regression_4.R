## Author: PGL  Porta Mana
## Created: 2021-03-20T10:07:17+0100
## Last-Updated: 2021-10-04T20:00:33+0200
################
## Script for direct regression, continuous RMSD
################

.libPaths(c("/cluster/home/user1/R",.libPaths()))
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
#library('LaplacesDemon') # used for Dirichlet generator
## library('ash')
## library('extraDistr')
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
meansccovs <- apply(alldata[1:ndata,..continuousCovs],2,mean)
varsccovs <- apply(alldata[1:ndata,..continuousCovs],2,function(x)var(x, na.rm=T))
## shape & scale parameters for the gamma distribution for tau
tauQccovs <- sapply(continuousCovs, function(acov){
    fn <- function(parms){
        (pinvgamma(varsccovs[acov]/sqrt(10), shape=parms[1], scale=parms[2]) - 0.005)^2 +
            (pinvgamma(varsccovs[acov]*sqrt(10), shape=parms[1], scale=parms[2]) - 0.995)^2
    }
    optim(c(1, 1), fn=fn,
              gr = function(x) pracma::grad(fn, x), 
              method = "L-BFGS-B",
              lower = 0, upper = Inf,
              control = list(factr = 1e-10, maxit = 100))$par
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
meansdcovs <- apply(alldata[1:ndata,..discreteCovs],2,mean)
varsdcovs <- apply(alldata[1:ndata,..discreteCovs],2,function(x)var(x, na.rm=T))
maxdcovs <- apply(alldata[1:ndata,..discreteCovs],2,max)
## prob and size parameters for the neg-binomial distribution for sizeD
sizeQdcovs <- sapply(discreteCovs, function(acov){
    fn <- function(parms){
        (pnbinom(round(maxdcovs[acov]/sqrt(10)), prob=parms[1], size=parms[2]) - 0.005)^2 +
            (pnbinom(round(maxdcovs[acov]*sqrt(10)), prob=parms[1], size=parms[2]) - 0.995)^2
    }
    optim(c(0.5, 1), fn=fn,
              gr = function(x) pracma::grad(fn, x), 
              method = "L-BFGS-B",
              lower = c(0.001,0), upper = c(1,Inf),
              control = list(factr = 1e-10, maxit = 1000))$par
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
## shape1 and shape2 parameters for the beta distribution for probD
shapesratio <- (maxdcovs-meansdcovs)/meansdcovs
alphadcovs <- sapply(discreteCovs, function(acov){
    fn <- function(parms){
        parms2 <- parms*shapesratio[acov]
        -(lbeta(parms,parms2) - (parms-1)*digamma(parms) - (parms2-1)*digamma(parms2) + (parms+parms2-2)*digamma(parms+parms2))
    }
    optim(1, fn=fn,
              gr = function(x) pracma::grad(fn, x), 
              method = "L-BFGS-B",
              lower = 0, upper = Inf,
              control = list(factr = 1e-10,
                             maxit = 100))$par
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
    C=rep(1,ndata) # rcat(n=ndata, prob=rep(1/nclusters,nclusters))
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
confmodel <- configureMCMC(Cmodel, monitors=c('q','meanC', 'tauC', 'probD', 'sizeD')) #, control=list(adaptive=FALSE))
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
version <- 'postHM10'
gc()
totalruntime <- Sys.time()
mcsamples <- runMCMC(Cmcmcsampler, nburnin=1, niter=102, thin=1, inits=initsFunction, setSeed=149)
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
## mcsamples <- readRDS(file='_mcsamples-Rpost9thinTpsC-V7-D6000-K100-I2000.rds')
parmList <- mcsamples2parmlist(mcsamples, parmNames=c('q', 'meanC', 'tauC', 'probD', 'sizeD'))
momentstraces <- moments12Samples(parmList)
allmomentstraces <- cbind(Dcov=plogis(momentstraces$Dcovars, scale=1/10), momentstraces$means, log(momentstraces$vars), momentstraces$covars)
##
diagnESS <- LaplacesDemon::ESS(allmomentstraces)
diagnBMK <- LaplacesDemon::BMK.Diagnostic(allmomentstraces, batches=2)[,1]
diagnMCSE <- 100*LaplacesDemon::MCSE(allmomentstraces)/apply(allmomentstraces, 2, sd)
diagnStat <- apply(allmomentstraces, 2, function(x){LaplacesDemon::is.stationary(as.matrix(x,ncol=1))})
diagnBurn <- apply(allmomentstraces, 2, function(x){LaplacesDemon::burnin(x)})
## print(summary(ess))
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
##
##
timecount <- Sys.time()
plan(sequential)
plan(multisession, workers = 6L)
samplesQuantiles <- foreach(asample=seq_len(nrow(parmList$q)), .combine=c)%:%foreach(acov=covNames, .combine=c)%dopar%{
    if(acov %in% continuousCovs){
        mixq <- function(x){sum(parmList$q[asample,] * pnorm(x, mean=parmList$meanC[asample,acov,], sd=1/sqrt(parmList$tauC[asample,acov,])))}
        fn <- function(z){(mixq(z[1]) - 0.005)^2 + (mixq(z[2]) - 0.995)^2}
        out <- optim(c(meansccovs[acov], meansccovs[acov]), fn,
              gr = function(x) pracma::grad(fn, x), 
              method = "L-BFGS-B",
              lower = -Inf, upper = Inf,
              control = list(factr = 1e-10,
                             maxit = 100))$par
    }else{
        searchgrid <- 0:max(parmList$sizeD[asample,acov,])
        dq <- colSums(c(parmList$q[asample,]) * pbinom(matrix(searchgrid, ncol=length(searchgrid), nrow=ncol(parmList$q), byrow=TRUE), prob=parmList$probD[asample,acov,], size=parmList$sizeD[asample,acov,]))
        out <- c(which.min(abs(dq - 0.005))-1, which.min(abs(dq - 0.995))-1)        
    }
    out
}
plan(sequential)
print(Sys.time()-timecount)
##
dim(samplesQuantiles) <- c(2, nccovs+ndcovs, nrow(parmList$q))
samplesQuantiles <- aperm(samplesQuantiles)
dimnames(samplesQuantiles) <- list(NULL, covNames, c('0.5%', '99.5%'))
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
alldataRanges <- dataQuantiles <- xlimits <- list()
for(acov in covNames){
    dataQuantiles[[acov]] <- quantile(alldata[[acov]], prob=c(0.005,0.995))
    alldataRanges[[acov]] <- range(alldata[[acov]])
    xlimits[[acov]] <- range(c(alldataRanges[[acov]], samplesQuantiles[,acov,]))
}
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
if(FALSE){
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


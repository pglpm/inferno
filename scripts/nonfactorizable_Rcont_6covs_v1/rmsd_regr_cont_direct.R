## Author: PGL  Porta Mana
## Created: 2021-03-20T10:07:17+0100
## Last-Updated: 2021-08-25T07:31:32+0200
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
## specify hyperparameters for hyperpriors of continuous variables
murmsd <- 1
musasa <- 1
mutani <- 0
mu0 <- c(murmsd,musasa,mutani)
names(mu0) <- c('RMSD','sasa','tanimoto')
varmurmsd <- 2^2
varmusasa <- 2^2
varmutani <- 2^2
sigmu0 <- c(varmurmsd,varmusasa,varmutani)
names(sigmu0) <- names(mu0)
expvarrmsd <- (1/2)^2
expvarsasa <- (1/2)^2
expvartani <- (1/2)^2
expvar0 <- c(expvarrmsd,expvarsasa,expvartani)
names(expvar0) <- names(mu0)
## diagvar0 <- (2)^2
## df0 <- (2*expvarsasa^2)/diagvar0 + 4
df0 <- 30
##
tanimoto2y <- function(x){
    -log(1/x-1)/2 ## faster than qlogis(x, scale=1/2)
}
##
## correct compared with study_prior_bayes_regr.nb
Dtanimoto2y <- function(x){
    1/(2*x*(1-x))
}
##
## y2tanimoto <- function(y){
##     1/2 + atan(y)/pi
## }
##
sasa2y <- function(x){
   log(x)/4
}
##
## correct compared with study_prior_bayes_regr.nb
Dsasa2y <- function(x){
    1/(4*x)
}
##
## Read and reorganize data
rm(data)
data <- fread('../processed_data_scaled.csv', sep=' ')
nameFeatures <- names(data)
nSamples <- nrow(data)
nFeatures <- ncol(data)
##
##
## GOOD values:
## sdsasa <- 1
## sdtani <- 1
## mymu0 <- c(5,0)
## mynu0 <- 20+1
## mykappa0 <- 0.05
## coefdiag <- (mykappa0+1)/(mykappa0*(mynu0-length(mymu0)-1))
## #diag(c(sdsasa,sdtani)^2/coefdiag)
## testhp <- setHyperparams(mu0=mymu0, kappa0=mynu0, R0=solve(diag(c(sdsasa,sdtani)^2/coefdiag)),nu0=mykappa0)
##
## sdsasa <- 1/2
## sdtani <- 1
## mymu0 <- c(5,0)
## mynu0 <- 20+1
## mykappa0 <- 0.05
## coefdiag <- (mykappa0+1)/(mykappa0*(mynu0-length(mymu0)-1))
## #diag(c(sdsasa,sdtani)^2/coefdiag)
## testhp <- setHyperparams(mu0=mymu0, kappa0=mynu0, R0=solve(diag(c(sdsasa,sdtani)^2/coefdiag)),nu0=mykappa0)
##
## sdsasa <- 1/2
## sdtani <- 1
## mymu0 <- c(5,0)
## mynu0 <- 50+1
## mykappa0 <- 0.05
## coefdiag <- (mykappa0+1)/(mykappa0*(mynu0-length(mymu0)-1))
## #diag(c(sdsasa,sdtani)^2/coefdiag)
## alpha <- -2
## testhp <- setHyperparams(mu0=mymu0, kappa0=mynu0, R0=solve(diag(c(sdsasa,sdtani)^2/coefdiag)),nu0=mykappa0)
##
## sdsasa <- 1/2
## sdtani <- 9/10
## mymu0 <- c(4,0)
## mynu0 <- 30+1
## mykappa0 <- 0.1
## coefdiag <- (mykappa0+1)/(mykappa0*(mynu0-length(mymu0)-1))
## #diag(c(sdsasa,sdtani)^2/coefdiag)
## testhp <- setHyperparams(mu0=mymu0, kappa0=mynu0, R0=solve(diag(c(sdsasa,sdtani)^2))*coefdiag,nu0=mykappa0)
## alpha <- 4
##
## NB: nu here = kappa in fmri paper
## kappa here = nu in fmri paper
##################################################
## Mixed-x, no y-model
## 5 covs, 5000 points: 14844 s
## 6 covs, 5000 pts, 1000e3+1000e3 its: 9.696259 hours
## 6 covs, 6000 pts, 2000e3+1000e3 its: 15.75 hours
## 6 covs, 6000 pts, 3000e3+2000e3 its: 1.138318 days

ndata <- 6000 # nSamples = 37969
#set.seed(222)
rmsdCol <- which(names(data)=='log_RMSD')
covNums <- which(colnames(data) %in%  c('log_RMSD',
                                        'scale_mcs_unbonded_polar_sasa',
                                        'scale_ec_tanimoto_similarity',
                                        'mcs_NumHeteroAtoms',
                                        ##'scale_fc_tanimoto_similarity'
                                        'docked_HeavyAtomCount',
                                        'mcs_RingCount',
                                        'docked_NumRotatableBonds'
                                        ))
covNames <- names(data)[covNums]
discreteCovs <- covNames[sapply(covNames, function(x){is.integer(data[[x]])})]
continuousCovs <- covNames[sapply(covNames, function(x){is.double(data[[x]])})]
allCovNums <- c(rmsdCol, covNums)
##
dimsC <- length(continuousCovs)
ilogrmsd <- continuousCovs[grepl('RMSD', continuousCovs)]
itanimoto <- continuousCovs[grepl('tanimoto', continuousCovs)]
isasa <- continuousCovs[grepl('sasa', continuousCovs)]
## Hyperparameters
hmu0 <- sapply(continuousCovs, function(x){mu0[sapply(names(mu0),function(y){grepl(y,x)})]})
hnu0 <- df0 + dimsC - 1
hkappa0 <- expvarsasa/varmusasa
## hkappa0 <- 0.01
hDelta0 <- solve(diag(sapply(continuousCovs, function(x){expvar0[sapply(names(expvar0),function(y){grepl(y,x)})]})))/(df0-2)
colnames(hDelta0) <- rownames(hDelta0) <- names(hmu0)
##
testhp <- setHyperparams(mu0=hmu0, kappa0=hnu0, R0=hDelta0, nu0=hkappa0)
##
rmsdVals <- 1:3
directcases <- c('fdata')
##
rm(mcmcrun)
gc()
starttime <- Sys.time()
plan(sequential)
##plan(multisession, workers = 3L)
mcmcrun <- foreach(case=directcases, .inorder=FALSE)%do%{
    outfile <- paste0('_mcoutput_',case)
    ##
    if(case==directcases[1]){
        seldata <- 1:ndata
        datamcr <- data[seldata, covNames, with=F]
    } else if(case==directcases[2]){
        print('THIS CASE SHOULD BE EXCLUDED')
        ## construct tranining set with equally occurring freqs of bin_RMSD
        testd <- data[, covNames, with=F]
        datamcr <- data.table()
        seldata <- c()
        for(bin in rmsdVals){
            whichbin <- which(testd$bin_RMSD==bin)[1:round(ndata/3)]
            datamcr <- rbind(datamcr, testd[whichbin])
            seldata <- c(seldata,whichbin)
        }
        rm(testd)
    } else {print('Error!')
    datamcr <- NULL}
    ## make sure all levels of discrete variables appear
    for(i in discreteCovs){
        datum <- datamcr[[i]]
        datum <- datum[!is.na(datum)]
        levels <- as.numeric(names(table(data[[i]])))
        for(level in setdiff(min(levels):max(levels), as.numeric(names(table(datum))))){ #print(paste0(case,' ',i,' ',level))
            dt <- data.table(level)
            names(dt) <- i
            datamcr <- rbind(datamcr, dt, fill=TRUE)
        }
    }
##
##
        c(case=case, profRegr(excludeY=TRUE, xModel='Mixed', nSweeps=2000e3, nBurn=3000e3, nFilter=2e3, nProgress=100e3, data=as.data.frame(datamcr), nClusInit=80, covNames=c(discreteCovs,continuousCovs), discreteCovs=discreteCovs, continuousCovs=continuousCovs, seed=147, output=outfile, useHyperpriorR1=FALSE, useNormInvWishPrior=TRUE, hyper=testhp, alpha=4))
}
plan(sequential)
names(mcmcrun) <- paste0('freqs',sapply(mcmcrun,function(i){i$case}))
elapsedtime <- Sys.time() - starttime
elapsedtime
## 2000: 7.64 min
## 2000e1: 49.2 min
## 2000e2: 6.41 h
## Save MCMC samples
MCMCdata <- as.list(rep(NA,length(mcmcrun)))
names(MCMCdata) <- names(mcmcrun)
##
for(case in directcases){
    outfile <- paste0('_mcoutput_',case)
    testmc <- mcmcrun[[paste0('freqs',case)]]
    ## log-likelihood and log-posteriors
    fc <- file(paste0(outfile,"_logPost.txt"))
    logPost <- sapply(strsplit(readLines(fc), " +"), as.numeric)
    rownames(logPost) <- c('log-post','log-likelihood','log-prior')
    close(fc)
    ## Samples of numbers of clusters
    fc <- file(paste0(outfile,'_nClusters.txt'))
    nList <- sapply(strsplit(readLines(fc), " +"), as.integer)
    close(fc)
    ## Samples of Dirichlet-process alpha
    fc <- file(paste0(outfile,'_alpha.txt'))
    alphaList <- sapply(strsplit(readLines(fc), " +"), as.numeric)
    close(fc)
    ## Samples of cluster weights
    fc <- file(paste0(outfile,'_psi.txt'))
    psiList <- lapply(strsplit(readLines(fc), " +"), function(x){x <- as.numeric(x); x[x>=0]})
    close(fc)
    ##  Samples of Dirichlet-distribution phis (discrete covariates)
    fc <- file(paste0(outfile,'_phi.txt'))
    phiList <- lapply(strsplit(readLines(fc), " +"), function(x){x <- as.numeric(x); x[x>=0]})
    close(fc)
    ##  Samples of normal-distribution means (continuous covariates)
    fc <- file(paste0(outfile,'_mu.txt'))
    muList <- lapply(strsplit(readLines(fc), " +"), as.numeric)
    close(fc)
    ##  Samples of normal-distribution covariances (continuous covariates)
    fc <- file(paste0(outfile,'_Sigma.txt'))
    sigmaList <- lapply(strsplit(readLines(fc), " +"), as.numeric)
    close(fc)
    ##
    nContCov <- testmc$nContinuousCovs
    nDiscrCov <- testmc$nDiscreteCovs
    nCat <- testmc$nCategories
    cumcats <- c(0,cumsum(nCat))
    for(i in 1:length(nList)){
        catgroups <- cumcats*nList[i]
        datum <- phiList[[i]]
        phiList[[i]] <- lapply(1:length(testmc$nCategories), function(j){
            y <- datum[(catgroups[j]+1):catgroups[j+1]]
            dim(y) <- c(nList[i], testmc$nCategories[j])
            aperm(y)
        })
        names(phiList[[i]]) <- testmc$discreteCovs
        dim(sigmaList[[i]]) <- c(nList[i], nContCov, nContCov)
        sigmaList[[i]] <- aperm(sigmaList[[i]])
        rownames(sigmaList[[i]]) <- testmc$continuousCovs
        colnames(sigmaList[[i]]) <- testmc$continuousCovs
        dim(muList[[i]]) <- c(nList[i], nContCov)
        muList[[i]] <- aperm(muList[[i]])
        rownames(muList[[i]]) <- testmc$continuousCovs
    }
    ## Save the samples above
    MCMCdata[[paste0('freqs',case)]] <- list(case=case, nList=nList, alphaList=alphaList, psiList=psiList, phiList=phiList, muList=muList, sigmaList=sigmaList, logPost=logPost)
}
##
save.image(file=paste0('_directmodel_contR_N',ndata,'_',length(covNums),'covs.RData'))
##
## Diagnostic plots
pdff('mcsummary')
for(j in 1:length(MCMCdata)){
    sd <- MCMCdata[[j]]
    matplot(sd$logPost[2,], type='l',ylim=range(sd$logPost[2,],na.rm=T,finite=T),ylab='log-likelihood',col=palette()[j], main=paste0('freqs: ',directcases[j]))
    }
for(j in 1:length(MCMCdata)){
    sd <- MCMCdata[[j]]
matplot(sd$nList,type='l',ylim=range(sd$nList,na.rm=T,finite=T),ylab='no. clusters',col=palette()[j], main=paste0('freqs: ',directcases[j]))
}
for(j in 1:length(MCMCdata)){
    sd <- MCMCdata[[j]]
    matplot(sd$alphaList,type='l',ylim=range(sd$alphaList,na.rm=T,finite=T),ylab='alpha',col=palette()[j], main=paste0('freqs: ',directcases[j]))
    }
for(j in 1:length(MCMCdata)){
    sd <- MCMCdata[[j]]
    for(i in c(1,3)){matplot(sd$logPost[i,],type='l',ylim=range(sd$logPost[i,],na.rm=T,finite=T),col=palette()[j], main=paste0('freqs: ',directcases[j]))}
    }
dev.off()
##

## Plots of sample distributions
isasa <- continuousCovs[grepl('sasa', continuousCovs)]
itani <- continuousCovs[grepl('tanimoto', continuousCovs)]
## Jacobian determinant for scaled quantities
jac <- function(X){
   Dtanimoto2y(X[itani])*Dsasa2y(X[isasa])
}
##
predictSasaTani <- function(dataobj, X){
    foreach(asample=seq_along(dataobj$nList), .combine='+', .inorder=FALSE)%dopar%{
        sum(dataobj$psiList[[asample]] *
            sapply(seq_len(dataobj$nList[asample]),function(cluster){
            ## prod(sapply(discreteCovs, function(covariate){
            ##     dataobj$phiList[[asample]][[covariate]][X[covariate], cluster]
            ## })) *
                dmvnorm(X, mean=dataobj$muList[[asample]][continuousCovs,cluster], sigma=as.matrix(dataobj$sigmaList[[asample]][continuousCovs,continuousCovs,cluster]))
        }))
       #drop(dataobj$phiList[[asample]]$bin_RMSD %*% weights)/sum(weights)
}/length(dataobj$nList)
}
predictSampleSasaTani <- function(dataobj,asample, X){
        sum(dataobj$psiList[[asample]] *
            sapply(seq_len(dataobj$nList[asample]),function(cluster){
            ## prod(sapply(discreteCovs, function(covariate){
            ##     dataobj$phiList[[asample]][[covariate]][X[covariate], cluster]
            ## })) *
                dmvnorm(X, mean=dataobj$muList[[asample]][continuousCovs,cluster], sigma=as.matrix(dataobj$sigmaList[[asample]][continuousCovs,continuousCovs,cluster]))
        }))
       #drop(dataobj$phiList[[asample]]$bin_RMSD %*% weights)/sum(weights)
}
##
##
## Plot sample densities
nplots <- 50
plan(sequential)
plan(multisession, workers = 3L)
xgrid <- seq(1, 200,length.out=20)
ygrid <- seq(0.001,0.999,length.out=20)
for(j in 1:length(MCMCdata)){
pdf(file=paste0('testsampledensities_seq_',directcases[j],'.pdf'),height=11.7,width=16.5*10)
## plot(NA,xlim=c(0,1),ylim=c(0,1),main='parameters')
## text(0,1,(paste0('transf = ','id','\ndim = ',dims,'\nmeanalpha = ',meanalpha,'\nsdalpha = ',sdalpha,'\nK = ',ka,'\nL = ',la,'\nM = ',mu,'\nresc = ',resc,'\nNclust = ',mean(clusters))),adj = c(0,1),cex=5)
##
for(asample in round(seq(1,1000,length.out=nplots))){
zgrid <- outer(xgrid,ygrid,function(x,y){
    mapply(function(xx,yy){
                    XX <- c(xx,yy)
                    names(XX) <- continuousCovs
                    YY <- c(sasa2y(xx),tanimoto2y(yy))
                    names(YY) <- continuousCovs
                    predictSampleSasaTani(MCMCdata[[j]], asample, YY)*jac(XX)
                    },
        x,y)}
        )
persp(xgrid,ygrid,zgrid,zlim=c(0,max(zgrid)),ticktype='detailed',theta = 45, phi = 15,xlab='sasa',ylab='tanimoto')
}
dev.off()
##
##
xgrid <- seq(-3, 3, length.out=20)
ygrid <- seq(-3, 3, length.out=20)
pdf(file=paste0('testsampledensities_scaled_seq_',directcases[j],'.pdf'),height=11.7,width=16.5*10)
## plot(NA,xlim=c(0,1),ylim=c(0,1),main='parameters')
## text(0,1,(paste0('transf = ','id','\ndim = ',dims,'\nmeanalpha = ',meanalpha,'\nsdalpha = ',sdalpha,'\nK = ',ka,'\nL = ',la,'\nM = ',mu,'\nresc = ',resc,'\nNclust = ',mean(clusters))),adj = c(0,1),cex=5)
##
for(asample in round(seq(1,1000,length.out=nplots))){
zgrid <- outer(xgrid,ygrid,function(x,y){
    mapply(function(xx,yy){
                    XX <- c(xx,yy)
                    names(XX) <- continuousCovs
                    predictSampleSasaTani(MCMCdata[[j]], asample, XX)
                    },
        x,y)}
        )
persp(xgrid,ygrid,zgrid,zlim=c(0,max(zgrid)),ticktype='detailed',theta = 45, phi = 15,xlab='scaled sasa',ylab='scaled tanimoto')
}
dev.off()
}

##################################################################
##################################################################
#### Valuation on test set
##
## Utility matrices
dgain <- diag(1,3)
cgain <- 1-sapply(1:3,function(x){abs(x-1:3)})/2
## function for randomly choosing ties
resample <- function(x, ...) x[sample.int(length(x), ...)]
## function to compute utilities in test set
utilities <- function(truevalues,priorY,utilitym,probX=matrix(rep(1L,length(truevalues)*length(priorY)),length(truevalues))){
    y <- apply(t(tcrossprod(utilitym, normalizem(t(t(probX) * priorY)))), 1,
               function(z){resample(which(z==max(z)))})
    diag(utilitym[truevalues,y])
}
## function to compute various prob-averages in test set
probscores <- function(truevalues,priorY,meanfunction,probX=matrix(rep(1L,length(truevalues)*length(priorY)),length(truevalues))){
    y <- normalizem(t(t(probX) * priorY))
    meanfunction(diag(y[,truevalues]))
}
## function to construct data table with results
metrics <- function(truevalues, probX, priorY, chanceprior=priorY){
    data.table(
        model=c('model', 'chance','min','max'),
        delta_gain=c(
            mean(utilities(truevalues,priorY,dgain,probX)),
            mean(utilities(truevalues,chanceprior,dgain)),
            0,1
        ),
        ##
        contig_gain=c(
            mean(utilities(truevalues,priorY,cgain,probX)),
            mean(utilities(truevalues,chanceprior,cgain)),
            0,1
        ),
        ##
        log_score=c(
            mean(probscores(truevalues,priorY,log,probX)),
            mean(probscores(truevalues,chanceprior,log)),
            -Inf,0
        ),
        ##
        mean_score=c(
            mean(probscores(truevalues,priorY,identity,probX)),
            mean(probscores(truevalues,chanceprior,identity)),
            0,1
        )
    )}
## old function
oldmetrics <- function(testres, priorP){
    data.table(
        model=c('model', 'chance','min','max'),
        delta_gain=c(
            mean(apply(testres,1,function(x){
                y <- c(dgain %*% normalize(x[-1] * priorP))
                dgain[x[1], resample(which(y==max(y)))]
            })),
            mean(apply(testres,1,function(x){
                y <- c(dgain %*% normalize(priorP))
                dgain[x[1], resample(which(y==max(y)))]
            })),
            0,1
        ),
        ##
        contig_gain=c(
            mean(apply(testres,1,function(x){
                y <- c(cgain %*% normalize(x[-1] * priorP))
                cgain[x[1], resample(which(y==max(y)))]
            })),
            mean(apply(testres,1,function(x){
                y <- c(cgain %*% normalize(priorP))
                cgain[x[1], resample(which(y==max(y)))]
            })),
            0,1
        ),
        ##
        log_score=c(
            mean(apply(testres,1,function(x){
                y <- normalize(x[-1] * priorP)
                log(y[x[1]])
            })),
            mean(apply(testres,1,function(x){
                y <- normalize(priorP)
                log(y[x[1]])
            })),
            -Inf,0
        ),
        ##
        mean_score=c(
            mean(apply(testres,1,function(x){
                y <- normalize(x[-1] * priorP)
                y[x[1]]
            })),
            mean(apply(testres,1,function(x){
                y <- normalize(priorP)
                y[x[1]]
            })),
            0,1
        )
    )}
##
source(file='calibration_plots.R')

#### Calculation of utilities and scores for test set
##
## Case with equally occurring RMSD categories
priorP <- rep(1,3)/3
unseldata <- setdiff(1:nrow(data), seldata)
nTest <- 2000
## construct test set
testd <- data[unseldata, c('bin_RMSD',covNames), with=F]
testdata <- data.table()
for(val in rmsdVals){
    testdata <- rbind(testdata, tail(testd[bin_RMSD==val],n=nTest))
}
rm(testd)
##
gc()
plan(sequential)
plan(multisession, workers = 6L)
condfreqs <- predictYX(MCMCdata, testdata, rmsdVals)
plan(sequential)
##
set.seed(247)
scores1 <- metrics(testdata[,bin_RMSD], condfreqs[,,1], priorP)
scores1
## set.seed(247)
## oldmetrics(cbind(testdata[,bin_RMSD], condfreqs[,,1]), priorP)
##
##
calibrationplots(condfreqs, priorP, divs=1, title='score1')
calibrationplots(condfreqs, priorP, divs=2, title='score1')

## Case with unequally occurring RMSD categories
priorP2 <- normalize(as.vector(table(data$bin_RMSD)))
unseldata <- setdiff(1:nrow(data), seldata)
nTest <- 2000
## construct test set
testd <- data[unseldata, c('bin_RMSD',covNames), with=F]
testdata <- data.table()
for(val in rmsdVals){
    testdata <- rbind(testdata, tail(testd[bin_RMSD==val],n=round(3*nTest*priorP2[val])))
}
rm(testd)
##
gc()
plan(sequential)
plan(multisession, workers = 6L)
condfreqs <- predictYX(MCMCdata, testdata, rmsdVals)
plan(sequential)
##
set.seed(247)
scores2 <- metrics(testdata[,bin_RMSD], condfreqs[,,1], priorP2)
scores2
## set.seed(247)
## oldmetrics(cbind(testdata[,bin_RMSD], condfreqs[,,1]), priorP2)
calibrationplots(condfreqs, priorP2, divs=1, title='score2')
calibrationplots(condfreqs, priorP2, divs=2, title='score2')
##
##
##
## Case with unequally occurring RMSD categories, uniform prior
##
set.seed(247)
scores3 <- metrics(testdata[,bin_RMSD], condfreqs[,,1], priorP, priorP2)
scores3
## set.seed(247)
## oldmetrics(cbind(testdata[,bin_RMSD], condfreqs[,,1]), priorP2)
calibrationplots(condfreqs, priorP, divs=1, title='score3')
calibrationplots(condfreqs, priorP, divs=2, title='score3')
##
## Save scores
save(list=c('scores1','scores2','scores3'), file=paste0('scores_T',nTest*3,'_direct_test_N',ndata,'_',length(covNums),'covs.RData'))


##################################################################
##################################################################
#### Calibration
## distfunction <- function(freq1,freq0){
##     sum(freq1*log(freq1/freq0),na.rm=TRUE)
## }
distfunction <- function(freq0,freq1){
    sum(freq1*log(freq1/freq0),na.rm=TRUE)
}
##
divs <- 1
##
nodes <- seq(0,1,1/divs)
centres1 <- nodes[-(divs+1)]+1/3/divs
centres2 <- nodes[-(divs+(0:1))]+2/3/divs
##
pbins <- rbind(
    foreach(i1=1:(divs+1), .combine=rbind)%:%foreach(i3=1:(divs-i1+2), .combine=rbind)%do%{
    p1 <- nodes[i1]
    p3 <- nodes[i3]
    c(p1, 1-p1-p3, p3)
    })
rownames(pbins) <- paste0('bin',1:nrow(pbins))
#pbins <- normalizem(pbins + 1e-6)
##
## pbins <- rbind(
##     foreach(i1=1:divs, .combine=rbind)%:%foreach(i3=1:(divs-i1+1), .combine=rbind)%do%{
##     p1 <- centres1[i1]
##     p3 <- centres1[i3]
##     c(p1, 1-p1-p3, p3)
## },
##     foreach(i3=1:(divs-1), .combine=rbind)%:%foreach(i1=1:(divs-i3), .combine=rbind)%do%{
##     p1 <- centres2[i1]
##     p3 <- centres2[i3]
##     c(p1, 1-p1-p3, p3)
##     })
## rownames(pbins) <- paste0('bin',1:nrow(pbins))
##
##
freqbins <- 0 * pbins
psample <- normalizem(t(t(condfreqs[,,1]) * priorP))
for(asample in 1:nrow(condfreqs)){
    distances <- apply(pbins,1,function(x){distfunction(psample[asample,],x)})
    bin <- which.min(distances)
    outcome <- testdata[asample,bin_RMSD]
    freqbins[bin,outcome] <- freqbins[bin,outcome] + 1
}
##
wfreqbins <- rowSums(freqbins)
rfreqbins <- freqbins/wfreqbins
##
##
pdf(file=paste0('calibration_plots_1_',nrow(pbins),'bins.pdf'),height=11.7,width=11.7)
matplot(x=rbind(nodes[1],nodes[divs+1],nodes[1],nodes[1]), y=rbind(nodes[1],nodes[1],nodes[divs+1],nodes[1]), type='l', lty=1, xlim=c(0,1), ylim=c(0,1), col=mygrey, xlab='p(1)', ylab='p(3)', cex.axis=2, cex.lab=2, cex.main=2, main='distribution of posterior probabilities')
psample <- normalizem(t(t(condfreqs[,,1]) * priorP))
matpoints(x=psample[,1],y=psample[,3],type='p',pch=1,col=mygrey)
##
matplot(x=rbind(nodes[1],nodes[divs+1],nodes[1],nodes[1]), y=rbind(nodes[1],nodes[1],nodes[divs+1],nodes[1]), type='l', lty=1, xlim=c(0,1), ylim=c(0,1), col=mygrey, xlab='p(1)', ylab='p(3)', cex.axis=2, cex.lab=2, cex.main=2, main='binning of probability simplex')
psampled <- rdirichlet(n=10000, alpha=rep(1,3))
for(i in 1:10000){
    dists <- apply(pbins,1,function(x){distfunction(psampled[i,],x)})
    cent <- which.min(dists)
    matpoints(x=psampled[i,1],y=psampled[i,3],type='p',pch=20,col=palette()[(cent%%7)+1])
}
matpoints(x=pbins[,1],y=pbins[,3],type='p',pch=18,cex=3,col='black')
##
for(bin in order(wfreqbins, decreasing=TRUE)){
    matplot(x=rbind(nodes[1],nodes[divs+1],nodes[1],nodes[1]), y=rbind(nodes[1],nodes[1],nodes[divs+1],nodes[1]), type='l', lty=1, xlim=c(0,1), ylim=c(0,1), col=mygrey, xlab='p(1)', ylab='p(3)', cex.axis=2, cex.lab=2, cex.main=2, main=paste0('W = ',wfreqbins[bin]))
    matpoints(x=pbins[,1], y=pbins[,3], type='p', pch=15, col='black')
    for(i in 1:10000){
        dists <- apply(pbins,1,function(x){distfunction(psampled[i,],x)})
        if(bin==which.min(dists)){
            matpoints(x=psampled[i,1],y=psampled[i,3],type='p',pch=20,col=mygrey)
        }
    }
#    matlines(x=rbind(pbins[bin,1],rfreqbins[bin,1]),y=rbind(pbins[bin,3],rfreqbins[bin,3]),type='l',lty=1,col=myblue)
    matpoints(x=pbins[bin,1],y=pbins[bin,3],type='p',pch=18,cex=3,col='black')
    matpoints(x=rfreqbins[bin,1],y=rfreqbins[bin,3],type='p',pch=20,cex=4,col=myred)
}
dev.off()





##################################################
##################################################
##################################################
##################################################
#### OLD STUFF
##################################################
##################################################
##################################################
##################################################










##
#### Evaluation
##
##discrMin <- sapply(data[,discreteCovs,with=F], min)-1
##
priorP <- rep(1,3)/3
unseldata <- setdiff(1:nrow(data), seldata)
##
## 500, 6 threads: 1269.53  954.36 9520.95  
nTest <- 1000
testd <- data[unseldata, c('bin_RMSD',covNames), with=F]
testdata <- data.table()
for(val in rmsdVals){
    testdata <- rbind(testdata, tail(testd[bin_RMSD==val],n=nTest))
}
rm(testd)
##
gc()
plan(sequential)
plan(multisession, workers = 6L)
condfreqs0 <- lapply(rmsdVals,function(val){predictYXcross(MCMCdata[[val]], testdata)})
plan(sequential)
##

plan(sequential)
plan(multisession, workers = 6L)
condfreqs <- foreach(val=rmsdVals, .combine=cbind)%dorng%{
    me <- condfreqs0[[val]]$means
    si <- diag(condfreqs0[[val]]$covariance)/100
    c(rgamma(n=length(me), shape=me^2/si, scale=si/me))}
plan(sequential)
evals1 <- metrics(testdata[,bin_RMSD], condfreqs, priorP)
evals1
#oldmetrics(cbind(testdata[,bin_RMSD], condfreqs), priorP)
##


plan(sequential)
plan(multisession, workers = 6L)
condfreqs <- foreach(val=rmsdVals, .combine=cbind)%dorng%{
    c(rmvnorm(n=1, mean=condfreqs0[[val]]$means, sigma=condfreqs0[[val]]$covariance/100))}
plan(sequential)
evals1 <- metrics(testdata[,bin_RMSD], condfreqs, priorP)
evals1
oldmetrics(cbind(testdata[,bin_RMSD], condfreqs), priorP)
##




evals1 <- metrics(testres, priorP)
save.image(file=paste0('_direct_test_N',ndata,'_',length(covNums),'covs.RData'))
evals1
##
predictYpar <- function(dataobj, X){
    freqs <- foreach(asample=seq_along(dataobj$nList), .combine=c, .inorder=FALSE)%dopar%{
        sum(dataobj$psiList[[asample]] *
            sapply(seq_len(dataobj$nList[asample]),function(cluster){
            prod(sapply(discreteCovs, function(covariate){
                dataobj$phiList[[asample]][[covariate]][X[covariate], cluster]
            })) *
                dmvnorm(X[continuousCovs], mean=dataobj$muList[[asample]][continuousCovs,cluster], sigma=as.matrix(dataobj$sigmaList[[asample]][continuousCovs,continuousCovs,cluster]))
        }))
       #drop(dataobj$phiList[[asample]]$bin_RMSD %*% weights)/sum(weights)
    }
    me <- mean(freqs)
    c(me, sqrt(mean(freqs^2) - me^2))
}
##

ntdata <- 3000

plan(multisession, workers = 6L)
tstart <- Sys.time()
memtest0 <- t(sapply(rmsdVals,function(val){predictYtest(MCMCdata[[val]], data[30000+(1:ntdata)])}))
dim(memtest0) <- c(3,ntdata,2)
Sys.time()-tstart
plan(sequential)

plan(multisession, workers = 6L)
tstart <- Sys.time()
memtest1 <- t(sapply(rmsdVals,function(val){predictYXtestsingle(MCMCdata[[val]], data[30000+(1:ntdata)])}))
dim(memtest1) <- c(3,ntdata,2)
Sys.time()-tstart
plan(sequential)

plan(multisession, workers = 6L)
tstart <- Sys.time()
memtest2 <- aperm(predictYXtest(MCMCdata, data[30000+(1:ntdata)], rmsdVals),c(2,1,3))
Sys.time()-tstart
plan(sequential)

plan(multisession, workers = 6L)
tstart <- Sys.time()
memtest3 <- aperm(predictYXtests(MCMCdata, data[30000+(1:ntdata)], rmsdVals), c(2,1,3))
Sys.time()-tstart
plan(sequential)



plan(multisession, workers = 6L)
tstart <- Sys.time()
memtest <- predictYtest(MCMCdata[[1]], data[36000:36100])
Sys.time()-tstart
plan(sequential)





unseldata <- setdiff(1:nrow(data), seldata)
##
## 500, 6 threads: 1269.53  954.36 9520.95  
nTest <- 500
testdata <- data[unseldata, c('bin_RMSD',covNames), with=F]
testd <- data.table()
for(val in rmsdVals){
    testd <- rbind(testd, tail(testdata[bin_RMSD==val],n=nTest))
}
testdata <- testd
rm(testd)
##
rm(testres)
gc()
system.time(testres <- t(apply(testdata, 1, function(datum){
    c(datum['bin_RMSD'], sapply(rmsdVals,function(val){predictYpar(MCMCdata[[val]],datum)}))
})))
plan(sequential)
##
evals1 <- metrics(testres, priorP)
save.image(file=paste0('_direct_test_N',ndata,'_',length(covNums),'covs.RData'))
evals1
## 6 covs, 5000 pts, 9735 s
##     model delta_gain contig_gain log_score mean_score
## 1:  model  0.5080000   0.7036667 -1.106211  0.4532316
## 2: chance  0.3333333   0.6666667 -1.098612  0.3333333
## 3:    min  0.0000000   0.0000000      -Inf  0.0000000
## 4:    max  1.0000000   1.0000000  0.000000  1.0000000
##
## 5 covs, 5000 points, 7348 s
##     model delta_gain contig_gain log_score mean_score
## 1:  model  0.5046667   0.6993333 -1.065359  0.4354056
## 2: chance  0.3333333   0.6666667 -1.098612  0.3333333
## 3:    min  0.0000000   0.0000000      -Inf  0.0000000
## 4:    max  1.0000000   1.0000000  0.000000  1.0000000
## > > > >     model delta_gain contig_gain log_score mean_score
## 1:  model  0.5053333   0.7020000 -1.067252  0.4357623
## 2: chance  0.3333333   0.6666667 -1.098612  0.3333333
## 3:    min  0.0000000   0.0000000      -Inf  0.0000000
## 4:    max  1.0000000   1.0000000  0.000000  1.0000000
##
## 5 covs, 3500 points, 7305 s
##     model delta_gain contig_gain log_score mean_score
## 1:  model  0.4906667   0.6933333 -1.087008  0.4355049
## 2: chance  0.3333333   0.6666667 -1.098612  0.3333333
## 3:    min  0.0000000   0.0000000      -Inf  0.0000000
## 4:    max  1.0000000   1.0000000  0.000000  1.0000000
##
## 5 covs, 2500 points, 7016 s
##     model delta_gain contig_gain log_score mean_score
## 1:  model  0.4926667   0.6846667 -1.133043  0.4361120
## 2: chance  0.3333333   0.6666667 -1.098612  0.3333333
## 3:    min  0.0000000   0.0000000      -Inf  0.0000000
## 4:    max  1.0000000   1.0000000  0.000000  1.0000000
##
## 4 covs, 3500 points, 7484 s
##     model delta_gain contig_gain log_score mean_score
## 1:  model  0.4766667   0.6816667 -1.071031  0.4194373
## 2: chance  0.3333333   0.6666667 -1.098612  0.3333333
## 3:    min  0.0000000   0.0000000      -Inf  0.0000000
## 4:    max  1.0000000   1.0000000  0.000000  1.0000000
##
## 4 covs, 2500 points, 6483 s
##     model delta_gain contig_gain log_score mean_score
## 1:  model  0.4720000   0.6806667 -1.089924  0.4195193
## 2: chance  0.3333333   0.6666667 -1.098612  0.3333333
## 3:    min  0.0000000   0.0000000      -Inf  0.0000000
## 4:    max  1.0000000   1.0000000  0.000000  1.0000000
##
## 4 covs, 1500 points
##     model delta_gain contig_gain log_score mean_score
## 1:  model  0.4666667   0.6696667 -1.161806  0.4261656
## 2: chance  0.3333333   0.6666667 -1.098612  0.3333333
## 3:    min  0.0000000   0.0000000      -Inf  0.0000000
## 4:    max  1.0000000   1.0000000  0.000000  1.0000000
##
##     model delta_gain contig_gain log_score mean_score
## 1:  model  0.4586667   0.6756667 -1.631715  0.4247677
## 2: chance  0.3333333   0.6666667 -1.098612  0.3333333
## 3:    min  0.0000000   0.0000000      -Inf  0.0000000
## 4:    max  1.0000000   1.0000000  0.000000  1.0000000
##
##
#### Evaluation 2
##
##discrMin <- sapply(data[,discreteCovs,with=F], min)-1
discreteCovs <- discreteCovs
continuousCovs <- continuousCovs
plan(sequential)
plan(multisession, workers = 6L)
priorP2 <- normalize(as.vector(table(data$bin_RMSD)))
##
predictYpar2 <- function(dataobj, X){
    freqs <- foreach(asample=seq_along(dataobj$nList), .combine='c', .inorder=FALSE)%dopar%{
        sum(dataobj$psiList[[asample]] *
            sapply(seq_len(dataobj$nList[asample]),function(cluster){
            prod(sapply(discreteCovs, function(covariate){
                dataobj$phiList[[asample]][[covariate]][X[covariate], cluster]
            })) *
                dmvnorm(X[continuousCovs], mean=dataobj$muList[[asample]][continuousCovs,cluster], sigma=as.matrix(dataobj$sigmaList[[asample]][continuousCovs,continuousCovs,cluster]))
        }))
       #drop(dataobj$phiList[[asample]]$bin_RMSD %*% weights)/sum(weights)
    }
    me <- mean(freqs)
    c(me, sqrt(mean(freqs^2)-me^2))
}
##
unseldata <- setdiff(1:nrow(data), seldata)
##
## 500: 1291.88  961.86 9558.62 
nTest <- 500
testdata <- data[unseldata, c('bin_RMSD',covNames), with=F]
testd <- data.table()
for(val in rmsdVals){
    testd <- rbind(testd, tail(testdata[bin_RMSD==val],n=round(3*nTest*priorP2[val])))
}
testdata <- testd
rm(testd)
##
## user  system elapsed  26.41   18.95  190.12
rm(testres2)
gc()
system.time(testres2 <- t(apply(testdata, 1, function(datum){
    c(datum['bin_RMSD'], normalize(sapply(rmsdVals,function(val){predictYpar2(MCMCdata[[val]],datum)})))
})))
##
evals2 <- metrics(testres2, priorP2)
save.image(file=paste0('_direct_test_N',ndata,'_',length(covNums),'covs.RData'))
evals2
## 6 covs, 5000 pts, 8411 s
## > > >     model delta_gain contig_gain  log_score mean_score
## 1:  model  0.6826667   0.7546667 -0.8769486  0.6062381
## 2: chance  0.6006667   0.6683333 -0.9284347  0.4487041
## 3:    min  0.0000000   0.0000000       -Inf  0.0000000
## 4:    max  1.0000000   1.0000000  0.0000000  1.0000000
##
## 5 covs, 5000 points, 7308 s
##     model delta_gain contig_gain log_score mean_score
## 1:  model  0.6806667   0.7500000 -0.8471781  0.5860935
## 2: chance  0.6006667   0.6683333 -0.9284347  0.4487041
## 3:    min  0.0000000   0.0000000       -Inf  0.0000000
## 4:    max  1.0000000   1.0000000  0.0000000  1.0000000
## > > > >     model delta_gain contig_gain  log_score mean_score
## 1:  model  0.6820000   0.7503333 -0.8500523  0.5854054
## 2: chance  0.6006667   0.6683333 -0.9284347  0.4487041
## 3:    min  0.0000000   0.0000000       -Inf  0.0000000
## 4:    max  1.0000000   1.0000000  0.0000000  1.0000000
##
## 5 covs, 3500 pts, 7348 s
##     model delta_gain contig_gain  log_score mean_score
## 1:  model  0.6700000   0.7480000 -0.8622665  0.5880038
## 2: chance  0.6006667   0.6683333 -0.9284347  0.4487041
## 3:    min  0.0000000   0.0000000       -Inf  0.0000000
## 4:    max  1.0000000   1.0000000  0.0000000  1.0000000
##
## 5 covs, 2500 points, 7032 s
##     model delta_gain contig_gain  log_score mean_score
## 1:  model  0.6733333   0.7450000 -0.8941766  0.5928397
## 2: chance  0.6006667   0.6683333 -0.9284347  0.4487041
## 3:    min  0.0000000   0.0000000       -Inf  0.0000000
## 4:    max  1.0000000   1.0000000  0.0000000  1.0000000
##
## 3500 points, 7590 s
##     model delta_gain contig_gain  log_score mean_score
## 1:  model  0.6653333   0.7406667 -0.8515224  0.5708112
## 2: chance  0.6006667   0.6683333 -0.9284347  0.4487041
## 3:    min  0.0000000   0.0000000       -Inf  0.0000000
## 4:    max  1.0000000   1.0000000  0.0000000  1.0000000
##
## 2500 datapoints, 6435 s
##     model delta_gain contig_gain  log_score mean_score
## 1:  model  0.6653333   0.7426667 -0.8641327  0.5757285
## 2: chance  0.6006667   0.6683333 -0.9284347  0.4487041
## 3:    min  0.0000000   0.0000000       -Inf  0.0000000
## 4:    max  1.0000000   1.0000000  0.0000000  1.0000000

##     model delta_gain contig_gain  log_score mean_score
## 1:  model  0.6606667   0.7363333 -0.9217125  0.5877849
## 2: chance  0.6006667   0.6683333 -0.9284347  0.4487041
## 3:    min  0.0000000   0.0000000       -Inf  0.0000000
## 4:    max  1.0000000   1.0000000  0.0000000  1.0000000

##     model delta_gain contig_gain  log_score mean_score
## 1:  model  0.6520000   0.7356667 -1.0144791  0.5784422
## 2: chance  0.6006667   0.6683333 -0.9284347  0.4487041
## 3:    min  0.0000000   0.0000000       -Inf  0.0000000
## 4:    max  1.0000000   1.0000000  0.0000000  1.0000000





## Plot predictive density
plan(sequential)
plan(multisession, workers = 6L)
xgrid <- seq(0.5, 200,length.out=20)
ygrid <- seq(0.001,0.999,length.out=20)
pdf(file=paste0('testdensities.pdf'),height=11.7,width=16.5*10)
## plot(NA,xlim=c(0,1),ylim=c(0,1),main='parameters')
## text(0,1,(paste0('transf = ','id','\ndim = ',dims,'\nmeanalpha = ',meanalpha,'\nsdalpha = ',sdalpha,'\nK = ',ka,'\nL = ',la,'\nM = ',mu,'\nresc = ',resc,'\nNclust = ',mean(clusters))),adj = c(0,1),cex=5)
##
zgrid <- outer(xgrid,ygrid,function(x,y){
    mapply(function(xx,yy){
                    XX <- c(xx,yy)
                    names(XX) <- continuousCovs
                    YY <- c(sasa2y(xx),tanimoto2y(yy))
                    names(YY) <- continuousCovs
                    predictSasaTani(MCMCdata[[1]], YY)*jac(XX)
                    },
        x,y)}
        )
persp(xgrid,ygrid,zgrid,zlim=c(0,max(zgrid)),ticktype='detailed',theta = 45, phi = 15)
dev.off()


for(i in 1:100){
    zgrid <- outer(xgrid,ygrid,function(x,y){
        XX <- c(x,y)
        names(XX) <- continuousCovs
        predictSasaTani(MCMCdata[[1]], XX)*jac(X)}
        )
persp(xgrid,ygrid,zgrid,zlim=c(0,max(zgrid)),ticktype='detailed')
}
dev.off()




xgrid <- seq(0.5, 200,length.out=20)
ygrid <- seq(0.001,0.999,length.out=20)
pdf(file=paste0('testsampledensities.pdf'),height=11.7,width=16.5*10)
layout(matrix(1:10,1,10,byrow=TRUE))
## plot(NA,xlim=c(0,1),ylim=c(0,1),main='parameters')
## text(0,1,(paste0('transf = ','id','\ndim = ',dims,'\nmeanalpha = ',meanalpha,'\nsdalpha = ',sdalpha,'\nK = ',ka,'\nL = ',la,'\nM = ',mu,'\nresc = ',resc,'\nNclust = ',mean(clusters))),adj = c(0,1),cex=5)
##
for(asample in round(seq(1,1000,length.out=10^2))){
zgrid <- outer(xgrid,ygrid,function(x,y){
    mapply(function(xx,yy){
                    XX <- c(xx,yy)
                    names(XX) <- continuousCovs
                    YY <- c(log(xx),qnorm(yy))
                    names(YY) <- continuousCovs
                    predictSampleSasaTani(MCMCdata[[1]], asample,YY)*jac(XX)
                    },
        x,y)}
        )
persp(xgrid,ygrid,zgrid,zlim=c(0,max(zgrid)),ticktype='detailed')
}
dev.off()

for(i in 1:100){
    zgrid <- outer(xgrid,ygrid,function(x,y){
        XX <- c(x,y)
        names(XX) <- continuousCovs
        predictSasaTani(MCMCdata[[1]], XX)*jac(X)}
        )
persp(xgrid,ygrid,zgrid,zlim=c(0,max(zgrid)),ticktype='detailed')
}
dev.off()






plan(sequential)
xgrid <- seq(0.1, 2,length.out=20)
ygrid <- seq(0.001,0.999,length.out=20)
pdf(file=paste0('testsampledensities.pdf'),height=11.7,width=16.5*10)
## plot(NA,xlim=c(0,1),ylim=c(0,1),main='parameters')
## text(0,1,(paste0('transf = ','id','\ndim = ',dims,'\nmeanalpha = ',meanalpha,'\nsdalpha = ',sdalpha,'\nK = ',ka,'\nL = ',la,'\nM = ',mu,'\nresc = ',resc,'\nNclust = ',mean(clusters))),adj = c(0,1),cex=5)
##
zgrid <- outer(xgrid,ygrid,function(x,y){
    mapply(function(xx,yy){
                    XX <- c(xx,yy)
                    names(XX) <- continuousCovs
                    YY <- c(log(xx),qnorm(yy))
                    names(YY) <- continuousCovs
                    predictSampleSasaTani(MCMCdata[[1]], 3,YY)*jac(XX)
                    },
        x,y)}
        )
persp(xgrid,ygrid,zgrid,zlim=c(0,max(zgrid)),ticktype='detailed')
dev.off()






###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################
### OLD ###


## this is much slower
predictYtest2 <- function(dataobj, testset){
    t(apply(testset[,c(discreteCovs,continuousCovs),with=FALSE], 1, function(X){
        freqs <- foreach(asample=seq_along(dataobj$nList), .combine=c, .inorder=FALSE)%dopar%{
        sum(dataobj$psiList[[asample]] *
            sapply(seq_len(dataobj$nList[asample]),function(cluster){
            prod(sapply(discreteCovs, function(covariate){
                dataobj$phiList[[asample]][[covariate]][X[covariate], cluster]
            })) *
                dmvnorm(X[continuousCovs], mean=dataobj$muList[[asample]][continuousCovs,cluster], sigma=as.matrix(dataobj$sigmaList[[asample]][continuousCovs,continuousCovs,cluster]))
        }))
       #drop(dataobj$phiList[[asample]]$bin_RMSD %*% weights)/sum(weights)
    }
    me <- mean(freqs)
    c(me, sqrt(mean(freqs^2) - me^2))
}))}


predictYpar <- function(dataobj, X){
    foreach(asample=seq_along(dataobj$nList), .combine='+', .inorder=FALSE)%dopar%{
        sum(dataobj$psiList[[asample]] *
            sapply(seq_len(dataobj$nList[asample]),function(cluster){
            prod(sapply(discreteCovs, function(covariate){
                dataobj$phiList[[asample]][[covariate]][X[covariate], cluster]
            })) *
                dmvnorm(X[continuousCovs], mean=dataobj$muList[[asample]][continuousCovs,cluster], sigma=as.matrix(dataobj$sigmaList[[asample]][continuousCovs,continuousCovs,cluster]))
        }))
       #drop(dataobj$phiList[[asample]]$bin_RMSD %*% weights)/sum(weights)
}/length(dataobj$nList)
}







##
unseldata <- setdiff(1:nrow(data), seldata)
##
## 500, 6 threads: 1269.53  954.36 9520.95  
nTest <- 500
testdata <- data[unseldata, c('bin_RMSD',covNames), with=F]
testd <- data.table()
for(val in rmsdVals){
    testd <- rbind(testd, tail(testdata[bin_RMSD==val],n=nTest))
}
testdata <- testd
rm(testd)
##
rm(testres)
gc()
system.time(testres <- t(apply(testdata, 1, function(datum){
    c(datum['bin_RMSD'], sapply(rmsdVals,function(val){predictYpar(MCMCdata[[val]],datum)}))
})))
save.image(file='_reverse_test.RData')
plan(sequential)

##
evals1 <- metrics(testres, priorP)
evals1
save.image(file='_reverse_test.RData')
##     model delta_gain contig_gain log_score mean_score
## 1:  model  0.4586667   0.6756667 -1.631715  0.4247677
## 2: chance  0.3333333   0.6666667 -1.098612  0.3333333
## 3:    min  0.0000000   0.0000000      -Inf  0.0000000
## 4:    max  1.0000000   1.0000000  0.000000  1.0000000




#### Evaluation 2
##
##discrMin <- sapply(data[,discreteCovs,with=F], min)-1
discreteCovs <- discreteCovs
continuousCovs <- continuousCovs
plan(sequential)
plan(multisession, workers = 6L)
priorP2 <- normalize(as.vector(table(data$bin_RMSD)))
##
predictYpar2 <- function(dataobj, X){
    foreach(asample=seq_along(dataobj$nList), .combine='+', .inorder=FALSE)%dopar%{
        sum(dataobj$psiList[[asample]] *
            sapply(seq_len(dataobj$nList[asample]),function(cluster){
            prod(sapply(discreteCovs, function(covariate){
                dataobj$phiList[[asample]][[covariate]][X[covariate], cluster]
            })) *
                dmvnorm(X[continuousCovs], mean=dataobj$muList[[asample]][continuousCovs,cluster], sigma=as.matrix(dataobj$sigmaList[[asample]][continuousCovs,continuousCovs,cluster]))
        }))
       #drop(dataobj$phiList[[asample]]$bin_RMSD %*% weights)/sum(weights)
}/length(dataobj$nList)
}
##
unseldata <- setdiff(1:nrow(data), seldata)
##
## 500: 1291.88  961.86 9558.62 
nTest <- 500
testdata <- data[unseldata, c('bin_RMSD',covNames), with=F]
testd <- data.table()
for(val in rmsdVals){
    testd <- rbind(testd, tail(testdata[bin_RMSD==val],n=round(3*nTest*priorP2[val])))
}
testdata <- testd
rm(testd)
##
## user  system elapsed  26.41   18.95  190.12
rm(testres2)
gc()
system.time(testres2 <- t(apply(testdata, 1, function(datum){
    c(datum['bin_RMSD'], normalize(sapply(rmsdVals,function(val){predictYpar2(MCMCdata[[val]],datum)})))
})))
save.image(file='_reverse_test.RData')
##
evals2 <- metrics(testres2, priorP2)
evals2
save.image(file='_reverse_test.RData')
##     model delta_gain contig_gain  log_score mean_score
## 1:  model  0.6520000   0.7356667 -1.0144791  0.5784422
## 2: chance  0.6006667   0.6683333 -0.9284347  0.4487041
## 3:    min  0.0000000   0.0000000       -Inf  0.0000000
## 4:    max  1.0000000   1.0000000  0.0000000  1.0000000














#####################################################################
#####################################################################
#####################################################################
#####################################################################
## Old pieces

#### Evaluation - slower
##
##discrMin <- sapply(data[,discreteCovs,with=F], min)-1
discreteCovs <- discreteCovs
continuousCovs <- continuousCovs
plan(sequential)
plan(multisession, workers = 6L)
priorP <- rep(1,3)
##
predictYpar <- function(dataobj, x){
normalize(foreach(val=rmsdVals, .combine=c)%:%foreach(asample=seq_along(dataobj[[val]]$nList), .combine='+', .inorder=FALSE)%dopar%{
        sum(exp(log(dataobj[[val]]$psiList[[asample]]) + sapply(seq_len(dataobj[[val]]$nList[asample]),function(j){
            sum(log(sapply(discreteCovs, function(elem){
                dataobj[[val]]$phiList[[asample]][[elem]][x[elem], j]
            }))) +
                dmvnorm(x[continuousCovs], mean=dataobj[[val]]$muList[[asample]][continuousCovs,j], sigma=as.matrix(dataobj[[val]]$sigmaList[[asample]][continuousCovs,continuousCovs,j]), log=TRUE)
        }))) * priorP[val]
       #drop(dataobj[[val]]$phiList[[asample]]$bin_RMSD %*% weights)/sum(weights)
}/length(dataobj[[val]]$nList)
)}
##
unseldata <- setdiff(1:nrow(data), seldata)
##
nTest <- 10
testdata <- data[unseldata, c('bin_RMSD',covNames), with=F]
testd <- data.table()
for(val in rmsdVals){
    testd <- rbind(testd, tail(testdata[bin_RMSD==val],n=nTest))
}
testdata <- testd
rm(testd)
gc()
##
## user  system elapsed  24.22   20.91  225.42 
system.time(testres <- t(apply(testdata, 1, function(datum){
    c(datum['bin_RMSD'], (predictYpar(MCMCdata,datum)))
})))
save.image(file='_reverse_test.RData')
##
mean(abs(testres[,1]==apply(testres[,-1],1,which.max)))
1/3
##
-mean(abs(testres[,1]-apply(testres[,-1],1,which.max)))
-mean(c(c(0,1,2), c(1,0,1), c(2,1,0)))
##
mean(log(diag(testres[,testres[,1]+1])))
log(1/3)
##
mean((diag(testres[,testres[,1]+1])))
(1/3)
plan(sequential)
## [1] 0.4313333
## > [1] 0.3333333
## > > [1] -0.7753333
## > [1] -0.8888889
## > > [1] -1.233647
## > [1] -1.098612
## [1] 0.3846496
## > [1] 0.3333333




##########
########## old stuff
##########



















predictYpar2 <- function(mcobj, x){
    foreach(asample=1:length(nList), .combine=rbind, .inorder=FALSE)%dopar%{
        discreteCovs <- mcobj$discreteCovs[mcobj$discreteCovs != 'bin_RMSD']
        continuousCovs <- mcobj$continuousCovs
        weights <- exp(log(psiList[[asample]]) + sapply(seq_len(nList[asample]),function(j){
            sum(log(mapply(function(xx,yy){xx[yy, j]}, phiList[[asample]][discreteCovs], x[discreteCovs]-discrMin[discreteCovs]))) +
                dmvnorm(x[continuousCovs], mean=muList[[asample]][continuousCovs,j], sigma=as.matrix(sigmaList[[asample]][continuousCovs,continuousCovs,j]), log=TRUE)
        }))
       drop(phiList[[asample]]$bin_RMSD %*% weights)/sum(weights)
}
}
##
unseldata <- setdiff(1:nrow(data), seldata)
##
nTest <- 100
testdata <- data[unseldata,testmc$covNames, with=F]
testdata <- rbind(head(testdata[bin_RMSD==0],n=nTest), head(testdata[bin_RMSD==1],n=nTest), head(testdata[bin_RMSD==2],n=nTest))
##
## user  system elapsed    1.71    0.14    4.31 
system.time(testres <- t(sapply(seq_len(nrow(testdata)), function(datum){
    datum <- unlist(testdata[datum])
    c(datum['bin_RMSD']-discrMin['bin_RMSD'], colMeans(predictYpar2(testmc,datum)))
})))
##
sum(sapply(1:nrow(testres), function(i){testres[i,1]==which.max(testres[i,-1])}))/nrow(testres)



##   user  system elapsed   1.73    0.13    4.73 
system.time(testres2 <- apply(data[unseldata[1:10], testmc$covNames, with=FALSE], 1, function(x){
    c(x[1], colMeans(predictYpar(x[-1])))
}))

## user  system elapsed    7.37    0.05    7.43 
system.time(testres2 <- foreach(i=1:10, .combine=c)%do%{
    testyx <-  unlist(data[unseldata[i], testmc$covNames, with=FALSE])
    colMeans(predictY2(testyx[-1]))[testyx[1]+1]
})




##################################################
## Continuous-x, no y-model
## Transform int variates to double
set.seed <- 147
ndata <- round(nSamples/5)
covNums <- c(1,3:8)
seldata <- sort(sample(nrow(data), ndata, replace=FALSE))
rm(datamc)
datamc <- as.data.frame(data[seldata, ..covNums])
discreteCovs <- names(datamc)[which(sapply(datamc, is.integer)==TRUE)]
for(i in discreteCovs[-1]){
    datamc[i] <- datamc[i] + runif(nrow(datamc)) - 0.5
}
discreteCovs <- names(datamc)[which(sapply(datamc, is.integer)==TRUE)]
continuousCovs <- names(datamc)[which(sapply(datamc, is.double)==TRUE)]
##
#
system.time(testmc <- profRegr(excludeY=TRUE, xModel='Mixed', nSweeps=10000, nBurn=2000, nFilter=10, data=datamc, nClusInit=10, covNames=c(discreteCovs,continuousCovs), discreteCovs=discreteCovs, continuousCovs=continuousCovs, nProgress=100, seed=222, output='_mc2output'))
##
fc <- file("_mc2output_logPost.txt")
logPost <- sapply(strsplit(readLines(fc), " +"), as.numeric)
close(fc)
##
pdff('cmsum_cont')
print(globalParsTrace(testmc, parameters = "nClusters",plotBurnIn=F))
print(globalParsTrace(testmc, parameters = "alpha",plotBurnIn=F))
for(i in 1:3){matplot(logPost[i,],type='l')}
dev.off()
##
fc <- file("_mc2output_nClusters.txt")
nList <- sapply(strsplit(readLines(fc), " +"), as.numeric)
close(fc)
#
fc <- file("_mc2output_psi.txt")
psiList <- lapply(strsplit(readLines(fc), " +"), function(x){y <- as.numeric(x);y[y>=0]})
close(fc)
##                                        #
fc <- file("_mc2output_phi.txt")
phiList <- lapply(strsplit(readLines(fc), " +"), function(x){y <- as.numeric(x);y[y>=0]})
close(fc)
##
fc <- file("_mc2output_mu.txt")
muList <- lapply(strsplit(readLines(fc), " +"), as.numeric)
close(fc)
##
fc <- file("_mc2output_Sigma.txt")
sigmaList <- lapply(strsplit(readLines(fc), " +"), as.numeric)
close(fc)
##
nCatY <- 3
nContCov <- testmc$nContinuousCovs
for(i in 1:length(nList)){
    dim(phiList[[i]]) <- c(nList[i],nCatY)
    phiList[[i]] <- aperm(phiList[[i]])
    dim(sigmaList[[i]]) <- c(nList[i], nContCov, nContCov)
    sigmaList[[i]] <- aperm(sigmaList[[i]])
    dim(muList[[i]]) <- c(nList[i], nContCov)
    muList[[i]] <- aperm(muList[[i]])
}

registerDoFuture()  ## tell foreach to use the future framework
plan(multisession, workers = 4L)

plan(sequential)



predictYpar <- function(x){
    foreach(i=1:length(nList), .combine=rbind, .inorder=FALSE)%dopar%{
        weights <- psiList[[i]] * sapply(1:nList[i],function(j){dmvnorm(x, mean=muList[[i]][,j], sigma=sigmaList[[i]][,,j])})
        drop(phiList[[i]] %*% weights)/sum(weights)
}
}

predictY <- function(x){
    foreach(i=1:length(nList), .combine=rbind)%do%{
        weights <- psiList[[i]] * sapply(1:nList[i],function(j){dmvnorm(x, mean=muList[[i]][,j], sigma=sigmaList[[i]][,,j])})
        drop(phiList[[i]] %*% weights)/sum(weights)
}
}

predictY2 <- function(x){
    t(sapply(1:length(nList), function(i){
        weights <- psiList[[i]] * sapply(1:nList[i],function(j){dmvnorm(x, mean=muList[[i]][,j], sigma=sigmaList[[i]][,,j])})
        drop(phiList[[i]] %*% weights)/sum(weights)
}))
}

unseldata <- setdiff(1:nrow(data), seldata)

## user  system elapsed    1.71    0.14    4.31 
system.time(testres <- t(sapply(1:100, function(i){
    testyx <-  unlist(data[unseldata[i], testmc$covNames, with=FALSE])
    c(testyx[1], colMeans(predictYpar(testyx[-1])))
})))

#   user  system elapsed   1.73    0.13    4.73 
system.time(testres2 <- apply(data[unseldata[1:10], testmc$covNames, with=FALSE], 1, function(x){
    c(x[1], colMeans(predictYpar(x[-1])))
}))

sum(sapply(1:nrow(testres), function(i){testres[i,1]+1==which.max(testres[i,-1])}))/nrow(testres)

## user  system elapsed    7.37    0.05    7.43 
system.time(testres2 <- foreach(i=1:10, .combine=c)%do%{
    testyx <-  unlist(data[unseldata[i], testmc$covNames, with=FALSE])
    colMeans(predictY2(testyx[-1]))[testyx[1]+1]
})






## Continuous-x, y-model
## Transform int variates to double
set.seed <- 147
ndata <- round(nSamples/5)
covNums <- c(1,3:8)
seldata <- sort(sample(nrow(data), ndata, replace=FALSE))
rm(datamc)
datamc <- as.data.frame(data[seldata, ..covNums])
for(i in discreteCovs[-1]){
    datamc[i] <- datamc[i] + runif(nrow(datamc)) - 0.5
}
yCov <- 'bin_RMSD'
discreteCovs <- names(datamc)[which(sapply(datamc, is.integer)==TRUE)]
continuousCovs <- names(datamc)[which(sapply(datamc, is.double)==TRUE)]
##
system.time(testmc <- profRegr(excludeY=FALSE, yModel='Categorical', xModel='Normal', nSweeps=2000, nBurn=10, nFilter=2, data=datamc, nClusInit=2, outcome=yCov, covNames=c(discreteCovs,continuousCovs), discreteCovs=discreteCovs, continuousCovs=continuousCovs, nProgress=100, seed=222, output='_mc2output'))
##
pdff('cmsum_cont_y')
print(globalParsTrace(testmc, parameters = "nClusters",plotBurnIn=F))
print(globalParsTrace(testmc, parameters = "alpha",plotBurnIn=F))
dev.off()

fc <- file("_mccoutput_phi.txt")
phiList <- lapply(strsplit(readLines(fc), " +"), function(x){y <- as.numeric(x);y[y>=0]})
close(fc)
#
fc <- file("_mccoutput_nClusters.txt")
nList <- sapply(strsplit(readLines(fc), " +"), as.numeric)
close(fc)
#
fc <- file("_mccoutput_psi.txt")
psiList <- lapply(strsplit(readLines(fc), " +"), function(x){y <- as.numeric(x);y[y>=0]})
close(fc)


## discreteCovs <- names(datamc[,which(names(datamc)!='bin_RMSD' & sapply(datamc, is.integer)==TRUE)])
## continuousCovs <- names(datamc[,which(names(datamc)!='bin_RMSD' & sapply(datamc, is.double)==TRUE)])
## #
## system.time(testmc <- profRegr(excludeY=FALSE, yModel='Categorical', xModel='Mixed', nSweeps=500, nBurn=100, nFilter=5, data=datamc, nClusInit=10, outcome='bin_RMSD', covNames=c(discreteCovs,continuousCovs), discreteCovs=discreteCovs, continuousCovs=continuousCovs, nProgress=100, seed=222, output='_mc2output'))
## #
## pdff('cmsum_y_discr')
## print(globalParsTrace(testmc, parameters = "nClusters",plotBurnIn=F))
## print(globalParsTrace(testmc, parameters = "alpha",plotBurnIn=F))
## dev.off()


### Updated up to above ###









misubdata <- mipairs[which(mipairs[,2]==(nFeatures+1)),]
misubdata <- misubdata[order(misubdata[,7], decreasing=FALSE),]
rmsdData <- data[[rmsdCol]]
pdff('jointplotslog_RMSD_data')
for(k in 1:nrow(misubdata)){
    i <- misubdata[k,1]
if(i!=rmsdCol){    datum <- data[[i]]
    matplot(x=datum, y=rmsdData, type='p', pch='.', col=paste0('#000000','88'),
            ylab='log-rmsd', xlab=nameFeatures[i]
                                        #, log='y'
            )
matlines(x=matrix(fivenum(datum)[c(1,5)],2,1), y=matrix(rep(logRmsdThreshold,2),2,3,byrow=T),
         lty=c(2,3,2), lwd=3, col=c(myredpurple,myyellow,myredpurple))
    title(paste0(nameFeatures[i],
                 ', conditional entropy = ',signif(misubdata[k,7],4),
                 ' bit  (from ',signif(misubdata[k,6],4),' bit)'))
    ##
    ## summa <- fivenum(datum)
    ## width <- diff(summa[c(2,4)])/nbinsq
    ## nbins <- round(diff(summa[c(1,5)])/width)
    ## if(is.integer(datum)){nbins <- min(nbins, diff(summa[c(1,5)])+1)}
    ## freqs <- bin2(x=cbind(log10(rmsdData), datum),
    ##               ab=rbind(summa[1,5]*c(1,1.01),
    ##                        ))
}}
dev.off()

library('rmarkdown')

##infoFeatures <- misubdata[2:20,1] 
bestfeatures <- c(33,23,34,27,28,30,26,29,25,24,22,21,13,5,20,12,32)
##
feats <- bestfeatures[c(2,1,3)]
set.seed(149)
subsamplesize <- 2^9
samplescale <- max(foreach(i=(-1):1, .combine=c)%do%{length(which(data$bin_RMSD==i))})/subsamplesize

plotpoints <- foreach(i=(-1):1)%do%{
    matr <- as.matrix(data[which(data$bin_RMSD==i),..feats])
    matr[sample(1:nrow(matr),round(nrow(matr)/samplescale)),]
}
names(plotpoints) <- c('yes','maybe','no')

##
pointsize <- 0.4
outputfile <- paste0('feats3.htm')
rmarkdown::render('3dgenerator.Rmd', output_file=outputfile,
                  params=list(plotpoints=plotpoints,
                              subsamplesize=1000,
                              pointsize=pointsize##, theta=NA, phi=NA
                              ))


library('rgl')

plot3d(x=matrix(0,4,3),type='p',size=0.1,
                   xlab=featurenames[1],ylab=featurenames[2],zlab=featurenames[3],
       xlim=xlim,ylim=ylim,zlim=zlim,col='white')



rmsdData <- data[[totFeatures+1]]
pdff('jointplotsbin_RMSD_fabiodata')
for(k in 1:nrow(misubdata)){
    i <- misubdata[k,1]
    datum <- data[[i]]
    matplot(x=datum, y=rmsdData, type='p', pch='.', col=paste0('#000000','88'),
            ylab='binned rmsd', xlab=nameFeatures[i]
                                        #, log='y'
            )
    title(paste0(nameFeatures[i],
                 ', MI = ',signif(misubdata[k,3],4),
                 ' bit, normMI = ',signif(misubdata[k,4],4)))
    ##
    ## summa <- fivenum(datum)
    ## width <- diff(summa[c(2,4)])/nbinsq
    ## nbins <- round(diff(summa[c(1,5)])/width)
    ## if(is.integer(datum)){nbins <- min(nbins, diff(summa[c(1,5)])+1)}
    ## freqs <- bin2(x=cbind(log10(rmsdData), datum),
    ##               ab=rbind(summa[1,5]*c(1,1.01),
    ##                        ))
}
dev.off()












##
## pdff('jointplotsRMSD_fabiodata2')
## rmsdData <- unlist(data[[rmsdCol]])
## for(i in okFeatures[-rmsdCol]){
##     datum <- unlist(data[[i]])
##     matplot(x=rmsdData, y=datum, type='p', pch=20, cex=0.5, col=paste0('#000000','22'),
##             xlab='rmsd', ylab=nameFeatures[i])
## }
## dev.off()


data2 <- cbind(data.table(logrmsd=log(data$rmsd)),
               data[,5:ncol(data)])
write.table(data2,'dp_data.csv',sep=',',row.names=FALSE, col.names=TRUE)



preci <- sapply(1:1000,function(sample){
a <-         sum(exp(
            log(dataobj$psiList[[sample]]) +
            sapply(seq_len(dataobj$nList[sample]),function(cluster){
            sum(log(sapply(discreteCovs, function(covariate){
                dataobj$phiList[[sample]][[covariate]][X[covariate], cluster]
            }))) +
                dmvnorm(X[continuousCovs], mean=dataobj$muList[[sample]][continuousCovs,cluster], sigma=as.matrix(dataobj$sigmaList[[sample]][continuousCovs,continuousCovs,cluster]), log=TRUE)
        })))
b <-         sum(
            (dataobj$psiList[[sample]]) *
            sapply(seq_len(dataobj$nList[sample]),function(cluster){
            prod((sapply(discreteCovs, function(covariate){
                dataobj$phiList[[sample]][[covariate]][X[covariate], cluster]
            }))) *
                dmvnorm(X[continuousCovs], mean=dataobj$muList[[sample]][continuousCovs,cluster], sigma=as.matrix(dataobj$sigmaList[[sample]][continuousCovs,continuousCovs,cluster]))
            }))
abs((a-b)/(a+b))})
summary(preci)


rm(preci)
gc()
system.time(preci <- sapply(rep(1:1000,20),function(sample){
 sum(exp(
            log(dataobj$psiList[[sample]]) +
            sapply(seq_len(dataobj$nList[sample]),function(cluster){
            sum(log(sapply(discreteCovs, function(covariate){
                dataobj$phiList[[sample]][[covariate]][X[covariate], cluster]
            }))) +
                dmvnorm(X[continuousCovs], mean=dataobj$muList[[sample]][continuousCovs,cluster], sigma=as.matrix(dataobj$sigmaList[[sample]][continuousCovs,continuousCovs,cluster]), log=TRUE)
            })))
 }))


rm(preci)
gc()
system.time(preci <- sapply(rep(1:1000,20),function(sample){
sum(
            (dataobj$psiList[[sample]]) *
            sapply(seq_len(dataobj$nList[sample]),function(cluster){
            prod((sapply(discreteCovs, function(covariate){
                dataobj$phiList[[sample]][[covariate]][X[covariate], cluster]
            }))) *
                dmvnorm(X[continuousCovs], mean=dataobj$muList[[sample]][continuousCovs,cluster], sigma=as.matrix(dataobj$sigmaList[[sample]][continuousCovs,continuousCovs,cluster]))
            }))
}))

abs((a-b)/(a+b))})
summary(preci)



##
predictYtest <- function(dataobj, testset){
    freqs <- foreach(X=t(testset[,c(discreteCovs,continuousCovs),with=FALSE]), .combine=rbind)%:%foreach(sample=seq_along(dataobj$nList), .combine=c, .inorder=FALSE)%dopar%{
        sum(dataobj$psiList[[sample]] *
            sapply(seq_len(dataobj$nList[sample]),function(cluster){
            prod(sapply(discreteCovs, function(covariate){
                dataobj$phiList[[sample]][[covariate]][X[covariate,], cluster]
            })) *
                dmvnorm(X[continuousCovs,], mean=dataobj$muList[[sample]][continuousCovs,cluster], sigma=as.matrix(dataobj$sigmaList[[sample]][continuousCovs,continuousCovs,cluster]))
        }))
       #drop(dataobj$phiList[[sample]]$bin_RMSD %*% weights)/sum(weights)
    }
    me <- rowMeans(freqs)
    cbind(means=me, stds=sqrt(rowMeans(freqs^2) - me^2))
}


##
predictYXtesto <- function(dataobj, testset, rvals){
    testset <- testset[, c(discreteCovs,continuousCovs), with=FALSE]
    freqs <- foreach(X=t(testset[,c(discreteCovs,continuousCovs),with=FALSE]), .combine=rbind)%:%foreach(val=rvals, .combine=rbind)%:%foreach(sample=seq_along(dataobj[[val]]$nList), .combine=rbind, .inorder=FALSE)%dopar%{
     colSums(exp(
        log(dataobj[[val]]$psiList[[sample]]) +
        t(vapply(seq_len(dataobj[[val]]$nList[sample]), function(cluster){
            rowSums(log(
                vapply(discreteCovs, function(covariate){
                    dataobj[[val]]$phiList[[sample]][[covariate]][X[,covariate], cluster]
                }, numeric(3))
            )) +
                dmvnorm(X[,continuousCovs], mean=dataobj[[val]]$muList[[sample]][continuousCovs,cluster], sigma=as.matrix(dataobj[[val]]$sigmaList[[sample]][continuousCovs,continuousCovs,cluster]), log=TRUE)
        }, numeric(nrow(X))))
    ))
        ## colSums(dataobj[[val]]$psiList[[sample]] *
        ##     t(sapply(seq_len(dataobj[[val]]$nList[sample]),function(cluster){
        ##     exp(rowSums(log(sapply(discreteCovs, function(covariate){
        ##         dataobj[[val]]$phiList[[sample]][[covariate]][X[,covariate], cluster]
        ##     })))) *
        ##         dmvnorm(X[,continuousCovs], mean=dataobj[[val]]$muList[[sample]][continuousCovs,cluster], sigma=as.matrix(dataobj[[val]]$sigmaList[[sample]][continuousCovs,continuousCovs,cluster]))
        ##     }))
        ##     )
    }
    me <- rowMeans(freqs)
    freqs <- cbind(means=me, stds=sqrt(rowMeans(freqs^2) - me^2))
    dim(freqs) <- c(length(rvals),nrow(testset),2)
    freqs
}
    ##



##    
predictYXtestsingle <- function(dataobj, X){
    X <- as.matrix(X[, c(discreteCovs,continuousCovs), with=FALSE])
    freqs <- foreach(sample=seq_along(dataobj$nList), .combine=cbind, .inorder=FALSE)%dopar%{
        colSums(exp(
        log(dataobj$psiList[[sample]]) +
        t(vapply(seq_len(dataobj$nList[sample]), function(cluster){
            rowSums(log(
                vapply(discreteCovs, function(covariate){
                    dataobj$phiList[[sample]][[covariate]][X[,covariate], cluster]
                }, numeric(nrow(X)))
            )) +
                dmvnorm(X[,continuousCovs], mean=dataobj$muList[[sample]][continuousCovs,cluster], sigma=as.matrix(dataobj$sigmaList[[sample]][continuousCovs,continuousCovs,cluster]), log=TRUE)
        }, numeric(nrow(X))))
    ))
        ## colSums(dataobj$psiList[[sample]] *
        ##     t(sapply(seq_len(dataobj$nList[sample]),function(cluster){
        ##     exp(rowSums(log(sapply(discreteCovs, function(covariate){
        ##         dataobj$phiList[[sample]][[covariate]][X[,covariate], cluster]
        ##     })))) *
        ##         dmvnorm(X[,continuousCovs], mean=dataobj$muList[[sample]][continuousCovs,cluster], sigma=as.matrix(dataobj$sigmaList[[sample]][continuousCovs,continuousCovs,cluster]))
        ##     }))
        ##     )
    }
    me <- rowMeans(freqs)
    cbind(means=me, stds=sqrt(rowMeans(freqs^2) - me^2))
}


testmax <- function(x,y,z){max(x,y,z)}
testmax <- Vectorize(testmax)


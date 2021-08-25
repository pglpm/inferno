## Author: PGL  Porta Mana
## Created: 2021-03-20T10:07:17+0100
## Last-Updated: 2021-08-24T07:48:18+0200
################
## Script for evaluation of regression:
## Unfactorizable prior
## r|x exchangeable
## x not exchangeable
## r continuous
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

##################################################################
##################################################################
## Function to calculate \sum_t F^{(t)}_{r|x}
##
boundaries <- log(c(2,3))
rangeR <- fivenum(data$log_RMSD, na.rm=T)
rangeR <- c(diff(rangeR[c(3,1)]), diff(rangeR[c(3,5)]))*1.1+rangeR[3]
rgrid <- seq(rangeR[1],rangeR[2],length.out=101)
lrgrid <- length(rgrid)
cC <- setdiff(continuousCovs,'log_RMSD')
predictYX <- function(dataobj, X){
    X <- as.matrix(X[, c(discreteCovs,cC), with=FALSE])
    ntestdata <- nrow(X)
    ##
    XX <- rep(t(X),lrgrid)
    dim(XX) <- c(ncol(X), nrow(X)*lrgrid)
    XX <- cbind( rep(rgrid, each=ntestdata), t(XX))
    colnames(XX) <- c('log_RMSD',colnames(X))
    ##
    freqs <- foreach(asample=seq_along(dataobj$nList), .combine=cbind, .inorder=FALSE)%dopar%{
        numer <- colSums(
            exp(
            log(dataobj$psiList[[asample]]) +
            t(vapply(seq_len(dataobj$nList[asample]), function(cluster){
                rowSums(log(
                    vapply(discreteCovs, function(covariate){
                        dataobj$phiList[[asample]][[covariate]][XX[,covariate], cluster]
                    }, numeric(ntestdata*lrgrid))
                )) +
                    dmvnorm(XX[,continuousCovs], mean=dataobj$muList[[asample]][continuousCovs,cluster], sigma=as.matrix(dataobj$sigmaList[[asample]][continuousCovs,continuousCovs,cluster]), log=TRUE)
            }, numeric(ntestdata*lrgrid)))
            )
            )

        denom <- colSums(
            exp(
            log(dataobj$psiList[[asample]]) +
            t(vapply(seq_len(dataobj$nList[asample]), function(cluster){
                ## rowSums(log(
                ##     vapply(discreteCovs, function(covariate){
                ##         dataobj$phiList[[asample]][[covariate]][XX[,covariate], cluster]
                ##     }, numeric(ntestdata))
                ## )) +
                    dnorm(XX[,'log_RMSD'], mean=dataobj$muList[[asample]]['log_RMSD',cluster], sd=sqrt(dataobj$sigmaList[[asample]]['log_RMSD','log_RMSD',cluster]), log=TRUE)
            }, numeric(ntestdata*lrgrid)))
            )
            )

        numer/denom
    }
    ## me <- rowMeans(freqs)
    ## freqs <- cbind(means=me, stds=sqrt(rowMeans(freqs^2) - me^2))
    ## dim(freqs) <- c(ntestdata, lrgrid, 2)
    ## dimnames(freqs) <- list(NULL, NULL, c('means','stds'))
    ## freqs
    freqs <- rowMeans(freqs)
    dim(freqs) <- c(ntestdata, lrgrid)
    freqs
}


## Utility matrices
dgain <- diag(1,3)
cgain <- 1-sapply(1:3,function(x){abs(x-1:3)})/2
## function for randomly choosing ties
resample <- function(x, ...) x[sample.int(length(x), ...)]
## function to compute utilities in test set
utilities <- function(truevalues,probX,utilitym){
    if(is.null(dim(probX))){probX <- t(matrix(probX,nrow=length(probX),ncol=length(truevalues)))}
    y <- apply(t(tcrossprod(utilitym, normalizem(probX))), 1,
               function(z){resample(which(z==max(z)))})
    diag(utilitym[truevalues,y])
}
## function to compute various prob-averages in test set
probscores <- function(truevalues,probX,meanfunction){
    if(is.null(dim(probX))){probX <- t(matrix(probX,nrow=length(probX),ncol=length(truevalues)))}
    y <- normalizem(probX)
    meanfunction(diag(y[,truevalues]))
}
## function to construct data table with results
metrics <- function(truevalues, probX, chanceprior){
    list(
        delta_gain = data.table(
                 model=utilities(truevalues,probX,dgain),
                 chance=utilities(truevalues,chanceprior,dgain),
                 range=range(dgain)
             ),
         ##
        contig_gain = data.table(
                 model=utilities(truevalues,probX,cgain),
                 chance=utilities(truevalues,chanceprior,cgain),
                 range=range(cgain)
             ),
         ##
        log_score = data.table(
                 model=probscores(truevalues,probX,log),
                 chance=probscores(truevalues,chanceprior,log),
                 range=c(-Inf,0)
             ),
         ##
        mean_score = data.table(
                 model=probscores(truevalues,probX,identity),
                 chance=probscores(truevalues,chanceprior,identity),
                 range=c(0,1)
             )
        )}
##
##
#### Calculation of utilities and scores for test set
## source(file='calibration_plots.R')
##

starttime <- Sys.time()
case <- 1 # 'fdata'
## Frequencies of RMSD categories in test set are same as full dataset.
## It does not make sense to change them, because the conditional frequencies r|x
## then also change - we can then make the model have a score as bad as we please
##
indicesRMSD <- list(rgrid < boundaries[1],
                    rgrid >= boundaries[1] & rgrid <= boundaries[2],
                    rgrid > boundaries[2])
##
scores <- list()
basepriors <- rbind(asdata=normalize(as.vector(table(data$bin_RMSD))), uniform=rep(1,3)/3)
for(dset in c('tail')){
    ## Compare with a model based only on prior frequencies of r
    unseldata <- setdiff(1:nrow(data), seldata)
    nTest <- 2000
    ## construct test set
    testdata <- data[do.call(dset,list(x=unseldata,n=3*nTest)), c('bin_RMSD',covNames), with=F]
    ## testd <- data[unseldata, covNames, with=F]
    ## testdata <- data.table()
    ## for(val in rmsdVals){
    ##     testdata <- rbind(testdata, tail(testd[bin_RMSD==val],n=round(3*nTest*priorR[val])))
    ## }
    ## rm(testd)
    ##
gc()
plan(sequential)
plan(multisession, workers = 6L)
precondfreqs <- t(predictYX(MCMCdata[[1]], testdata))
plan(sequential)
precondfreqs2 <- sapply(1:length(indicesRMSD), function(x){
                colSums(precondfreqs[indicesRMSD[[x]],])})

load(file='_freqs_for_tests_cont_reverse.RData')
    for(ii in 1:nrow(basepriors)){
        priorR <- basepriors[ii,]
        condfreqs <- normalizem(t(t(precondfreqs2) * priorR))
    ##
    set.seed(247)
    scores[[dset]][[rownames(basepriors)[ii]]] <- metrics(truevalues=testdata[,bin_RMSD], probX=condfreqs, chanceprior=priorR)
    print(dset)
    pdff(paste0('scores_cont_reverse_unfact_P',rownames(basepriors)[ii],'_',dset))
    basescores <- scores[[dset]][[ii]]
    ##
print(    ggplot(dt <- melt(basescores$delta_gain[,1:2])) +
        geom_histogram(aes(x=value, fill=variable), color='white', bins=length(unique(c(dgain))), alpha=0.5, position='identity') + scale_fill_bright() + ylim(c(0,nTest*3))+
        xlab('delta utility') +
        geom_vline(data=dt[,.(value=mean(value)), by=variable], aes(xintercept=value, color=variable, linetype=variable), alpha=0.5, size=2))
    ##
print(    ggplot(dt <- melt(basescores$contig_gain[,1:2])) +
        geom_histogram(aes(x=value, fill=variable), color='white', bins=length(unique(c(cgain))), alpha=0.5, position='identity') + scale_fill_bright() + ylim(c(0,nTest*3))+
        xlab('contiguous utility') +
        geom_vline(data=dt[,.(value=mean(value)), by=variable], aes(xintercept=value, color=variable, linetype=variable), alpha=0.5, size=2))
    ##
print(    ggplot(dt <- melt(basescores$log_score[,1:2])) +
        geom_histogram(aes(x=value, y=..density.., fill=variable), color='white', bins=10, alpha=0.5, position='identity') + scale_fill_bright() + ylim(c(0,2.5))+
        xlab('log-probability') + xlim(c(-4,0)) +
        geom_vline(data=dt[,.(value=mean(value)), by=variable], aes(xintercept=value, color=variable, linetype=variable), alpha=0.5, size=2))
    ##
print(    ggplot(dt <- melt(basescores$mean_score[,1:2])) +
        geom_histogram(aes(x=value, y=..density.., fill=variable), color='white', bins=10, alpha=0.5, position='identity') + scale_fill_bright() + ylim(c(0,7.5+2.5*3/4)) +
        xlab('probability') + xlim(c(0,1)) +
        geom_vline(data=dt[,.(value=mean(value)), by=variable], aes(xintercept=value, color=variable, linetype=variable), alpha=0.5, size=2))
    ##
    dev.off()
    }
}
save(list=c('scores'), file=paste0('scores_cont_reverse_unfact_T',nTest*3,'_N',ndata,'_',length(covNums),'covs.RData'))
elapsedtime <- Sys.time() - starttime
elapsedtime

save.image(file='_freqs_for_tests_cont_reverse.RData')


stop('End')





scores <- list()
basepriors <- rbind(asdata=normalize(as.vector(table(data$bin_RMSD))), uniform=rep(1,3)/3)
for(dset in c('tail', 'head')){
    for(ii in 1:nrow(basepriors)){
        priorR <- basepriors[ii,]
    ## Compare with a model based only on prior frequencies of r
    unseldata <- setdiff(1:nrow(data), seldata)
    nTest <- 2000
    ## construct test set
    testdata <- data[do.call(dset,list(x=unseldata,n=3*nTest)), c('bin_RMSD',covNames), with=F]
    ## testd <- data[unseldata, covNames, with=F]
    ## testdata <- data.table()
    ## for(val in rmsdVals){
    ##     testdata <- rbind(testdata, tail(testd[bin_RMSD==val],n=round(3*nTest*priorR[val])))
    ## }
    ## rm(testd)
    ##
    gc()
    plan(sequential)
    plan(multisession, workers = 6L)
    condfreqs <- predictYX(MCMCdata[[case]], testdata)
        plan(sequential)



        ##
    set.seed(247)
    scores[[dset]][[rownames(basepriors)[ii]]] <- metrics(truevalues=testdata[,bin_RMSD], probX=condfreqs[,,1], chanceprior=priorR)
    print(dset)
    pdff(paste0('scores_direct_nonfact_conditional_P',rownames(basepriors)[ii],'_',dset))
    basescores <- scores[[dset]][[ii]]
    ##
print(    ggplot(dt <- melt(basescores$delta_gain[,1:2])) +
        geom_histogram(aes(x=value, fill=variable), color='white', bins=length(unique(c(dgain))), alpha=0.5, position='identity') + scale_fill_bright() +
        xlab('delta utility') +
        geom_vline(data=dt[,.(value=mean(value)), by=variable], aes(xintercept=value, color=variable, linetype=variable), alpha=0.5, size=2))
    ##
print(    ggplot(dt <- melt(basescores$contig_gain[,1:2])) +
        geom_histogram(aes(x=value, fill=variable), color='white', bins=length(unique(c(cgain))), alpha=0.5, position='identity') + scale_fill_bright() +
        xlab('contiguous utility') +
        geom_vline(data=dt[,.(value=mean(value)), by=variable], aes(xintercept=value, color=variable, linetype=variable), alpha=0.5, size=2))
    ##
print(    ggplot(dt <- melt(basescores$log_score[,1:2])) +
        geom_histogram(aes(x=value, y=..density.., fill=variable), color='white', bins=10, alpha=0.5, position='identity') + scale_fill_bright() +
        xlab('log-probability') + xlim(c(-4,0)) +
        geom_vline(data=dt[,.(value=mean(value)), by=variable], aes(xintercept=value, color=variable, linetype=variable), alpha=0.5, size=2))
    ##
print(    ggplot(dt <- melt(basescores$mean_score[,1:2])) +
        geom_histogram(aes(x=value, y=..density.., fill=variable), color='white', bins=10, alpha=0.5, position='identity') + scale_fill_bright() +
        xlab('probability') + xlim(c(0,1)) +
        geom_vline(data=dt[,.(value=mean(value)), by=variable], aes(xintercept=value, color=variable, linetype=variable), alpha=0.5, size=2))
    ##
    dev.off()
    }}
save(list=c('scores'), file=paste0('scores_direct_continuousr_conditional_nonfact_T',nTest*3,'_N',ndata,'_',length(covNums),'covs.RData'))






stop('End')









if(false){
## TEST
plan(multisession, workers = 6L)
testX <- as.matrix(data[2])
cdata <- 2
testrx <- foreach(asample=1:length(MCMCdata[[cdata]]$nList), .combine=cbind, .inorder=F)%dopar%{
    dC <- setdiff(discreteCovs,'bin_RMSD')
    psis <- MCMCdata[[cdata]]$psiList[[asample]]
    mus <- MCMCdata[[cdata]]$muList[[asample]][continuousCovs,]
    sigmas <- MCMCdata[[cdata]]$sigmaList[[asample]][continuousCovs,continuousCovs,]
    phis <- MCMCdata[[cdata]]$phiList[[asample]]
    ##
    numer <- foreach(cluster=1:length(psis), .combine=cbind)%do%{
        psis[cluster] *
            dmvnorm(x=testX[,continuousCovs], mean=mus[,cluster], sigma=sigmas[,,cluster], log=F) *
            prod(foreach(dv=dC, .combine=c)%do%{ phis[[dv]][testX[,dv], cluster] }) *
            phis[['bin_RMSD']][,cluster]
    }
    ##
    denom <- foreach(cluster=1:length(psis), .combine=c)%do%{
        psis[cluster] *
            dmvnorm(x=testX[,continuousCovs], mean=mus[,cluster], sigma=sigmas[,,cluster], log=F) *
            prod(foreach(dv=dC, .combine=c)%do%{ phis[[dv]][testX[,dv], cluster] })
    }
    ##
    rowSums(numer)/sum(denom)
}
plan(sequential)

}

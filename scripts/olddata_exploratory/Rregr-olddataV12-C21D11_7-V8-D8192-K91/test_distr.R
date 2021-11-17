## Author: PGL  Porta Mana
## Created: 2021-03-20T10:07:17+0100
## Last-Updated: 2021-11-17T18:49:13+0100
################
## Script for evaluation of regression:
## Unfactorizable prior
## r|x exchangeable
## x not exchangeable
## r continuous
################
if(file.exists("/cluster/home/pglpm/R")){
    .libPaths(c("/cluster/home/pglpm/R",.libPaths()))
}


seed <- 149
mcsamplesfile <- '_mcsamples-Rregr-olddataV12-C21D11_7-V8-D8192-K91-I2048.rds'
utilityfile <- 'utilitiestestdata-Rregr-olddataV12-C21D11_7-V8-D8192-K91-I2048.csv'
plotsfile <- 'plotstestdata-Rregr-olddataV12-C21D11_7-V8-D8192-K91-I2048'
##
utilitymatrix <- rbind(c(0, -5, -10), c(-1, 0, -5), c(-2, -1, 0))+10
colnames(utilitymatrix) <- paste0('true_', 1:3)
rownames(utilitymatrix) <- paste0('chosen_', 1:3)
#ndata <- 48
nplots <- 128
nsubsamples <- 64
## baseversion <- 'regr-olddataV12-C21D11_'
## nclusters <- as.integer(91) #as.integer(2^6)
## ndata <- as.integer(2^13) #as.integer(2^13) # nSamples = 37969
## niter <- as.integer(2^11) #as.integer(2^11)
## niter0 <- as.integer(2^10) #as.integer(2^10)
## nstages <- as.integer(15)
## ncheckpoints <- as.integer(8)
covNames <-  c('log_RMSD'
              ,'Xtransf_ec_tanimoto_similarity'
               ## ,'Xtransf_fc_tanimoto_similarity'
               ,'docked_HeavyAtomCount'
              ,'mcs_RingCount'
              ,'template_HeavyAtomCount'
              ,'mcs_docked_NumHAcceptors'
              ,'mcs_template_NumHAcceptors'
              ,'mcs_NumHeteroAtoms'
               )

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
## library('nimble')
## library('coda')
#### End custom setup ####

##load(file='_directmodel_contR_N6000_7covs.RData')


## initial.options <- commandArgs(trailingOnly = FALSE)
## thisscriptname <- sub('--file=', "", initial.options[grep('--file=', initial.options)])
## file.copy(from=thisscriptname, to=paste0('evaluationscript-R',baseversion,'-V',length(covNames),'-D',ndata,'-K',nclusters,'.Rscript'))

if(exists('alldata')){rm(alldata)}
datafile <- '../test_oldsmiledata_id_processed_transformed_shuffled.csv'
if(!file.exists(datafile)){
    datafile <- paste0('../', datafile)
}
alldata <- fread(datafile, sep=',')
discreteCovs <- covNames[sapply(covNames, function(x){is.integer(alldata[[x]])})]
continuousCovs <- covNames[sapply(covNames, function(x){is.double(alldata[[x]])})]
covNames <- c(continuousCovs, discreteCovs)
nccovs <- length(continuousCovs)
ndcovs <- length(discreteCovs)
##
bincovNames <- c(covNames, 'bin_RMSD')
alldata <- alldata[,..bincovNames]
if(!exists('ndata')){ndata <- nrow(alldata)}

source('../functions_rmsdregr_nimble_binom.R')
resample <- function(x, ...) x[sample.int(length(x), ...)]
##
boundaries <- log(c(2,3))
rangeR <- extendrange(alldata$log_RMSD)
lRgrid <- as.integer(2^7) # as.integer(2^11)
Rgrid <- seq(rangeR[1], rangeR[2], length.out=lRgrid)
costgrid <- utilitymatrix %*% (1*sapply(Rgrid,function(aRvalue){
    c( aRvalue <= boundaries[1],
        aRvalue > boundaries[1] && aRvalue <= boundaries[2],
        aRvalue > boundaries[2] )
}))

mcsamples <- readRDS(file=mcsamplesfile)
parmList <- mcsamples2parmlist(mcsamples)
nclusters <- ncol(parmList$q)
niter <- nrow(parmList$q)

gc()

actualutilities <- foreach(adatum=1:ndata, .combine=c)%dopar%{
    datum <- alldata[adatum,]
    distrfR <- samplesfRgivenX1(datum, parmList=parmList, RMSDgrid=Rgrid)
    exputil <- (costgrid %*% normalize(rowMeans(distrfR)))[,1]
    choice <- resample(which(exputil==max(exputil)),1)
    utilitymatrix[datum$bin_RMSD, choice]
}
write.table(matrix(actualutilities,nrow=1), utilityfile, row.names=F, col.names=F, sep=',')
print(paste0('average utility:'))
print(mean(actualutilities))

transparency <- '44'
##
plotdistributions <- foreach(adatum=1:nplots, .combine=cbind)%dopar%{
    samplesfRgivenX1(alldata[adatum,], parmList=parmList, RMSDgrid=Rgrid, nfsamples=nsubsamples)
}
dim(plotdistributions) <- c(lRgrid, nsubsamples, nplots)

ufreqs <- tabulate(actualutilities+1)
names(ufreqs) <- urange <- min(actualutilities):max(actualutilities)
pdff('utilities_testset_hist')
par(mar=c(5,5,1,0))
bp <- barplot(ufreqs/ndata, xlab='utility', ylab='relative frequency', col=1, cex.axis=2, cex=2, cex.lab=2)
grid(nx=NULL,ny=NULL,lwd=3)
## bp <- barplot(ufreqs, xlab='utility', ylab='frequency', col=1, cex.axis=2, cex=2, cex.lab=2)
meanline <- diff(bp)[1]/diff(urange)[1] * mean(actualutilities) + bp[1]- diff(bp)[1]/diff(urange)[1]*urange[1]
abline(v=meanline, lty=2, lwd=6, col=3)
legend('top', legend=paste0('mean utility: ', signif(mean(actualutilities),3)), bty='n', cex=2, lty=2, lwd=3, col=3)
dev.off()


##
cC <- setdiff(covNames, 'log_RMSD')
nplots <- min(nplots, ndata)
pdff(plotsfile)
for(atest in 1:nplots){
    datum <- alldata[atest,]
    dat <- plotdistributions[,,atest]
    meandat <- rowMeans(dat)
    exputil <- (costgrid %*% normalize(meandat))[,1]
    choice <- resample(which(exputil==max(exputil)),1)
    actualutility <- utilitymatrix[datum$bin_RMSD, choice]
    ##
    matplot(x=Rgrid, y=dat, type='l', col=paste0(palette()[7], transparency), lty=1, lwd=2, xlab='log-RMSD', ylab='probability density',cex.axis=1.5, cex.lab=1.5)
    matlines(x=Rgrid, y=meandat, col=1, lty=1, lwd=4)
    matlines(x=rep(datum$log_RMSD, 2), y=c(0,max(dat)),lty=2,lwd=5,col=palette()[2])
##    matlines(x=rep(c(Rgrid %*% normalize(meandat)),2),y=c(0,max(dat)),lty=4,lwd=3,col=palette()[3])
    matlines(x=rep(boundaries[1],2),y=c(0,max(dat)),lty=3,lwd=4,col=palette()[4])
    matlines(x=rep(boundaries[2],2),y=c(0,max(dat)),lty=3,lwd=4,col=palette()[4])
    grid(lwd=1,lty=1)
    legend('topleft',legend=c( '2-3 \u00c5',
                         'predicted distributions',
                         'mean (final predictive)',
                         '',
                         paste0('true value: ', signif(datum$log_RMSD,3)),
                         paste0('true category: ', datum$bin_RMSD),
                         paste0('considered as: ', choice)
                         ),
           lty=c(3,1,1, NA, 2,NA,NA),
           col=palette()[c(4,7,1, 8, 2,8)],
           lwd=c(4,2,4, 1, 4,0),
           bty='n',cex=1.25)
        legend('left', legend=c('FEATURES:', paste0(cC,' = ',signif(datum[1,..cC],3))), bty='n', horiz=F, cex=0.75, inset=-0.01)
        legend('topright', legend=paste0('test datum #', atest), bty='n', cex=1, inset=-0.005)
}
dev.off()





## distrf <- samplesfRgivenX(alldata, parmList=parmList, RMSDgrid=Rgrid)
## saveRDS(distrf, file=paste0('testdistr-R',baseversion,'-V',length(covNames),'-D',ndata,'-K',nclusters,'-I',nrow(mcsamples),'-T',nTest,'.rds'))


## ##
## scores <- foreach(autility=1:nrow(utilities), .combine=rbind)%do%{
##     maxloss <- utilities[autility, 1]
##     maxgain <- utilities[autility, 2]
##     costgrid <- plogis(Rgrid,location=mean(boundaries), scale=sd(boundaries)/3) * (maxloss-maxgain) + maxgain
##     choices <- sapply(1:nTest, function(atest){normalize(rowMeans(distrf[atest,,])) %*% costgrid})
##     truegains <- maxloss * (choices > 0 & testdata[,'log_RMSD'] > mean(boundaries)) +
##         maxgain * (choices > 0 & testdata[,'log_RMSD'] < mean(boundaries))
##     ##
##     maxloss50 <- -maxgain
##     costgrid50 <- plogis(Rgrid,location=mean(boundaries), scale=sd(boundaries)/3) * (maxloss50-maxgain) + maxgain
##     choices50 <- sapply(1:nTest, function(atest){normalize(rowMeans(distrf[atest,,])) %*% costgrid50})
##     truegains50 <- maxloss * (choices50 > 0 & testdata[,'log_RMSD'] > mean(boundaries)) +
##         maxgain * (choices50 > 0 & testdata[,'log_RMSD'] < mean(boundaries))
##     ##
##     c(mean(truegains), mean(truegains50))
## }

## for(autility in 1:nrow(utilities)){
##     print(paste0('false pos: ', utilities[autility,1], ' | true pos: ', utilities[autility,2]))
##     print(paste0('gain: ', scores[autility,1]))
##     print(paste0('gain with 50% threshold: ', scores[autility,2]))
## }


## maxloss <- -5
## maxgain <- 1
## costgrid <- plogis(Rgrid,location=mean(boundaries), scale=sd(boundaries)/3) * (maxloss-maxgain) + maxgain
## ##
## nplots <- 128
## ndistrsamples <- 32
## transparency <- '44'
## ##
## pdff(paste0('evaluation-R',baseversion,'-V',length(covNames),'-D',ndata,'-K',nclusters,'-I',nrow(mcsamples),'-T',nTest,'.rds'))
## for(atest in 1:nplots){
##     dat <- distrf[atest,,]
##     mdat <- rowMeans(dat)
##     expgain <- normalize(mdat) %*% costgrid
##     dat <- dat[,seq(1, ncol(dat), length.out=ndistrsamples)]
##     matplot(x=Rgrid, y=dat, type='l', col=paste0(palette()[7], transparency), lty=1, lwd=2, xlab='log-RMSD', ylab='probability density',cex.axis=1.5, cex.lab=1.5)
##     matlines(x=Rgrid, y=mdat, col=1, lty=1, lwd=4)
##     matlines(x=rep(testdata[atest,'log_RMSD'], 2), y=c(0,max(dat)),lty=2,lwd=5,col=palette()[2])
## ##    matlines(x=rep(c(Rgrid %*% normalize(mdat)),2),y=c(0,max(dat)),lty=4,lwd=3,col=palette()[3])
##     matlines(x=rep(boundaries[1],2),y=c(0,max(dat)),lty=3,lwd=4,col=palette()[4])
##     matlines(x=rep(boundaries[2],2),y=c(0,max(dat)),lty=3,lwd=4,col=palette()[4])
##     grid(lwd=1,lty=1)
##     legend('topleft',legend=c( '2-3 \u00c5',
##                          'predicted distributions',
##                          'mean (final predictive)',
##                          '',
##                          paste0('true value: ', signif(testdata[atest,'log_RMSD'],3)),
##                          paste0('EXPECTED GAIN = ', signif(expgain,2))
##                          ),
##            lty=c(3,1,1, NA, 2,NA),
##            col=palette()[c(4,7,1, 8, 2,8)],
##            lwd=c(4,2,4, 1, 4,0),
##            bty='n',cex=1.25)
##         legend('topright', legend=c('FEATURES:', paste0(colnames(testdata)[-1],' = ',signif(testdata[atest,-1],3))), bty='n', horiz=F, inset=-0.005)
## }
## dev.off()


## ## #######################################################################

## ## rm(alldata)
## ## alldata <- fread('../data_processed_transformed_rescaled_shuffled.csv', sep=' ')
## ## nameFeatures <- names(alldata)
## ## nSamples <- nrow(alldata)
## ## nFeatures <- ncol(alldata)
## ## ##
## ## ##
## ## #set.seed(222)
## ## covNames <-  c('log_RMSD',
## ##                'log_mcs_unbonded_polar_sasa',
## ##                'logit_ec_tanimoto_similarity',
## ##                'mcs_NumHeteroAtoms',
## ##                ##'scale_fc_tanimoto_similarity'
## ##                'docked_HeavyAtomCount',
## ##                'mcs_RingCount',
## ##                'docked_NumRotatableBonds'
## ##                )
## ## discreteCovs <- covNames[sapply(covNames, function(x){is.integer(alldata[[x]])})]
## ## continuousCovs <- covNames[sapply(covNames, function(x){is.double(alldata[[x]])})]
## ## covNames <- c(continuousCovs, discreteCovs)
## ## nclusters <- 100
## ## ndata <- 6000 # nSamples = 37969
## ## nccovs <- length(continuousCovs)
## ## ndcovs <- length(discreteCovs)
## ## ##
## ## mcsamples <- readRDS(file='_seminar_mcsamples.rds')
## ## parmNames <- c('q', 'meanC', 'tauC', 'probD', 'sizeD')
## ## parmList <- foreach(var=parmNames)%dopar%{
## ##     out <- mcsamples[,grepl(paste0(var,'\\['), colnames(mcsamples))]
## ##     if(grepl('C', var)){
## ##         dim(out) <- c(nrow(mcsamples), nccovs, nclusters)
## ##         dimnames(out) <- list(NULL, continuousCovs, NULL)
## ##     } else if(grepl('D', var)){
## ##         dim(out) <- c(nrow(mcsamples), ndcovs, nclusters)
## ##         dimnames(out) <- list(NULL, discreteCovs, NULL)
## ##     } else {dim(out) <- c(nrow(mcsamples), nclusters) }
## ##     out
## ## }
## ## names(parmList) <- parmNames

## ## seldata <- 1:ndata
## ## unseldata <- setdiff(1:nrow(alldata), seldata)
## ## nTest <- 200
## ## ## construct test set
## ## testdata <- as.matrix(alldata[do.call('tail',list(x=unseldata,n=nTest)), c(covNames), with=F])
## ## ##
## ## boundaries <- log(c(2,3))
## ## rangeR <- range(alldata$log_RMSD, na.rm=T)
## ## rangeR <- (rangeR - mean(rangeR)) * 1.1 + mean(rangeR)
## ## lRgrid <- 201
## ## Rgrid <- seq(rangeR[1],rangeR[2],length.out=lRgrid)

## ## source('functions_rmsdregr_nimble_negbinom.R')
## ## normalize <- function(freqs){freqs/sum(freqs)}
## ## gc()
## ## plan(sequential)
## ## plan(multisession, workers = 6L)
## ## distrf <- samplesfRgivenX(testdata, parmList=parmList, RMSDgrid=Rgrid)
## ## plan(sequential)
## ## ##dimnames(distrf) <- list(NULL, paste0('testid',1:nrow(testdata)), NULL)


## ## maxloss50 <- -1
## ## maxgain <- 1
## ## costgrid <- plogis(Rgrid,location=mean(boundaries), scale=sd(boundaries)/3) * (maxloss-maxgain) + maxgain
## ## ## matplot(Rgrid, costgrid, type='l')
## ## ## matlines(x=rep(boundaries[1],2),y=c(maxloss,maxgain),lty=3,lwd=4,col=palette()[4])
## ## ## matlines(x=rep(boundaries[2],2),y=c(maxloss,maxgain),lty=3,lwd=4,col=palette()[4])
## ## choices <- sapply(1:nTest,function(atest){normalize(rowMeans(distrf[atest,,])) %*% costgrid})
## ## truegains <- maxloss * (choices > 0 & testdata[,'log_RMSD'] > mean(boundaries)) +
## ##         maxgain * (choices > 0 & testdata[,'log_RMSD'] < mean(boundaries))
## ## mean(truegains)


## ## maxloss <- -5
## ## maxgain <- 1
## ## costgrid <- plogis(Rgrid,location=mean(boundaries), scale=sd(boundaries)/3) * (maxloss-maxgain) + maxgain
## ## ## matplot(Rgrid, costgrid, type='l')
## ## ## matlines(x=rep(boundaries[1],2),y=c(maxloss,maxgain),lty=3,lwd=4,col=palette()[4])
## ## ## matlines(x=rep(boundaries[2],2),y=c(maxloss,maxgain),lty=3,lwd=4,col=palette()[4])
## ## choices <- sapply(1:nTest,function(atest){normalize(rowMeans(distrf[atest,,])) %*% costgrid})
## ## truegains <- maxloss * (choices > 0 & testdata[,'log_RMSD'] > mean(boundaries)) +
## ##         maxgain * (choices > 0 & testdata[,'log_RMSD'] < mean(boundaries))
## ## ##
## ## maxloss50 <- -1
## ## costgrid50 <- plogis(Rgrid,location=mean(boundaries), scale=sd(boundaries)/3) * (maxloss50-maxgain) + maxgain
## ## choices50 <- sapply(1:nTest,function(atest){normalize(rowMeans(distrf[atest,,])) %*% costgrid50})
## ## truegains50 <- maxloss * (choices50 > 0 & testdata[,'log_RMSD'] > mean(boundaries)) +
## ##         maxgain * (choices50 > 0 & testdata[,'log_RMSD'] < mean(boundaries))
## ## c(mean(truegains), mean(truegains50))

## ## ##
## ## nplotsamples <- nTest
## ## transparency <- '44'
## ## ##
## ## pdff(paste0('_seminar_distribution_examples'))
## ## for(atest in 1:nrow(testdata)){
## ##     dat <- distrf[atest,,]
## ##     mdat <- rowMeans(dat)
## ##     expgain <- normalize(mdat) %*% costgrid
## ##     dat <- dat[,seq(1, ncol(dat), length.out=ndistrsamples)]
## ##     matplot(x=Rgrid, y=dat, type='l', col=paste0(palette()[7], transparency), lty=1, lwd=2, xlab='log-RMSD', ylab='probability density',cex.axis=1.5, cex.lab=1.5)
## ##     matlines(x=Rgrid, y=mdat, col=palette()[1], lty=1, lwd=4)
## ##     matlines(x=rep(testdata[atest,'log_RMSD'], 2), y=c(0,max(dat)),lty=2,lwd=5,col=palette()[2])
## ## ##    matlines(x=rep(c(Rgrid %*% normalize(mdat)),2),y=c(0,max(dat)),lty=4,lwd=3,col=palette()[3])
## ##     matlines(x=rep(boundaries[1],2),y=c(0,max(dat)),lty=3,lwd=4,col=palette()[4])
## ##     matlines(x=rep(boundaries[2],2),y=c(0,max(dat)),lty=3,lwd=4,col=palette()[4])
## ##     grid(lwd=1,lty=1)
## ##     legend('topleft',legend=c( '2-3 \u00c5',
## ##                          'predicted distributions',
## ##                          'mean (final predictive)',
## ##                          '',
## ##                          'true value',
## ##                          paste0('EXPECTED GAIN = ', signif(expgain,2))
## ##                          ),
## ##            lty=c(3,1,1, NA, 2,NA),
## ##            col=palette()[c(4,7,1, 8, 2,8)],
## ##            lwd=c(4,2,4, 1, 4,0),
## ##            bty='n',cex=1.25)
## ##         legend('topright',legend=paste0(colnames(testdata),' = ',signif(testdata[atest,],3)), bty='n',horiz=F,inset=-0.005)
## ## }
## ## dev.off()


## ## datum1 <- 3
## ## grid1 <- 5
## ## fsample1 <- 7
## ## testq <- parmList$q[fsample1,]
## ## testmeanC <- parmList$meanC[fsample1,,]
## ## testtauC <- parmList$tauC[fsample1,,]
## ## testprobD <- parmList$probD[fsample1,,]
## ## testsizeD <- parmList$sizeD[fsample1,,]
## ## ##
## ## testdatum <- testdata[datum1,]
## ## ## testdatum['log_RMSD'] <- Rgrid[grid1]
## ## testcC <- setdiff(continuousCovs, 'log_RMSD')
## ## ##
## ## testW <- foreach(acluster=1:length(testq), .combine=c)%do%{
## ##     testq[acluster] *
## ##         prod(dnorm(testdatum[testcC], mean=testmeanC[testcC,acluster], sd=1/sqrt(testtauC[testcC,acluster]))) *
## ##         prod(dnbinom(testdatum[discreteCovs], prob=testprobD[,acluster], size=testsizeD[,acluster]))
## ## }
## ## testcond <- sum(testW * dnorm(Rgrid[grid1], mean=testmeanC['log_RMSD',], sd=1/sqrt(testtauC['log_RMSD',])))/sum(testW)
## ## testcond
## ## distrf[datum1,grid1,fsample1]




## ## ## pdff(paste0('distributions_samples_','tail'))
## ## ## for(atest in 1:nrow(testdata)){
## ## ##     print(
## ## ##         dat <- as.data.table(t(distrf[,atest,]))
## ## ##         colnames(dat) <- paste0('log_RMSD')
## ## ##         dat$sample <- paste0('sample',1:dim(distrf)[3])
## ## ##         dat <- melt(dat, id.vars='sample')
## ## ##         mdat$log_RMSD = as.numeric(gsub("time", "", mdat$variable))        
## ## ##         trueRvalue <- testdata[atest,log_RMSD]
## ## ##         ggplot(dat) +
## ## ##         geom_line(aes(x=log_RMSD, y=probability_density), size=1.5, colour=palette()[1]) +
## ## ##         geom_vline(aes(xintercept=trueR, colour='true value'), linetype=2, alpha=0.75, size=2) + #scale_color_identity() +
## ## ##         geom_vline(aes(xintercept=I(boundaries[1])), colour=palette()[4], linetype=3, alpha=0.75, size=1.5) + 
## ## ##         geom_vline(aes(xintercept=I(boundaries[2])), colour=palette()[4], linetype=3, alpha=0.75, size=1.5)  +   #scale_color_identity(name = "Model fit", breaks = c("black", "red", "blue"), labels = c("Linear", "Quadratic", "Cubic"), guide = "legend") +
## ## ##         scale_color_manual(name='', breaks=c('true value','2-3 A'), values=c('true value'=palette()[2], '2-3 A'=palette()[4]) ) +
## ## ##         scale_linetype_manual(name='', breaks=c('true value','2-3 A'), values=c('true value'='twodash', '2-3 A'='dotted') ) +
## ## ##         theme(legend.pos='top')
## ## ##     )
## ## ## }
## ## ## dev.off()




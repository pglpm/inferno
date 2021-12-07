## Author: PGL  Porta Mana
## Created: 2021-03-20T10:07:17+0100
## Last-Updated: 2021-09-18T15:14:42+0200
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
scale_colour_continuous <- scale_colour_bright
scale_colour_identity <- scale_colour_bright
#library('cowplot')
library('png')
library('foreach')
library('doFuture')
library('doRNG')
registerDoFuture()
#library('LaplacesDemon') # used for Dirichlet generator
library('ash')
## library('extraDistr')
## library('PReMiuM')
## library('mvtnorm')
options(bitmapType='cairo')

load(file='_directmodel_contR_N6000_7covs.RData')

ndata <- 6000
seldata <- 1:ndata
unseldata <- setdiff(1:nrow(alldata), seldata)
nTest <- 30
## construct test set
testdata <- as.matrix(alldata[do.call('tail',list(x=unseldata,n=3*nTest)), c(covNames), with=F])
##
boundaries <- log(c(2,3))
rangeR <- range(alldata$log_RMSD, na.rm=T)
rangeR <- (rangeR - mean(rangeR)) * 1.1 + mean(rangeR)
lRgrid <- 201
Rgrid <- seq(rangeR[1],rangeR[2],length.out=lRgrid)

gc()
plan(sequential)
plan(multisession, workers = 6L)
distrf <- samplesfRgivenX(testdata, parmList=parmList, RMSDgrid=Rgrid)
plan(sequential)
##dimnames(distrf) <- list(NULL, paste0('testid',1:nrow(testdata)), NULL)

##
nplotsamples <- 100
transparency <- '44'
pdff(paste0('post_distributionsamples_direct_','tail','-run',version))
for(atest in 1:nrow(testdata)){
    dat <- distrf[atest,,]
    mdat <- rowMeans(dat)
    dat <- dat[,seq(1, ncol(dat), length.out=nplotsamples)]
    matplot(x=Rgrid, y=dat, type='l', col=paste0(palette()[7], transparency), lty=1, lwd=2, xlab='log-RMSD', ylab='probability density',cex.axis=1.5, cex.lab=1.5)
    matlines(x=Rgrid, y=mdat, col=palette()[1], lty=1, lwd=4)
    matlines(x=rep(testdata[atest,'log_RMSD'], 2), y=c(0,max(dat)),lty=2,lwd=5,col=palette()[2])
    matlines(x=rep(c(Rgrid %*% normalize(mdat)),2),y=c(0,max(dat)),lty=4,lwd=3,col=palette()[3])
    matlines(x=rep(boundaries[1],2),y=c(0,max(dat)),lty=3,lwd=4,col=palette()[4])
    matlines(x=rep(boundaries[2],2),y=c(0,max(dat)),lty=3,lwd=4,col=palette()[4])
    grid(lwd=1,lty=1)
    legend('topleft',legend=c('predicted distributions','mean (final predictive)','true value','quadratic-gain choice', '2-3 \u00c5'),lty=c(1,1,2,4,3),col=palette()[c(7,1,2,3,4)], lwd=c(2,4,4,3,4),bty='n',cex=1.25)
        legend('topright',legend=paste0(colnames(testdata),' = ',signif(testdata[atest,],3)), bty='n',horiz=F,inset=-0.005)
}
dev.off()


datum1 <- 3
grid1 <- 5
fsample1 <- 7
testq <- parmList$q[fsample1,]
testmeanC <- parmList$meanC[fsample1,,]
testtauC <- parmList$tauC[fsample1,,]
testprobD <- parmList$probD[fsample1,,]
testsizeD <- parmList$sizeD[fsample1,,]
##
testdatum <- testdata[datum1,]
## testdatum['log_RMSD'] <- Rgrid[grid1]
testcC <- setdiff(continuousCovs, 'log_RMSD')
##
testW <- foreach(acluster=1:length(testq), .combine=c)%do%{
    testq[acluster] *
        prod(dnorm(testdatum[testcC], mean=testmeanC[testcC,acluster], sd=1/sqrt(testtauC[testcC,acluster]))) *
        prod(dnbinom(testdatum[discreteCovs], prob=testprobD[,acluster], size=testsizeD[,acluster]))
}
testcond <- sum(testW * dnorm(Rgrid[grid1], mean=testmeanC['log_RMSD',], sd=1/sqrt(testtauC['log_RMSD',])))/sum(testW)
testcond
distrf[datum1,grid1,fsample1]




## pdff(paste0('distributions_samples_','tail'))
## for(atest in 1:nrow(testdata)){
##     print(
##         dat <- as.data.table(t(distrf[,atest,]))
##         colnames(dat) <- paste0('log_RMSD')
##         dat$sample <- paste0('sample',1:dim(distrf)[3])
##         dat <- melt(dat, id.vars='sample')
##         mdat$log_RMSD = as.numeric(gsub("time", "", mdat$variable))        
##         trueRvalue <- testdata[atest,log_RMSD]
##         ggplot(dat) +
##         geom_line(aes(x=log_RMSD, y=probability_density), size=1.5, colour=palette()[1]) +
##         geom_vline(aes(xintercept=trueR, colour='true value'), linetype=2, alpha=0.75, size=2) + #scale_color_identity() +
##         geom_vline(aes(xintercept=I(boundaries[1])), colour=palette()[4], linetype=3, alpha=0.75, size=1.5) + 
##         geom_vline(aes(xintercept=I(boundaries[2])), colour=palette()[4], linetype=3, alpha=0.75, size=1.5)  +   #scale_color_identity(name = "Model fit", breaks = c("black", "red", "blue"), labels = c("Linear", "Quadratic", "Cubic"), guide = "legend") +
##         scale_color_manual(name='', breaks=c('true value','2-3 A'), values=c('true value'=palette()[2], '2-3 A'=palette()[4]) ) +
##         scale_linetype_manual(name='', breaks=c('true value','2-3 A'), values=c('true value'='twodash', '2-3 A'='dotted') ) +
##         theme(legend.pos='top')
##     )
## }
## dev.off()




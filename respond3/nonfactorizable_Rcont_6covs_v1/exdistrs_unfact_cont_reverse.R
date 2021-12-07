## Author: PGL  Porta Mana
## Created: 2021-03-20T10:07:17+0100
## Last-Updated: 2021-08-25T08:05:16+0200
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
library('extraDistr')
library('PReMiuM')
library('mvtnorm')
options(bitmapType='cairo')

load(file='_directmodel_contR_N6000_7covs.RData')
##################################################################
##################################################################
## Function to calculate \sum_t F^{(t)}_{r|x}
##
rangeR <- fivenum(data$log_RMSD, na.rm=T)
rangeR <- c(diff(rangeR[c(3,1)]), diff(rangeR[c(3,5)]))*1.1+rangeR[3]
cC <- setdiff(continuousCovs,'log_RMSD')
## predictfXall <- function(dataobj, X){
##     X <- as.matrix(X[, c(discreteCovs,continuousCovs), with=FALSE])
##     ntestdata <- nrow(X)
##     freqs <- foreach(asample=seq_along(dataobj$nList), .combine=cbind, .inorder=FALSE)%dopar%{
##         numer <- colSums(
##             exp(
##             log(dataobj$psiList[[asample]]) +
##             t(vapply(seq_len(dataobj$nList[asample]), function(cluster){
##                 rowSums(log(
##                     vapply(discreteCovs, function(covariate){
##                         dataobj$phiList[[asample]][[covariate]][X[,covariate], cluster]
##                     }, numeric(ntestdata*lrgrid))
##                 )) +
##                     dmvnorm(X[,continuousCovs], mean=dataobj$muList[[asample]][continuousCovs,cluster], sigma=as.matrix(dataobj$sigmaList[[asample]][continuousCovs,continuousCovs,cluster]), log=TRUE)
##             }, numeric(ntestdata*lrgrid)))
##             )
##             )

##         denom <- colSums(
##             exp(
##             log(dataobj$psiList[[asample]]) +
##             t(vapply(seq_len(dataobj$nList[asample]), function(cluster){
##                 ## rowSums(log(
##                 ##     vapply(discreteCovs, function(covariate){
##                 ##         dataobj$phiList[[asample]][[covariate]][X[,covariate], cluster]
##                 ##     }, numeric(ntestdata))
##                 ## )) +
##                     dnorm(X[,'log_RMSD'], mean=dataobj$muList[[asample]]['log_RMSD',cluster], sd=sqrt(dataobj$sigmaList[[asample]]['log_RMSD','log_RMSD',cluster]), log=TRUE)
##             }, numeric(ntestdata*lrgrid)))
##             )
##             )

##         numer/denom
##     }
##         ## me <- rowMeans(freqs)
##     ## freqs <- cbind(means=me, stds=sqrt(rowMeans(freqs^2) - me^2))
##     ## dim(freqs) <- c(ntestdata, lrgrid, 2)
##     ## dimnames(freqs) <- list(NULL, NULL, c('means','stds'))
##     ## freqs
##     freqs <- rowMeans(freqs)
##     dim(freqs) <- c(ntestdata, lrgrid)
##     freqs
## }
##
predictfXall <- function(dataobj, X){
    X <- as.matrix(X[, c(discreteCovs,continuousCovs), with=FALSE])
    ntestdata <- nrow(X)
    freqs <- foreach(asample=seq_along(dataobj$nList), .combine=cbind, .inorder=FALSE)%dopar%{
        numer <- colSums(
            exp(
            log(dataobj$psiList[[asample]]) +
            t(vapply(seq_len(dataobj$nList[asample]), function(cluster){
                rowSums(log(
                    vapply(discreteCovs, function(covariate){
                        dataobj$phiList[[asample]][[covariate]][X[,covariate], cluster]
                    }, numeric(ntestdata))
                )) +
                    dmvnorm(X[,continuousCovs], mean=dataobj$muList[[asample]][continuousCovs,cluster], sigma=as.matrix(dataobj$sigmaList[[asample]][continuousCovs,continuousCovs,cluster]), log=TRUE)
            }, numeric(ntestdata)))
            )
            )

        denom <- colSums(
            exp(
            log(dataobj$psiList[[asample]]) +
            t(vapply(seq_len(dataobj$nList[asample]), function(cluster){
                ## rowSums(log(
                ##     vapply(discreteCovs, function(covariate){
                ##         dataobj$phiList[[asample]][[covariate]][X[,covariate], cluster]
                ##     }, numeric(ntestdata))
                ## )) +
                    dnorm(X[,'log_RMSD'], mean=dataobj$muList[[asample]]['log_RMSD',cluster], sd=sqrt(dataobj$sigmaList[[asample]]['log_RMSD','log_RMSD',cluster]), log=TRUE)
            }, numeric(ntestdata)))
            )
            )

        numer/denom
    }
    freqs
}

    unseldata <- setdiff(1:nrow(data), seldata)
    nTest <- 30
    ## construct test set
    testdata <- data[do.call('tail',list(x=unseldata,n=3*nTest)), c(covNames), with=F]
    ## testd <- data[unseldata, covNames, with=F]
    ## testdata <- data.table()
    ## for(val in rmsdVals){
    ##     testdata <- rbind(testdata, tail(testd[bin_RMSD==val],n=round(3*nTest*priorR[val])))
    ## }
    ## rm(testd)
##
rgrid <- seq(rangeR[1],rangeR[2],length.out=101)
XX <- data.table()
for(atest in 1:nrow(testdata)){
    XX <- rbind(XX,
                cbind(testid=atest, log_RMSD=rgrid, testdata[atest,!'log_RMSD'])
                )
}
##
gc()
plan(sequential)
plan(multisession, workers = 6L)
distrf <- predictfXall(MCMCdata[[1]], XX)
plan(sequential)
dim(distrf) <- c(length(rgrid), nrow(testdata), ncol(distrf))
dimnames(distrf) <- list(NULL, paste0('testid',1:nrow(testdata)), NULL)
##

nsamples <- 100
pdff(paste0('distributions_samples_reverse_','tail'))
for(atest in 1:nrow(testdata)){
    dat <- t(normalizem(t(distrf[,atest,round(seq(1,dim(distrf)[3],length.out=nsamples))])))
    matplot(x=rgrid, y=dat, type='l', col=paste0(palette()[5],'44'), lty=1, lwd=2, xlab='log-RMSD', ylab='probability density',cex.axis=1.5, cex.lab=1.5)
    matlines(x=rgrid, y=rowMeans(dat), col=palette()[1], lty=1, lwd=3)
    matlines(x=rep(testdata[atest,log_RMSD],2),y=c(0,max(dat)),lty=2,lwd=4,col=palette()[2])
    matlines(x=rep(c(rgrid %*% normalize(rowMeans(dat))),2),y=c(0,max(dat)),lty=4,lwd=3,col=palette()[3])
    matlines(x=rep(boundaries[1],2),y=c(0,max(dat)),lty=3,lwd=4,col=palette()[4])
    matlines(x=rep(boundaries[2],2),y=c(0,max(dat)),lty=3,lwd=4,col=palette()[4])
    grid(lwd=1,lty=1)
    legend('topleft',legend=c('predicted distributions','mean (final predictive)','true value','quadratic-gain choice', '2-3 \u00c5'),lty=c(1,1,2,4,3),col=palette()[c(5,1,2,3,4)], lwd=c(2,3,4,3,4),bty='n',cex=1.25)
        legend('topright',legend=paste0(names(testdata[atest]),' = ',signif(testdata[atest],3)), bty='n',horiz=F,inset=-0.005)
}
dev.off()
    
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




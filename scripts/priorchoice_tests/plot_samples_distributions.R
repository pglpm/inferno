## Author: PGL  Porta Mana
## Created: 2021-03-20T10:07:17+0100
## Last-Updated: 2021-09-02T13:04:02+0200
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
pdff <- function(filename){pdf(file=paste0(filename,'.pdf'),paper='a4r',height=11.7,width=16.5)} # to output in pdf format
pngf <- function(filename,res=300){png(file=paste0(filename,'.png'),height=11.7*1.2,width=16.5,units='in',res=res,pointsize=36)} # to output in pdf format
#### End custom setup ####

indir <- 'output4valpha-prior_N0_Ncov7/'
##
## setwd(indir)
mcmcoutput <- readRDS(paste0(indir, '_mcmcoutput.rds'))
mcmcrun <- readRDS(paste0(indir, '_mcmcrun.rds'))
continuousCovs <- mcmcrun$continuousCovs
discreteCovs <- mcmcrun$discreteCovs
covNames <- mcmcrun$covNames
##################################################################
##################################################################
## Function to calculate \sum_t F^{(t)}_{r,x}
##
##
predictfXall <- function(dataobj, X){
    cC <- colnames(X)[colnames(X) %in% continuousCovs]
    dC <- colnames(X)[colnames(X) %in% discreteCovs]
    if(is.data.table(X)){ X <- as.matrix(X[, c(dC,cC), with=FALSE]) }
    ntestdata <- nrow(X)
    freqs <- foreach(asample=seq_along(dataobj$nList), .combine=cbind, .inorder=FALSE)%dopar%{
        colSums(
            exp(
                log(dataobj$psiList[[asample]]) +
                t(vapply(seq_len(dataobj$nList[asample]), function(cluster){
                    ## discrete covariates
                    if(length(dC) > 0){
                        rowSums(log(
                            vapply(dC, function(covariate){
                                dataobj$phiList[[asample]][[covariate]][X[,covariate], cluster]
                            }, numeric(ntestdata))
                        ))
                    }else{0} +
                        ## continuous covariates
                        if(length(cC) > 0){
                            dmvnorm(X[,cC,drop=F], mean=dataobj$muList[[asample]][cC,cluster], sigma=as.matrix(dataobj$sigmaList[[asample]][cC,cC,cluster]), log=TRUE)
                        }else{0}
                    ##            }))
                }, numeric(ntestdata)))
            )
        )
    }
    freqs
}
##
## rm(alldata)
## alldata <- fread('../processed_data_scaled.csv', sep=' ')
plotVars <- c('log_RMSD'
  , 'scale_mcs_unbonded_polar_sasa'
 ## , 'scale_ec_tanimoto_similarity'
  ## , 'mcs_NumHeteroAtoms'
  ## , 'docked_HeavyAtomCount'
  ## , 'mcs_RingCount'
  ## , 'docked_NumRotatableBonds'
  )
##
ngridpoints <- 100
plotgrids <- lapply(plotVars,function(var){
    varrange <- fivenum(alldata[[var]], na.rm=T)
    varrange <- c(diff(varrange[c(3,1)]), diff(varrange[c(3,5)]))*1.1+varrange[3]
    seq(varrange[1], varrange[2], length.out=ngridpoints)
})
names(plotgrids) <- plotVars
##
##gridpoints2 <- mesh(plotgrids[[1]], plotgrids[[2]])
gridpoints <- cbind(rep(plotgrids[[1]],each=length(plotgrids[[2]])),
                    rep(plotgrids[[2]],length(plotgrids[[1]])))
colnames(gridpoints) <- plotVars
##
plan(sequential)
plan(multisession, workers = 6L)
zvalues <- predictfXall(mcmcoutput, gridpoints)
plan(sequential)
dim(zvalues) <- c(ngridpoints, ngridpoints, ncol(zvalues))
##
pdff(paste0(indir, 'samplesvars'))
for(asample in round(seq(1, ncol(zvalues), length.out=100))){
persp(plotgrids[[1]],plotgrids[[2]],zvalues[,,asample],zlim=c(0,max(zvalues)),ticktype='detailed',theta = 45, phi = 15,xlab=plotVars[1],ylab=plotVars[2],zlab='freq. density')
}
dev.off()


require('plotly')
##
fig <- plot_ly(x=plotgrids[[1]], y=plotgrids[[2]], z=zvalues[,,2]) %>% layout(xaxis = list(title = 'X Axis Title'))  %>% add_surface()
##
fig

Sys.setenv(PATH = paste0(Sys.getenv('PATH'),"C:\\Program Files\\orca\\orca.exe;"))

orca(p=fig, file='test3d2.svg')

pdff('test3d')
print(fig)
dev.off()


##
dvalues <- mapply(function(xx,yy){
    predictfXall(mcmcoutput, alldata[1:5,plotVars,with=F])
}, gridpoints$x, gridpoints$y)
##
plan(sequential)



dim(dvalues) <- sapply(plotgrids,length)



    
dvalues <- outer(plotgrids[[1]],plotgrids[[2]],function(x,y){
    mapply(function(xx,yy){
                    XX <- c(xx,yy)
                    names(XX) <- continuousCovs
                    YY <- c(sasa2y(xx),tanimoto2y(yy))
                    names(YY) <- continuousCovs
                    predictSampleSasaTani(MCMCdata[[j]], asample, YY)*jac(XX)
                    },
        c(grid),y)}
        )


require('plotly')

for(j in 1:length(MCMCdata)){
pdf(file=paste0('testsampledensities_seq_',cases[j],'.pdf'),height=11.7,width=16.5*10)
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


fig <- plot_ly(x=1:nrow(volcano), y=1:ncol(volcano), z=volcano) %>% add_surface()


mesh()

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
boundaries <- log(c(2,3))
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
pdff(paste0('distributions_samples_direct_','tail'))
for(atest in 1:nrow(testdata)){
    dat <- distrf[,atest,round(seq(1,dim(distrf)[3],length.out=nsamples))]
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




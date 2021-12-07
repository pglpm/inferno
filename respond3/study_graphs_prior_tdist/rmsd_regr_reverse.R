## Author: PGL  Porta Mana
## Created: 2021-03-20T10:07:17+0100
## Last-Updated: 2021-07-15T08:46:19+0200
################
## Script for reverse regression
################

#### Custom setup ####
## Colour-blind friendly palettes, from https://personal.sron.nl/~pault/
## (consider using khroma package instead)
library('RColorBrewer')
mypurpleblue <- '#4477AA'
myblue <- '#66CCEE'
mygreen <- '#228833'
myyellow <- '#CCBB44'
myred <- '#EE6677'
myredpurple <- '#AA3377'
mygrey <- '#BBBBBB'
mypalette <- c(myblue, myred, mygreen, myyellow, myredpurple, mypurpleblue, mygrey, 'black')
palette(mypalette)
barpalette <- colorRampPalette(c(mypurpleblue,'white',myredpurple),space='Lab')
barpalettepos <- colorRampPalette(c('white','black'),space='Lab')
#dev.off()
####
library('data.table')
library('khroma')
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

## Good values:
## musasa <- 1
## mutani <- 0
## mu0 <- c(musasa,mutani)
## names(mu0) <- c('sasa','tanimoto')
## varmusasa <- 2^2
## varmutani <- 2^2
## sigmu0 <- c(varmusasa,varmutani)
## names(sigmu0) <- c('sasa','tanimoto')
## expvarsasa <- (1/2)^2
## expvartani <- (1/2)^2
## expvar0 <- c(expvarsasa,expvartani)
## names(expvar0) <- c('sasa','tanimoto')
## df0 <- 30
## alpha <- 4 or 5

doplots <- FALSE
## specify hyperparams
musasa <- 1
mutani <- 0
mu0 <- c(musasa,mutani)
names(mu0) <- c('sasa','tanimoto')
varmusasa <- 2^2
varmutani <- 2^2
sigmu0 <- c(varmusasa,varmutani)
names(sigmu0) <- c('sasa','tanimoto')
expvarsasa <- (1/2)^2
expvartani <- (1/2)^2
expvar0 <- c(expvarsasa,expvartani)
names(expvar0) <- c('sasa','tanimoto')
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
##
## Skip all preprocessing below if data is already saved
##if(TRUE){
## Read and reorganize data
rm(data)
data <- as.data.table(read.csv(file='../dataset_template_based_docking_predict_rmsd.csv', header=T, sep=','))
data <- data[!is.na(rmsd)]
data <- data[, which(sapply(data, is.numeric)==TRUE), with=FALSE]
origdata <- data
## doublecols <- which(sapply(data, function(x){all(is.double(x))}))
## posdoublecols <- which(sapply(data, function(x){all(is.double(x))&&all(x>0)}))
data$rmsd <- log(data$rmsd)
names(data)[which(names(data)=='rmsd')] <- 'log_RMSD'
rmsdThreshold <- c(2, 2.5, 3)
logRmsdThreshold <- log(rmsdThreshold)
## transform 'sasa' features to log-scale
indx <- grepl('sasa', colnames(data))
for(elem in colnames(data)[indx]){
    datum <- sasa2y(data[[elem]])
    eps <- max(diff(sort(unique(datum[abs(datum)!=Inf]))))
    datum[datum==-Inf] <- min(datum[abs(datum)!=Inf])-2*eps
    data[, elem] <- datum
}
names(data)[indx] <- paste0('scale_',names(data)[indx])
## transform 'tanimoto' features to logit scale
indx <- grepl('tanimoto', colnames(data))
for(elem in colnames(data)[indx]){
    datum <- tanimoto2y(data[[elem]])
    eps <- max(diff(sort(unique(datum[abs(datum)!=Inf]))))
    datum[datum==Inf] <- max(datum[abs(datum)!=Inf])+2*eps
    datum[datum==-Inf] <- min(datum[abs(datum)!=Inf])-2*eps
    data[, elem] <- datum
}
names(data)[indx] <- paste0('scale_',names(data)[indx])
## shift integer features to start from value 1
indx <- sapply(1:ncol(data), function(x){is.integer(data[[x]])})
for(elem in colnames(data)[indx]){
    data[, elem] <- data[, ..elem] - min(data[, ..elem],na.rm=TRUE) +1L
}
##
data <- data.table(bin_RMSD=as.integer(1+(data$log_RMSD>logRmsdThreshold[1])+(data$log_RMSD>logRmsdThreshold[3])), data)
nameFeatures <- names(data)
nSamples <- nrow(data)
nFeatures <- ncol(data)
##
## Format bins to calculate mutual info
nbinsq <- 6
##
breakFeatures <- list()
for(i in 1:ncol(data)){
    datum <- data[[i]]
    summa <- fivenum(datum)
    drange <- diff(range(datum))
    #print(paste0('i',i));print(drange)
    if(is.integer(datum)){
        breaks <- (summa[1]:(summa[5]+1))-0.5
    } else {
        width <- diff(summa[c(2,4)])/nbinsq
        nbins <- round(drange/width)
        breaks <- seq(summa[1]-drange/(nbins*100), summa[5]+drange/(nbins*100), length.out=nbins)
    }
    breakFeatures[[i]] <- breaks
}
names(breakFeatures) <- names(data)
##
##
## Mutual infos with RMSD
minfos <- matrix(NA,4,nFeatures)
rownames(minfos) <- c('MI','norm_MI','entropy','cond_entropy')
yVar <- 'log_RMSD'
##yVar <- 'bin_RMSD'
colRmsd <- which(names(data)==yVar)
rangeRmsd <- range(breakFeatures[[yVar]])
colnames(minfos) <- names(data)
for(i in names(data)){
        freqs <- normalize(bin2(x=cbind(data[[yVar]],data[[i]]),
                                ab=rbind(rangeRmsd,range(breakFeatures[[i]])),
                                nbin=c(length(breakFeatures[[yVar]]), length(breakFeatures[[i]]))-1 )$nc)
        mi <- mutualinfo(freqs)
        en <- entropy(colSums(freqs))
        enrmsd <- entropy(rowSums(freqs))
        conden21 <- condentropy21(t(freqs)) 
        minfos[,i] <- c(mi, mi/min(en,enrmsd),en, conden21)
}
##
reorder <- order(minfos[1,], decreasing=TRUE)
minfos <- minfos[,reorder]
data <- data[, ..reorder]
origdata <- origdata[, setdiff(reorder-1,0), with=FALSE]
breakFeatures <- breakFeatures[reorder]
##
##
## Plots
if(doplots==TRUE){
pdff('histograms_scaled_data')
for(i in 1:ncol(data)){
    datum <- data[[i]]
    breaks <- breakFeatures[[i]]
    print(ggplot(data[,..i], aes_(x=as.name(names(data)[i]))) + geom_histogram(breaks=breaks))
}
dev.off()
##
pdff('plotslogMI')
for(k in 1:ncol(minfos)){
    mi <- signif(minfos[1,k],4)
    nmi <- signif(minfos[2,k],4)
    en <- signif(minfos[3,k],4)
    conden <- signif(minfos[4,k],4)
    matplot(x=data[[k]], y=data$log_RMSD, type='p', pch='.', col=paste0('#000000','88'),
            xlab=paste0(colnames(minfos)[k], ', H = ',en,' bit'),
            ylab=paste0('log-RMSD')
            )
    title(paste0(colnames(minfos)[k],
                         ', MI = ',mi,' bit, norm = ',nmi,', cond entr = ',conden, ' bit'))
}
dev.off()
##
pdff('histograms_data')
for(i in 1:ncol(origdata)){
    datum <- origdata[[i]]
    summa <- fivenum(datum)
    drange <- diff(range(datum))
    if(is.integer(datum)){
        breaks <- (summa[1]:(summa[5]+1))-0.5
    } else {
        width <- diff(summa[c(2,4)])/nbinsq
        nbins <- round(drange/width)
        breaks <- seq(summa[1]-drange/(nbins*100), summa[5]+drange/(nbins*100), length.out=nbins)
    }
    print(ggplot(origdata[,..i], aes_(x=as.name(names(origdata)[i]))) + geom_histogram(breaks=breaks))
}
dev.off()}
rm(origdata)
gc()
##
## Shuffle the data for training and test
set.seed(222)
data <- data[sample(1:nrow(data))]
##
fwrite(data,'../processed_data_scaled.csv', sep=' ')
##} else {
##
##data <- fread('../processed_data_scaled.csv', sep=' ')
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
##   user  system elapsed    0.08    0.03  371.17 
ndata <- 1 # nSamples = 37969
#set.seed(222)
seldata <- 1:ndata
rmsdCol <- which(names(data)=='bin_RMSD')
covNums <- which(colnames(data) %in%  c('scale_mcs_unbonded_polar_sasa', 'scale_ec_tanimoto_similarity', 'mcs_NumHeteroAtoms', # 'scale_fc_tanimoto_similarity'
                                        # 'docked_HeavyAtomCount', 
                                        'mcs_RingCount'))
covNames <- names(data)[covNums]
discreteCovs <- covNames[sapply(covNames, function(x){is.integer(data[[x]])})]
continuousCovs <- covNames[sapply(covNames, function(x){is.double(data[[x]])})]
allCovNums <- c(rmsdCol, covNums)
##
dimsC <- length(continuousCovs)
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
rmsdVals <- 1
##
rm(testmcall)
gc()
plan(sequential)
val <- 1
outfile <- paste0('_mcoutput',val)
        datamcr <- data[seldata]
        datamcr <- datamcr[bin_RMSD==val, covNames, with=F]
        for(i in discreteCovs){
            datum <- datamcr[[i]]
            datum <- datum[!is.na(datum)]
            levels <- as.numeric(names(table(data[[i]])))
            for(level in setdiff(min(levels):max(levels), as.numeric(names(table(datum))))){
                #print(paste0(val,' ',i,' ',level))
                dt <- data.table(level)
                names(dt) <- i
                datamcr <- rbind(datamcr, dt, fill=TRUE)
            }
        }
##
datamcr <- datamcr[-1]
##
regr <- profRegr(excludeY=TRUE, xModel='Mixed', nSweeps=50e3, nBurn=20e3, nFilter=50, data=as.data.frame(datamcr), nClusInit=80, covNames=c(discreteCovs,continuousCovs), discreteCovs=discreteCovs, continuousCovs=continuousCovs, nProgress=1000, seed=147, output=outfile, useHyperpriorR1=FALSE, useNormInvWishPrior=TRUE, hyper=testhp, alpha=4)
testmcall <- list()
testmcall[[1]] <- c(val=val,regr)
names(testmcall) <- paste0('bin',sapply(testmcall,function(i){i$val}))
##
##
## registerDoFuture()
## plan(multisession, workers = 4L)
##
sampledata <- as.list(rep(NA,length(testmcall)))
names(sampledata) <- names(testmcall)
##
for(val in rmsdVals){
    outfile <- paste0('_mcoutput',val)
    testmc <- testmcall[[paste0('bin',val)]]
##
fc <- file(paste0(outfile,"_logPost.txt"))
logPost <- sapply(strsplit(readLines(fc), " +"), as.numeric)
rownames(logPost) <- c('log-post','log-likelihood','log-prior')
close(fc)
fc <- file(paste0(outfile,'_nClusters.txt'))
nList <- sapply(strsplit(readLines(fc), " +"), as.integer)
close(fc)
fc <- file(paste0(outfile,'_alpha.txt'))
alphaList <- sapply(strsplit(readLines(fc), " +"), as.numeric)
close(fc)
##
fc <- file(paste0(outfile,'_psi.txt'))
psiList <- lapply(strsplit(readLines(fc), " +"), function(x){x <- as.numeric(x); x[x>=0]})
close(fc)
##                                        #
fc <- file(paste0(outfile,'_phi.txt'))
phiList <- lapply(strsplit(readLines(fc), " +"), function(x){x <- as.numeric(x); x[x>=0]})
close(fc)
##
fc <- file(paste0(outfile,'_mu.txt'))
muList <- lapply(strsplit(readLines(fc), " +"), as.numeric)
close(fc)
##
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
    ##
sampledata[[paste0('bin',val)]] <- list(val=val, nList=nList, alphaList=alphaList, psiList=psiList, phiList=phiList, muList=muList, sigmaList=sigmaList, logPost=logPost)
}
##
pdff('mcsummary')
for(j in 1:length(sampledata)){
    sd <- sampledata[[j]]
    matplot(sd$logPost[2,], type='l',ylim=range(sd$logPost[2,],na.rm=T,finite=T),ylab='log-likelihood',col=mypalette[j], main=paste0('bin = ',j))
    }
for(j in 1:length(sampledata)){
    sd <- sampledata[[j]]
matplot(sd$nList,type='l',ylim=range(sd$nList,na.rm=T,finite=T),ylab='no. clusters',col=mypalette[j], main=paste0('bin = ',j))
}
for(j in 1:length(sampledata)){
    sd <- sampledata[[j]]
    matplot(sd$alphaList,type='l',ylim=range(sd$alphaList,na.rm=T,finite=T),ylab='alpha',col=mypalette[j], main=paste0('bin = ',j))
    }
for(j in 1:length(sampledata)){
    sd <- sampledata[[j]]
    for(i in c(1,3)){matplot(sd$logPost[i,],type='l',ylim=range(sd$logPost[i,],na.rm=T,finite=T),col=mypalette[j], main=paste0('bin = ',j))}
    }
dev.off()
#plan(sequential)
save.image(file='_reverse_test.RData')
##
## Plots
##
##discrMin <- sapply(data[,discreteCovs,with=F], min)-1
dC <- discreteCovs
cC <- continuousCovs
plan(sequential)
#plan(multisession, workers = 6L)
priorP <- rep(1,3)/3
##
isasa <- cC[grepl('sasa', cC)]
itani <- cC[grepl('tanimoto', cC)]
##
jac <- function(X){
   Dtanimoto2y(X[itani])*Dsasa2y(X[isasa])
}
##
predictSasaTani <- function(dataobj, X){
    foreach(sample=seq_along(dataobj$nList), .combine='+', .inorder=FALSE)%dopar%{
        sum(dataobj$psiList[[sample]] *
            sapply(seq_len(dataobj$nList[sample]),function(cluster){
            ## prod(sapply(dC, function(covariate){
            ##     dataobj$phiList[[sample]][[covariate]][X[covariate], cluster]
            ## })) *
                dmvnorm(X, mean=dataobj$muList[[sample]][cC,cluster], sigma=as.matrix(dataobj$sigmaList[[sample]][cC,cC,cluster]))
        }))
       #drop(dataobj$phiList[[sample]]$bin_RMSD %*% weights)/sum(weights)
}/length(dataobj$nList)
}
predictSampleSasaTani <- function(dataobj,sample, X){
        sum(dataobj$psiList[[sample]] *
            sapply(seq_len(dataobj$nList[sample]),function(cluster){
            ## prod(sapply(dC, function(covariate){
            ##     dataobj$phiList[[sample]][[covariate]][X[covariate], cluster]
            ## })) *
                dmvnorm(X, mean=dataobj$muList[[sample]][cC,cluster], sigma=as.matrix(dataobj$sigmaList[[sample]][cC,cC,cluster]))
        }))
       #drop(dataobj$phiList[[sample]]$bin_RMSD %*% weights)/sum(weights)
}
##
##
## Plot sample densities
nplots <- 50
plan(sequential)
##plan(multisession, workers = 6L)
xgrid <- seq(1, 200,length.out=20)
ygrid <- seq(0.001,0.999,length.out=20)
pdf(file=paste0('testsampledensities_seq.pdf'),height=11.7,width=16.5*10)
## plot(NA,xlim=c(0,1),ylim=c(0,1),main='parameters')
## text(0,1,(paste0('transf = ','id','\ndim = ',dims,'\nmeanalpha = ',meanalpha,'\nsdalpha = ',sdalpha,'\nK = ',ka,'\nL = ',la,'\nM = ',mu,'\nresc = ',resc,'\nNclust = ',mean(clusters))),adj = c(0,1),cex=5)
##
for(sample in round(seq(1,1000,length.out=nplots))){
zgrid <- outer(xgrid,ygrid,function(x,y){
    mapply(function(xx,yy){
                    XX <- c(xx,yy)
                    names(XX) <- cC
                    YY <- c(sasa2y(xx),tanimoto2y(yy))
                    names(YY) <- cC
                    predictSampleSasaTani(sampledata[[1]], sample, YY)*jac(XX)
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
pdf(file=paste0('testsampledensities_scaled_seq.pdf'),height=11.7,width=16.5*10)
## plot(NA,xlim=c(0,1),ylim=c(0,1),main='parameters')
## text(0,1,(paste0('transf = ','id','\ndim = ',dims,'\nmeanalpha = ',meanalpha,'\nsdalpha = ',sdalpha,'\nK = ',ka,'\nL = ',la,'\nM = ',mu,'\nresc = ',resc,'\nNclust = ',mean(clusters))),adj = c(0,1),cex=5)
##
for(sample in round(seq(1,1000,length.out=nplots))){
zgrid <- outer(xgrid,ygrid,function(x,y){
    mapply(function(xx,yy){
                    XX <- c(xx,yy)
                    names(XX) <- cC
                    predictSampleSasaTani(sampledata[[1]], sample, XX)
                    },
        x,y)}
        )
persp(xgrid,ygrid,zgrid,zlim=c(0,max(zgrid)),ticktype='detailed',theta = 45, phi = 15,xlab='scaled sasa',ylab='scaled tanimoto')
}
dev.off()

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
                    names(XX) <- cC
                    YY <- c(sasa2y(xx),tanimoto2y(yy))
                    names(YY) <- cC
                    predictSasaTani(sampledata[[1]], YY)*jac(XX)
                    },
        x,y)}
        )
persp(xgrid,ygrid,zgrid,zlim=c(0,max(zgrid)),ticktype='detailed',theta = 45, phi = 15)
dev.off()


for(i in 1:100){
    zgrid <- outer(xgrid,ygrid,function(x,y){
        XX <- c(x,y)
        names(XX) <- cC
        predictSasaTani(sampledata[[1]], XX)*jac(X)}
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
for(sample in round(seq(1,1000,length.out=10^2))){
zgrid <- outer(xgrid,ygrid,function(x,y){
    mapply(function(xx,yy){
                    XX <- c(xx,yy)
                    names(XX) <- cC
                    YY <- c(log(xx),qnorm(yy))
                    names(YY) <- cC
                    predictSampleSasaTani(sampledata[[1]], sample,YY)*jac(XX)
                    },
        x,y)}
        )
persp(xgrid,ygrid,zgrid,zlim=c(0,max(zgrid)),ticktype='detailed')
}
dev.off()

for(i in 1:100){
    zgrid <- outer(xgrid,ygrid,function(x,y){
        XX <- c(x,y)
        names(XX) <- cC
        predictSasaTani(sampledata[[1]], XX)*jac(X)}
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
                    names(XX) <- cC
                    YY <- c(log(xx),qnorm(yy))
                    names(YY) <- cC
                    predictSampleSasaTani(sampledata[[1]], 3,YY)*jac(XX)
                    },
        x,y)}
        )
persp(xgrid,ygrid,zgrid,zlim=c(0,max(zgrid)),ticktype='detailed')
dev.off()





#
file.copy('testdensities.pdf',paste0('priordensities-','id_','D',dims,'_ma',meanalpha,'_sa',sdalpha,'_K',ka,'_L',la,'_M',mu,'_sc',resc,'_C',mean(clusters),'.pdf'))



predictYpar <- function(dataobj, X){
    foreach(sample=seq_along(dataobj$nList), .combine='+', .inorder=FALSE)%dopar%{
        sum(dataobj$psiList[[sample]] *
            sapply(seq_len(dataobj$nList[sample]),function(cluster){
            prod(sapply(dC, function(covariate){
                dataobj$phiList[[sample]][[covariate]][X[covariate], cluster]
            })) *
                dmvnorm(X[cC], mean=dataobj$muList[[sample]][cC,cluster], sigma=as.matrix(dataobj$sigmaList[[sample]][cC,cC,cluster]))
        }))
       #drop(dataobj$phiList[[sample]]$bin_RMSD %*% weights)/sum(weights)
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
    c(datum['bin_RMSD'], sapply(rmsdVals,function(val){predictYpar(sampledata[[val]],datum)}))
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
dC <- discreteCovs
cC <- continuousCovs
plan(sequential)
plan(multisession, workers = 6L)
priorP2 <- normalize(as.vector(table(data$bin_RMSD)))
##
predictYpar2 <- function(dataobj, X){
    foreach(sample=seq_along(dataobj$nList), .combine='+', .inorder=FALSE)%dopar%{
        sum(dataobj$psiList[[sample]] *
            sapply(seq_len(dataobj$nList[sample]),function(cluster){
            prod(sapply(dC, function(covariate){
                dataobj$phiList[[sample]][[covariate]][X[covariate], cluster]
            })) *
                dmvnorm(X[cC], mean=dataobj$muList[[sample]][cC,cluster], sigma=as.matrix(dataobj$sigmaList[[sample]][cC,cC,cluster]))
        }))
       #drop(dataobj$phiList[[sample]]$bin_RMSD %*% weights)/sum(weights)
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
    c(datum['bin_RMSD'], normalize(sapply(rmsdVals,function(val){predictYpar2(sampledata[[val]],datum)})))
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
dC <- discreteCovs
cC <- continuousCovs
plan(sequential)
plan(multisession, workers = 6L)
priorP <- rep(1,3)
##
predictYpar <- function(dataobj, x){
normalize(foreach(val=rmsdVals, .combine=c)%:%foreach(sample=seq_along(dataobj[[val]]$nList), .combine='+', .inorder=FALSE)%dopar%{
        sum(exp(log(dataobj[[val]]$psiList[[sample]]) + sapply(seq_len(dataobj[[val]]$nList[sample]),function(j){
            sum(log(sapply(dC, function(elem){
                dataobj[[val]]$phiList[[sample]][[elem]][x[elem], j]
            }))) +
                dmvnorm(x[cC], mean=dataobj[[val]]$muList[[sample]][cC,j], sigma=as.matrix(dataobj[[val]]$sigmaList[[sample]][cC,cC,j]), log=TRUE)
        }))) * priorP[val]
       #drop(dataobj[[val]]$phiList[[sample]]$bin_RMSD %*% weights)/sum(weights)
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
    c(datum['bin_RMSD'], (predictYpar(sampledata,datum)))
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
    foreach(sample=1:length(nList), .combine=rbind, .inorder=FALSE)%dopar%{
        dC <- mcobj$discreteCovs[mcobj$discreteCovs != 'bin_RMSD']
        cC <- mcobj$continuousCovs
        weights <- exp(log(psiList[[sample]]) + sapply(seq_len(nList[sample]),function(j){
            sum(log(mapply(function(xx,yy){xx[yy, j]}, phiList[[sample]][dC], x[dC]-discrMin[dC]))) +
                dmvnorm(x[cC], mean=muList[[sample]][cC,j], sigma=as.matrix(sigmaList[[sample]][cC,cC,j]), log=TRUE)
        }))
       drop(phiList[[sample]]$bin_RMSD %*% weights)/sum(weights)
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
            sum(log(sapply(dC, function(covariate){
                dataobj$phiList[[sample]][[covariate]][X[covariate], cluster]
            }))) +
                dmvnorm(X[cC], mean=dataobj$muList[[sample]][cC,cluster], sigma=as.matrix(dataobj$sigmaList[[sample]][cC,cC,cluster]), log=TRUE)
        })))
b <-         sum(
            (dataobj$psiList[[sample]]) *
            sapply(seq_len(dataobj$nList[sample]),function(cluster){
            prod((sapply(dC, function(covariate){
                dataobj$phiList[[sample]][[covariate]][X[covariate], cluster]
            }))) *
                dmvnorm(X[cC], mean=dataobj$muList[[sample]][cC,cluster], sigma=as.matrix(dataobj$sigmaList[[sample]][cC,cC,cluster]))
            }))
abs((a-b)/(a+b))})
summary(preci)


rm(preci)
gc()
system.time(preci <- sapply(rep(1:1000,20),function(sample){
 sum(exp(
            log(dataobj$psiList[[sample]]) +
            sapply(seq_len(dataobj$nList[sample]),function(cluster){
            sum(log(sapply(dC, function(covariate){
                dataobj$phiList[[sample]][[covariate]][X[covariate], cluster]
            }))) +
                dmvnorm(X[cC], mean=dataobj$muList[[sample]][cC,cluster], sigma=as.matrix(dataobj$sigmaList[[sample]][cC,cC,cluster]), log=TRUE)
            })))
 }))


rm(preci)
gc()
system.time(preci <- sapply(rep(1:1000,20),function(sample){
sum(
            (dataobj$psiList[[sample]]) *
            sapply(seq_len(dataobj$nList[sample]),function(cluster){
            prod((sapply(dC, function(covariate){
                dataobj$phiList[[sample]][[covariate]][X[covariate], cluster]
            }))) *
                dmvnorm(X[cC], mean=dataobj$muList[[sample]][cC,cluster], sigma=as.matrix(dataobj$sigmaList[[sample]][cC,cC,cluster]))
            }))
}))

abs((a-b)/(a+b))})
summary(preci)

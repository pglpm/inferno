## Author: PGL  Porta Mana
## Created: 2021-03-20T10:07:17+0100
## Last-Updated: 2021-07-03T22:54:07+0200
################
## Script for exploring Fabio's dataset
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

## Read data
rm(data)
data <- as.data.table(read.csv(file='dataset_template_based_docking_predict_rmsd.csv', header=T, sep=','))
data <- data[complete.cases(data),]
data <- data[, which(sapply(data, is.numeric)==TRUE), with=FALSE]
## doublecols <- which(sapply(data, function(x){all(is.double(x))}))
## posdoublecols <- which(sapply(data, function(x){all(is.double(x))&&all(x>0)}))
data$rmsd <- log(data$rmsd)
names(data)[which(names(data)=='rmsd')] <- 'logRMSD'
rmsdThreshold <- c(2, 2.5, 3)
logRmsdThreshold <- log(rmsdThreshold)
data <- data.table(binrmsd=as.integer(1-(data$logRMSD<logRmsdThreshold[1])+(data$logRMSD>logRmsdThreshold[3])), data)
nameFeatures <- names(data)
nSamples <- nrow(data)
nFeatures <- ncol(data)
colBinRmsd <- which(nameFeatures=='binrmsd')
colLogRmsd <- which(nameFeatures=='logrmsd')

##
nbinsq <- 5
##binsCheck <- numeric(length=length(okFeatures))*NA
pdff('histograms_fabiodata')
## i <- rmsdCol
## datum <- log(unlist(data[[i]]))
## summa <- fivenum(datum)
## ##binsCheck[i] <- nbins
## hist(x=datum, breaks=nbins, main=paste0('log-',nameFeatures[i]), xlab=paste0('log-',nameFeatures[i]))
breakFeatures <- list()
for(i in 1:length(data)){
    datum <- unlist(data[[i]])
    summa <- fivenum(datum)
    if(is.integer(datum)){
        breaks <- (summa[1]:(summa[5]+1))-0.5
    } else {
        width <- diff(summa[c(2,4)])/nbinsq
        nbins <- round(diff(summa[c(1,5)])/width)
        breaks <- seq(summa[1], summa[5]*1.01, length.out=nbins)
    }
    breakFeatures[[i]] <- breaks
    print(ggplot(data[, ..i], aes_(x=as.name(names(data)[i]))) + geom_histogram(breaks=breaks))
    #hist(x=datum, breaks=breaks, main=nameFeatures[i], xlab=nameFeatures[i])
    ## if(nameFeatures[i] == 'logrmsd'){
    ##     matlines(y=matrix(c(-1,5000),2,1), x=matrix(rep(logRmsdThreshold,2),2,3,byrow=T),
    ##              lty=c(2,3,2), lwd=3, col=c(myredpurple,myyellow,myredpurple))
    ## }
}
dev.off()

## ##
## nbinsq <- 5
## ##binsCheck <- numeric(length=length(okFeatures))*NA
## pdff('loghistograms_fabiodata')
## ## i <- rmsdCol
## ## datum <- log(unlist(data[[i]]))
## ## summa <- fivenum(datum)
## ## ##binsCheck[i] <- nbins
## ## hist(x=datum, breaks=nbins, main=paste0('log-',nameFeatures[i]), xlab=paste0('log-',nameFeatures[i]))
## breakFeatures <- list()
## for(i in 1:length(data)){
##     datum <- unlist(data[[i]])
##     if(is.integer(datum)){
##         datum <- log(datum+2)
##         summa <- fivenum(datum)
##                 width <- abs(diff(summa[c(2,4)]))/nbinsq
##         nbins <- round(abs(diff(summa[c(1,5)]))/width)
##         breaks <- seq(summa[1], summa[5], length.out=nbins)
## #        breaks <- (summa[1]:(summa[5]+1))-0.5
##     } else if(i==2){
##         summa <- fivenum(datum)
##         width <- diff(summa[c(2,4)])/nbinsq
##         nbins <- round(diff(summa[c(1,5)])/width)
##         breaks <- seq(summa[1], summa[5]*1.01, length.out=nbins)
##     } else {
##         datum <- log(datum+1)
##         summa <- fivenum(datum)
##         width <- abs(diff(summa[c(2,4)]))/nbinsq
##         nbins <- round(abs(diff(summa[c(1,5)]))/width)
##         breaks <- seq(summa[1], summa[5], length.out=nbins)
##     }
##     df <- data.frame(x=datum)
##     names(df) <- paste0("log-",names(data)[i])
##         print(ggplot(df, aes_(x=as.name(names(df)))) + geom_histogram(breaks=breaks))
##     rm(df)
##     ## if(nameFeatures[i] == 'logrmsd'){
##     ##     matlines(y=matrix(c(-1,5000),2,1), x=matrix(rep(logRmsdThreshold,2),2,3,byrow=T),
##     ##              lty=c(2,3,2), lwd=3, col=c(myredpurple,myyellow,myredpurple))
##     ## }
## }
## dev.off()



## ## Pairwise mutual infos
## mipairs <- matrix(NA,choose(length(okFeatures)+1,2),7)
## colnames(mipairs) <- c('i','j','I(i,j)','norm MI','H(i)','H(j)','H(j|i)')
## pair <- 0
## for(i in okFeatures){
##     for(j in (i+1):(totFeatures+1)){
##         pair <- pair+1
##         freqs <- normalize(bin2(x=cbind(data[[i]],data[[j]]),
##                                 ab=rbind(range(breakFeatures[[i]]),range(breakFeatures[[j]])),
##                                 nbin=c(length(breakFeatures[[i]]), length(breakFeatures[[j]]))-1 )$nc)
##         mi <- mutualinfo(freqs)
##         en1 <- entropy(rowSums(freqs))
##         en2 <- entropy(colSums(freqs))
##         conden21 <- condentropy21(freqs) 
##         mipairs[pair,] <- c(i,j,mi,mi/min(en1,en2),en1,en2,conden21)
##     }
## }
## mipairs <- mipairs[order(mipairs[,4], decreasing=TRUE),]
## ##
## ##

## Mutual infos with logRmsd
minfos <- matrix(NA,4,nFeatures)
rownames(minfos) <- c('MI','norm_MI','entropy','cond_entropy')
colnames(minfos) <- nameFeatures
rangeBinRmsd <- range(breakFeatures[[colBinRmsd]])
lengthBinRmsd <- length(breakFeatures[[colBinRmsd]])
for(i in 1:ncol(minfos)){
        freqs <- normalize(bin2(x=cbind(data$binrmsd,data[[i]]),
                                ab=rbind(rangeBinRmsd,range(breakFeatures[[i]])),
                                nbin=c(lengthBinRmsd, length(breakFeatures[[i]]))-1 )$nc)
        mi <- mutualinfo(freqs)
        en <- entropy(colSums(freqs))
        enrmsd <- entropy(rowSums(freqs))
        conden21 <- condentropy21(t(freqs)) 
        minfos[,i] <- c(mi, mi/min(en,enrmsd),en, conden21)
}
reorder <- order(minfos[1,], decreasing=TRUE)
minfos <- minfos[,reorder]
data <- data[, ..reorder]

##
pdff('plotsbinMI')
for(k in 1:ncol(minfos)){
    mi <- signif(minfos[1,k],4)
    nmi <- signif(minfos[2,k],4)
    en <- signif(minfos[3,k],4)
    conden <- signif(minfos[4,k],4)
    matplot(x=data[[k]], y=data$binrmsd, type='p', pch='.', col=paste0('#000000','88'),
            xlab=paste0(colnames(minfos)[k], ', H = ',en,' bit'),
            ylab=paste0('binned RMSD')
            )
    title(paste0(colnames(minfos)[k],
                         ', MI = ',mi,' bit, norm = ',nmi,', cond entr = ',conden, ' bit'))
}
dev.off()

##
pdff('plotslogMI')
for(k in 1:ncol(minfos)){
    mi <- signif(minfos[1,k],4)
    nmi <- signif(minfos[2,k],4)
    en <- signif(minfos[3,k],4)
    conden <- signif(minfos[4,k],4)
    matplot(x=data[[k]], y=data$logrmsd, type='p', pch='.', col=paste0('#000000','88'),
            xlab=paste0(colnames(minfos)[k], ', H = ',en,' bit'),
            ylab=paste0('log-RMSD')
            )
    title(paste0(colnames(minfos)[k],
                         ', MI = ',mi,' bit, norm = ',nmi,', cond entr = ',conden, ' bit'))
}
dev.off()


## ## With y-model
## ##   user  system elapsed  27.93    0.22   28.18 
## set.seed <- 147
## ndata <- round(nSamples/10)
## covNums <- c(1,3:8)
## seldata <- sort(sample(nrow(data), ndata, replace=FALSE))
## rm(datamc)
## datamc <- as.data.frame(data[seldata, ..covNums])
## yCov <- 'binrmsd'
## discreteCovs <- names(datamc[,which(names(datamc)!='binrmsd' & sapply(datamc, is.integer)==TRUE)])
## continuousCovs <- names(datamc[,which(names(datamc)!='binrmsd' & sapply(datamc, is.double)==TRUE)])
## #
## system.time(testmc <- profRegr(excludeY=FALSE, yModel='Categorical', xModel='Mixed', nSweeps=500, nBurn=100, nFilter=5, data=datamc, nClusInit=10, outcome='binrmsd', covNames=c(discreteCovs,continuousCovs), discreteCovs=discreteCovs, continuousCovs=continuousCovs, nProgress=100, seed=222, output='_mcoutput'))
## #
## pdff('cmsum_y_discr')
## print(globalParsTrace(testmc, parameters = "nClusters",plotBurnIn=F))
## print(globalParsTrace(testmc, parameters = "alpha",plotBurnIn=F))
## dev.off()


##################################################
## Mixed-x, no y-model
##    user  system elapsed    7.86    0.33    8.25 
ndata <- round(nSamples/5)
covNums <- c(1,3:6)
set.seed(147)
seldata <- sort(sample(nrow(data), ndata, replace=FALSE))
rm(datamc)
datamc <- data[seldata, ..covNums]
discreteCovs <- names(datamc)[which(sapply(datamc, is.integer)==TRUE)]
continuousCovs <- names(datamc)[which(sapply(datamc, is.double)==TRUE)]
for(i in discreteCovs){
    datum <- unlist(datamc[,i,with=FALSE])
    if(length(table(datum)) < diff(range(datum, na.rm=TRUE))+1){
        print(i)
        for(j in min(datum,na.rm=TRUE):max(datum,na.rm=TRUE)){
            print(j)
            dt <- data.table(j)
            names(dt) <- i
            datamc <- rbind(datamc, dt, fill=TRUE)
        }
    }
}


##
system.time(testmc <- profRegr(excludeY=TRUE, xModel='Mixed', nSweeps=1e4, nBurn=2e4, nFilter=10, data=as.data.frame(datamc), nClusInit=10, covNames=c(discreteCovs,continuousCovs), discreteCovs=discreteCovs, continuousCovs=continuousCovs, nProgress=100, seed=222, output='_mcoutput'))
##
fc <- file("_mcoutput_logPost.txt")
logPost <- sapply(strsplit(readLines(fc), " +"), as.numeric)
rownames(logPost) <- c('log-post','log-likelihood','log-prior')
close(fc)
fc <- file("_mcoutput_nClusters.txt")
nList <- sapply(strsplit(readLines(fc), " +"), as.integer)
close(fc)
fc <- file("_mcoutput_alpha.txt")
alphaList <- sapply(strsplit(readLines(fc), " +"), as.numeric)
close(fc)
##
pdff('mcsummary')
matplot(logPost[2,],type='l',ylim=range(logPost[2,],na.rm=T,finite=T),ylab='log-likelihood')
matplot(nList,type='l',ylim=range(nList,na.rm=T,finite=T),ylab='no. clusters')
matplot(alphaList,type='l',ylim=range(alphaList,na.rm=T,finite=T),ylab='alpha')
for(i in c(1,3)){matplot(logPost[i,],type='l',ylim=range(logPost[i,],na.rm=T,finite=T))}
dev.off()
##
fc <- file("_mcoutput_psi.txt")
psiList <- lapply(strsplit(readLines(fc), " +"), function(x){x <- as.numeric(x); x[x>=0]})
close(fc)
##                                        #
fc <- file("_mcoutput_phi.txt")
phiList <- lapply(strsplit(readLines(fc), " +"), function(x){x <- as.numeric(x); x[x>=0]})
close(fc)
##
fc <- file("_mcoutput_mu.txt")
muList <- lapply(strsplit(readLines(fc), " +"), as.numeric)
close(fc)
##
fc <- file("_mcoutput_Sigma.txt")
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

registerDoFuture()
plan(multisession, workers = 4L)

plan(sequential)

discrMin <- sapply(data[,testmc$discreteCovs,with=F], min)-1
##
predictYpar <- function(mcobj, x){
    foreach(sample=1:length(nList), .combine=rbind, .inorder=FALSE)%dopar%{
        dC <- mcobj$discreteCovs[mcobj$discreteCovs != 'binrmsd']
        cC <- mcobj$continuousCovs
        weights <- exp(log(psiList[[sample]]) + sapply(seq_len(nList[sample]),function(j){
            sum(log(sapply(dC, function(elem){
                phiList[[sample]][[elem]][x[elem]-discrMin[elem], j]
            }))) +
                dmvnorm(x[cC], mean=muList[[sample]][cC,j], sigma=as.matrix(sigmaList[[sample]][cC,cC,j]), log=TRUE)
        }))
       drop(phiList[[sample]]$binrmsd %*% weights)/sum(weights)
}
}
##
unseldata <- setdiff(1:nrow(data), seldata)
##
nTest <- 100
testdata <- data[unseldata,testmc$covNames, with=F]
testdata <- rbind(head(testdata[binrmsd==0],n=nTest), head(testdata[binrmsd==1],n=nTest), head(testdata[binrmsd==2],n=nTest))
##
## user  system elapsed    1.71    0.14    4.31 
system.time(testres <- t(sapply(seq_len(nrow(testdata)), function(datum){
    datum <- unlist(testdata[datum])
    c(datum['binrmsd']-discrMin['binrmsd'], colMeans(predictYpar(testmc,datum)))
})))
##
sum(sapply(1:nrow(testres), function(i){testres[i,1]==which.max(testres[i,-1])}))/nrow(testres)


predictYpar2 <- function(mcobj, x){
    foreach(sample=1:length(nList), .combine=rbind, .inorder=FALSE)%dopar%{
        dC <- mcobj$discreteCovs[mcobj$discreteCovs != 'binrmsd']
        cC <- mcobj$continuousCovs
        weights <- exp(log(psiList[[sample]]) + sapply(seq_len(nList[sample]),function(j){
            sum(log(mapply(function(xx,yy){xx[yy, j]}, phiList[[sample]][dC], x[dC]-discrMin[dC]))) +
                dmvnorm(x[cC], mean=muList[[sample]][cC,j], sigma=as.matrix(sigmaList[[sample]][cC,cC,j]), log=TRUE)
        }))
       drop(phiList[[sample]]$binrmsd %*% weights)/sum(weights)
}
}
##
unseldata <- setdiff(1:nrow(data), seldata)
##
nTest <- 100
testdata <- data[unseldata,testmc$covNames, with=F]
testdata <- rbind(head(testdata[binrmsd==0],n=nTest), head(testdata[binrmsd==1],n=nTest), head(testdata[binrmsd==2],n=nTest))
##
## user  system elapsed    1.71    0.14    4.31 
system.time(testres <- t(sapply(seq_len(nrow(testdata)), function(datum){
    datum <- unlist(testdata[datum])
    c(datum['binrmsd']-discrMin['binrmsd'], colMeans(predictYpar2(testmc,datum)))
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
system.time(testmc <- profRegr(excludeY=TRUE, xModel='Mixed', nSweeps=10000, nBurn=2000, nFilter=10, data=datamc, nClusInit=10, covNames=c(discreteCovs,continuousCovs), discreteCovs=discreteCovs, continuousCovs=continuousCovs, nProgress=100, seed=222, output='_mcoutput'))
##
fc <- file("_mcoutput_logPost.txt")
logPost <- sapply(strsplit(readLines(fc), " +"), as.numeric)
close(fc)
##
pdff('cmsum_cont')
print(globalParsTrace(testmc, parameters = "nClusters",plotBurnIn=F))
print(globalParsTrace(testmc, parameters = "alpha",plotBurnIn=F))
for(i in 1:3){matplot(logPost[i,],type='l')}
dev.off()
##
fc <- file("_mcoutput_nClusters.txt")
nList <- sapply(strsplit(readLines(fc), " +"), as.numeric)
close(fc)
#
fc <- file("_mcoutput_psi.txt")
psiList <- lapply(strsplit(readLines(fc), " +"), function(x){y <- as.numeric(x);y[y>=0]})
close(fc)
##                                        #
fc <- file("_mcoutput_phi.txt")
phiList <- lapply(strsplit(readLines(fc), " +"), function(x){y <- as.numeric(x);y[y>=0]})
close(fc)
##
fc <- file("_mcoutput_mu.txt")
muList <- lapply(strsplit(readLines(fc), " +"), as.numeric)
close(fc)
##
fc <- file("_mcoutput_Sigma.txt")
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
yCov <- 'binrmsd'
discreteCovs <- names(datamc)[which(sapply(datamc, is.integer)==TRUE)]
continuousCovs <- names(datamc)[which(sapply(datamc, is.double)==TRUE)]
##
system.time(testmc <- profRegr(excludeY=FALSE, yModel='Categorical', xModel='Normal', nSweeps=2000, nBurn=10, nFilter=2, data=datamc, nClusInit=2, outcome=yCov, covNames=c(discreteCovs,continuousCovs), discreteCovs=discreteCovs, continuousCovs=continuousCovs, nProgress=100, seed=222, output='_mcoutput'))
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


## discreteCovs <- names(datamc[,which(names(datamc)!='binrmsd' & sapply(datamc, is.integer)==TRUE)])
## continuousCovs <- names(datamc[,which(names(datamc)!='binrmsd' & sapply(datamc, is.double)==TRUE)])
## #
## system.time(testmc <- profRegr(excludeY=FALSE, yModel='Categorical', xModel='Mixed', nSweeps=500, nBurn=100, nFilter=5, data=datamc, nClusInit=10, outcome='binrmsd', covNames=c(discreteCovs,continuousCovs), discreteCovs=discreteCovs, continuousCovs=continuousCovs, nProgress=100, seed=222, output='_mcoutput'))
## #
## pdff('cmsum_y_discr')
## print(globalParsTrace(testmc, parameters = "nClusters",plotBurnIn=F))
## print(globalParsTrace(testmc, parameters = "alpha",plotBurnIn=F))
## dev.off()


### Updated up to above ###









misubdata <- mipairs[which(mipairs[,2]==(nFeatures+1)),]
misubdata <- misubdata[order(misubdata[,7], decreasing=FALSE),]
rmsdData <- data[[rmsdCol]]
pdff('jointplotslogRMSD_fabiodata')
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
samplescale <- max(foreach(i=(-1):1, .combine=c)%do%{length(which(data$binrmsd==i))})/subsamplesize

plotpoints <- foreach(i=(-1):1)%do%{
    matr <- as.matrix(data[which(data$binrmsd==i),..feats])
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
pdff('jointplotsbinRMSD_fabiodata')
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

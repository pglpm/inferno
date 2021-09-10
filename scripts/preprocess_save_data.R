## Author: PGL  Porta Mana
## Created: 2021-03-20T10:07:17+0100
## Last-Updated: 2021-09-10T21:57:19+0200
################
## Script for reverse regression
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

doplots <- FALSE#TRUE
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
##
tanimoto2y <- function(x){
    -log(1/x-1) ## faster than qlogis(x, scale=1/2)
}
##
## correct compared with study_prior_bayes_regr.nb
Dtanimoto2y <- function(x){
    1/(x*(1-x))
}
##
## y2tanimoto <- function(y){
##     1/2 + atan(y)/pi
## }
##
sasa2y <- function(x){
   log(x)
}
##
## correct compared with study_prior_bayes_regr.nb
Dsasa2y <- function(x){
    1/(x)
}
##
## Read and reorganize data
rm(data)
data <- as.data.table(read.csv(file='dataset_template_based_docking_predict_rmsd.csv', header=T, sep=','))
## remove datapoints with missing or erroneous features
data <- data[!is.na(rmsd)]
data <- data[, which(sapply(data, is.numeric)==TRUE), with=FALSE]
minusfeatures <- which(apply(data,2,min)==-1)
for(feat in minusfeatures){
    data <- data[!(data[[feat]]==-1)]
}
## copy of data before rescaling, to plot histograms in the orig. scale
origdata <- data
## add log-RMSD
data$log_RMSD <- log(data$rmsd)
rmsdThreshold <- c(2, 2.5, 3)
logRmsdThreshold <- log(rmsdThreshold)
##
## Add column with three, binned RMSD values
data$bin_RMSD <- as.integer(1+(data$log_RMSD>logRmsdThreshold[1])+(data$log_RMSD>logRmsdThreshold[3]))
##
## find log-'sasa' features
indxs <- grepl('sasa', colnames(data))
## add rescaled 'sasa' features
for(elem in colnames(data)[indxs]){
    datum <- data[[elem]]
    datum <- datum/signif(mean(datum),1)
    data[[paste0('scale_',elem)]] <- datum
## add 'sasa' feature in standardized log-scale
    datum <- sasa2y(datum)
    datum <- datum/signif(sd(datum[abs(datum)!=Inf]),1)
    eps <- max(diff(sort(unique(datum[abs(datum)!=Inf]))))
    datum[datum==-Inf] <- min(datum[abs(datum)!=Inf]) - 2 * eps
    data[[paste0('log_',elem)]] <- datum
}
## add 'tanimoto' features in logit-scale
indxt <- grepl('tanimoto', colnames(data))
for(elem in colnames(data)[indxt]){
    datum <- tanimoto2y(data[[elem]])
    datum <- datum/signif(sd(datum[abs(datum)!=Inf]),1)
    eps <- max(diff(sort(unique(datum[abs(datum)!=Inf]))))
    datum[datum==Inf] <- max(datum[abs(datum)!=Inf])+2*eps
    datum[datum==-Inf] <- min(datum[abs(datum)!=Inf])-2*eps
    data[[paste0('logit_',elem)]] <- datum
}
## add shifted integer features to start from value 1
indx <- sapply(1:ncol(data), function(x){is.integer(data[[x]])})
for(elem in colnames(data)[indx]){
    datum <- data[[elem]]
    datum <- datum - min(datum, na.rm=TRUE) + 1L
    data[[paste0('shift_',elem)]] <- datum
}
##

nameFeatures <- names(data)
nSamples <- nrow(data)
nFeatures <- ncol(data)
##
## Format bins to calculate mutual info
nbinsq <- 10
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
## Reorder features according to mutual info with RMSD
reorder <- order(minfos[1,], decreasing=TRUE)
minfos <- minfos[,reorder]
data <- data[, ..reorder]
breakFeatures <- breakFeatures[reorder]
##origdata <- origdata[, setdiff(reorder-1,0), with=FALSE]
##
##
## Plots
if(doplots==TRUE){
    pdff(paste0('histograms_data_transf_scaled_',nbinsq,'bins'))
    for(i in 1:ncol(data)){
        datum <- data[[i]]
        breaks <- breakFeatures[[i]]
        print(ggplot(data[,..i], aes_(x=as.name(names(data)[i]))) + geom_histogram(breaks=breaks))
    }
    dev.off()
    ##
    pdff(paste0('plotslogMI_transf_scaled_',nbinsq,'bins'))
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
}
##

if(doplots==TRUE){
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

##
fwrite(data,'data_processed_transformed_rescaled.csv', sep=' ')
## Shuffle the data for training and test
set.seed(222)
data <- data[sample(1:nrow(data))]
fwrite(data,'data_processed_transformed_rescaled_shuffled.csv', sep=' ')

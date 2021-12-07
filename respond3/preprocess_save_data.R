## Author: PGL  Porta Mana
## Created: 2021-03-20T10:07:17+0100
## Last-Updated: 2021-10-28T07:44:09+0200
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
## library('extraDistr')
## library('mvtnorm')
options(bitmapType='cairo')
pdff <- function(filename){pdf(file=paste0(filename,'.pdf'),paper='a4r',height=11.7,width=16.5)} # to output in pdf format
pngf <- function(filename,res=300){png(file=paste0(filename,'.png'),height=11.7*1.2,width=16.5,units='in',res=res,pointsize=36)} # to output in pdf format
#### End custom setup ####

doplots <- TRUE
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
## normalizem <- function(freqs){freqs/rowSums(freqs)}
##
##
tanimoto2y <- function(x, lambda=1){
    (x^lambda - (1-x)^lambda) * 2^lambda/(4*lambda)
}
##
## correct compared with study_prior_bayes_regr.nb
## Dtanimoto2y <- function(x){
##     1/(x*(1-x))
## }
##
## y2tanimoto <- function(y){
##     1/2 + atan(y)/pi
## }
##
sasa2y <- function(x, lambda=1){
   (x^lambda - 1)/lambda
}
##
## correct compared with study_prior_bayes_regr.nb
## Dsasa2y <- function(x){
##     1/(x)
## }
##
## Read and reorganize data
rm(alldata)
alldata <- as.data.table(read.csv(file='dataset_template_based_docking_predict_rmsd.csv', header=T, sep=','))
## remove datapoints with missing or erroneous features
alldata <- alldata[!is.na(rmsd)]
## alldata <- alldata[, which(sapply(alldata, is.numeric)==TRUE), with=FALSE]
minusfeatures <- which(apply(alldata,2,min)==-1)
for(feat in minusfeatures){
    alldata <- alldata[!(alldata[[feat]]==-1)]
}
## copy of data before rescaling, to plot histograms in the orig. scale
origdata <- alldata
## add log-RMSD
alldata$log_RMSD <- log(alldata$rmsd)
rmsdThreshold <- c(2, 2.5, 3)
logRmsdThreshold <- log(rmsdThreshold)
##
## Add column with three, binned RMSD values
alldata$bin_RMSD <- as.integer(1+(alldata$log_RMSD>logRmsdThreshold[1])+(alldata$log_RMSD>logRmsdThreshold[3]))
##
source('nimbleMCMC_hpc/functions_rmsdregr_nimble_binom.R')
outfile <- file('data_transformation_parameters.txt', 'wb')
write(x='Transformation parameters', file=outfile)
## find 'sasa' features
indxs <- grepl('sasa', colnames(alldata))
## add 'sasa' features in transformed scale
for(elem in colnames(alldata)[indxs]){
    write(x='', file=outfile, append=TRUE)
    write(x=paste0(elem, ':'), file=outfile, append=TRUE)
    datum <- alldata[[elem]]
    optimf <- function(par){
        lambda <- exp(par[1])
        mu <- par[2]
        si <- exp(par[3])
        -sum(dnorm(x=sasa2y(datum, lambda), mean=mu, sd=si, log=TRUE))
    }
    optout <- myoptim(c(0, 0, 0), optimf)
    if(optout$convergence!=0){print(optout)}
    loglambda <- round(log2(exp(optout$par[1])))
    write(x=paste0('loglambda: ', loglambda), file=outfile, append=TRUE)
    ## transform
    datum <- sasa2y(datum, lambda=2^loglambda)
    ## centre and standardize for numerical efficiency
    datummean <- signif(mean(datum), 1)
    write(x=paste0('mean: ', datummean), file=outfile, append=TRUE)
    datumsd <- signif(sd(datum), 1)
    write(x=paste0('SD: ', datumsd), file=outfile, append=TRUE)
    datum <- (datum - datummean)/datumsd
    alldata[[paste0('Xtransf_', elem)]] <- datum
}
## find 'tanimoto' features
indxt <- grepl('tanimoto', colnames(alldata))
## add 'tanimoto' features in transformed scale
for(elem in colnames(alldata)[indxt]){
    write(x='', file=outfile, append=TRUE)
    write(x=paste0(elem, ':'), file=outfile, append=TRUE)
    datum <- alldata[[elem]]
    optimf <- function(par){
        lambda <- exp(par[1])
        mu <- par[2]
        si <- exp(par[3])
        -sum(dnorm(x=tanimoto2y(datum, lambda), mean=mu, sd=si, log=TRUE))
    }
    optout <- myoptim(c(0, 0, 0), optimf)
    if(optout$convergence!=0){print(optout)}
    loglambda <- round(log2(exp(optout$par[1])))
    write(x=paste0('loglambda: ', loglambda), file=outfile, append=TRUE)
    ## transform
    datum <- tanimoto2y(datum, lambda=2^loglambda)
    ## centre and standardize for numerical efficiency
    datummean <- signif(mean(datum), 1)
    write(x=paste0('mean: ', datummean), file=outfile, append=TRUE)
    datumsd <- signif(sd(datum), 1)
    write(x=paste0('SD: ', datumsd), file=outfile, append=TRUE)
    datum <- (datum - datummean)/datumsd
    alldata[[paste0('Xtransf_', elem)]] <- datum
}
## ## add shifted integer features to start from value 1
## indx <- sapply(1:ncol(alldata), function(x){is.integer(alldata[[x]])})
## for(elem in colnames(alldata)[indx]){
##     datum <- alldata[[elem]]
##     datum <- datum - min(datum, na.rm=TRUE) + 1L
##     alldata[[paste0('shift_',elem)]] <- datum
## }
## ##
close(outfile)

## alldata <- alldata[, which(sapply(alldata, is.numeric)==TRUE), with=FALSE]
nameFeatures <- names(which(sapply(alldata, is.numeric)==TRUE))
nSamples <- nrow(alldata)
nFeatures <- length(nameFeatures)
##
## Format bins to calculate mutual info
nbinsq <- 16
##
breakFeatures <- list()
for(i in nameFeatures){
    datum <- alldata[[i]]
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
names(breakFeatures) <- nameFeatures
##
##
## Mutual infos with RMSD
minfos <- matrix(NA,4,nFeatures)
rownames(minfos) <- c('MI','norm_MI','entropy','cond_entropy')
yVar <- 'log_RMSD'
##yVar <- 'bin_RMSD'
colRmsd <- which(names(alldata)==yVar)
rangeRmsd <- range(breakFeatures[[yVar]])
colnames(minfos) <- nameFeatures
for(i in nameFeatures){
        freqs <- normalize(bin2(x=cbind(alldata[[yVar]],alldata[[i]]),
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
##origdata <- origdata[, setdiff(reorder-1,0), with=FALSE]
##
##
## Plots
if(doplots==TRUE){
    pdff(paste0('histograms_data_transf_scaled_',nbinsq,'bins'))
    for(i in nameFeatures[reorder]){
        datum <- alldata[[i]]
        breaks <- breakFeatures[[i]]
        print(ggplot(alldata[,..i], aes_(x=as.name(i))) + geom_histogram(breaks=breaks))
    }
    dev.off()
    ##
    pdff(paste0('plotslogMI_transf_scaled_',nbinsq,'bins'))
    for(k in nameFeatures[reorder]){
        mi <- signif(minfos[1,k],4)
        nmi <- signif(minfos[2,k],4)
        en <- signif(minfos[3,k],4)
        conden <- signif(minfos[4,k],4)
        matplot(x=alldata[[k]], y=alldata$log_RMSD, type='p', pch='.', col=paste0('#000000','88'),
                xlab=paste0(k, ', H = ',en,' bit'),
                ylab=paste0('log-RMSD')
                )
        title(paste0(k,
                     ', MI = ',mi,' bit, norm = ',nmi,', cond entr = ',conden, ' bit'))
    }
    dev.off()
}
##

if(doplots==TRUE){
    pdff('histograms_data')
    for(i in intersect(nameFeatures, colnames(origdata))){
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
        print(ggplot(origdata[,..i], aes_(x=as.name(i))) + geom_histogram(breaks=breaks))
    }
    dev.off()}
## rm(origdata)
gc()
##

##
fwrite(alldata,'data_id_processed_transformed.csv', sep=',')
## Shuffle the data for training and test
set.seed(222)
alldata <- alldata[sample(1:nrow(alldata))]
fwrite(alldata,'data_id_processed_transformed_shuffled.csv', sep=',')


## > nameFeatures[reorder]
##  [1] "log_RMSD"                         "rmsd"                            
##  [3] "bin_RMSD"                         "ec_tanimoto_similarity"          
##  [5] "Xtransf_ec_tanimoto_similarity"   "Xtransf_mcs_unbonded_polar_sasa" 
##  [7] "mcs_unbonded_polar_sasa"          "fc_tanimoto_similarity"          
##  [9] "Xtransf_fc_tanimoto_similarity"   "Xtransf_mcs_unbonded_apolar_sasa"
## [11] "Xtransf_sasa_bonded_apolar"       "Xtransf_sasa_unbonded_apolar"    
## [13] "docked_HeavyAtomCount"            "mcs_bonded_polar_sasa"           
## [15] "mcs_unbonded_apolar_sasa"         "sasa_unbonded_apolar"            
## [17] "mcs_NumHeteroAtoms"               "docked_NumRotatableBonds"        
## [19] "sasa_bonded_apolar"               "Xtransf_mcs_bonded_apolar_sasa"  
## [21] "mcs_RingCount"                    "mcs_NOCount"                     
## [23] "Xtransf_sasa_unbonded_polar"      "mcs_HeavyAtomCount"              
## [25] "sasa_unbonded_polar"              "mcs_bonded_apolar_sasa"          
## [27] "sasa_bonded_polar"                "Xtransf_sasa_bonded_polar"       
## [29] "mcs_template_NumHAcceptors"       "mcs_docked_NumHAcceptors"        
## [31] "template_HeavyAtomCount"          "mcs_docked_NumHDonors"           
## [33] "mcs_docked_NHOHCount"             "mcs_template_NHOHCount"          
## [35] "mcs_template_NumHDonors"          "template_NumRotatableBonds"      
## [37] "Xtransf_mcs_bonded_polar_sasa"    "docked_NumHeteroAtoms"           
## [39] "docked_RingCount"                 "docked_NumHAcceptors"            
## [41] "docked_NOCount"                   "mcs_template_NumRotatableBonds"  
## [43] "mcs_docked_NumRotatableBonds"     "template_NumHeteroAtoms"         
## [45] "template_NOCount"                 "docked_NHOHCount"                
## [47] "template_RingCount"               "template_NHOHCount"              
## [49] "template_NumHAcceptors"           "docked_NumHDonors"               
## [51] "template_NumHDonors"             

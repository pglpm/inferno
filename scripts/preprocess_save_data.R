## Author: PGL  Porta Mana
## Created: 2021-03-20T10:07:17+0100
## Last-Updated: 2021-09-29T11:09:42+0200
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
## find log-'sasa' features
indxs <- grepl('sasa', colnames(alldata))
## add rescaled 'sasa' features
for(elem in colnames(alldata)[indxs]){
    datum <- alldata[[elem]]
    datum <- datum/signif(mean(datum),1)
    if(min(datum) == 0){
        datum <- datum + diff(sort(unique(datum))[1:2])/100
    }
    alldata[[paste0('scale_',elem)]] <- datum
## add 'sasa' feature in standardized log-scale
    datum <- sasa2y(datum)
    ## eps <- max(diff(sort(unique(datum[abs(datum)!=Inf]))))
    ## datum[datum==-Inf] <- min(datum[abs(datum)!=Inf]) - 2 * eps
    datum <- datum/signif(sd(datum),1)
    alldata[[paste0('log_',elem)]] <- datum
}
## add 'tanimoto' features in logit-scale
indxt <- grepl('tanimoto', colnames(alldata))
for(elem in colnames(alldata)[indxt]){
    datum <- alldata[[elem]]
    eps <- min(c(1-datum[datum < 1], datum[datum > 0]))/100
    datum <- tanimoto2y((2*datum - 1)*(0.5-eps) + 0.5)
    ## eps <- max(diff(sort(unique(datum[abs(datum)!=Inf]))))
    ## datum[datum==Inf] <- max(datum[abs(datum)!=Inf])+2*eps
    ## datum[datum==-Inf] <- min(datum[abs(datum)!=Inf])-2*eps
    datum <- datum/signif(sd(datum[abs(datum)!=Inf]),1)
    alldata[[paste0('logit_',elem)]] <- datum
}
## add shifted integer features to start from value 1
indx <- sapply(1:ncol(alldata), function(x){is.integer(alldata[[x]])})
for(elem in colnames(alldata)[indx]){
    datum <- alldata[[elem]]
    datum <- datum - min(datum, na.rm=TRUE) + 1L
    alldata[[paste0('shift_',elem)]] <- datum
}
##


## alldata <- alldata[, which(sapply(alldata, is.numeric)==TRUE), with=FALSE]
nameFeatures <- names(which(sapply(alldata, is.numeric)==TRUE))
nSamples <- nrow(alldata)
nFeatures <- length(nameFeatures)
##
## Format bins to calculate mutual info
nbinsq <- 10
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
    for(i in nameFeatures){
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
rm(origdata)
gc()
##

##
fwrite(alldata,'data_id_processed_transformed_rescaled.csv', sep=',')
## Shuffle the data for training and test
set.seed(222)
alldata <- alldata[sample(1:nrow(alldata))]
fwrite(alldata,'data_id_processed_transformed_rescaled_shuffled.csv', sep=',')


## > nameFeatures[reorder]
##  [1] "log_RMSD"                            
##  [2] "rmsd"                                
##  [3] "bin_RMSD"                            
##  [4] "shift_bin_RMSD"                      
##  [5] "log_mcs_unbonded_polar_sasa"         
##  [6] "mcs_unbonded_polar_sasa"             
##  [7] "scale_mcs_unbonded_polar_sasa"       
##  [8] "ec_tanimoto_similarity"              
##  [9] "logit_ec_tanimoto_similarity"        
## [10] "logit_fc_tanimoto_similarity"        
## [11] "fc_tanimoto_similarity"              
## [12] "docked_HeavyAtomCount"               
## [13] "shift_docked_HeavyAtomCount"         
## [14] "mcs_NumHeteroAtoms"                  
## [15] "shift_mcs_NumHeteroAtoms"            
## [16] "docked_NumRotatableBonds"            
## [17] "shift_docked_NumRotatableBonds"      
## [18] "mcs_RingCount"                       
## [19] "shift_mcs_RingCount"                 
## [20] "mcs_NOCount"                         
## [21] "shift_mcs_NOCount"                   
## [22] "log_mcs_unbonded_apolar_sasa"        
## [23] "mcs_HeavyAtomCount"                  
## [24] "shift_mcs_HeavyAtomCount"            
## [25] "log_sasa_bonded_apolar"              
## [26] "log_sasa_unbonded_apolar"            
## [27] "mcs_unbonded_apolar_sasa"            
## [28] "scale_mcs_unbonded_apolar_sasa"      
## [29] "sasa_unbonded_apolar"                
## [30] "scale_sasa_unbonded_apolar"          
## [31] "mcs_template_NumHAcceptors"          
## [32] "shift_mcs_template_NumHAcceptors"    
## [33] "mcs_docked_NumHAcceptors"            
## [34] "shift_mcs_docked_NumHAcceptors"      
## [35] "log_mcs_bonded_apolar_sasa"          
## [36] "mcs_docked_NumHDonors"               
## [37] "shift_mcs_docked_NumHDonors"         
## [38] "mcs_docked_NHOHCount"                
## [39] "shift_mcs_docked_NHOHCount"          
## [40] "sasa_bonded_apolar"                  
## [41] "scale_sasa_bonded_apolar"            
## [42] "mcs_bonded_polar_sasa"               
## [43] "scale_mcs_bonded_polar_sasa"         
## [44] "mcs_template_NHOHCount"              
## [45] "shift_mcs_template_NHOHCount"        
## [46] "mcs_template_NumHDonors"             
## [47] "shift_mcs_template_NumHDonors"       
## [48] "template_HeavyAtomCount"             
## [49] "shift_template_HeavyAtomCount"       
## [50] "log_sasa_unbonded_polar"             
## [51] "mcs_bonded_apolar_sasa"              
## [52] "scale_mcs_bonded_apolar_sasa"        
## [53] "template_NumRotatableBonds"          
## [54] "shift_template_NumRotatableBonds"    
## [55] "sasa_unbonded_polar"                 
## [56] "scale_sasa_unbonded_polar"           
## [57] "log_sasa_bonded_polar"               
## [58] "sasa_bonded_polar"                   
## [59] "scale_sasa_bonded_polar"             
## [60] "docked_RingCount"                    
## [61] "shift_docked_RingCount"              
## [62] "docked_NumHAcceptors"                
## [63] "shift_docked_NumHAcceptors"          
## [64] "docked_NumHeteroAtoms"               
## [65] "shift_docked_NumHeteroAtoms"         
## [66] "mcs_template_NumRotatableBonds"      
## [67] "shift_mcs_template_NumRotatableBonds"
## [68] "mcs_docked_NumRotatableBonds"        
## [69] "shift_mcs_docked_NumRotatableBonds"  
## [70] "docked_NOCount"                      
## [71] "shift_docked_NOCount"                
## [72] "template_RingCount"                  
## [73] "shift_template_RingCount"            
## [74] "template_NumHeteroAtoms"             
## [75] "shift_template_NumHeteroAtoms"       
## [76] "docked_NHOHCount"                    
## [77] "shift_docked_NHOHCount"              
## [78] "template_NOCount"                    
## [79] "shift_template_NOCount"              
## [80] "template_NHOHCount"                  
## [81] "shift_template_NHOHCount"            
## [82] "template_NumHAcceptors"              
## [83] "shift_template_NumHAcceptors"        
## [84] "docked_NumHDonors"                   
## [85] "shift_docked_NumHDonors"             
## [86] "log_mcs_bonded_polar_sasa"           
## [87] "template_NumHDonors"                 
## [88] "shift_template_NumHDonors"           

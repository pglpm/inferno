## Author: PGL  Porta Mana
## Created: 2022-10-07T12:13:20+0200
## Last-Updated: 2023-01-18T14:53:44+0100
################
## Combine multiple Monte Carlo chains
################
if(!exists('tplot')){source('~/work/pglpm_plotfunctions.R')}

#rm(list=ls())

outputdir <- 'inferencep4'
extratitle <- 'longrunprobs'
## totsamples <- 4096L
datafile <- 'ingrid_data_nogds6.csv'
predictorfile <- 'predictors.csv'
predictandfile <- NULL # 'predictors.csv'
functionsfile <- 'functionsmcmc_2212120902.R'
showdata <- TRUE # 'histogram' 'scatter' FALSE
plotmeans <- TRUE
quantilestoshow <- c(1,7)/8# c(1,31)/32
ndata <- NULL # set this if you want to use fewer data
shuffledata <- FALSE # useful if subsetting data
family <- 'Palatino'


## load customized plot functions
if(!exists('tplot')){source('~/work/pglpm_plotfunctions.R')}
library('data.table')
library('png')
library('foreach')
library('doFuture')
library('doRNG')
registerDoFuture()
print('availableCores:')
print(availableCores())
print('availableCores-multicore:')
print(availableCores('multicore'))
if(Sys.info()['nodename']=='luca-HP-Z2-G9'){
    ncores <- 20}else{
    ncores <- 6}
print(paste0('using ',ncores,' cores'))
if(ncores>1){
    if(.Platform$OS.type=='unix'){
        plan(multicore, workers=ncores)
    }else{
        plan(multisession, workers=ncores)
    }
}else{
    plan(sequential)
}

#### Load files with variate information and data ####
if(!exists('outputdir') || is.na(outputdir) || is.null(outputdir)){
    outputdir <- as.character(commandArgs(trailingOnly=TRUE))[1]
    if(is.na(outputdir)){
        outputdir <- './'
    }
}
outputdir <- paste0(sub('(.+)/','\\1',outputdir),'/')
setwd(outputdir)
origdir <- '../'
source(paste0(origdir,functionsfile)) # load functions for post-MCMC

varinfofile <- list.files(pattern='^_varinfo.*\\.rds')
if(length(varinfofile) !=1){stop('Problems with varinfo file')}
varinfo <- readRDS(varinfofile)
##
variate <- lapply(variatetypes, function(x)names(varinfo[['type']])[varinfo[['type']]==x])
len <- lapply(variate,length)
names(variate) <- names(len) <- variatetypes

variatenames <- unlist(variate)
if(!is.null(predictorfile)){predictorfile <- paste0(origdir,predictorfile) }
if(!is.null(predictandfile)){predictandfile <- paste0(origdir,predictandfile) }
if(!is.null(predictorfile) && !is.null(predictandfile)){
    predictors <- as.vector(unlist(read.csv(predictorfile, header=F)))
    predictands <- as.vector(unlist(read.csv(predictandfile, header=F)))
}else if(!is.null(predictorfile) && is.null(predictandfile)){
    predictors <- as.vector(unlist(read.csv(predictorfile, header=F)))
    predictands <- setdiff(unlist(variate), predictors)
}else if(is.null(predictorfile) && !is.null(predictandfile)){
    predictands <- as.vector(unlist(read.csv(predictandfile, header=F)))
    predictors <- setdiff(unlist(variate), predictands)
}else{warning('predictors and predictands both missing')}

mcsamplesfile <- list.files(pattern='^_jointmcsamples.*\\.rds')
if(length(mcsamplesfile) !=1){stop('Problems with mcsamples file')}
mcsamples <- readRDS(paste0(mcsamplesfile))

data0 <- fread(paste0(origdir,datafile), sep=',')
if(!all(unlist(variate) %in% colnames(data0))){cat('\nERROR: variates missing from datafile')}
data0 <- data0[, unlist(variate), with=F]
## shuffle data
if(exists('shuffledata') && shuffledata){data0 <- data0[sample(1:nrow(data0))]}
if(!exists('ndata') || is.null(ndata) || is.na(ndata)){ndata <- nrow(data0)}
data0 <- data0[1:ndata]

cogvars <- c('ANARTERR_neuro', 'AVDEL30MIN_neuro', 'AVDELTOT_neuro', 'CATANIMSC_neuro', 'GDTOTAL_gds', 'RAVLT_immediate', 'TRAASCOR_neuro', 'TRABSCOR_neuro')
genvars <- c('Apoe4_', 'LRHHC_n_long')
demovars <- c('Gender_num_', 'AGE')

predictandvalues <- cbind(seq(varinfo[['min']][predictands],varinfo[['max']][predictands],length.out=varinfo[['n']][predictands]))
colnames(predictandvalues) <- predictands
rownames(predictandvalues) <- c('sMCI', 'cAD')

um <- matrix(c(1,0.9,0.8,0,
               0,0.3,0.5,1), 4,2)
rownames(um) <- c(1:4)

npsamples <- 1024
nlrsamples <- round(seq(1, nrow(mcsamples), length.out=1024))
meanrm <- function(x){mean(x, na.rm=T)}
##system.time(
longrunprobsamples <- vapply(nlrsamples, function(asample){
    print(asample)
        points <- t(generateVariates(Ynames=c(predictors,predictands),
                                     X=NULL,
                                     mcsamples=mcsamples[asample,,drop=F],
                                     varinfo=varinfo,
                                     n=npsamples)[,,1])
        ## test2 <- samplesFDistribution(Y=points[,predictands,drop=F],
        ##                         X=points[,predictors,drop=F],
        ##                         mcsamples=mcsamples[1:9,], varinfo=varinfo,
        ##                         jacobian=TRUE, fn=identity)
        ##
        samplesFDistribution(Y=points[,predictands,drop=F],
                                X=points[,predictors,drop=F],
                                mcsamples=mcsamples, varinfo=varinfo,
                             jacobian=TRUE, fn=meanrm)
    }, numeric(npsamples))
    ## rows: within longrunsample, cols: across longrunsamples
##)
###saveRDS(longrunprobsamples, paste0('longrunprobsamples',npsamples,'.rds'))


npsamples <- 1024
longrunprobsamples <- readRDS(paste0('_longrunprobsamples',npsamples,'.rds'))



print('done')

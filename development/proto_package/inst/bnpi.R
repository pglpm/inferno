#### This script call all libraries and functions needed by the protopackage
## Needed libraries
library('data.table')
library('foreach')
library('doParallel')
library('doRNG')
# library('khroma')
loadNamespace('LaplacesDemon')
loadNamespace('nimble')

## Main functions
source('buildmetadata.R')
source('inferpopulation.R')
## source('__orig_inferpopulation.R') # original version, for debugging
source('samplesFDistribution.R')
source('plotFsamples.R')
source('mutualinfo.R')

## Utility functions (should be invisible to the user in the package)
source('tplotfunctions.R')
source('util_buildauxmetadata.R')
## source('__orig_buildauxmetadata.R')
source('util_vtransform.R')
source('util_mcsubset.R')
source('util_proposethinning.R')
source('util_proposeburnin.R')
source('util_mcmclength.R')

## These are not needed in normal use,
## but may be needed for further development
if(FALSE){
  source('util_createQfunction.R')
}


## ## To join MC samples in list form
## joinmc <- function(mc1, mc2){
##     if(is.null(mc1)){
##         mc2
##     }else{
##         mapply(function(xx,yy){
##             temp <- c(xx,yy)
##             dx <- dim(xx)[-length(dim(xx))]
##             dim(temp) <- c(dx, length(temp)/prod(dx))
##             temp
##         },
##         mc1, mc2)
##     }
## }
## ## To remove iterations with non-finite values
## toremove <- sort(unique(unlist(lapply(mcsamplesl,function(xx){temp <- which(is.na(xx),arr.ind=T);temp[,ncol(temp)]}))))

## cleanmc <- function(mcx, tokeep){
##     lapply(mcx,function(xx){
##         do.call('[',c(list(xx),rep(TRUE,length(dim(xx))-1), list(toremove), list(drop=FALSE)) )
##     })
## }

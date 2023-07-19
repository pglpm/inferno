library('data.table')
library('foreach')
library('doParallel')
library('doRNG')
loadNamespace('LaplacesDemon')
loadNamespace('nimble')
##
source('buildmetadata.R')
source('buildauxmetadata.R')
source('samplesFDistribution.R')
source('inferpopulation.R')
source('plotFsamples.R')

## ## To join MC samples in list form
## joinmc <- function(mc1, mc2){
##     if(is.null(mc1)){
##         mc2
##     }else{
##         mapply(function(xx,yy){
##             temp <- c(xx,yy)
##             dx <- dim(xx)[-length(dim(xx))]
##             dim(temp) <- c(dx, length(temp)/dx)
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

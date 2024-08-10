#############################################################
## Run through elements
#############################################################

#### custom
dataset <- read.csv('/home/pglpm/repositories/bayes_nonparametric_inference/development/tests/custom/dataset_custom30.csv', na.strings='')
mcoutput <- readRDS('/home/pglpm/repositories/bayes_nonparametric_inference/development/tests/custom/_newdeletepackagetest-vrt10_dat30_smp120/Fdistribution.rds')
##
source('util_vtransform.R', local=T)
source('samplesFDistribution.R', local=T)
source('../development/functions_tests/testSFD.R')
##
out <- sapply(1:ncol(dataset), function(j){ sapply(1:nrow(dataset), function(i){
    v0 <- samplesFDistribution(
        Y=dataset[i,j,drop=F],
        X=NULL,
        mcoutput = mcoutput, jacobian = T, silent=T, parallel=4)
    ##
    v1 <- testSFD(
        Y=dataset[i,j,drop=F],
        X=NULL,
        mcoutput = mcoutput, jacobian = T, reduceprobx = T)
    ##
    max(abs((v1-v0)/(v1+v0)*2))
    })
})
colnames(out) <- mcoutput$auxmetadata$mcmctype
apply(abs(out),2,max)
max(abs(out))


#### nadpark
source('util_vtransform.R', local=T)
source('samplesFDistribution.R', local=T)
dataset <- read.csv('/home/pglpm/repositories/parkinsonbayes/data/nadpark/toydata.csv', na.strings='')[,-1]
mcoutput <- readRDS('/home/pglpm/repositories/parkinsonbayes/data/nadpark/_testtoy_analysis_2_3A4-1-vrt7_dat30_smp1024/Fdistribution.rds')
out <- t(sapply(1:5, function(ipoint){
    sapply(1:ncol(dataset), function(i){
        v0 <- samplesFDistribution(Y=dataset[ipoint,i,drop=F], X=NULL, mcoutput = mcoutput, jacobian = T, silent=T)
        v1 <- testSFD(Y=dataset[ipoint,i,drop=F], X=NULL, mcoutput = mcoutput)
        max(abs((v1-v0)/(v1+v0)*2))
    })
}))
colnames(out) <- mcoutput$auxmetadata$name
out
max(abs(out))




#### BN_functdependence
## B1='yes' is correlated with N1=c('A','B')
## B1='no' is correlated with N1=c('C')
dataset <- data.frame(
    B1 = c('yes','yes','no','no','yes'),
    N1 = c('A','C','A','C','B')
)
mcoutput <- readRDS('/home/pglpm/repositories/bayes_nonparametric_inference/development/tests/BN_functdependence/Fdistribution.rds')
##
source('util_vtransform.R', local=T)
source('samplesFDistribution.R', local=T)
source('../development/functions_tests/testSFD.R')
##
samplesFDistribution(
        Y=dataset[,'B1',drop=F],
        X=NULL,
        mcoutput = mcoutput, jacobian = T, silent=T, parallel=4)


samplesFDistribution(
        Y=dataset[,'N1',drop=F],
        X=NULL,
        mcoutput = mcoutput, jacobian = T, silent=T, parallel=4)

samplesFDistribution(
        Y=dataset[,'N1',drop=F],
        X=dataset[,'B1',drop=F],
        mcoutput = mcoutput, jacobian = T, silent=T, parallel=4)

samplesFDistribution(
        Y=dataset[,'B1',drop=F],
        X=dataset[,'N1',drop=F],
        mcoutput = mcoutput, jacobian = T, silent=T, parallel=4)



out <- sapply(1:ncol(dataset), function(j){ sapply(1:nrow(dataset), function(i){
    v0 <- samplesFDistribution(
        Y=dataset[i,j,drop=F],
        X=NULL,
        mcoutput = mcoutput, jacobian = T, silent=T, parallel=4)
    ##
    v1 <- testSFD(
        Y=dataset[i,j,drop=F],
        X=NULL,
        mcoutput = mcoutput, jacobian = T, reduceprobx = T)
    ##
    max(abs((v1-v0)/(v1+v0)*2))
    })
})
colnames(out) <- mcoutput$auxmetadata$mcmctype
apply(abs(out),2,max)
max(abs(out))



#############################################################
## Combinations YX
#############################################################

#### custom
dataset <- read.csv('/home/pglpm/repositories/bayes_nonparametric_inference/development/test_custom/datanew_test_custom_30.csv', na.strings='')
mcoutput <- readRDS('/home/pglpm/repositories/bayes_nonparametric_inference/development/test_custom/_newdeletepackagetest-vrt10_dat30_smp120/Fdistribution.rds')
jac <- T
##
nY <- sample(1:(ncol(dataset)-1), 1)
Ys <- sample(colnames(dataset), nY)
Xs <- sample(setdiff(colnames(dataset), Ys), ncol(dataset)-nY)
npoints <- nrow(dataset)
iY <- sample(1:nrow(dataset),npoints)
iX <- sample(1:nrow(dataset),npoints)
##
source('util_vtransform.R', local=T)
source('samplesFDistribution.R', local=T)
source('/home/pglpm/repositories/bayes_nonparametric_inference/development/scripts_tests/testSFD.R', local=T)
out <- sapply(1:length(iY), function(i){
    v0 <- samplesFDistribution(
        Y=dataset[iY[i],Ys,drop=F],
        X=dataset[iX[i],Xs,drop=F],
        mcoutput = mcoutput, jacobian = T, silent=T, parallel=1)
    ##
    v1 <- testSFD(
        Y=dataset[iY[i],Ys,drop=F],
        X=dataset[iX[i],Xs,drop=F],
        mcoutput = mcoutput, jacobian = T)
    ##
    max(abs((v1-v0)/(v1+v0)*2))
})
mcoutput$auxmetadata$mcmctype[mcoutput$auxmetadata$name %in% Ys]
mcoutput$auxmetadata$mcmctype[mcoutput$auxmetadata$name %in% Xs]
max(abs(out))


#### mtcars
data(mtcars)
dataset <- mtcars
mcoutput <- readRDS('/home/pglpm/repositories/bayes_nonparametric_inference/development/tests/mtcars/_newD_test_mtcars-3-vrt11_dat32_smp512/Fdistribution.rds')
jac <- T
##
nY <- sample(1:(ncol(dataset)-1), 1)
Ys <- sample(colnames(dataset), nY)
Xs <- sample(setdiff(colnames(dataset), Ys), ncol(dataset)-nY)
npoints <- nrow(dataset)
iY <- sample(1:nrow(dataset),npoints)
iX <- sample(1:nrow(dataset),npoints)
##
source('util_vtransform.R', local=T)
source('samplesFDistribution.R', local=T)
source('/home/pglpm/repositories/bayes_nonparametric_inference/development/scripts_tests/testSFD.R', local=T)
out <- sapply(1:length(iY), function(i){
    v0 <- samplesFDistribution(
        Y=dataset[iY[i],Ys,drop=F],
        X=dataset[iX[i],Xs,drop=F],
        mcoutput = mcoutput, jacobian = jac, silent=T, parallel=4)
    ##
    v1 <- testSFD(
        Y=dataset[iY[i],Ys,drop=F],
        X=dataset[iX[i],Xs,drop=F],
        mcoutput = mcoutput, jacobian = jac)
    ##
    max(abs((v1-v0)/(v1+v0)*2))
})
mcoutput$auxmetadata$mcmctype[mcoutput$auxmetadata$name %in% Ys]
mcoutput$auxmetadata$mcmctype[mcoutput$auxmetadata$name %in% Xs]
max(abs(out))

#### nadpark
dataset <- read.csv('/home/pglpm/repositories/parkinsonbayes/data/nadpark/toydata.csv', na.strings='')[,-1]
mcoutput <- readRDS('/home/pglpm/repositories/parkinsonbayes/data/nadpark/_testtoy_analysis_2_3A4-1-vrt7_dat30_smp1024/Fdistribution.rds')
jac <- T
##
nY <- sample(1:(ncol(dataset)-1), 1)
Ys <- sample(colnames(dataset), nY)
Xs <- sample(setdiff(colnames(dataset), Ys), ncol(dataset)-nY)
npoints <- nrow(dataset)
iY <- sample(1:nrow(dataset),npoints)
iX <- sample(1:nrow(dataset),npoints)
##
source('util_vtransform.R', local=T)
source('samplesFDistribution.R', local=T)
source('/home/pglpm/repositories/bayes_nonparametric_inference/development/scripts_tests/testSFD.R', local=T)
out <- sapply(1:length(iY), function(i){
    v0 <- samplesFDistribution(
        Y=dataset[iY[i],Ys,drop=F],
        X=dataset[iX[i],Xs,drop=F],
        mcoutput = mcoutput, jacobian = jac, silent=T, parallel=4)
    ##
    v1 <- testSFD(
        Y=dataset[iY[i],Ys,drop=F],
        X=dataset[iX[i],Xs,drop=F],
        mcoutput = mcoutput, jacobian = jac)
    ##
    max(abs((v1-v0)/(v1+v0)*2))
})
mcoutput$auxmetadata$mcmctype[mcoutput$auxmetadata$name %in% Ys]
mcoutput$auxmetadata$mcmctype[mcoutput$auxmetadata$name %in% Xs]
max(abs(out))



#############################################################
## Run through elements
#############################################################

#### custom – test time use of util_lprob.R
dataset <- read.csv('~/repos/bayes_nonparametric_inference/development/tests/custom/dataset_custom30.csv', na.strings='')
learnt <- readRDS('~/repos/bayes_nonparametric_inference/development/tests/custom/test_custom-240816T094427-vrt10_dat30_smp120/learnt.rds')
##
source('util_vtransform.R', local=T)
source('samplesFDistribution.R', local=T)
source('util_lprob.R', local=T)
source('samplesFDistribution2.R', local=T)

system.time(
    out0 <- sapply(1:ncol(dataset), function(j){ sapply(1:nrow(dataset), function(i){
        samplesFDistribution(
            Y=dataset[i,j,drop=F],
            X=NULL,
            learnt = learnt, jacobian = T, silent=T, parallel=1)
    })})
)

system.time(
    out2 <- sapply(1:ncol(dataset), function(j){ sapply(1:nrow(dataset), function(i){
        samplesFDistribution2(
            Y=dataset[i,j,drop=F],
            X=NULL,
            learnt = learnt, jacobian = T, silent=T, parallel=1)
    })})
)

max(abs(2*(out0-out2)/(out0+out2)))






#### custom
dataset <- read.csv('~/repositories/bayes_nonparametric_inference/development/tests/custom/dataset_custom30.csv', na.strings='')
learnt <- readRDS('~/repositories/bayes_nonparametric_inference/development/tests/custom/_newdeletepackagetest-vrt10_dat30_smp120/learnt.rds')
##
source('util_vtransform.R', local=T)
source('samplesFDistribution.R', local=T)
source('../development/functions_tests/testSFD.R')
##
out <- sapply(1:ncol(dataset), function(j){ sapply(1:nrow(dataset), function(i){
    v0 <- samplesFDistribution(
        Y=dataset[i,j,drop=F],
        X=NULL,
        learnt = learnt, jacobian = T, silent=T, parallel=4)
    ##
    v1 <- testSFD(
        Y=dataset[i,j,drop=F],
        X=NULL,
        learnt = learnt, jacobian = T, reduceprobx = T)
    ##
    max(abs((v1-v0)/(v1+v0)*2))
    })
})
colnames(out) <- learnt$auxmetadata$mcmctype
apply(abs(out),2,max)
max(abs(out))


#### nadpark
source('util_vtransform.R', local=T)
source('samplesFDistribution.R', local=T)
dataset <- read.csv('~/repositories/parkinsonbayes/data/nadpark/toydata.csv', na.strings='')[,-1]
learnt <- readRDS('~/repositories/parkinsonbayes/data/nadpark/_testtoy_analysis_2_3A4-1-vrt7_dat30_smp1024/learnt.rds')
out <- t(sapply(1:5, function(ipoint){
    sapply(1:ncol(dataset), function(i){
        v0 <- samplesFDistribution(Y=dataset[ipoint,i,drop=F], X=NULL, learnt = learnt, jacobian = T, silent=T)
        v1 <- testSFD(Y=dataset[ipoint,i,drop=F], X=NULL, learnt = learnt)
        max(abs((v1-v0)/(v1+v0)*2))
    })
}))
colnames(out) <- learnt$auxmetadata$name
out
max(abs(out))




#### BN_functdependence
## B1='yes' is correlated with N1=c('A','B')
## B1='no' is correlated with N1=c('C')
dataset <- data.frame(
    B1 = c('yes','yes','no','no','yes'),
    N1 = c('A','C','A','C','B')
)
learnt <- readRDS('~/repositories/bayes_nonparametric_inference/development/tests/BN_functdependence/learnt.rds')
##
source('util_vtransform.R', local=T)
source('samplesFDistribution.R', local=T)
source('../development/functions_tests/testSFD.R')
##
out <- samplesFDistribution(
        Y=dataset[,'B1',drop=F],
        X=NULL,
        learnt = learnt, jacobian = T, silent=T, parallel=4)
out
rowMeans(out)


out <- samplesFDistribution(
        Y=dataset[,'N1',drop=F],
        X=NULL,
        learnt = learnt, jacobian = T, silent=T, parallel=4)
out
rowMeans(out)

samplesFDistribution(
        Y=dataset[,'N1',drop=F],
        X=dataset[,'B1',drop=F],
        learnt = learnt, jacobian = T, silent=T, parallel=4)

samplesFDistribution(
        Y=dataset[,'B1',drop=F],
        X=dataset[,'N1',drop=F],
        learnt = learnt, jacobian = T, silent=T, parallel=4)

samplesFDistribution(
        Y=cbind(B1='yes'),
        X=cbind(N1=LETTERS[1:3]),
        learnt = learnt, jacobian = T, silent=T, parallel=4)

samplesFDistribution(
        Y=cbind(N1=LETTERS[1:3]),
        X=cbind(B1='yes'),
        learnt = learnt, jacobian = T, silent=T, parallel=4)

samplesFDistribution(
        Y=cbind(N1=LETTERS[1:3]),
        X=cbind(B1='no'),
        learnt = learnt, jacobian = T, silent=T, parallel=4)



out <- sapply(1:ncol(dataset), function(j){ sapply(1:nrow(dataset), function(i){
    v0 <- samplesFDistribution(
        Y=dataset[i,j,drop=F],
        X=NULL,
        learnt = learnt, jacobian = T, silent=T, parallel=4)
    ##
    v1 <- testSFD(
        Y=dataset[i,j,drop=F],
        X=NULL,
        learnt = learnt, jacobian = T, reduceprobx = T)
    ##
    max(abs((v1-v0)/(v1+v0)*2))
    })
})
colnames(out) <- learnt$auxmetadata$mcmctype
apply(abs(out),2,max)
max(abs(out))



#############################################################
## Combinations YX
#############################################################

#### custom - test new samplesFDistribution
dataset <- read.csv('~/repositories/bayes_nonparametric_inference/development/tests/custom/dataset_custom150.csv', na.strings='')
learnt <- readRDS('~/repositories/bayes_nonparametric_inference/development/tests/custom/_newdeletepackagetest-vrt10_dat30_smp120/learnt.rds')
##

library('inferno')

nY <- sample(1:(ncol(dataset)-1), 1)
Ys <- sample(colnames(dataset), nY)
Xs <- sample(setdiff(colnames(dataset), Ys), ncol(dataset)-nY)
npoints <- nrow(dataset)
iY <- sample(1:nrow(dataset),npoints)
iX <- sample(1:nrow(dataset),npoints)
##
Ys
Xs

system.time(
out <- sapply(1:length(iY), function(i){
    samplesFDistribution(
        Y=dataset[iY[i],Ys,drop=F],
        X=dataset[iX[i],Xs,drop=F],
        learnt = learnt, jacobian = T, silent=T, parallel=1)
})
)

system.time(
out2 <- sapply(1:length(iY), function(i){
    samplesFDistribution2(
        Y=dataset[iY[i],Ys,drop=F],
        X=dataset[iX[i],Xs,drop=F],
        learnt = learnt, jacobian = T, silent=T, parallel=1)
})
)

max(abs(2*(out-out2)/(out+out2)))





#### custom
dataset <- read.csv('~/repositories/bayes_nonparametric_inference/development/test_custom/datanew_test_custom_30.csv', na.strings='')
learnt <- readRDS('~/repositories/bayes_nonparametric_inference/development/test_custom/_newdeletepackagetest-vrt10_dat30_smp120/learnt.rds')
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
source('~/repositories/bayes_nonparametric_inference/development/scripts_tests/testSFD.R', local=T)
out <- sapply(1:length(iY), function(i){
    v0 <- samplesFDistribution(
        Y=dataset[iY[i],Ys,drop=F],
        X=dataset[iX[i],Xs,drop=F],
        learnt = learnt, jacobian = T, silent=T, parallel=1)
    ##
    v1 <- testSFD(
        Y=dataset[iY[i],Ys,drop=F],
        X=dataset[iX[i],Xs,drop=F],
        learnt = learnt, jacobian = T)
    ##
    max(abs((v1-v0)/(v1+v0)*2))
})
learnt$auxmetadata$mcmctype[learnt$auxmetadata$name %in% Ys]
learnt$auxmetadata$mcmctype[learnt$auxmetadata$name %in% Xs]
max(abs(out))


#### mtcars
data(mtcars)
dataset <- mtcars
learnt <- readRDS('~/repositories/bayes_nonparametric_inference/development/tests/mtcars/_newD_test_mtcars-3-vrt11_dat32_smp512/learnt.rds')
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
source('~/repositories/bayes_nonparametric_inference/development/scripts_tests/testSFD.R', local=T)
out <- sapply(1:length(iY), function(i){
    v0 <- samplesFDistribution(
        Y=dataset[iY[i],Ys,drop=F],
        X=dataset[iX[i],Xs,drop=F],
        learnt = learnt, jacobian = jac, silent=T, parallel=4)
    ##
    v1 <- testSFD(
        Y=dataset[iY[i],Ys,drop=F],
        X=dataset[iX[i],Xs,drop=F],
        learnt = learnt, jacobian = jac)
    ##
    max(abs((v1-v0)/(v1+v0)*2))
})
learnt$auxmetadata$mcmctype[learnt$auxmetadata$name %in% Ys]
learnt$auxmetadata$mcmctype[learnt$auxmetadata$name %in% Xs]
max(abs(out))


dataset <- mtcars
learnt <- readRDS('~/repositories/bayes_nonparametric_inference/development/tests/mtcars/_newD_test_mtcars-3-vrt11_dat32_smp512/learnt.rds')


library('inferno')

nY <- sample(1:(ncol(dataset)-1), 1)
Ys <- sample(colnames(dataset), nY)
Xs <- sample(setdiff(colnames(dataset), Ys), ncol(dataset)-nY)
npoints <- nrow(dataset)
iY <- sample(1:nrow(dataset),npoints)
iX <- sample(1:nrow(dataset),npoints)
##
Ys
Xs

system.time(
out <- sapply(1:length(iY), function(i){
    samplesFDistribution(
        Y=dataset[iY[i],Ys,drop=F],
        X=dataset[iX[i],Xs,drop=F],
        learnt = learnt, jacobian = T, silent=T, parallel=1)
})
)

system.time(
out2 <- sapply(1:length(iY), function(i){
    samplesFDistribution2(
        Y=dataset[iY[i],Ys,drop=F],
        X=dataset[iX[i],Xs,drop=F],
        learnt = learnt, jacobian = T, silent=T, parallel=1)
})
)

max(abs(2*(out-out2)/(out+out2)))




#### nadpark
dataset <- read.csv('~/repositories/parkinsonbayes/data/nadpark/toydata.csv', na.strings='')[,-1]
learnt <- readRDS('~/repositories/parkinsonbayes/data/nadpark/_testtoy_analysis_2_3A4-1-vrt7_dat30_smp1024/learnt.rds')
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
source('~/repositories/bayes_nonparametric_inference/development/scripts_tests/testSFD.R', local=T)
out <- sapply(1:length(iY), function(i){
    v0 <- samplesFDistribution(
        Y=dataset[iY[i],Ys,drop=F],
        X=dataset[iX[i],Xs,drop=F],
        learnt = learnt, jacobian = jac, silent=T, parallel=4)
    ##
    v1 <- testSFD(
        Y=dataset[iY[i],Ys,drop=F],
        X=dataset[iX[i],Xs,drop=F],
        learnt = learnt, jacobian = jac)
    ##
    max(abs((v1-v0)/(v1+v0)*2))
})
learnt$auxmetadata$mcmctype[learnt$auxmetadata$name %in% Ys]
learnt$auxmetadata$mcmctype[learnt$auxmetadata$name %in% Xs]
max(abs(out))


#####################################################
#### test Pr()
#####################################################


testrd <- function(x,y){2*max(abs((x-y)/(x+y)), na.rm=T)}
#### custom – test time use of util_lprob.R
dataset <- read.csv('~/repos/bayes_nonparametric_inference/development/tests/custom/dataset_custom30.csv', na.strings='')
learnt <- readRDS('~/repos/bayes_nonparametric_inference/development/tests/custom/test_custom-240816T094427-vrt10_dat30_smp120/learnt.rds')

nY <- sample(1:(ncol(dataset)-1), 1)
Ys <- sample(colnames(dataset), nY)
Xs <- sample(setdiff(colnames(dataset), Ys), ncol(dataset)-nY)
##
YY <- dataset[sample(1:nrow(dataset),2),Ys, drop=F]
XX <- NULL#dataset[sample(1:nrow(dataset),20),Xs, drop=F]
##
## load('~/repos/bayes_nonparametric_inference/R/sysdata.rda')
## source('~/repos/bayes_nonparametric_inference/R/util_lprob.R', local=T)
## source('~/repos/bayes_nonparametric_inference/R/util_vtransform.R', local=T)
## source('~/repos/bayes_nonparametric_inference/R/samplesFDistribution.R', local=T)

test0 <- rbind(samplesFDistribution(
    Y = YY,
    X = XX,
    learnt = learnt, parallel=2))

##
## rowMeans(test0, na.rm = T)
## apply(test0, 1, quantile, prob=c(5,95)/100, type = 6, na.rm=T)
##
## source('~/repos/bayes_nonparametric_inference/R/util_vtransform.R', local=T)
## source('~/repos/bayes_nonparametric_inference/R/util_lprob.R', local=T)
## load('~/repos/bayes_nonparametric_inference/R/sysdata.rda')
## source('~/repos/bayes_nonparametric_inference/R/Pr.R', local=T)
test1 <- Pr(
    Y = YY,
    X = XX,
    learnt = learnt, parallel=2)

##
testrd(diag(test1$values) , rowMeans(test0, na.rm=T))
## test1$quantiles
testrd(t(apply(test1$quantiles,3,diag)),
    apply(test0, 1, quantile, prob=c(5,95)/100, type = 6, na.rm=T))

test1$values
rowMeans(test0, na.rm = T)

## > Ys
## [1] "Nvrt"  "N2vrt" "Rvrt"  "Cvrt"  "Ruvrt"
## > Xs
## [1] "Bvrt"  "Duvrt" "Ovrt"  "Dvrt"  "Rpvrt"


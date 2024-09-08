source('util_vtransform.R', local=T)

######################################################################
#### mtcars
data(mtcars)
dataset <- mtcars
learnt <- readRDS('~/repos/inferno/development/tests/mtcars/_newD_test_mtcars-2-vrt11_dat32_smp512/learnt.rds')
jac <- T

library('foreach')
source('util_tplotfunctions.R', local=T)
source('util_vtransform.R', local=T)
source('samplesFDistribution.R', local=T)
source('plotFsamples.R', local=T)
plotFsamples(
    file='_justatest',
    learnt = learnt,
    data = dataset,
    plotprobability = TRUE,
    plotvariability = 'samples',
    nFsamples = NULL,
    #datahistogram = !(missing(data) || is.null(data)),
    #datascatter = !(missing(data) || is.null(data)),
    parallel = 4,
    silent = FALSE)


######################################################################
#### nadpark
dataset <- read.csv('~/repositories/parkinsonbayes/data/nadpark/toydata.csv', na.strings='')
learnt <- readRDS('~/repositories/parkinsonbayes/data/nadpark/_testtoy_analysis_2_3A4-1-vrt7_dat30_smp1024/learnt.rds')
jac <- T

library('foreach')
source('util_tplotfunctions.R', local=T)
source('util_vtransform.R', local=T)
source('samplesFDistribution.R', local=T)
source('plotFsamples.R', local=T)
for(pv in c('samples', 'quantiles')) {
plotFsamples(
    file= paste0('_justatest_', pv),
    learnt = learnt,
    data = dataset,
    plotprobability = TRUE,
    plotvariability = pv,
    nFsamples = NULL,
    #datahistogram = !(missing(data) || is.null(data)),
    #datascatter = !(missing(data) || is.null(data)),
    parallel = 4,
    silent = FALSE)
}


######################################################################
#### custom
dataset <- read.csv('~/repos/inferno/development/tests/custom/dataset_custom30.csv', na.strings='')
learnt <- readRDS('~/repos/inferno/development/tests/custom/test_custom-vrt10_dat30_smp120/learnt.rds')

unloadNamespace('inferno')
library('inferno')
##
for(pv in c('samples', 'quantiles')) {
plotFsamples(
    file= paste0('_justatest_', pv),
    learnt = learnt,
    data = dataset,
    plotprobability = TRUE,
    plotvariability = pv,
    nFsamples = NULL,
    #datahistogram = !(missing(data) || is.null(data)),
    #datascatter = !(missing(data) || is.null(data)),
    parallel = 1,
    silent = FALSE)
}




dataset <- read.csv('~/repositories/parkinsonbayes/data/nadpark/toydata.csv', na.strings='')[,-1]
learnt <- readRDS('~/repos/inferno/development/tests/custom/test_custom-240901T130028-vrt10_dat30_smp120/learnt.rds')
out <- t(sapply(1:5, function(ipoint){
    sapply(1:ncol(dataset), function(i){
        v0 <- samplesFDistribution(Y=dataset[ipoint,i,drop=F], X=NULL, learnt = learnt, jacobian = T, silent=T)
        v1 <- testSFD(Y=dataset[ipoint,i,drop=F], X=NULL, learnt = learnt)
        max(abs((v1-v0)/(v1+v0)*2))
    })
}))
colnames(out) <- learnt$auxmetadata$name
out



source('util_vtransform.R', local=T)
source('samplesFDistribution.R', local=T)
data(mtcars)
dataset <- mtcars
learnt <- readRDS('~/repos/inferno/development/tests/mtcars/test_mtcars_learnt.rds')
jac <- T
##
out <- t(sapply(1:32, function(ipoint){
    sapply(1:ncol(dataset), function(i){
        v0 <- samplesFDistribution(Y=dataset[ipoint,i,drop=F], X=NULL, learnt = learnt, jacobian = jac, silent=T)
        v1 <- testSFD(Y=dataset[ipoint,i,drop=F], X=NULL, learnt = learnt, jacobian = jac)
        max(abs((v1-v0)/(v1+v0)*2))
    })
}))
colnames(out) <- learnt$auxmetadata$name
out
max(out)





data <- read.csv('~/repos/inferno/development/tests/custom/dataset_custom30.csv', na.strings='')
learnt <- readRDS('~/repos/inferno/development/tests/custom/test_custom-240822T190544-vrt10_dat30_smp120/learnt.rds')

unloadNamespace('inferno')
library('inferno')
##
plotFsamples(file='justatest',learnt=learnt,data=data,plotvariability = 'quantiles',nFsamples = c(5,95)/100,plotprobability = T,datahistogram = T,datascatter = T,parallel=4)


plotFsamples(file='justatest',learnt=learnt,data=data,plotvariability = 'samples',nFsamples = 100,plotprobability = T,datahistogram = T,datascatter = T,parallel=4)

Y <- data.frame(Rvrt=seq(-1,1,length.out=10))
X <- NULL
Pr(Y=Y, X=X, learnt=learnt,nsamples=3)




Y <- cbind(Age=65)

X <- NULL

learnt <- readRDS('~/repos/parkinsonbayes/data/nadpark/_output_learn_toydata-1-240810T161910-vrt7_dat30_smp3600/learnt.rds')


prediction <- Pr(Y=Y, X=X, learnt=learnt)

prediction$value

Y <- cbind(diff.MDS.UPRS.III = (-10):10)
X <- cbind(Sex='Male', TreatmentGroup='NR')
predictionM <- Pr(Y=Y, X=X, learnt=learnt)
predictionM$value

tplot(x=(-10):10, y=predictionM$values)


Y <- cbind(diff.MDS.UPRS.III = (-10):10)
X <- cbind(Sex='Female', TreatmentGroup='NR')
predictionF <- Pr(Y=Y, X=X, learnt='~/repos/parkinsonbayes/data/nadpark/_output_learn_toydata-1-240810T161910-vrt7_dat30_smp3600')
predictionF$value

plot(x=(-10):10, y=predictionF$value, type='l',col='blue')

plotquantiles(x=(-10):10, y=predictionM$quantiles[,1,],col=1)
tplot(x=(-10):10, y=predictionM$value, add=T,col=1)
plotquantiles(x=(-10):10, y=predictionF$quantiles[,1,],add=T,col=2)
tplot(x=(-10):10, y=predictionF$value, add=T,col=2)


Y <- cbind(diff.MDS.UPRS.III = (-132):132)
X <- cbind(Sex='Male', TreatmentGroup='NR')
predictionM <- Pr(Y=Y, X=X, learnt=learnt)
predictionM$value

Y <- cbind(diff.MDS.UPRS.III = (-132):132)
X <- cbind(Sex='Female', TreatmentGroup='NR')
predictionF <- Pr(Y=Y, X=X, learnt=learnt)
predictionF$value


sum(predictionM$values * ((-132):132))
sum(predictionF$values * ((-132):132))

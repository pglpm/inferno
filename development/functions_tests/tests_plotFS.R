source('util_vtransform.R', local=T)

######################################################################
#### mtcars
data(mtcars)
dataset <- mtcars
learnt <- readRDS('~/repositories/bayes_nonparametric_inference/development/tests/mtcars/_newD_test_mtcars-2-vrt11_dat32_smp512/learnt.rds')
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
dataset <- read.csv('~/repositories/bayes_nonparametric_inference/development/test_custom/datanew_test_custom_30.csv', na.strings='')
learnt <- readRDS('~/repositories/bayes_nonparametric_inference/development/test_custom/_newdeletepackagetest-vrt10_dat30_smp120/learnt.rds')
jac <- T
##
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




dataset <- read.csv('~/repositories/parkinsonbayes/data/nadpark/toydata.csv', na.strings='')[,-1]
learnt <- readRDS('~/repositories/parkinsonbayes/data/nadpark/_testtoy_analysis_3A4-1-vrt7_dat30_smp1024/learnt.rds')
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
learnt <- readRDS('~/repositories/bayes_nonparametric_inference/development/tests/mtcars/test_mtcars_learnt.rds')
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

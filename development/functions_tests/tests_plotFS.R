source('util_vtransform.R', local=T)

######################################################################
#### mtcars
data(mtcars)
dataset <- mtcars
agent <- readRDS('/home/pglpm/repositories/bayes_nonparametric_inference/development/tests/mtcars/_newD_test_mtcars-2-vrt11_dat32_smp512/agent.rds')
jac <- T

library('foreach')
source('util_tplotfunctions.R', local=T)
source('util_vtransform.R', local=T)
source('samplesFDistribution.R', local=T)
source('plotFsamples.R', local=T)
plotFsamples(
    file='_justatest',
    agent = agent,
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
dataset <- read.csv('/home/pglpm/repositories/parkinsonbayes/data/nadpark/toydata.csv', na.strings='')
agent <- readRDS('/home/pglpm/repositories/parkinsonbayes/data/nadpark/_testtoy_analysis_2_3A4-1-vrt7_dat30_smp1024/agent.rds')
jac <- T

library('foreach')
source('util_tplotfunctions.R', local=T)
source('util_vtransform.R', local=T)
source('samplesFDistribution.R', local=T)
source('plotFsamples.R', local=T)
for(pv in c('samples', 'quantiles')) {
plotFsamples(
    file= paste0('_justatest_', pv),
    agent = agent,
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
dataset <- read.csv('/home/pglpm/repositories/bayes_nonparametric_inference/development/test_custom/datanew_test_custom_30.csv', na.strings='')
agent <- readRDS('/home/pglpm/repositories/bayes_nonparametric_inference/development/test_custom/_newdeletepackagetest-vrt10_dat30_smp120/agent.rds')
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
    agent = agent,
    data = dataset,
    plotprobability = TRUE,
    plotvariability = pv,
    nFsamples = NULL,
    #datahistogram = !(missing(data) || is.null(data)),
    #datascatter = !(missing(data) || is.null(data)),
    parallel = 4,
    silent = FALSE)
}




dataset <- read.csv('/home/pglpm/repositories/parkinsonbayes/data/nadpark/toydata.csv', na.strings='')[,-1]
agent <- readRDS('/home/pglpm/repositories/parkinsonbayes/data/nadpark/_testtoy_analysis_3A4-1-vrt7_dat30_smp1024/agent.rds')
out <- t(sapply(1:5, function(ipoint){
    sapply(1:ncol(dataset), function(i){
        v0 <- samplesFDistribution(Y=dataset[ipoint,i,drop=F], X=NULL, agent = agent, jacobian = T, silent=T)
        v1 <- testSFD(Y=dataset[ipoint,i,drop=F], X=NULL, agent = agent)
        max(abs((v1-v0)/(v1+v0)*2))
    })
}))
colnames(out) <- agent$auxmetadata$name
out



source('util_vtransform.R', local=T)
source('samplesFDistribution.R', local=T)
data(mtcars)
dataset <- mtcars
agent <- readRDS('/home/pglpm/repositories/bayes_nonparametric_inference/development/tests/mtcars/test_mtcars_agent.rds')
jac <- T
##
out <- t(sapply(1:32, function(ipoint){
    sapply(1:ncol(dataset), function(i){
        v0 <- samplesFDistribution(Y=dataset[ipoint,i,drop=F], X=NULL, agent = agent, jacobian = jac, silent=T)
        v1 <- testSFD(Y=dataset[ipoint,i,drop=F], X=NULL, agent = agent, jacobian = jac)
        max(abs((v1-v0)/(v1+v0)*2))
    })
}))
colnames(out) <- agent$auxmetadata$name
out
max(out)

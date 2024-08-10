#### custom
dataset <- read.csv('/home/pglpm/repositories/bayes_nonparametric_inference/development/tests/custom/dataset_custom30.csv', na.strings='')
mcoutput <- readRDS('/home/pglpm/repositories/bayes_nonparametric_inference/development/tests/custom/_newdeletepackagetest-vrt10_dat30_smp120/Fdistribution.rds')
##
Y1s <- sample(colnames(dataset),2)
Y2s <- sample(setdiff(colnames(dataset), Y1s),2)
Xs <- sample(setdiff(colnames(dataset), c(Y1s,Y2s)),2)
##
source('util_vtransform.R', local=T)
source('samplesFDistribution.R', local=T)
source('mutualinfo.R', local=T)
set.seed(16)
##
out <- mutualinfo(
    Y1names = Y1s,
    Y2names = Y2s,
    X = NULL,#mtcars[1, Xs, drop=F],
    mcoutput = mcoutput,
    nsamples = 3600,
    parallel = 4
)
out

out <- mutualinfo(
    Y1names = 'N2vrt',
    Y2names = c('Rvrt', 'Bvrt'),
    X = NULL,
    mcoutput = mcoutput,
    nsamples = 3600,
    parallel = 4
)
out



#### mtcars
data(mtcars)
dataset <- mtcars
mcoutput <- readRDS('/home/pglpm/repositories/bayes_nonparametric_inference/development/tests/mtcars/_newD_test_mtcars-3-vrt11_dat32_smp512/Fdistribution.rds')
##
Y1s <- sample(colnames(mtcars),2)
Y2s <- sample(setdiff(colnames(mtcars), Y1s),2)
Xs <- sample(setdiff(colnames(mtcars), c(Y1s,Y2s)),2)
##
source('util_vtransform.R', local=T)
source('samplesFDistribution.R', local=T)
source('mutualinfo.R', local=T)
## set.seed(16)
##
out <- mutualinfo(
    Y1names = Y1s,
    Y2names = Y2s,
    X = NULL,#mtcars[1, Xs, drop=F],
    mcoutput = mcoutput,
    nsamples = 3600,
    parallel = 4
)
out


#### BN_functdependence
## B1='yes' is correlated with N1=c('A','B')
## B1='no' is correlated with N1=c('C')
mcoutput <- readRDS('/home/pglpm/repositories/bayes_nonparametric_inference/development/tests/BN_functdependence/Fdistribution.rds')
##
source('util_vtransform.R', local=T)
source('samplesFDistribution.R', local=T)
source('mutualinfo.R', local=T)
set.seed(16)
out <- mutualinfo(
    Y1names = 'B1',
    Y2names = 'N1',
    X = NULL,
    mcoutput = mcoutput,
    nsamples = 3600,
    parallel = 4
)
out

source('mutualinfo.R', local=T)
out <- mutualinfo(
    Y1names = 'N1',
    Y2names = NULL,
    X = cbind(B1='no'),
    mcoutput = mcoutput,
    nsamples = 3600,
    parallel = 4
)
out

source('mutualinfo.R', local=T)
out <- mutualinfo(
    Y1names = 'B1',
    Y2names = NULL,
    X = cbind(N1='C'),
    mcoutput = mcoutput,
    nsamples = 3600,
    parallel = 4
)
out

source('mutualinfo.R', local=T)
out <- mutualinfo(
    Y1names = 'N1',
    Y2names = 'B1',
    X = NULL,
    mcoutput = mcoutput,
    nsamples = 3600,
    parallel = 4
)
out

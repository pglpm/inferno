#### mtcars
data(mtcars)
dataset <- mtcars
mcoutput <- readRDS('/home/pglpm/repositories/bayes_nonparametric_inference/development/tests/mtcars/_newD_test_mtcars-3-vrt11_dat32_smp512/Fdistribution.rds')
##
Y1s <- sample(colnames(mtcars),2)
Y2s <- sample(setdiff(colnames(mtcars), Y1s),2)
Xs <- sample(setdiff(colnames(mtcars), c(Y1s,Y2s)),2)

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


#### BN_functdependence
## B1='yes' is correlated with N1=c('A','B')
## B1='no' is correlated with N1=c('C')
mcoutput <- readRDS('/home/pglpm/repositories/bayes_nonparametric_inference/development/tests/BN_functdependence/Fdistribution.rds')
##
Y1s <- 'B1'
Y2s <- 'N1'
Xs <- NULL
##
source('util_vtransform.R', local=T)
source('samplesFDistribution.R', local=T)
source('mutualinfo.R', local=T)
set.seed(16)
##
out <- mutualinfo(
    Y1names = Y1s,
    Y2names = Y2s,
    X = Xs,
    mcoutput = mcoutput,
    nsamples = 3600,
    parallel = 4
)
out

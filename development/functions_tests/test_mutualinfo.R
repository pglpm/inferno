#### custom
dataset <- read.csv('~/repositories/bayes_nonparametric_inference/development/tests/custom/dataset_custom30.csv', na.strings='')
learned <- readRDS('~/repositories/bayes_nonparametric_inference/development/tests/custom/_newdeletepackagetest-vrt10_dat30_smp120/learned.rds')
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
    learned = learned,
    nsamples = 3600,
    parallel = 4
)
out

out <- mutualinfo(
    Y1names = 'N2vrt',
    Y2names = c('Rvrt', 'Bvrt'),
    X = NULL,
    learned = learned,
    nsamples = 3600,
    parallel = 4
)
out


library('predict')



dataset <- read.csv('~/repositories/bayes_nonparametric_inference/development/tests/custom/dataset_custom30.csv', na.strings='')
learned <- readRDS('~/repositories/bayes_nonparametric_inference/development/tests/custom/_newdeletepackagetest-vrt10_dat30_smp120/learned.rds')
##
Y1s <- sample(colnames(dataset),3)
Y2s <- sample(setdiff(colnames(dataset), Y1s),3)
Xs <- sample(setdiff(colnames(dataset), c(Y1s,Y2s)),3)

set.seed(16)
##
system.time(
    out <- mutualinfo(
    Y1names = Y1s,
    Y2names = Y2s,
    X = dataset[1,Xs],
    learned = learned,
    nsamples = 3600,
    parallel = 4
    )
)


set.seed(16)
##
system.time(
out2 <- mutualinfo2(
    Y1names = Y1s,
    Y2names = Y2s,
    X = dataset[1,Xs],
    learned = learned,
    nsamples = 3600,
    parallel = 1
)
)


identical(out,out2)







#### mtcars
data(mtcars)
dataset <- mtcars
learned <- readRDS('~/repositories/bayes_nonparametric_inference/development/tests/mtcars/_newD_test_mtcars-3-vrt11_dat32_smp512/learned.rds')
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
    learned = learned,
    nsamples = 3600,
    parallel = 4
)
out


#### BN_functdependence
## B1='yes' is correlated with N1=c('A','B')
## B1='no' is correlated with N1=c('C')
learned <- readRDS('~/repositories/bayes_nonparametric_inference/development/tests/BN_functdependence/learned.rds')

source('util_vtransform.R', local=T)
source('samplesFDistribution.R', local=T)
source('mutualinfo.R', local=T)
set.seed(16)
out <- mutualinfo(
    Y1names = 'B1',
    Y2names = 'N1',
    X = NULL,
    learned = learned,
    nsamples = 3600,
    parallel = 4
)
out

source('mutualinfo.R', local=T)
out <- mutualinfo(
    Y1names = 'N1',
    Y2names = NULL,
    X = cbind(B1='no'),
    learned = learned,
    nsamples = 3600,
    parallel = 4
)
out

source('mutualinfo.R', local=T)
out <- mutualinfo(
    Y1names = 'B1',
    Y2names = NULL,
    X = cbind(N1='C'),
    learned = learned,
    nsamples = 3600,
    parallel = 4
)
out

source('mutualinfo.R', local=T)
out <- mutualinfo(
    Y1names = 'N1',
    Y2names = 'B1',
    X = NULL,
    learned = learned,
    nsamples = 3600,
    parallel = 4
)
out


set.seed(16)
system.time(
out <- mutualinfo(
    Y1names = 'N1',
    Y2names = 'B1',
    X = NULL,
    learned = learned,
    nsamples = 3600,
    parallel = 1
)
)

set.seed(16)
system.time(
out2 <- mutualinfo2(
    Y1names = 'N1',
    Y2names = 'B1',
    X = NULL,
    learned = learned,
    nsamples = 3600,
    parallel = 1
)
)







### test indentation

sthsth <- function(asnthsnth,
    sthsnth,
    sthsnh) {astnhsnth
    tnsnthsnth
    snthsth
}

a23456 <- if(a23456 > a23456 +
                 a23456 < a23456) {
    a23456
}

a23456 <- function(a23456 +
                       a23456) {
    a23456
}

a23456 <- function(a23456,
    a23456) {a23456,
    a23456
}

a23456 <- function(a23456 +
                       a23456) {a23456 +
                                    a23456
}







object <- ifelse(condition1, out1,
    ifelse(condition2, out2, out3))

fun <- function(argument1,
    argument2
    argument3) {
    body
}

(se)

testa <- matrix(rnorm(64*3600),64,3600)
tests <- abs(testa)
##
testd <- rnorm(128)


system.time(
    x <- vapply(testa, dnorm, mean=testa, sd=testd, FUN.VALUE=testa)
)


system.time(
    x <- sapply(testa, dnorm, mean=testa, sd=testd)
)

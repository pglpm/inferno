#### custom
dataset <- read.csv('~/repos/inferno/development/tests/custom/dataset_custom30.csv', na.strings='')
learnt <- readRDS('~/repos/inferno/development/tests/custom/_newdeletepackagetest-vrt10_dat30_smp120/learnt.rds')
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
    learnt = learnt,
    nsamples = 3600,
    parallel = 4
)
out

out <- mutualinfo(
    Y1names = 'N2vrt',
    Y2names = c('Rvrt', 'Bvrt'),
    X = NULL,
    learnt = learnt,
    nsamples = 3600,
    parallel = 4
)
out


library('inferno')



dataset <- read.csv('~/repos/inferno/development/tests/custom/dataset_custom30.csv', na.strings='')
learnt <- readRDS('~/repos/inferno/development/tests/custom/_newdeletepackagetest-vrt10_dat30_smp120/learnt.rds')
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
    learnt = learnt,
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
    learnt = learnt,
    nsamples = 3600,
    parallel = 1
)
)


identical(out,out2)







#### mtcars
data(mtcars)
dataset <- mtcars
learnt <- readRDS('~/repos/inferno/development/tests/mtcars/_newD_test_mtcars-3-vrt11_dat32_smp512/learnt.rds')
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
    learnt = learnt,
    nsamples = 3600,
    parallel = 4
)
out


#### BN_functdependence
## B1='yes' is correlated with N1=c('A','B')
## B1='no' is correlated with N1=c('C')
learnt <- '~/repos/inferno/development/tests/BN_functdependence/learnt.rds'

load('sysdata.rda')
source('util_mcsubset.R')
source('util_lprob.R')
##
source('util_vtransform.R')
source('_util_vtransform.R')
source('Pr.R')
source('_Pr.R')

source('mutualinfo.R')
source('../development/proto_package/_mutualinfo.R')
set.seed(16)
test0 <- mutualinfo2(
    Y1names = 'B1',
    Y2names = 'N1',
    X = NULL,
    learnt = learnt,
    nsamples = 3600,
    parallel = 1
)
test0
##
set.seed(16)
test1 <- mutualinfo(
    Y1names = 'B1',
    Y2names = 'N1',
    X = NULL,
    learnt = learnt,
    nsamples = 3600,
    parallel = 1
)
test1
for(i in names(test0)[1:5]){
    print(max(2*abs(test0[[i]]-test1[[i]])/abs(test0[[i]]+test1[[i]]), na.rm=T))
}
##
Pr(Y=data.frame(B1=c('yes','no')), X=data.frame(N1=LETTERS[1:3]), learnt=learnt)$values
Pr(X=data.frame(B1=c('yes','no')), Y=data.frame(N1=LETTERS[1:3]), learnt=learnt)$values

source('mutualinfo.R')
source('../development/proto_package/_mutualinfo.R')
set.seed(16)
test0 <- mutualinfo2(
    Y1names = 'N1',
    Y2names = NULL,
    X = data.frame(B1='no'),
    learnt = learnt,
    nsamples = 3600,
    parallel = 1
)
test0
##
set.seed(16)
test1 <- mutualinfo(
    Y1names = 'N1',
    Y2names = NULL,
    X = data.frame(B1='no'),
    learnt = learnt,
    nsamples = 3600,
    parallel = 1
)
test1
for(i in names(test0)[1:5]){
    print(max(2*abs(test0[[i]]-test1[[i]])/abs(test0[[i]]+test1[[i]]), na.rm=T))
}
##
Pr(X=data.frame(B1=c('yes','no')), Y=data.frame(N1=LETTERS[1:3]), learnt=learnt)$values

source('mutualinfo.R')
source('../development/proto_package/_mutualinfo.R')
set.seed(16)
test0 <- mutualinfo2(
    Y1names = 'B1',
    Y2names = NULL,
    X = data.frame(N1='A'),
    learnt = learnt,
    nsamples = 3600,
    parallel = 1
)
test0
##
set.seed(16)
test1 <- mutualinfo(
    Y1names = 'B1',
    Y2names = NULL,
    X = data.frame(N1='A'),
    learnt = learnt,
    nsamples = 3600,
    parallel = 1
)
test1
for(i in names(test0)[1:5]){
    print(max(2*abs(test0[[i]]-test1[[i]])/abs(test0[[i]]+test1[[i]]), na.rm=T))
}


## Knowing B1 almost selects one of the two gaussians making up R1

temp <- Pr(Y=data.frame(R1=(ygrid <- seq(-20,20,length.out=1024))), X=data.frame(B1=c('yes','no')),learnt=learnt)$values
tplot(x=ygrid, y=temp)


source('mutualinfo.R')
source('../development/proto_package/_mutualinfo.R')
set.seed(100)
test0 <- mutualinfo2(
    Y1names = 'R1',
    Y2names = 'B1',
    X = NULL,
    learnt = learnt,
    nsamples = 3600,
    parallel = 1
)
test0
##
set.seed(100)
test1 <- mutualinfo(
    Y1names = 'R1',
    Y2names = 'B1',
    X = NULL,
    learnt = learnt,
    nsamples = 3600,
    parallel = 1
)
test1
for(i in names(test0)[1:5]){
    print(max(2*abs(test0[[i]]-test1[[i]])/abs(test0[[i]]+test1[[i]]), na.rm=T))
}
print('true:')
print(0.5*log(1*2*sqrt(2*pi*exp(1)),2)+0.5*log(2*2*sqrt(2*pi*exp(1)),2)+1)
print(log(1*2*sqrt(2*pi*exp(1)),2))
print(log(2*2*sqrt(2*pi*exp(1)),2)+1)


source('mutualinfo.R')
source('../development/proto_package/_mutualinfo.R')
set.seed(100)
test1 <- mutualinfo(
    Y1names = 'R1',
    Y2names = NULL,
    X = data.frame(B1='yes'),
    learnt = learnt,
    nsamples = 3600,
    parallel = 1
)
test1
## for(i in names(test0)[1:5]){
##     print(max(2*abs(test0[[i]]-test1[[i]])/abs(test0[[i]]+test1[[i]]), na.rm=T))
## }
print('true:')
print(0.5*log(1*2*sqrt(2*pi*exp(1)),2)+0.5*log(2*2*sqrt(2*pi*exp(1)),2)+1)
print(log(1*2*sqrt(2*pi*exp(1)),2))
print(log(2*2*sqrt(2*pi*exp(1)),2))





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

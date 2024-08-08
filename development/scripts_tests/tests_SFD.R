source('samplesFDistribution.R', local=T)
source('util_vtransform.R', local=T)

testSFD <- function(
    Y,
    X,
    mcoutput,
    reduceprobx = TRUE,
    jacobian = TRUE
) {

    auxmetadata <- mcoutput$auxmetadata

    probX <- log(mcoutput$W)
    if(!all(is.na(X)) && !is.null(X)){
        for(i in 1:ncol(X)){
            auxm <- auxmetadata[auxmetadata$name == colnames(X)[i],]
            id <- auxm$id

            if(auxm$mcmctype == 'R'){
                thisx <- c(as.matrix(vtransform(X[,i,drop=F], auxmetadata = auxm,
                    Rout='normalized')))

                probX <- probX +
                        colSums(dnorm(x=thisx,
                            mean=mcoutput$Rmean[id,,,drop=F],
                            sd=sqrt(mcoutput$Rvar[id,,,drop=F]),
                            log=T))
            }

            if(auxm$mcmctype == 'C') {
                thisx <- c(as.matrix(vtransform(X[,i,drop=F],
                    auxmetadata = auxm,
                    Cout='normalized')))
                thisxo <-X[,i]
                tscale <- auxm$tscale
                probX <- probX + colSums(
                    (
                        if(thisxo <= auxm$leftbound){
                            pnorm(q=auxm$tleftbound,
                            mean=mcoutput$Cmean[id,,,drop=F],
                            sd=sqrt(mcoutput$Cvar[id,,,drop=F]), log.p=T)
                        } else if(thisxo >= auxm$rightbound){
                            pnorm(q= -auxm$trightbound,
                            mean= -mcoutput$Cmean[id,,,drop=F],
                            sd= -sqrt(mcoutput$Cvar[id,,,drop=F]), log.p=T)
                        } else {
                        dnorm(x=thisx,
                            mean=mcoutput$Cmean[id,,,drop=F],
                            sd=sqrt(mcoutput$Cvar[id,,,drop=F]), log=T)
                        }
                    )
                )
            }

            if(auxm$mcmctype == 'D') {
                thisx <- c(as.matrix(vtransform(X[,i,drop=F],
                    auxmetadata = auxm,
                    Dout='normalized')))
                thisxo <-X[,i]
                hstep <- auxm$halfstep
                tscale <- auxm$tscale
                probX <- probX + colSums(
                    (
                        if(thisxo + hstep <= auxm$leftbound){
                            pnorm(q=auxm$tleftbound,
                            mean=mcoutput$Dmean[id,,,drop=F],
                            sd=sqrt(mcoutput$Dvar[id,,,drop=F]), log.p=T)
                        } else if(thisxo - hstep >= auxm$rightbound){
                            pnorm(q= -auxm$trightbound,
                            mean= -mcoutput$Dmean[id,,,drop=F],
                            sd= -sqrt(mcoutput$Dvar[id,,,drop=F]), log.p=T)
                        } else if(thisxo + hstep >= auxm$rightbound){
                        pnorm(q= -thisx-hstep/tscale,
                            mean= -mcoutput$Dmean[id,,,drop=F],
                            sd= -sqrt(mcoutput$Dvar[id,,,drop=F]), log.p=T)
                        } else if(thisxo - hstep <= auxm$leftbound){
                        pnorm(q=thisx+hstep/tscale,
                            mean=mcoutput$Dmean[id,,,drop=F],
                            sd=sqrt(mcoutput$Dvar[id,,,drop=F]), log.p=T)
                        } else {
                        log(pnorm(q=thisx+hstep/tscale,
                            mean=mcoutput$Dmean[id,,,drop=F],
                            sd=sqrt(mcoutput$Dvar[id,,,drop=F])) -
                        pnorm(q=thisx-hstep/tscale,
                            mean=mcoutput$Dmean[id,,,drop=F],
                            sd=sqrt(mcoutput$Dvar[id,,,drop=F])))
                        }
                    )
                )
            }

            if(auxm$mcmctype == 'N'){
                thisx <- c(as.matrix(vtransform(X[,i,drop=F], auxmetadata = auxm,
                    Nout='numeric')))
                probX <- probX +
                    colSums(log(
                        mcoutput$Nprob[id,,thisx,,drop=F]
                    ))[,1,]
            }

            if(auxm$mcmctype == 'O'){
                thisx <- c(as.matrix(vtransform(X[,i,drop=F], auxmetadata = auxm,
                    Oout='numeric')))
                probX <- probX +
                    colSums(log(
                        mcoutput$Oprob[id,,thisx,,drop=F]
                    ))[,1,]
            }

            if(auxm$mcmctype == 'B'){
                thisx <- c(as.matrix(vtransform(X[,i,drop=F], auxmetadata = auxm,
                    Bout='numeric')))
                probX <- probX +
                    colSums(dbinom(
                        x=thisx,
                        prob=mcoutput$Bprob[id,,,drop=F],
                        size=1, log=T
                    )
                    )
            }

        } # end loop through X columns
    } # end X calculation

    jaco <- 1
    probY <- 0
        for(i in 1:ncol(Y)){
            auxm <- auxmetadata[auxmetadata$name == colnames(Y)[i],]
            id <- auxm$id

            if(auxm$mcmctype == 'R'){
                thisx <- c(as.matrix(vtransform(Y[,i,drop=F], auxmetadata = auxm,
                    Rout='normalized')))

                probY <- probY +
                        colSums(dnorm(x=thisx,
                            mean=mcoutput$Rmean[id,,,drop=F],
                            sd=sqrt(mcoutput$Rvar[id,,,drop=F]),
                            log=T))

                jaco <- jaco *
                    c(as.matrix(vtransform(Y[,i,drop=F], auxmetadata = auxm,
                    invjacobian=T)))
            }

            if(auxm$mcmctype == 'C') {
                thisx <- c(as.matrix(vtransform(Y[,i,drop=F],
                    auxmetadata = auxm,
                    Cout='normalized')))
                thisxo <-Y[,i]
                tscale <- auxm$tscale
                probY <- probY + colSums(
                    (
                        if(thisxo <= auxm$leftbound){
                            pnorm(q=auxm$tleftbound,
                            mean=mcoutput$Cmean[id,,,drop=F],
                            sd=sqrt(mcoutput$Cvar[id,,,drop=F]), log.p=T)
                        } else if(thisxo >= auxm$rightbound){
                            pnorm(q= -auxm$trightbound,
                            mean= -mcoutput$Cmean[id,,,drop=F],
                            sd= -sqrt(mcoutput$Cvar[id,,,drop=F]), log.p=T)
                        } else {
                        dnorm(x=thisx,
                            mean=mcoutput$Cmean[id,,,drop=F],
                            sd=sqrt(mcoutput$Cvar[id,,,drop=F]), log=T)
                        }
                    )
                )

                jaco <- jaco *
                    (if(thisxo <= auxm$leftbound ||
                        thisxo >= auxm$rightbound) {
                        1
                    } else {
                    c(as.matrix(vtransform(Y[,i,drop=F], auxmetadata = auxm,
                        invjacobian=T)))
                    })
            }

            if(auxm$mcmctype == 'D') {
                thisx <- c(as.matrix(vtransform(Y[,i,drop=F],
                    auxmetadata = auxm,
                    Dout='normalized')))
                thisxo <-Y[,i]
                hstep <- auxm$halfstep
                tscale <- auxm$tscale
                probY <- probY + colSums(
                    (
                        if(thisxo + hstep <= auxm$leftbound){
                            pnorm(q=auxm$tleftbound,
                            mean=mcoutput$Dmean[id,,,drop=F],
                            sd=sqrt(mcoutput$Dvar[id,,,drop=F]), log.p=T)
                        } else if(thisxo - hstep >= auxm$rightbound){
                            pnorm(q= -auxm$trightbound,
                            mean= -mcoutput$Dmean[id,,,drop=F],
                            sd= -sqrt(mcoutput$Dvar[id,,,drop=F]), log.p=T)
                        } else if(thisxo + hstep >= auxm$rightbound){
                        pnorm(q= -thisx-hstep/tscale,
                            mean= -mcoutput$Dmean[id,,,drop=F],
                            sd= -sqrt(mcoutput$Dvar[id,,,drop=F]), log.p=T)
                        } else if(thisxo - hstep <= auxm$leftbound){
                        pnorm(q=thisx+hstep/tscale,
                            mean=mcoutput$Dmean[id,,,drop=F],
                            sd=sqrt(mcoutput$Dvar[id,,,drop=F]), log.p=T)
                        } else {
                        log(pnorm(q=thisx+hstep/tscale,
                            mean=mcoutput$Dmean[id,,,drop=F],
                            sd=sqrt(mcoutput$Dvar[id,,,drop=F])) -
                        pnorm(q=thisx-hstep/tscale,
                            mean=mcoutput$Dmean[id,,,drop=F],
                            sd=sqrt(mcoutput$Dvar[id,,,drop=F])))
                        }
                    )
                )
            }

            if(auxm$mcmctype == 'N'){
                thisx <- c(as.matrix(vtransform(Y[,i,drop=F], auxmetadata = auxm,
                    Nout='numeric')))
                probY <- probY +
                    colSums(log(
                        mcoutput$Nprob[id,,thisx,,drop=F]
                    ))[,1,]
            }

            if(auxm$mcmctype == 'O'){
                thisx <- c(as.matrix(vtransform(Y[,i,drop=F], auxmetadata = auxm,
                    Oout='numeric')))
                probY <- probY +
                    colSums(log(
                        mcoutput$Oprob[id,,thisx,,drop=F]
                    ))[,1,]
            }

            if(auxm$mcmctype == 'B'){
                thisx <- c(as.matrix(vtransform(Y[,i,drop=F], auxmetadata = auxm,
                    Bout='numeric')))
                probY <- probY +
                    colSums(dbinom(
                        x=thisx,
                        prob=mcoutput$Bprob[id,,,drop=F],
                        size=1, log=T
                    )
                    )
            }

        } # end loop through Y columns

    if(reduceprobx){
        probX <- apply(probX, 2, function(xx) {
            xx - max(xx[is.finite(xx)])
        })
    }
    colSums(exp(probX + probY)) / colSums(exp(probX)) /
        (if(jacobian) { jaco } else { 1 })
}


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

#############################################################
## Combs YX
#############################################################

## nadpark
dataset <- read.csv('/home/pglpm/repositories/parkinsonbayes/data/nadpark/toydata.csv', na.strings='')[,-1]
mcoutput <- readRDS('/home/pglpm/repositories/parkinsonbayes/data/nadpark/_testtoy_analysis_2_3A4-1-vrt7_dat30_smp1024/Fdistribution.rds')
jac <- F
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


## mtcars
data(mtcars)
dataset <- mtcars
mcoutput <- readRDS('/home/pglpm/repositories/bayes_nonparametric_inference/development/tests/mtcars/_newD_test_mtcars-3-vrt11_dat32_smp512/Fdistribution.rds')
jac <- F
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



out <- t(sapply(1:nrow(dataset), function(ipoint){
    sapply(1:ncol(dataset), function(i){
        v0 <- samplesFDistribution(Y=dataset[ipoint,i,drop=F], X=NULL, mcoutput = mcoutput, jacobian = jac, silent=T, parallel=4)
        v1 <- testSFD(Y=dataset[ipoint,i,drop=F], X=NULL, mcoutput = mcoutput, jacobian = jac)
        max(abs((v1-v0)/(v1+v0)*2))
    })
}))
colnames(out) <- mcoutput$auxmetadata$mcmctype
apply(abs(out),2,max)
max(abs(out))




source('util_vtransform.R', local=T)
source('samplesFDistribution.R', local=T)
source('mutualinfo.R', local=T)
data(mtcars)
dataset <- mtcars
mcoutput <- readRDS('/home/pglpm/repositories/bayes_nonparametric_inference/development/tests/mtcars/_newD_test_mtcars-3-vrt11_dat32_smp512/Fdistribution.rds')
##
Y1s <- sample(colnames(mtcars),2)
Y2s <- sample(setdiff(colnames(mtcars), Y1s),2)
Xs <- sample(setdiff(colnames(mtcars), c(Y1s,Y2s)),2)
mutualinfo(
    Y1names = Y1s,
    Y2names = Y2s,
    X = mtcars[1, Xs, drop=F],
    mcoutput = mcoutput,
    nsamples = 3600,
    parallel = 4
)

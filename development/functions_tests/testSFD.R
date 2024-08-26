source('samplesFDistribution.R', local=T)
source('util_vtransform.R', local=T)
testSFD <- function(
    Y,
    X,
    learnt,
    reduceprobx = TRUE,
    jacobian = TRUE
) {

    auxmetadata <- learnt$auxmetadata

    probX <- log(learnt$W)
    if(!all(is.na(X)) && !is.null(X)){
        for(i in 1:ncol(X)){
            auxm <- auxmetadata[auxmetadata$name == colnames(X)[i],]
            id <- auxm$id

            if(auxm$mcmctype == 'R'){
                thisx <- c(as.matrix(vtransform(X[,i,drop=F], auxmetadata = auxm,
                    Rout='normalized')))

                probX <- probX +
                        colSums(dnorm(x=thisx,
                            mean=learnt$Rmean[id,,,drop=F],
                            sd=sqrt(learnt$Rvar[id,,,drop=F]),
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
                        if(thisxo <= auxm$domainminplushs){
                            pnorm(q=auxm$tdomainminplushs,
                            mean=learnt$Cmean[id,,,drop=F],
                            sd=sqrt(learnt$Cvar[id,,,drop=F]), log.p=T)
                        } else if(thisxo >= auxm$domainmaxminushs){
                            pnorm(q= -auxm$tdomainmaxminushs,
                            mean= -learnt$Cmean[id,,,drop=F],
                            sd= -sqrt(learnt$Cvar[id,,,drop=F]), log.p=T)
                        } else {
                        dnorm(x=thisx,
                            mean=learnt$Cmean[id,,,drop=F],
                            sd=sqrt(learnt$Cvar[id,,,drop=F]), log=T)
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
                        if(thisxo + hstep <= auxm$domainminplushs){
                            pnorm(q=auxm$tdomainminplushs,
                            mean=learnt$Dmean[id,,,drop=F],
                            sd=sqrt(learnt$Dvar[id,,,drop=F]), log.p=T)
                        } else if(thisxo - hstep >= auxm$domainmaxminushs){
                            pnorm(q= -auxm$tdomainmaxminushs,
                            mean= -learnt$Dmean[id,,,drop=F],
                            sd= -sqrt(learnt$Dvar[id,,,drop=F]), log.p=T)
                        } else if(thisxo + hstep >= auxm$domainmaxminushs){
                        pnorm(q= -thisx-hstep/tscale,
                            mean= -learnt$Dmean[id,,,drop=F],
                            sd= -sqrt(learnt$Dvar[id,,,drop=F]), log.p=T)
                        } else if(thisxo - hstep <= auxm$domainminplushs){
                        pnorm(q=thisx+hstep/tscale,
                            mean=learnt$Dmean[id,,,drop=F],
                            sd=sqrt(learnt$Dvar[id,,,drop=F]), log.p=T)
                        } else {
                        log(pnorm(q=thisx+hstep/tscale,
                            mean=learnt$Dmean[id,,,drop=F],
                            sd=sqrt(learnt$Dvar[id,,,drop=F])) -
                        pnorm(q=thisx-hstep/tscale,
                            mean=learnt$Dmean[id,,,drop=F],
                            sd=sqrt(learnt$Dvar[id,,,drop=F])))
                        }
                    )
                )
            }

            if(auxm$mcmctype == 'N'){
                thisx <- c(as.matrix(vtransform(X[,i,drop=F], auxmetadata = auxm,
                    Nout='numeric')))
                probX <- probX +
                    colSums(log(
                        learnt$Nprob[id,,thisx,,drop=F]
                    ))[,1,]
            }

            if(auxm$mcmctype == 'O'){
                thisx <- c(as.matrix(vtransform(X[,i,drop=F], auxmetadata = auxm,
                    Oout='numeric')))
                probX <- probX +
                    colSums(log(
                        learnt$Oprob[id,,thisx,,drop=F]
                    ))[,1,]
            }

            if(auxm$mcmctype == 'B'){
                thisx <- c(as.matrix(vtransform(X[,i,drop=F], auxmetadata = auxm,
                    Bout='numeric')))
                probX <- probX +
                    colSums(dbinom(
                        x=thisx,
                        prob=learnt$Bprob[id,,,drop=F],
                        size=1, log=T
                    )
                    )
            }

        } # end loop through X columns
    } # end X calculation

    ljaco <- 0
    probY <- 0
        for(i in 1:ncol(Y)){
            auxm <- auxmetadata[auxmetadata$name == colnames(Y)[i],]
            id <- auxm$id

            if(auxm$mcmctype == 'R'){
                thisx <- c(as.matrix(vtransform(Y[,i,drop=F], auxmetadata = auxm,
                    Rout='normalized')))

                probY <- probY +
                        colSums(dnorm(x=thisx,
                            mean=learnt$Rmean[id,,,drop=F],
                            sd=sqrt(learnt$Rvar[id,,,drop=F]),
                            log=T))

                ljaco <- ljaco +
                    log(c(as.matrix(vtransform(Y[,i,drop=F], auxmetadata = auxm,
                    invjacobian=TRUE))))
            }

            if(auxm$mcmctype == 'C') {
                thisx <- c(as.matrix(vtransform(Y[,i,drop=F],
                    auxmetadata = auxm,
                    Cout='normalized')))
                thisxo <-Y[,i]
                tscale <- auxm$tscale
                probY <- probY + colSums(
                    (
                        if(thisxo <= auxm$domainminplushs){
                            pnorm(q=auxm$tdomainminplushs,
                            mean=learnt$Cmean[id,,,drop=F],
                            sd=sqrt(learnt$Cvar[id,,,drop=F]), log.p=T)
                        } else if(thisxo >= auxm$domainmaxminushs){
                            pnorm(q= -auxm$tdomainmaxminushs,
                            mean= -learnt$Cmean[id,,,drop=F],
                            sd= -sqrt(learnt$Cvar[id,,,drop=F]), log.p=T)
                        } else {
                        dnorm(x=thisx,
                            mean=learnt$Cmean[id,,,drop=F],
                            sd=sqrt(learnt$Cvar[id,,,drop=F]), log=T)
                        }
                    )
                )

                ljaco <- ljaco + log(
                    (if(thisxo <= auxm$domainminplushs ||
                        thisxo >= auxm$domainmaxminushs) {
                        1
                    } else {
                    c(as.matrix(vtransform(Y[,i,drop=F], auxmetadata = auxm,
                        invjacobian=TRUE)))
                    }))
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
                        if(thisxo + hstep <= auxm$domainminplushs){
                            pnorm(q=auxm$tdomainminplushs,
                            mean=learnt$Dmean[id,,,drop=F],
                            sd=sqrt(learnt$Dvar[id,,,drop=F]), log.p=T)
                        } else if(thisxo - hstep >= auxm$domainmaxminushs){
                            pnorm(q= -auxm$tdomainmaxminushs,
                            mean= -learnt$Dmean[id,,,drop=F],
                            sd= -sqrt(learnt$Dvar[id,,,drop=F]), log.p=T)
                        } else if(thisxo + hstep >= auxm$domainmaxminushs){
                        pnorm(q= -thisx-hstep/tscale,
                            mean= -learnt$Dmean[id,,,drop=F],
                            sd= -sqrt(learnt$Dvar[id,,,drop=F]), log.p=T)
                        } else if(thisxo - hstep <= auxm$domainminplushs){
                        pnorm(q=thisx+hstep/tscale,
                            mean=learnt$Dmean[id,,,drop=F],
                            sd=sqrt(learnt$Dvar[id,,,drop=F]), log.p=T)
                        } else {
                        log(pnorm(q=thisx+hstep/tscale,
                            mean=learnt$Dmean[id,,,drop=F],
                            sd=sqrt(learnt$Dvar[id,,,drop=F])) -
                        pnorm(q=thisx-hstep/tscale,
                            mean=learnt$Dmean[id,,,drop=F],
                            sd=sqrt(learnt$Dvar[id,,,drop=F])))
                        }
                    )
                )
            }

            if(auxm$mcmctype == 'N'){
                thisx <- c(as.matrix(vtransform(Y[,i,drop=F], auxmetadata = auxm,
                    Nout='numeric')))
                probY <- probY +
                    colSums(log(
                        learnt$Nprob[id,,thisx,,drop=F]
                    ))[,1,]
            }

            if(auxm$mcmctype == 'O'){
                thisx <- c(as.matrix(vtransform(Y[,i,drop=F], auxmetadata = auxm,
                    Oout='numeric')))
                probY <- probY +
                    colSums(log(
                        learnt$Oprob[id,,thisx,,drop=F]
                    ))[,1,]
            }

            if(auxm$mcmctype == 'B'){
                thisx <- c(as.matrix(vtransform(Y[,i,drop=F], auxmetadata = auxm,
                    Bout='numeric')))
                probY <- probY +
                    colSums(dbinom(
                        x=thisx,
                        prob=learnt$Bprob[id,,,drop=F],
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
    colSums(exp(probX + probY - (if(jacobian){ljaco}else{0}))) / colSums(exp(probX))
}

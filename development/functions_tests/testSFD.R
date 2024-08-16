source('samplesFDistribution.R', local=T)
source('util_vtransform.R', local=T)
testSFD <- function(
    Y,
    X,
    agent,
    reduceprobx = TRUE,
    jacobian = TRUE
) {

    auxmetadata <- agent$auxmetadata

    probX <- log(agent$W)
    if(!all(is.na(X)) && !is.null(X)){
        for(i in 1:ncol(X)){
            auxm <- auxmetadata[auxmetadata$name == colnames(X)[i],]
            id <- auxm$id

            if(auxm$mcmctype == 'R'){
                thisx <- c(as.matrix(vtransform(X[,i,drop=F], auxmetadata = auxm,
                    Rout='normalized')))

                probX <- probX +
                        colSums(dnorm(x=thisx,
                            mean=agent$Rmean[id,,,drop=F],
                            sd=sqrt(agent$Rvar[id,,,drop=F]),
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
                            mean=agent$Cmean[id,,,drop=F],
                            sd=sqrt(agent$Cvar[id,,,drop=F]), log.p=T)
                        } else if(thisxo >= auxm$rightbound){
                            pnorm(q= -auxm$trightbound,
                            mean= -agent$Cmean[id,,,drop=F],
                            sd= -sqrt(agent$Cvar[id,,,drop=F]), log.p=T)
                        } else {
                        dnorm(x=thisx,
                            mean=agent$Cmean[id,,,drop=F],
                            sd=sqrt(agent$Cvar[id,,,drop=F]), log=T)
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
                            mean=agent$Dmean[id,,,drop=F],
                            sd=sqrt(agent$Dvar[id,,,drop=F]), log.p=T)
                        } else if(thisxo - hstep >= auxm$rightbound){
                            pnorm(q= -auxm$trightbound,
                            mean= -agent$Dmean[id,,,drop=F],
                            sd= -sqrt(agent$Dvar[id,,,drop=F]), log.p=T)
                        } else if(thisxo + hstep >= auxm$rightbound){
                        pnorm(q= -thisx-hstep/tscale,
                            mean= -agent$Dmean[id,,,drop=F],
                            sd= -sqrt(agent$Dvar[id,,,drop=F]), log.p=T)
                        } else if(thisxo - hstep <= auxm$leftbound){
                        pnorm(q=thisx+hstep/tscale,
                            mean=agent$Dmean[id,,,drop=F],
                            sd=sqrt(agent$Dvar[id,,,drop=F]), log.p=T)
                        } else {
                        log(pnorm(q=thisx+hstep/tscale,
                            mean=agent$Dmean[id,,,drop=F],
                            sd=sqrt(agent$Dvar[id,,,drop=F])) -
                        pnorm(q=thisx-hstep/tscale,
                            mean=agent$Dmean[id,,,drop=F],
                            sd=sqrt(agent$Dvar[id,,,drop=F])))
                        }
                    )
                )
            }

            if(auxm$mcmctype == 'N'){
                thisx <- c(as.matrix(vtransform(X[,i,drop=F], auxmetadata = auxm,
                    Nout='numeric')))
                probX <- probX +
                    colSums(log(
                        agent$Nprob[id,,thisx,,drop=F]
                    ))[,1,]
            }

            if(auxm$mcmctype == 'O'){
                thisx <- c(as.matrix(vtransform(X[,i,drop=F], auxmetadata = auxm,
                    Oout='numeric')))
                probX <- probX +
                    colSums(log(
                        agent$Oprob[id,,thisx,,drop=F]
                    ))[,1,]
            }

            if(auxm$mcmctype == 'B'){
                thisx <- c(as.matrix(vtransform(X[,i,drop=F], auxmetadata = auxm,
                    Bout='numeric')))
                probX <- probX +
                    colSums(dbinom(
                        x=thisx,
                        prob=agent$Bprob[id,,,drop=F],
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
                            mean=agent$Rmean[id,,,drop=F],
                            sd=sqrt(agent$Rvar[id,,,drop=F]),
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
                        if(thisxo <= auxm$leftbound){
                            pnorm(q=auxm$tleftbound,
                            mean=agent$Cmean[id,,,drop=F],
                            sd=sqrt(agent$Cvar[id,,,drop=F]), log.p=T)
                        } else if(thisxo >= auxm$rightbound){
                            pnorm(q= -auxm$trightbound,
                            mean= -agent$Cmean[id,,,drop=F],
                            sd= -sqrt(agent$Cvar[id,,,drop=F]), log.p=T)
                        } else {
                        dnorm(x=thisx,
                            mean=agent$Cmean[id,,,drop=F],
                            sd=sqrt(agent$Cvar[id,,,drop=F]), log=T)
                        }
                    )
                )

                ljaco <- ljaco + log(
                    (if(thisxo <= auxm$leftbound ||
                        thisxo >= auxm$rightbound) {
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
                        if(thisxo + hstep <= auxm$leftbound){
                            pnorm(q=auxm$tleftbound,
                            mean=agent$Dmean[id,,,drop=F],
                            sd=sqrt(agent$Dvar[id,,,drop=F]), log.p=T)
                        } else if(thisxo - hstep >= auxm$rightbound){
                            pnorm(q= -auxm$trightbound,
                            mean= -agent$Dmean[id,,,drop=F],
                            sd= -sqrt(agent$Dvar[id,,,drop=F]), log.p=T)
                        } else if(thisxo + hstep >= auxm$rightbound){
                        pnorm(q= -thisx-hstep/tscale,
                            mean= -agent$Dmean[id,,,drop=F],
                            sd= -sqrt(agent$Dvar[id,,,drop=F]), log.p=T)
                        } else if(thisxo - hstep <= auxm$leftbound){
                        pnorm(q=thisx+hstep/tscale,
                            mean=agent$Dmean[id,,,drop=F],
                            sd=sqrt(agent$Dvar[id,,,drop=F]), log.p=T)
                        } else {
                        log(pnorm(q=thisx+hstep/tscale,
                            mean=agent$Dmean[id,,,drop=F],
                            sd=sqrt(agent$Dvar[id,,,drop=F])) -
                        pnorm(q=thisx-hstep/tscale,
                            mean=agent$Dmean[id,,,drop=F],
                            sd=sqrt(agent$Dvar[id,,,drop=F])))
                        }
                    )
                )
            }

            if(auxm$mcmctype == 'N'){
                thisx <- c(as.matrix(vtransform(Y[,i,drop=F], auxmetadata = auxm,
                    Nout='numeric')))
                probY <- probY +
                    colSums(log(
                        agent$Nprob[id,,thisx,,drop=F]
                    ))[,1,]
            }

            if(auxm$mcmctype == 'O'){
                thisx <- c(as.matrix(vtransform(Y[,i,drop=F], auxmetadata = auxm,
                    Oout='numeric')))
                probY <- probY +
                    colSums(log(
                        agent$Oprob[id,,thisx,,drop=F]
                    ))[,1,]
            }

            if(auxm$mcmctype == 'B'){
                thisx <- c(as.matrix(vtransform(Y[,i,drop=F], auxmetadata = auxm,
                    Bout='numeric')))
                probY <- probY +
                    colSums(dbinom(
                        x=thisx,
                        prob=agent$Bprob[id,,,drop=F],
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

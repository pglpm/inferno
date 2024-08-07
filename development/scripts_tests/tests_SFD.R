source('samplesFDistribution.R', local=T)
source('util_vtransform.R', local=T)
testSFD <- function(
    Y,
    X,
    mcoutput,
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
                            log=T
                        ))
            }

            if(auxm$mcmctype == 'D') {
                thisx <- c(as.matrix(vtransform(X[,i,drop=F],
                    auxmetadata = auxm,
                    Dout='normalized')))
                thisxo <-X[,i]
                hstep <- auxm$halfstep
                tscale <- auxm$tscale
                probX <- probX + colSums(log(
                    (
                        if(thisxo + hstep <= auxm$leftbound){
                            pnorm(q=auxm$tleftbound,
                            mean=mcoutput$Dmean[id,,,drop=F],
                            sd=sqrt(mcoutput$Dvar[id,,,drop=F]))
                        } else if(thisxo - hstep >= auxm$rightbound){
                            pnorm(q=-auxm$trightbound,
                            mean=-mcoutput$Dmean[id,,,drop=F],
                            sd=-sqrt(mcoutput$Dvar[id,,,drop=F]))
                        } else if(thisxo + hstep >= auxm$rightbound){
                        pnorm(q=-thisx-hstep/tscale,
                            mean=-mcoutput$Dmean[id,,,drop=F],
                            sd=-sqrt(mcoutput$Dvar[id,,,drop=F]))
                        } else if(thisxo - hstep <= auxm$leftbound){
                        pnorm(q=thisx+hstep/tscale,
                            mean=mcoutput$Dmean[id,,,drop=F],
                            sd=sqrt(mcoutput$Dvar[id,,,drop=F]))
                        } else {
                        pnorm(q=thisx+hstep/tscale,
                            mean=mcoutput$Dmean[id,,,drop=F],
                            sd=sqrt(mcoutput$Dvar[id,,,drop=F])) -
                        pnorm(q=thisx-hstep/tscale,
                            mean=mcoutput$Dmean[id,,,drop=F],
                            sd=sqrt(mcoutput$Dvar[id,,,drop=F]))
                        }
                    )
                ))
            }

            if(auxm$mcmctype == 'N'){
                thisx <- c(as.matrix(vtransform(X[,i,drop=F], auxmetadata = auxm,
                    Nout='numeric')))
                probX <- probX +
                    colSums(log(
                        mcoutput$Nprob[id,,thisx,,drop=F]
                    ))[,1,]
            }
        }
    }
    probX
}




ipoint <- 100
sapply(1:ncol(iris), function(i){
    v0 <- samplesFDistribution(Y=iris[ipoint,i,drop=F], X=NULL, mcoutput = testf, jacobian = T, silent=T)
    v1 <- colSums(exp(testSFD(X=iris[ipoint,i,drop=F], Y=NULL, mcoutput = testf)))
    max(abs((v1-v0)/(v1+v0)*2))
    })

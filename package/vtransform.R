## Transformation from variate to internal variable
vtransform <- function(x, varinfoaux, Cout='init', Dout='data', Oout='data', Bout='numeric', Nout='numeric', Ofunction='Ofunction512', variates=NULL){
    x <- cbind(x)
    if(!is.null(variates)){colnames(x) <- variates}
    matrix(sapply(colnames(x), function(v){
        ##
        datum <- unlist(x[,v,with=F])
        info <- as.list(varinfoaux[name == v])
        ##
        if(info$mcmctype == 'R'){
            if (info$transform == 'log'){
                datum <- log(datum-info$domainmin)
            }else if (info$transform == 'logminus'){
                datum <- log(info$domainmax-datum)
            }else if (info$transform == 'probit'){
                datum <- qnorm((datum-info$domainmin)/(info$domainmax-info$domainmin))
            }
            datum <- (datum-info$tlocation)/info$tscale
            ##
        } else if(info$mcmctype == 'O'){ # ordinal
            if(is.character(Ofunction)){
                Ofunction <- readRDS(paste0(Ofunction,'.rds'))
            }
            datum <- round((datum-info$tlocation)/info$tscale) # output is in range 0 to n-1
            if(Oout == 'init'){ # in sampling functions or init MCMC
                datum[is.na(datum)] <- info$Nvalues/2+0.5
                datum <- Ofunction((datum-0.5)/info$Nvalues)
            } else if(Oout == 'left'){ # as left for MCMC
                datum <- Ofunction(pmax(0,datum-1)/info$Nvalues)
                datum[is.na(datum)] <- -Inf
            } else if(Oout == 'right'){ # as right for MCMC
                datum <- Ofunction(pmin(info$Nvalues,datum)/info$Nvalues)
                datum[is.na(datum)] <- +Inf
            } else if(Oout == 'aux'){ # aux variable in MCMC
                sel <- is.na(datum)
                datum[sel] <- NA
                datum[!sel] <- 1L
            } else if(Oout == 'index'){ # in output functions
                datum <- datum-1L
            } else if(Oout == 'original'){ # in output functions
                datum <- datum * info$tscale + info$tlocation
            }
            ##
        } else if(info$mcmctype == 'D'){ # discretized
            leftbound <- datum - info$step
            rightbound <- datum + info$step
            if (info$transform == 'log'){
                datum <- log(datum-info$domainmin)
                leftbound <- log(leftbound-info$domainmin)
                rightbound <- log(rightbound-info$domainmin)
            } else if (info$transform == 'logminus'){
                datum <- log(info$domainmax-datum)
                leftbound <- log(info$domainmax-leftbound)
                rightbound <- log(info$domainmax-rightbound)
            } else if (info$transform == 'probit'){
                datum <- qnorm((datum-info$domainmin)/(info$domainmax-info$domainmin))
                leftbound <- qnorm((leftbound-info$domainmin)/(info$domainmax-info$domainmin))
                rightbound <- qnorm((rightbound-info$domainmin)/(info$domainmax-info$domainmin))
            }
            datum <- (datum-info$tlocation)/info$tscale
            rightbound <- (rightbound-info$tlocation)/info$tscale
            leftbound <- (leftbound-info$tlocation)/info$tscale
            if(Dout == 'left'){
                sel <- is.na(datum)
                datum[sel] <- -Inf
                datum[!sel] <- leftbound[!sel]
            } else if(Dout == 'right'){
                sel <- is.na(datum)
                datum[sel] <- +Inf
                datum[!sel] <- rightbound[!sel]
            } else if(Dout == 'init'){ #init in MCMC
                datum[is.na(datum)] <- 0L
            } else if(Dout == 'aux'){ # aux variable in MCMC
                sel <- is.na(datum)
                datum[sel] <- NA
                datum[!sel] <- 1L
            }
            ##
        } else if(info$mcmctype == 'C'){ # censored
            leftbound <- info$censormin
            rightbound <- info$censormax
            if (info$transform == 'log'){
                datum <- log(datum-info$domainmin)
                leftbound <- log(leftbound-info$domainmin)
                rightbound <- log(rightbound-info$domainmin)
            } else if (info$transform == 'logminus'){
                datum <- log(info$domainmax-datum)
                leftbound <- log(info$domainmax-leftbound)
                rightbound <- log(info$domainmax-rightbound)
            } else if (info$transform == 'probit'){
                datum <- qnorm((datum-info$domainmin)/(info$domainmax-info$domainmin))
                leftbound <- qnorm((leftbound-info$domainmin)/(info$domainmax-info$domainmin))
                rightbound <- qnorm((rightbound-info$domainmin)/(info$domainmax-info$domainmin))
            }
            datum <- (datum-info$tlocation)/info$tscale
            rightbound <- (rightbound-info$tlocation)/info$tscale
            leftbound <- (leftbound-info$tlocation)/info$tscale
            xv <- data.matrix(x[,..v])
            if(Cout == 'left'){ # in MCMC
                sel <- is.na(datum) | (xv < info$censormax)
                datum[sel] <- -Inf
                datum[!sel] <- rightbound
            } else if(Cout == 'right'){ # in MCMC
                sel <- is.na(datum) | (xv > info$censormin)
                datum[sel] <- +Inf
                datum[!sel] <- leftbound
            } else if(Cout == 'lat'){ # latent variable in MCMC
                sel <- is.na(datum) | (xv >= info$censormax) | (xv <= info$censormin)
                datum[sel] <- NA
            } else if(Cout == 'init'){ #init in MCMC
                datum[is.na(datum)] <- 0L
                datum[!is.na(datum)] <- NA
                datum[xv >= info$censormax] <- rightbound + 0.125
                datum[xv <= info$censormin] <- leftbound - 0.125
            } else if(Cout == 'aux'){ # aux variable in MCMC
                sel <- is.na(datum)
                datum[sel] <- NA
                datum[!sel] <- 1L
            } else if(Cout == 'index'){ #in sampling functions
                datum[xv >= info$censormax] <- +Inf
                datum[xv <= info$censormin] <- -Inf
            }
        } else if(info$mcmctype == 'B'){ # binary
            bvalues <- 0:1
            names(bvalues) <- unlist(info[c('V1','V2')])
            if(Bout == 'numeric'){
                datum <- bvalues[as.character(datum)]
            } else if(Bout == 'original'){
                datum <- names(bvalues[datum+1])
            }
            ##
        } else if(info$mcmctype == 'N'){ # nominal
            bvalues <- 1:info$Nvalues
            names(bvalues) <- unlist(info[paste0('V', bvalues)])
            if(Nout == 'numeric'){
                datum <- bvalues[as.character(datum)]
            } else if(Nout == 'original'){
                datum <- names(bvalues[datum])
            }
        }
            ##
        datum
    }),ncol=ncol(x),dimnames=list(NULL,colnames(x)))
}

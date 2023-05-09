## Transformation from variate to internal variable
vtransform <- function(x, varinfoaux, Cout='init', Dout='data', Oout='data', Bout='data', Ofunction='Ofunction512', variates=NULL){
    x <- cbind(data.matrix(x))
    if(!is.null(variates)){colnames(x) <- variates}
    matrix(sapply(colnames(x), function(v){
        ##
        datum <- x[,v,drop=F]
        info <- as.list(varinfoaux[name == v])
        ##
        if(info$mcmctype == 'R'){
            if (info$transform == 'log'){
                datum <- log(datum-info$domainmin)
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
            } else if(Oout == 'index'){ # in output functions
                datum <- datum-1L
            } else if(Oout == 'correct'){ # in output functions
                datum <- datum * info$tscale + info$tlocation
            }
            ##
        } else if(info$mcmctype == 'B'){ # binary
            datum <- round((datum-info$tlocation)/info$tscale) # output in range 0 to 1
            if(Bout == 'correct'){ # in output functions
                datum <- datum * info$tscale + info$tlocation
            }
            ##
        } else if(info$mcmctype == 'D'){ # censored
            if (info$transform == 'log'){
                datum <- log(datum-info$domainmin)
                rightbound <- log(info$tmax-info$domainmin)
                leftbound <- log(info$tmin-info$domainmin)
            }else if (info$transform == 'probit'){
                datum <- qnorm((datum-info$domainmin)/(info$domainmax-info$domainmin))
                rightbound <- qnorm((info$tmax-info$domainmin)/(info$domainmax-info$domainmin))
                leftbound <- qnorm((info$tmin-info$domainmin)/(info$domainmax-info$domainmin))
            }
            datum <- (datum-info$tlocation)/info$tscale
            rightbound <- (rightbound-info$tlocation)/info$tscale
            leftbound <- (leftbound-info$tlocation)/info$tscale
            if(Dout == 'left'){ # in MCMC
                sel <- is.na(datum) | (x[,v] < info$tmax)
                datum[sel] <- -Inf
                datum[!sel] <- rightbound
            } else if(Dout == 'right'){ # in MCMC
                sel <- is.na(datum) | (x[,v] > info$tmin)
                datum[sel] <- +Inf
                datum[!sel] <- leftbound
            } else if(Dout == 'data'){ # data in MCMC
                sel <- is.na(datum) | (x[,v] >= info$tmax) | (x[,v] <= info$tmin)
                datum[sel] <- NA
            } else if(Dout == 'init'){ #init in MCMC
                datum[is.na(datum)] <- 0L
                datum[x[,v] >= info$tmax] <- rightbound + 0.125
                datum[x[,v] <= info$tmin] <- leftbound - 0.125
            } else if(Dout == 'index'){ #in sampling functions
                datum[x[,v] >= info$tmax] <- +Inf
                datum[x[,v] <= info$tmin] <- -Inf
            }
            ##
        ## } else if(info$mcmctype == 'O'){ # one-censored
        ##     if (info$transform == 'log'){
        ##         datum <- log(datum-info$domainmin)
        ##         rightbound <- log(info$tmax-info$domainmin)
        ##     }else if (info$transform == 'probit'){
        ##         datum <- qnorm((datum-info$domainmin)/(info$domainmax-info$domainmin))
        ##         rightbound <- qnorm((info$tmax-info$domainmin)/(info$domainmax-info$domainmin))
        ##     }
        ##     datum <- (datum-info$tlocation)/info$tscale
        ##     rightbound <- (rightbound-info$tlocation)/info$tscale
        ##     if(Oout == 'left'){ # in MCMC
        ##         sel <- is.na(datum) | (x[,v] < info$tmax)
        ##         datum[sel] <- -Inf
        ##         datum[!sel] <- rightbound
        ##     } else if(Oout == 'data'){ # data in MCMC
        ##         sel <- is.na(datum) | (x[,v] >= info$tmax)
        ##         datum[sel] <- NA
        ##     } else if(Oout == 'init'){ #init in MCMC
        ##         datum[is.na(datum)] <- 0L
        ##         datum[x[,v] >= info$tmax] <- rightbound+0.125
        ##     } else if(Oout == 'index'){ #in sampling functions
        ##         datum[x[,v] >= info$tmax] <- +Inf
        ##     }
        ##     ##
        }
        datum
    }),ncol=ncol(x),dimnames=list(NULL,colnames(x)))
}

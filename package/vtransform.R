## Transformation from variate to internal variable
vtransform <- function(x, auxmetadata, Cout='init', Dout='data', Oout='data', Bout='numeric', Nout='numeric', Qfunction='Qfunction8192', variates=NULL, invjacobian=FALSE){

    if(is.character(Qfunction)){
        Qfunction <- readRDS(paste0(Qfunction,'.rds'))
    }

    x <- as.data.table(cbind(x))
    if(!is.null(variates)){colnames(x) <- variates}
    matrix(sapply(colnames(x), function(v){
        ##
        datum <- unlist(x[,v,with=F])
info <- as.list(auxmetadata[name == v])
        ##
        if(invjacobian){
#### Calculation of reciprocal Jacobian factors
            if(info$transform == 'log'){
                datum <- (datum - info$domainmin) * info$tscale
            }else if(info$transform == 'logminus'){
                datum <- (info$domainmax - datum) * info$tscale
            }else if(info$transform == 'probit'){
                datum <- qnorm((datum-info$domainmin)/(info$domainmax-info$domainmin))
                datum <- dnorm(datum) * info$tscale * (info$domainmax-info$domainmin) 
            }else{
                datum <- rep(info$tscale, length(datum))
            }
            xv <- data.matrix(x[,..v])
            if(info$mcmctype %in% c('C','D')){
                datum[(xv >= info$censormax) | (xv <= info$censormin)] <- 1L
            }
            datum[is.na(xv)] <- 1L
        }else{
#### Transformation to internal value for MCMC
            if(info$mcmctype == 'R'){# continuous
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
                datum <- round((datum-info$tlocation)/info$tscale) # output is in range 1 to Nvalues
                if(Oout == 'init'){ # in sampling functions or init MCMC
                    datum[is.na(datum)] <- info$Nvalues/2+0.5
                    datum <- Qfunction((datum-0.5)/info$Nvalues)
                    datum[datum==+Inf] <- 1e6
                    datum[datum==-Inf] <- -1e6
                } else if(Oout == 'left'){ # as left for MCMC
                    datum <- Qfunction(pmax(0,datum-1L)/info$Nvalues)
                    datum[is.na(datum)] <- -Inf
                } else if(Oout == 'right'){ # as right for MCMC
                    datum <- Qfunction(pmin(info$Nvalues,datum)/info$Nvalues)
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
                leftbound <- pmax(datum - info$step, info$domainmin, na.rm=T)
                leftbound[leftbound <= info$censormin] <- info$domainmin
                rightbound <- pmin(datum + info$step, info$domainmax, na.rm=T)
                rightbound[rightbound >= info$censormax] <- info$domainmax
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
                if(Dout == 'left'){
                    datum <- leftbound
                } else if(Dout == 'right'){
                    datum <- rightbound
                } else if(Dout == 'init'){ #init in MCMC
                    datum[is.na(datum)] <- 0L
                } else if(Dout == 'aux'){ # aux variable in MCMC
                    sel <- is.na(datum)
                    datum[sel] <- NA
                    datum[!sel] <- 1L
                } else if(Dout == 'index'){ #in sampling functions
                    datum[xv >= info$censormax] <- +Inf
                    datum[xv <= info$censormin] <- -Inf
                } else if(Dout == 'sleft'){ #in sampling functions
                    datum[!is.na(xv) & (xv <= info$censormin)] <- leftbound[!is.na(xv) & (xv <= info$censormin)]
                } else if(Dout == 'sright'){ #in sampling functions
                    datum[!is.na(xv) & (xv >= info$censormax)] <- rightbound[!is.na(xv) & (xv >= info$censormax)]
                }
                ##
            } else if(info$mcmctype == 'C'){ # censored
                leftbound <- rep(info$domainmin, length(datum))
                leftbound[datum >= info$censormax] <- info$censormax
                rightbound <- rep(info$domainmax, length(datum))
                rightbound[datum <= info$censormin] <- info$censormin
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
                    datum <- leftbound
                } else if(Cout == 'right'){ # in MCMC
                    datum <- rightbound
                } else if(Cout == 'lat'){ # latent variable in MCMC
                    sel <- is.na(datum) | (xv >= info$censormax) | (xv <= info$censormin)
                    datum[sel] <- NA
                } else if(Cout == 'init'){ #init in MCMC
                    datum[is.na(xv)] <- 0L
                    datum[!is.na(xv) & (xv <= info$censormin)] <- rightbound[!is.na(xv) & (xv <= info$censormin)] - 0.125
                    datum[!is.na(xv) & (xv >= info$censormax)] <- leftbound[!is.na(xv) & (xv >= info$censormax)] + 0.125
                    datum[!is.na(xv) & (xv < info$censormax) & (xv > info$censormin)] <- NA
                } else if(Cout == 'aux'){ # aux variable in MCMC
                    sel <- is.na(datum)
                    datum[sel] <- NA
                    datum[!sel] <- 1L
                } else if(Cout == 'index'){ #in sampling functions
                    datum[xv >= info$censormax] <- +Inf
                    datum[xv <= info$censormin] <- -Inf
                } else if(Cout == 'sleft'){ #in sampling functions
                    datum[!is.na(xv) & (xv <= info$censormin)] <- leftbound[!is.na(xv) & (xv <= info$censormin)]
                } else if(Cout == 'sright'){ #in sampling functions
                    datum[!is.na(xv) & (xv >= info$censormax)] <- rightbound[!is.na(xv) & (xv >= info$censormax)]
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
        }
        ##
        datum
    }),ncol=ncol(x),dimnames=list(NULL,colnames(x)))
}

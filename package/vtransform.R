## Transformation from variate to internal variable
vtransform <- function(x, auxmetadata, Cout='init', Dout='data', Oout='data', Bout='numeric', Nout='numeric', variates=NULL, invjacobian=FALSE, useOquantiles=TRUE){

        Qf <- readRDS('Qfunction8192.rds')
        DQf <- readRDS('DQfunction2048.rds')

    x <- as.data.table(cbind(x))
    if(!is.null(variates)){colnames(x) <- variates}
    matrix(sapply(colnames(x), function(v){
        ##
        datum <- unlist(x[,v,with=F])
info <- as.list(auxmetadata[name == v])
        ##
        if(invjacobian){
#### Calculation of reciprocal Jacobian factors
            if(info$mcmctype %in% c('B','N','O')){
                datum <- rep(1L, length(datum))
            }else{
                if(info$transform == 'log'){
                    datum <- (datum - info$domainmin) * info$tscale
                }else if(info$transform == 'logminus'){
                    datum <- (info$domainmax - datum) * info$tscale
                }else if(info$transform == 'Q'){
                    datum <- Qf((datum-info$domainmin)/(info$domainmax-info$domainmin))
                    datum <- DQf(datum) * info$tscale * (info$domainmax-info$domainmin) 
                }else{
                    datum <- rep(info$tscale, length(datum))
                }
                xv <- data.matrix(x[,..v])
                if(info$mcmctype %in% c('C','D')){
                    datum[(xv >= info$censormax) | (xv <= info$censormin)] <- 1L
                }
                datum[is.na(xv)] <- 1L
            }
        }else{
            
#### Transformation to internal value for MCMC
            if(info$mcmctype == 'R'){# continuous
                if (info$transform == 'log'){
                    datum <- log(datum-info$domainmin)
                }else if (info$transform == 'logminus'){
                    datum <- log(info$domainmax-datum)
                }else if (info$transform == 'Q'){
                    datum <- Qf((datum-info$domainmin)/(info$domainmax-info$domainmin))
                }
                datum <- (datum-info$tlocation)/info$tscale
                ##
            } else if(info$mcmctype == 'O'){ # ordinal
                olocation <- (info$Nvalues*info$domainmin - info$domainmax)/(info$Nvalues - 1)
                oscale <- (info$domainmax - info$domainmin)/(info$Nvalues - 1)
                ##
                datum <- round((datum-olocation)/oscale) # output is in range 1 to Nvalues
                if(Oout == 'init'){ # in sampling functions or init MCMC
                    datum[is.na(datum)] <- info$Nvalues/2+0.5
                    datum <- Qf((datum-0.5)/info$Nvalues)
                    datum[datum==+Inf] <- 1e6
                    datum[datum==-Inf] <- -1e6
                    if(useOquantiles){
                        datum <- (datum-info$tlocation)/info$tscale
                    }
                } else if(Oout == 'left'){ # as left for MCMC
                    datum <- Qf(pmax(0,datum-1L)/info$Nvalues)
                    if(useOquantiles){
                        datum <- (datum-info$tlocation)/info$tscale
                    }
                    datum[is.na(datum)] <- -Inf
                } else if(Oout == 'right'){ # as right for MCMC
                    datum <- Qf(pmin(info$Nvalues,datum)/info$Nvalues)
                    if(useOquantiles){
                        datum <- (datum-info$tlocation)/info$tscale
                    }
                    datum[is.na(datum)] <- +Inf
                } else if(Oout == 'aux'){ # aux variable in MCMC
                    sel <- is.na(datum)
                    datum[sel] <- NA
                    datum[!sel] <- 1L
                } else if(Oout == 'boundisinf'){ # in output functions
                    datum <- datum-1L
                } else if(Oout == 'original'){ # in output functions
                    datum <- datum * oscale + olocation
                }
                ##
            } else if(info$mcmctype == 'D'){ # discretized
                xv <- data.matrix(x[,..v])
                censmin <- info$censormin
                censmax <- info$censormax
                leftbound <- pmax(datum - info$step, info$domainmin, na.rm=T)
                leftbound[leftbound <= censmin] <- info$domainmin
                leftbound[leftbound >= censmax] <- censmax
                rightbound <- pmin(datum + info$step, info$domainmax, na.rm=T)
                rightbound[rightbound >= censmax] <- info$domainmax
                rightbound[rightbound <= censmin] <- censmin
                if (info$transform == 'log'){
                    datum <- log(datum-info$domainmin)
                    leftbound <- log(leftbound-info$domainmin)
                    rightbound <- log(rightbound-info$domainmin)
                    censmin <- log(censmin-info$domainmin)
                    censmax <- log(censmax-info$domainmin)
                } else if (info$transform == 'logminus'){
                    datum <- log(info$domainmax-datum)
                    leftbound <- log(info$domainmax-leftbound)
                    rightbound <- log(info$domainmax-rightbound)
                    censmin <- log(info$domainmax-censmin)
                    censmax <- log(info$domainmax-censmax)
                } else if (info$transform == 'Q'){
                    datum <- Qf((datum-info$domainmin)/(info$domainmax-info$domainmin))
                    leftbound <- Qf((leftbound-info$domainmin)/(info$domainmax-info$domainmin))
                    rightbound <- Qf((rightbound-info$domainmin)/(info$domainmax-info$domainmin))
                    censmin <- Qf((censmin-info$domainmin)/(info$domainmax-info$domainmin))
                    censmax <- Qf((censmax-info$domainmin)/(info$domainmax-info$domainmin))
                }
                ## datum <- (datum-info$tlocation)/info$tscale
                ## rightbound <- (rightbound-info$tlocation)/info$tscale
                ## leftbound <- (leftbound-info$tlocation)/info$tscale
                ## censmax <- (censmax-info$tlocation)/info$tscale
                ## censmin <- (censmin-info$tlocation)/info$tscale
                if(Dout == 'left'){
                    datum <- leftbound
                } else if(Dout == 'right'){
                    datum <- rightbound
                } else if(Dout == 'init'){ #init in MCMC
                    datum[is.na(xv)] <- 0L
                    datum[!is.na(xv) & (xv <= info$censormin)] <- censmin - 0.125*info$tscale
                    datum[!is.na(xv) & (xv >= info$censormax)] <- censmax + 0.125*info$tscale
                } else if(Dout == 'aux'){ # aux variable in MCMC
                    sel <- is.na(xv)
                    datum[sel] <- NA
                    datum[!sel] <- 1L
                } else if(Dout == 'boundisinf'){ #in sampling functions
                    datum[xv >= info$censormax] <- +Inf
                    datum[xv <= info$censormin] <- -Inf
                } else if(Dout == 'sleft'){ #in sampling functions
                    datum <- rep(censmin, length(datum))
                } else if(Dout == 'sright'){ #in sampling functions
                    datum <- rep(censmax, length(datum))
                }
                if(Dout != 'aux'){
                    datum <- (datum-info$tlocation)/info$tscale
                    }
                ##
            } else if(info$mcmctype == 'C'){ # censored
                xv <- data.matrix(x[,..v])
                censmin <- info$censormin
                censmax <- info$censormax
                if (info$transform == 'log'){
                    datum <- log(datum-info$domainmin)
                    censmin <- log(censmin-info$domainmin)
                    censmax <- log(censmax-info$domainmin)
                } else if (info$transform == 'logminus'){
                    datum <- log(info$domainmax-datum)
                    censmin <- log(info$domainmax-censmin)
                    censmax <- log(info$domainmax-censmax)
                } else if (info$transform == 'Q'){
                    datum <- Qf((datum-info$domainmin)/(info$domainmax-info$domainmin))
                    censmin <- Qf((censmin-info$domainmin)/(info$domainmax-info$domainmin))
                    censmax <- Qf((censmax-info$domainmin)/(info$domainmax-info$domainmin))
                }
                ## datum <- (datum-info$tlocation)/info$tscale
                ## censmax <- (censmax-info$tlocation)/info$tscale
                ## censmin <- (censmin-info$tlocation)/info$tscale
                if(Cout == 'left'){ # in MCMC
                    sel <- is.na(xv) | (xv < info$censormax)
                    datum[sel] <- -Inf
                    datum[!sel] <- censmax
                } else if(Cout == 'right'){ # in MCMC
                    sel <- is.na(xv) | (xv > info$censormin)
                    datum[sel] <- +Inf
                    datum[!sel] <- censmin
                } else if(Cout == 'lat'){ # latent variable in MCMC
                    sel <- is.na(xv) | (xv >= info$censormax) | (xv <= info$censormin)
                    datum[sel] <- NA
                } else if(Cout == 'init'){ #init in MCMC
                    datum[is.na(xv)] <- 0L
                    datum[!is.na(xv) & (xv <= info$censormin)] <- censmin - 0.125*info$tscale
                    datum[!is.na(xv) & (xv >= info$censormax)] <- censmax + 0.125*info$tscale
                    datum[!is.na(xv) & (xv < info$censormax) & (xv > info$censormin)] <- NA
                } else if(Cout == 'aux'){ # aux variable in MCMC
                    sel <- is.na(xv)
                    datum[sel] <- NA
                    datum[!sel] <- 1L
                } else if(Cout == 'boundisinf'){ #in sampling functions
                    datum[xv >= info$censormax] <- +Inf
                    datum[xv <= info$censormin] <- -Inf
                } else if(Cout == 'sleft'){ #in sampling functions
                    datum <- rep(censmin, length(datum))
                } else if(Cout == 'sright'){ #in sampling functions
                    datum <- rep(censmax, length(datum))
                }
                if(Cout != 'aux'){
                    datum <- (datum-info$tlocation)/info$tscale
                }
                ##                
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

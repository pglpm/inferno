buildauxmetadata <- function(data, metadata, file=TRUE){
    sdoveriqr <- 0.5/qnorm(0.75)

    datafile <- NULL
    if(is.character(data) && file.exists(data)){
        datafile <- data
        data <- fread(datafile, na.strings='')
    }
    data <- as.data.table(data)
    if(is.character(metadata) && file.exists(metadata)){
        metadata <- fread(metadata, na.strings='')
    }
    metadata <- as.data.table(metadata)
    ## consistency checks
    if(!identical(metadata$name, colnames(data))){
        stop('ERROR: mismatch in variate names or order')
    }
    ##
    ##Q <- readRDS('Qfunction512.rds')
    ##
    idR <- idC <- idD <- idO <- idB <- idN <- 1L
    auxmetadata <- data.table()
    for(xn in colnames(data)){
        x <- data[[xn]]
        x <- x[!is.na(x)]
        xinfo <- as.list(metadata[name == xn])
        xinfo$type <- tolower(xinfo$type)
        ordinal <- NA
        cens <- any(is.finite(c(xinfo$censormin, xinfo$censormax)))
        rounded <- NA
        transf <- 'identity' # temporary
        vval <- xinfo[grep('^V[0-9]+$', names(xinfo))]
        ## print(xn)
        ## str(vval)
        Q1 <- NA
        Q2 <- NA
        Q3 <- NA
        if(xinfo$type == 'binary'){# seems binary variate
            if(length(unique(x)) != 2){
                cat('Warning: inconsistencies with variate ', xn, '\n')
            }
            vtype <- 'B'
            vid <- idB
            idB <- idB+1L
            vn <- xinfo$Nvalues
            vd <- xinfo$rounding/2
            domainmin <- 0
            domainmax <- 1
            censormin <- -Inf
            censormax <- +Inf
            location <- 0
            scale <- 1
            plotmin <- NA
            plotmax <- NA
            mctest1 <- vval[1]
            mctest2 <- vval[round(vn/2)]
            mctest3 <- vval[vn]
        }else if(xinfo$type == 'nominal'){# nominal variate
            vtype <- 'N'
            vid <- idN
            idN <- idN+1L
            vn <- xinfo$Nvalues
            vd <- 0.5
            domainmin <- 1 # Nimble index categorical from 1
            domainmax <- vn
            censormin <- -Inf
            censormax <- +Inf
            location <- 0
            scale <- 1
            plotmin <- NA
            plotmax <- NA
            mctest1 <- vval[1]
            mctest2 <- vval[round(vn/2)]
            mctest3 <- vval[vn]
        }else if(xinfo$type == 'ordinal'){# ordinal variate
            vtype <- 'O'
            vid <- idO
            idO <- idO+1L
            transf <- 'identity'
            ordinal <- TRUE
            vn <- xinfo$Nvalues
            vd <- 0.5
            domainmin <- xinfo$domainmin
            domainmax <- xinfo$domainmax
            censormin <- -Inf
            censormax <- +Inf
            ##vval <- as.vector(xinfo[paste0('V',1:vn)], mode='character')
            location <- (vn*domainmin - domainmax)/(vn - 1)
            scale <- (domainmax - domainmin)/(vn - 1)
            plotmin <- xinfo$plotmin
            plotmax <- xinfo$plotmax
            Q1 <- mctest1 <- quantile(x, probs=0.25, type=6)
            Q2 <- mctest2 <- quantile(x, probs=0.5, type=6)
            Q3 <- mctest3 <- quantile(x, probs=0.75, type=6)
        }else if(xinfo$type == 'continuous'){# continuous variate (R,C,D)
            vn <- +Inf
            vd <- xinfo$rounding/2
            rounded <- (vd > 0)
            domainmin <- xinfo$domainmin
            domainmax <- xinfo$domainmax
            censormin <- max(domainmin, xinfo$censormin, na.rm=T)
            censormax <- min(domainmax, xinfo$censormax, na.rm=T)
            ## cens <- (censormin > domainmin) || (censormax < domainmax)
            location <- xinfo$centralvalue
            scale <- abs(xinfo$highvalue - xinfo$lowvalue)
            Q1 <- mctest1 <- quantile(x, probs=0.25, type=6)
            Q2 <- mctest2 <- quantile(x, probs=0.5, type=6)
            Q3 <- mctest3 <- quantile(x, probs=0.75, type=6)
            plotmin <- xinfo$plotmin
            plotmax <- xinfo$plotmax
            if(is.finite(xinfo$domainmin) && is.finite(xinfo$domainmax)){ # needs transformation
                Qf <- readRDS('Qfunction8192.rds')
                transf <- 'Q'
                location <- Qf((location-domainmin)/(domainmax-domainmin))
                scale <- abs(Qf((xinfo$highvalue-domainmin)/(domainmax-domainmin)) - Qf((xinfo$lowvalue-domainmin)/(domainmax-domainmin)))*sdoveriqr
            }else if(is.finite(xinfo$domainmin)){
                transf <- 'log'
                location <- log(location-domainmin)
                scale <- abs(log(xinfo$highvalue-domainmin) - log(xinfo$lowvalue-domainmin))*sdoveriqr
            }else if(is.finite(xinfo$domainmax)){
                transf <- 'logminus'
                location <- log(domainmax-location)
                scale <- abs(log(domainmax-xinfo$highvalue) - log(domainmax-xinfo$lowvalue))*sdoveriqr
            }
            if(xinfo$rounding > 0){# discretized
                ## if(diff(range(x))/xinfo$rounding > 256){
                ##     Message('\nVariate ',xn,' is reported as "rounded".\nConsider the possibility of treating it as continuous setting "rounding" to 0.\n')
                ## }
                vtype <- 'D'
                vid <- idD
                idD <- idD+1L
            }
            else if(cens){# censored
                vtype <- 'C'
                vid <- idC
                idC <- idC+1L
            }else{ # continuous
                vtype <- 'R'
                vid <- idR
                idR <- idR+1L
            }
        }else{
            stop(paste0('ERROR: unknown variate type for ', xn))
        }
        ##
        ## print(auxmetadata[nrow(auxmetadata)])
        ## print(                         as.data.table(c(list(name=xn, type=vtype, transform=transf, Nvalues=vn, step=vd, domainmin=domainmin, domainmax=domainmax, censormin=censormin, censormax=censormax, tlocation=location, tscale=scale, plotmin=plotmin, plotmax=plotmax, Q1=Q1, Q2=Q2, Q3=Q3),
        ##                    vval
        ##                    )))
        auxmetadata <- rbind(auxmetadata,
                             c(list(name=xn, mcmctype=vtype, id=vid, # censored=cens,
                                    rounded=rounded, transform=transf, Nvalues=vn, step=vd, domainmin=domainmin, domainmax=domainmax, censormin=censormin, censormax=censormax, tlocation=location, tscale=scale, plotmin=plotmin, plotmax=plotmax, Q1=Q1, Q2=Q2, Q3=Q3, mctest1=mctest1, mctest2=mctest2, mctest3=mctest3),
                           vval
                           ), fill=FALSE)
    }

    if(!missing(file) && file!=FALSE){# must save to file
        if(is.character(file)){
            file <- paste0(sub('.rds$', '', file), '.rds')
        }else{
            file <- paste0('auxmetadata_', datafile)
            file <- paste0(sub('.csv$', '', file), '.rds')
        }
        if(file.exists(file)){
            file.rename(from=file, to=paste0(sub('.rds$', '', file), '_bak',format(Sys.time(), '%y%m%dT%H%M%S'),'.rds'))
        }
        saveRDS(auxmetadata, file)
        cat(paste0('Saved proposal aux-metadata file as ', file, '\n'))
    }else{
        auxmetadata
    }
}

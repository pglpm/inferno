buildvarinfoaux <- function(data, varinfo, file=TRUE){
    require('data.table')
    
    if(is.character(data) && file.exists(data)){data <- fread(data, na.strings='')}
    data <- as.data.table(data)
    if(is.character(varinfo) && file.exists(varinfo)){
        varinfoname <- varinfo
        varinfo <- fread(varinfo)
    }
    varinfo <- as.data.table(varinfo)
    ## consistency checks
    if(!identical(varinfo$name, colnames(data))){
        stop('ERROR: mismatch in variate names or order')
    }
    ##
    ##Q <- readRDS('Qfunction512.rds')
    ##
    idR <- idC <- idD <- idO <- idB <- idN <- 1L
    varinfoaux <- data.table()
    for(xn in colnames(data)){
        x <- data[[xn]]
        x <- x[!is.na(x)]
        xinfo <- as.list(varinfo[name == xn])
        xinfo$type <- tolower(xinfo$type)
        ordinal <- NA
        cens <- FALSE
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
            vmin <- 0
            vmax <- 1
            tmin <- -Inf
            tmax <- +Inf
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
            vmin <- 1 # Nimble index categorical from 1
            vmax <- vn
            tmin <- -Inf
            tmax <- +Inf
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
            vmin <- xinfo$domainmin
            vmax <- xinfo$domainmax
            tmin <- -Inf
            tmax <- +Inf
            ##vval <- as.vector(xinfo[paste0('V',1:vn)], mode='character')
            location <- (vn*vmin - vmax)/(vn - 1)
            scale <- (vmax - vmin)/(vn - 1)
            plotmin <- xinfo$plotmin
            plotmax <- xinfo$plotmax
            Q1 <- mctest1 <- quantile(x, probs=0.25, type=6)
            Q2 <- mctest2 <- quantile(x, probs=0.5, type=6)
            Q3 <- mctest3 <- quantile(x, probs=0.75, type=6)
        }else if(xinfo$type == 'continuous'){# continuous variate (R,C,D)
            vn <- +Inf
            vd <- xinfo$rounding/2
            rounded <- (vd > 0)
            vmin <- xinfo$domainmin
            vmax <- xinfo$domainmax
            tmin <- xinfo$censormin # max(xinfo$censormin, vmin, na.rm=TRUE)
            tmax <- xinfo$censormax # min(xinfo$censormax, vmax, na.rm=TRUE)
            cens <- (tmin > vmin) || (tmax < vmax)
            location <- xinfo$centralvalue
            scale <- abs(xinfo$highvalue - xinfo$lowvalue)
            Q1 <- mctest1 <- quantile(x, probs=0.25, type=6)
            Q2 <- mctest2 <- quantile(x, probs=0.5, type=6)
            Q3 <- mctest3 <- quantile(x, probs=0.75, type=6)
            plotmin <- xinfo$plotmin
            plotmax <- xinfo$plotmax
            if(is.finite(xinfo$domainmin) && is.finite(xinfo$domainmax)){ # needs transformation
                transf <- 'probit'
                location <- qnorm((location-vmin)/(vmax-vmin))
                scale <- abs(qnorm((xinfo$highvalue-vmin)/(vmax-vmin)) - qnorm((xinfo$lowvalue-vmin)/(vmax-vmin)))
            }else if(is.finite(xinfo$domainmin)){
                transf <- 'log'
                location <- log(location-vmin)
                scale <- abs(log(xinfo$highvalue-vmin) - log(xinfo$lowvalue-vmin))
            }else if(is.finite(xinfo$domainmax)){
                transf <- 'logminus'
                location <- log(vmax-location)
                scale <- abs(log(vmax-xinfo$highvalue) - log(vmax-xinfo$lowvalue))
            }
            if(xinfo$rounding > 0){# discretized
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
        ## print(varinfoaux[nrow(varinfoaux)])
        ## print(                         as.data.table(c(list(name=xn, type=vtype, transform=transf, Nvalues=vn, step=vd, domainmin=vmin, domainmax=vmax, censormin=tmin, censormax=tmax, tlocation=location, tscale=scale, plotmin=plotmin, plotmax=plotmax, Q1=Q1, Q2=Q2, Q3=Q3),
        ##                    vval
        ##                    )))
        varinfoaux <- rbind(varinfoaux,
                         c(list(name=xn, mcmctype=vtype, id=vid, censored=cens, rounded=rounded, transform=transf, Nvalues=vn, step=vd, domainmin=vmin, domainmax=vmax, censormin=tmin, censormax=tmax, tlocation=location, tscale=scale, plotmin=plotmin, plotmax=plotmax, Q1=Q1, Q2=Q2, Q3=Q3, mctest1=mctest1, mctest2=mctest2, mctest3=mctest3),
                           vval
                           ), fill=FALSE)
    }
    if(is.character(file) || (is.logical(file) && file)){ # must save to file
        if(is.character(file)){
            file <- paste0(sub('.rds$', '', file), '.rds')
        }else{
            if(file.exists('varinfoaux.rds')){
                file <- paste0('varinfoaux_',format(Sys.time(), '%y%m%dT%H%M%S'),'.rds')
            }else{
                file <- 'varinfoaux.rds'
            }
        }
        saveRDS(varinfoaux, file)
        cat(paste0('Saved auxiliary variate-info file to ', file, '\n'))
    }else{
        varinfoaux
    }
}

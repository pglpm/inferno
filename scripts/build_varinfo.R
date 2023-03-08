library('data.table')
library('png')
library('foreach')
library('doFuture')
library('doRNG')
registerDoFuture()
cat('\navailableCores: ')
cat(availableCores())
cat('\navailableCores-multicore: ')
cat(availableCores('multicore'))
if(Sys.info()['nodename']=='luca-HP-Z2-G9'){
    ncores <- 20}else{
    ncores <- 6}
cat(paste0('\nusing ',ncores,' cores\n'))
if(ncores>1){
    if(.Platform$OS.type=='unix'){
        plan(multicore, workers=ncores)
    }else{
        plan(multisession, workers=ncores)
    }
}else{
    plan(sequential)
}
    
buildvarinfo <- function(data, file=NULL){
    gcd2 <- function(a, b){ if (b == 0) a else Recall(b, a %% b) }
    gcd <- function(...) Reduce(gcd2, c(...))
    ##
    if(is.character(data) && file.exists(data)){data <- fread(data)}
    data <- as.data.table(data)
    varinfo <- data.table()
    for(x in data){
        x <- x[!is.na(x)]
        transf <- 'identity' # temporary
        Q1 <- NA
        Q2 <- NA
        Q3 <- NA
        if(length(unique(x)) == 2){# seems binary variate
            vtype <- 'B'
            vn <- 2
            vd <- NA
            vmin <- NA
            vmax <- NA
            tmin <- NA
            tmax <- NA
            vval <- as.character(unique(x))
            names(vval) <- paste0('V',1:2)
            location <- NA
            scale <- NA
            plotmin <- NA
            plotmax <- NA
        }else if(!is.numeric(x)){# nominal variate
            vtype <- 'N'
            vn <- length(unique(x))
            vd <- NA
            vmin <- NA # Nimble index categorical from 1
            vmax <- NA
            tmin <- NA
            tmax <- NA
            vval <- as.character(unique(x))
            names(vval) <- paste0('V',1:vn)
            location <- NA
            scale <- NA
            plotmin <- NA
            plotmax <- NA
        }else{# discrete, continuous, or boundary-singular variate
            ud <- unique(signif(diff(sort(unique(x))),3)) # differences
            multi <- 10^(-min(floor(log10(ud))))
            dd <- round(gcd(ud*multi))/multi # greatest common difference
            ##
            Q1 <- quantile(x, probs=0.25, type=6)
            Q2 <- quantile(x, probs=0.5, type=6)
            Q3 <- quantile(x, probs=0.75, type=6)
            if(dd == 0){ # consider it as continuous
                ## temporary values
                vtype <- 'R'
                vn <- Inf
                vd <- 0
                vmin <- -Inf
                vmax <- +Inf
                tmin <- NA
                tmax <- NA
                location <- Q2
                scale <- (Q3-Q1)/2
                plotmin <- min(x) - (Q3-Q1)/2
                plotmax <- max(x) + (Q3-Q1)/2
                ##
                ix <- x[!(x %in% range(x))] # exclude boundary values
                repindex <- mean(table(ix)) # average of repeated inner values
                ## contindex <- length(unique(diff(sort(unique(ix)))))/length(ix) # check for repeated values
                if(sum(x == min(x)) > repindex){ # seems to be left-singular
                    tmin <- min(x)
                    plotmin <- tmin
                }
                if(sum(x == max(x)) > repindex){ # seems to be right-singular
                    tmax <- max(x)
                    plotmax <- tmax
                }
                if(all(x > 0)){ # seems to be strictly positive
                    transf <- 'log'
                    vmin <- 0
                    ## location <- log(Q2)
                    ## scale <- (log(Q3) - log(Q1))/2
                    plotmin <- max(vmin, plotmin)
                }
            }else{# ordinal
                vtype <- 'O'
                if(dd >= 1){ # seems originally integer
                    transf <- 'Q'
                    vmin <- min(1, x)
                    vmax <- max(x)
                    vn <- vmax - vmin + 1
                    vd <- 1
                    tmin <- NA
                    tmax <- NA
                    location <- NA # (vn*vmin-vmax)/(vn-1)
                    scale <- NA # (vmax-vmin)/(vn-1)
                    plotmin <- max(vmin, min(x) - IQR(x)/2)
                    plotmax <- max(x)
                }else{ # seems a rounded continuous variate
                    vtype <- 'R'
                    vn <- Inf
                    vd <- dd
                    vmin <- -Inf
                    vmax <- +Inf
                    tmin <- NA
                    tmax <- NA
                    location <- Q2
                    scale <- (Q3-Q1)/2
                    plotmin <- min(x) - (Q3-Q1)/2
                    plotmax <- max(x) + (Q3-Q1)/2
                    if(all(x > 0)){ # seems to be strictly positive
                        transf <- 'log'
                        vmin <- 0
                        ## location <- log(Q2)
                        ## scale <- (log(Q3) - log(Q1))/2
                        plotmin <- max(vmin, plotmin)
                    }
                }# end rounded
            }# end integer
            vval <- NULL
        }# end numeric
        ##
        varinfo <- rbind(varinfo,
                         c(list(type=vtype, Nvalues=vn, step=vd, domainmin=vmin, domainmax=vmax, truncmin=tmin, truncmax=tmax, location=location, scale=scale, plotmin=plotmin, plotmax=plotmax),
                           as.list(vval)
                         ), fill=TRUE)
    }
    varinfo <- cbind(name=names(data), varinfo)
    if(!is.null(file)){
        cat(paste0('Saved proposal variate-info file to ',file,'\n'))
        fwrite(varinfo, file)
    }else{
        varinfo
    }
}

testf <- function(a){
    if(a==1){
        stop('abort')
    }else{
        print('continuing')
    }
    34
}


buildvarinfoaux <- function(data, varinfo){
    if(is.character(data) && file.exists(data)){data <- fread(data)}
    data <- as.data.table(data)
    if(is.character(varinfo) && file.exists(varinfo)){varinfo <- fread(varinfo)}
    varinfo <- as.data.table(varinfo)
    ## consistency checks
    if(!identical(varinfo$name, colnames(data))){
        stop('ERROR: mismatch in variate names or order')
    }
    ##
    Q <- readRDS('Qfunction512.rds')
    ##
    varinfoaux <- data.table()
    for(xn in colnames(data)){
        x <- data[[xn]]
        x <- x[!is.na(x)]
        xinfo <- as.list(varinfo[name == xn])
        transf <- 'identity' # temporary
        Q1 <- NA
        Q2 <- NA
        Q3 <- NA
        if(varinfo$type == 'B'){# seems binary variate
            if(length(unique(x)) != 2){
                cat(paste0('Warning: inconsistencies with variate ',xn,'\n'))
            }
            vtype <- 'B'
            vn <- xinfo$Nvalues
            vd <- xinfo$step/2
            vmin <- 0
            vmax <- 1
            tmin <- -Inf
            tmax <- +Inf
            vval <- as.vector(xinfo[paste0('V',1:vn)], mode='character')
            location <- 0
            scale <- 1
            plotmin <- 0
            plotmax <- 1
        }else if(varinfo$type == 'N'){# nominal variate
            if(!is.numeric(x)){
                cat(paste0('Warning: inconsistencies with variate ',xn,'\n'))
            }
            vtype <- 'C'
            vn <- xinfo$Nvalues
            vd <- 0.5
            vmin <- 1 # Nimble index categorical from 1
            vmax <- vn
            tmin <- -Inf
            tmax <- +Inf
            vval <- as.vector(xinfo[paste0('V',1:vn)], mode='character')
            location <- 0
            scale <- 1
            plotmin <- 1
            plotmax <- vn
        }else{# discrete, continuous, or boundary-singular variate
            ## ud <- unique(signif(diff(sort(unique(x))),3)) # differences
            ## multi <- 10^(-min(floor(log10(ud))))
            dd <- xinfo$step/2 # greatest common difference
            ##
            Q1 <- quantile(x, probs=0.25, type=6)
            Q2 <- quantile(x, probs=0.5, type=6)
            Q3 <- quantile(x, probs=0.75, type=6)
            if(dd == 0){ # consider it as continuous
                ## temporary values
                vtype <- 'R'
                vn <- Inf
                vd <- 0
                vmin <- -Inf
                vmax <- +Inf
                tmin <- xinfo$truncmin
                tmax <- xinfo$truncmax
                location <- xinfo$location
                scale <- xinfo$scale
                plotmin <- xinfo$plotmin
                plotmax <- xinfo$plotmax
                ##
                ix <- x[!(x %in% range(x))] # exclude boundary values
                repindex <- mean(table(ix)) # average of repeated inner values
                ## contindex <- length(unique(diff(sort(unique(ix)))))/length(ix) # check for repeated values
                if(is.finite(xinfo$domainmin) && !is.finite(xinfo$domainmax)){ # seems to be strictly positive
                    transf <- 'log'
                    vmin <- xinfo$domainmin
                    location <- log(xinfo$location - xinfo$domainmin)
                    ## scale <- (log(xinfo$location + xinfo$scale - xinfo$domainmin) - log(xinfo$location - xinfo$scale + xinfo$domainmin))/2
                    scale <- xinfo$scale/(xinfo$location - xinfo$domainmin)
                    tmin <- log(tmin)
                    tmax <- log(tmax)
                }
                if(is.finite(xinfo$domainmin) && is.finite(xinfo$domainmax)){ # seems to be doubly bounded
                    transf <- 'Q'
                    vmin <- xinfo$domainmin
                    vmax <- xinfo$domainmax
                    location <- Q((xinfo$location - xinfo$domainmin)/(xinfo$domainmax - xinfo$domainmin))
                    scale <- (Q((xinfo$location + xinfo$scale - xinfo$domainmin)/(xinfo$domainmax - xinfo$domainmin)) - Q((xinfo$location - xinfo$scale - xinfo$domainmin)/(xinfo$domainmax - xinfo$domainmin)))/2
                    tmin <- Q((tmin - xinfo$domainmin)/(xinfo$domainmax - xinfo$domainmin))
                    tmax <- Q((tmax - xinfo$domainmin)/(xinfo$domainmax - xinfo$domainmin))
                }
                if(is.finite(tmin) || is.finite(tmax)){ # seems to be singular at boundary 
                    vtype <- 'D'
                }

            }else if(xinfo$type == 'O'){# ordinal
                vtype <- 'I'
                v
                
                if(dd >= 1){ # seems originally integer
                    transf <- 'Q'
                    vmin <- min(1, x)
                    vmax <- max(x)
                    vn <- vmax - vmin + 1
                    vd <- 1
                    tmin <- NA
                    tmax <- NA
                    location <- NA # (vn*vmin-vmax)/(vn-1)
                    scale <- NA # (vmax-vmin)/(vn-1)
                    plotmin <- max(vmin, min(x) - IQR(x)/2)
                    plotmax <- max(x)
                }else{ # seems a rounded continuous variate
                    vtype <- 'R'
                    vn <- Inf
                    vd <- dd
                    vmin <- -Inf
                    vmax <- +Inf
                    tmin <- NA
                    tmax <- NA
                    location <- Q2
                    scale <- (Q3-Q1)/2
                    plotmin <- min(x) - (Q3-Q1)/2
                    plotmax <- max(x) + (Q3-Q1)/2
                    if(all(x > 0)){ # seems to be strictly positive
                        transf <- 'log'
                        vmin <- 0
                        ## location <- log(Q2)
                        ## scale <- (log(Q3) - log(Q1))/2
                        plotmin <- max(vmin, plotmin)
                    }
                }# end rounded
            }# end integer
            vval <- NULL
        }# end numeric
        ##
        varinfo <- rbind(varinfo,
                         c(list(type=vtype, Nvalues=vn, step=vd, domainmin=vmin, domainmax=vmax, truncmin=tmin, truncmax=tmax, location=location, scale=scale, plotmin=plotmin, plotmax=plotmax),
                           as.list(vval)
                         ), fill=TRUE)
    }
    varinfo <- cbind(name=names(data), varinfo)
    if(!is.null(file)){
        cat(paste0('Saved proposal variate-info file to ',file,'\n'))
        fwrite(varinfo, file)
    }else{
        varinfo
    }
}


## gcd <- function(vect){Reduce(function(x,y) ifelse(y, Recall(y, x %% y), x), as.list(vect))}
## gcdm <- function(...){Reduce(function(x,y) ifelse(y, Recall(y, x %% y), x), list(...))}


gcd2 <- function(a, b) {
  if (b == 0) a else Recall(b, a %% b)
}
gcd <- function(...) Reduce(gcd2, c(...))

dt <- fread('ingrid_data_nogds6.csv')
data(iris)
iris <- as.data.table(iris)
iris2 <- iris
iris2$Species <- as.integer(iris2$Species)

dtx <- dt
dtx$extra <- rnorm(nrow(dtx))
dtx <- dtx[sample(1:nrow(dtx), min(10,nrow(dtx)))]
t(sapply(dtx, function(xx){
    xx <- xx[!is.na(xx)]
    testd <- unique(signif(diff(sort(unique(xx))),6))
    multi <- 10^(-min(floor(log10(testd))))
    round(gcd(testd*multi))/multi
}))





dtx <- iris
dtx$extra <- rnorm(nrow(dtx))
dtx <- dtx[sample(1:nrow(dtx), min(150,nrow(dtx)))]
t(sapply(dtx, function(xx){
    xx <- xx[!is.na(xx)]
    if(length(unique(xx)) > 2){ix <- xx[!(xx %in% range(xx))]}else{ix <- xx}
    dd0 <- diff(sort(xx))
    dd <- diff(sort(unique(xx)))
    q1 <- tquant(ix,0.25)
    q3 <- tquant(ix,0.75)
    c(
      left=sum(xx==min(xx)),
      right=sum(xx==max(xx)),
      diffratio=length(unique(dd))/length(dd),
      diffratio0=length(unique(dd0))/length(dd0),
      diffindex=length(unique(dd))/length(xx),
      diffindex0=length(unique(dd0))/length(xx),
      meanrep=mean(table(ix)),
      meanrepiqr=mean(table(ix[ix >= q1 & ix <=q3])),
      iqrrange=IQR(ix)/diff(range(xx)),
      min=min(xx),
      max=max(xx),
      int=is.integer(xx),
      unique=length(unique(xx)),
      uniqueratio=length(unique(xx))/length(xx),
      rg=diff(range(xx)),
      NULL
      )
}))
rm(dtx)

dtx <- dt
dtx$extra <- rnorm(nrow(dtx))
summary(t(sapply(1:100, function(xxx){set.seed(xxx)
    dtx <- dtx[sample(1:nrow(dtx), 10)]
    test <- t(sapply(dtx, function(xx){
        xx <- xx[!is.na(xx)]
        if(length(unique(xx)) > 2){ix <- xx[!(xx %in% range(xx))]}else{ix <- xx}
        dd <- diff(sort(unique(xx)))
        c(
            sum(xx==min(xx)),
            sum(xx==max(xx)),
            length(unique(dd))/length(dd),
            length(unique(dd))/length(xx),
            mean(table(ix))
        )
    }))
    c(sum(test['extra',5] < test[-c(13:14),5]),
      sum(test['extra',4] > test[-c(13:14),4]),
      sum(test['extra',3] > test[-c(13:14),3])
      )
})))

dtx <- dt
dtx$extra <- rnorm(nrow(dtx))
resu <- (t(sapply(1:100, function(xxx){set.seed(xxx)
    dtx <- dtx[sample(1:nrow(dtx), 150)]
    test <- t(sapply(dtx, function(xx){
        xx <- xx[!is.na(xx)]
        if(length(unique(xx)) > 2){ix <- xx[!(xx %in% range(xx))]}else{ix <- xx}
        dd <- diff(sort(unique(xx)))
        c(
            sum(xx==min(xx)),
            sum(xx==max(xx)),
            length(unique(dd))/length(dd),
            length(unique(dd))/length(xx),
            mean(table(ix))
        )
    }))
    c(
      min(test[13:14,4]),
      max(test[-c(13:14),4])
      )
})))
summary(resu)

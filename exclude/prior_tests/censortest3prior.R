library('data.table')
library('nimble')

varnames <- c('Apoe4_', 'Subgroup_num_',  'Gender_num_', 'GDTOTAL_gds',  'AVDEL30MIN_neuro', 'AVDELTOT_neuro',  'CATANIMSC_neuro', 'ANARTERR_neuro',  'RAVLT_immediate', 'TRAASCOR_neuro',  'TRABSCOR_neuro', 'AGE',  'LRHHC_n_long')
##
dt0 <- fread('~/repositories/ADBayes/scripts/bl.csv')
all(varnames %in% colnames(dt0))
dt <- dt0[,..varnames]
##
nas <- (1:nrow(dt))[apply(dt,1,function(x)any(is.na(x)))]
gds6 <- (1:nrow(dt))[dt[['GDTOTAL_gds']] >= 6]
## ##
## dt3 <- fread('~/repositories/ADBayes/scripts/3.13_pycharet_K50.csv')
## diffs <- sapply(varnames,function(x){
##     x1 <- sort(dt[[x]][-nas])
##     x2 <- sort(dt3[[x]])
##     abs(x1-x2)/(x1+x2+(x1+x2==0))/2
## })
## range(diffs)
## ## [1] 0.00000e+00 8.94359e-15
## ##
## dtred <- dt[-nas]
## ##
## fwrite(dt, 'ingrid_data_all.csv')
## fwrite(dt[-gds6], 'ingrid_data_nogds6.csv')
## fwrite(dt[-nas], 'ingrid_data_nonan.csv')
## fwrite(dt[-c(nas,gds6)], 'ingrid_data_nonangds6.csv')
## rm(dt0,dt3)

##
rtran <- function(x, location=0, scale=1){
    (log(x)-location)/scale
}
rtraninv <- function(y, location=0, scale=1){
    exp(y*scale+location)
}
##
itran0 <- function(x, n, min=0, max=n-1){
    round((x-min)/(max-min)*(n-1))
}
itraninv0 <- function(ind, n, min=0, max=n-1){
    ind/(n-1)*(max-min) + min
}
itran <- function(ind, n){
    qnorm((ind+0.5)/n)
}
itraninv <- function(y, n, min=0, max=n-1){
    pnorm(y)*n-0.5
}
##
varinfo <- cbind(
    min=c(0,0,0,0,0,0,0,0,0,0,0,0,0),
    max=c(1,1,1,6,15,15,63,50,75,150,300,Inf,Inf),
    n=c(2,2,2,7,16,16,64,51,76,Inf,Inf,Inf,Inf),
    location=signif(sapply(rtran(dt[,..varnames]),median,na.rm=T),3),
    scale=signif(sapply(rtran(dt[,..varnames]),IQR,na.rm=T),3),
    Q1=signif(sapply(rtran(dt[,..varnames]),tquant,probs=0.25,na.rm=T),3),
    Q3=signif(sapply(rtran(dt[,..varnames]),tquant,probs=0.75,na.rm=T),3)
)
rownames(varinfo) <- varnames

subsel <- 4:nrow(varinfo)
nclusters <- 64
##
set.seed(229)
idata <- cbind(
    dt[[4]],
    dt[[5]],
    dt[[6]],
    dt[[7]],
    dt[[8]],
    dt[[9]]
)
rdata <- cbind(
    log(dt[[12]]),
    log(dt[[13]])
)
nivars <- ncol(idata)
nrvars <- ncol(rdata)
ndata <- nrow(rdata)
nint <- ncol(idata)
iint <- 4:9
irea <- 12:13
##
initint <-  t((t(idata)-varinfo[iint,'min'])/(varinfo[iint,'max']-varinfo[iint,'min']))*2-1
initint[is.na(initint)] <- 0


##
nimax <- max(varinfo[iint,'n'])
mbreaks <- t(sapply(varinfo[iint,'n'], function(n){
    c(itran(0.5:(n-1.5), n), rep(+Inf, nimax-n))
}))

##
## censored <- foreach(datum=rdata[,1:ncens,drop=F], breaks=bounds, .combine=cbind)%do%{
##     sapply(datum, function(x){
##         which(x <= c(breaks,+Inf) & x > c(-Inf,breaks))-1
##     })
## }
## rdata2 <- rdata
## rdata[,1:ncens] <- NA
##

datapoints <- c(
    list()
)
ndata <- 2
##
constants <- list(ndata=ndata, nclusters=nclusters, nrvars=nrvars, nivars=nivars, mbreaks=mbreaks, maxnbreaks=nimax-1)
##
initsFunction <- function(){
    c(
        list(
            palphas = c(1/3,1/3,1/3),
            qalpha0 = rep(0.25/nclusters, nclusters), # cluster probabilities
            alpha0 = 3,
            q=rep(1/nclusters,nclusters),
            Labels = rep(1, ndata), # all data to first cluster
            ##
            Rrates = rinvgamma(nrvars, shape=0.5, scale=1),
            Rvariances = matrix(rinvgamma(nrvars*nclusters, shape=1, rate=1), nrow=nrvars, ncol=nclusters),
            Rmeans = matrix(rnorm(nrvars*nclusters, mean=0, sd=3), nrow=nrvars, ncol=nclusters),
            Irates = rinvgamma(nivars, shape=0.5, scale=1),
            Ivariances = matrix(rinvgamma(nivars*nclusters, shape=1, rate=1), nrow=nivars, ncol=nclusters),
            Imeans = matrix(rnorm(nivars*nclusters, mean=0, sd=3), nrow=nivars, ncol=nclusters),
            ##
            imean=0,
            ivariance=1*1,
            ishapemicro=1,
            ishapemacro=1,
            iscaleprec=1/8^-2,
            ##
            rmean=0,
            rvariance=3*3,
            rshapemicro=0.5,
            rshapemacro=0.5,
            rscaleprec=1^-2
        )
    )
}

infmixture <- nimbleCode({
    alpha0 ~ dcat(prob=palphas[1:3])
    qalpha[1:nclusters] <- 2^alpha0 * qalpha0[1:nclusters]
    q[1:nclusters] ~ ddirch(alpha=qalpha[1:nclusters])
    ##
    for(variate in 1:nrvars){
        Rrates[variate] ~ dinvgamma(shape=rshapemicro, scale=rscaleprec)
    }
    for(variate in 1:nivars){
        Irates[variate] ~ dinvgamma(shape=ishapemicro, scale=iscaleprec)
    }
    for(cluster in 1:nclusters){
        for(variate in 1:nrvars){# real variates
            Rmeans[variate, cluster] ~ dnorm(mean=rmean, var=rvariance)
            Rvariances[variate, cluster] ~ dinvgamma(shape=rshapemacro, rate=Rrates[variate])
        }
        for(variate in 1:nivars){# real variates
            Imeans[variate, cluster] ~ dnorm(mean=imean, var=ivariance)
            Ivariances[variate, cluster] ~ dinvgamma(shape=ishapemacro, rate=Irates[variate])
        }
    }
    ##
    for(datum in 1:ndata){
        Labels[datum] ~ dcat(prob=q[1:nclusters])
        ##
        for(variate in 1:nrvars){
            Rdata[datum, variate] ~ dnorm(mean=Rmeans[variate, Labels[datum]], var=Rvariances[variate, Labels[datum]])
        }
        for(variate in 1:nivars){# real variates
            Idata[datum, variate] ~ dinterval(Ihidden[datum, variate], mbreaks[variate,1:maxnbreaks])
            Ihidden[datum, variate] ~ dnorm(mean=Imeans[variate, Labels[datum]], var=Ivariances[variate, Labels[datum]])
        }
    }
})

model <- nimbleModel(code=infmixture, name='model',
                     constants=constants,
                     inits=initsFunction(),
                     data=datapoints,
                     dimensions=c(
                         list(
                              Labels=ndata,
                              q=nclusters,
                              Rrates=nrvars,
                              Irates=nivars,
                              Rmeans=c(nrvars, nclusters),
                              Imeans=c(nivars, nclusters),
                              Rvariances=c(nrvars, nclusters),
                             Ivariances=c(nivars, nclusters),
                             Ihidden=c(ndata,nivars),
                             Rdata=c(ndata,nrvars),
                             Idata=c(ndata,nivars)
                         )
                     )
                     )

Cmodel <- compileNimble(model, showCompilerOutput=FALSE)

confmodel <- configureMCMC(Cmodel,
                           monitors=c('q','alpha0',
                                      ## 'Labels',
                                      ## 'Rdata[2, 1]',
                                      ## 'Rdata[4, 2]',
                                      'Rrates', 'Rmeans', 'Rvariances',
                                      'Irates', 'Imeans', 'Ivariances'
                                      ),
                           monitors2=c('Rdata','Idata')
                           )

mcsampler <- buildMCMC(confmodel)
Cmcsampler <- compileNimble(mcsampler, resetFunctions = TRUE)

initsFunction <- function(){
    c(
        list(
            palphas = c(1/3,1/3,1/3),
            qalpha0 = rep(0.25/nclusters, nclusters), # cluster probabilities
            alpha0 = 3,
            q=rep(1/nclusters,nclusters),
            Labels = rep(1, ndata), # all data to first cluster
            ##
            Rrates = rinvgamma(nrvars, shape=0.5, scale=1),
            Rvariances = matrix(rinvgamma(nrvars*nclusters, shape=1, rate=1), nrow=nrvars, ncol=nclusters),
            Rmeans = matrix(rnorm(nrvars*nclusters, mean=0, sd=3), nrow=nrvars, ncol=nclusters),
            Irates = rinvgamma(nivars, shape=0.5, scale=1),
            Ivariances = matrix(rinvgamma(nivars*nclusters, shape=1, rate=1), nrow=nivars, ncol=nclusters),
            Imeans = matrix(rnorm(nivars*nclusters, mean=0, sd=3), nrow=nivars, ncol=nclusters),
            ##
            imean=0,
            ivariance=1^2,
            ishapemicro=1,
            ishapemacro=1,
            iscaleprec=1/16^-2,
            ##
            rmean=0,
            rvariance=2^2,
            rshapemicro=0.5,
            rshapemacro=0.5,
            rscaleprec=1/4^-2
        )
    )
}
##
timestart <- Sys.time()
Cmodel$setInits(initsFunction())
output <- Cmcsampler$run(niter=100, thin=1, reset=TRUE, resetMV=TRUE, nburnin=1, time=T)
output <- Cmcsampler$run(niter=1024*8, thin=1, reset=FALSE, resetMV=TRUE, nburnin=0, time=T)
timestop <- Sys.time() - timestart
timestop
##
dpsamples <- as.matrix(Cmcsampler$mvSamples)
dpsamples2 <- as.matrix(Cmcsampler$mvSamples2)
##
pdff('testpriors_cens')
for(i in 1:nivars){
    data <- itraninv0(c(dpsamples2[, paste0('Idata[1, ',i,']')],
                        dpsamples2[, paste0('Idata[2, ',i,']')]),
                      n=varinfo[iint,'n'][i],
                      min=varinfo[iint,'min'][i],
                      max=varinfo[iint,'max'][i])
    his <- thist(data, n='i')
    xgrid <- seq(varinfo[iint,][i,'min'], varinfo[iint,][i,'max'],
                 length.out=varinfo[iint,][i,'n'])
    tplot(x=xgrid, y=his$density, ylim=c(0,NA))
}
for(i in 1:nrvars){
    data <- rtraninv(c(dpsamples2[, paste0('Rdata[1, ',i,']')],
                       dpsamples2[, paste0('Rdata[2, ',i,']')]),
                      location=varinfo[irea,'location'][i],
                      scale=varinfo[irea,'scale'][i])
        data <- data[data<=1.5*max(dt[[rownames(varinfo)[irea][i]]],na.rm=T) &
                     data>=0.5*min(dt[[rownames(varinfo)[irea][i]]],na.rm=T)]
    his <- thist(data)
    ## xgrid <- seq(varinfo[iint,][i,'min'], varinfo[iint,][i,'max'],
    ##              length.out=varinfo[iint,][i,'n'])
    tplot(x=his$mids, y=his$density,ylim=c(0,NA))
}
dev.off()
    



times <- Cmcsampler$getTimes()
names(times) <- sapply(confmodel$getSamplers(),function(x)x$target)
## sum(times[grepl('^Rdata',names(times))])
##
prefs <- unique(sub('^([^[]+)(\\[.*\\])', '\\1', names(times)))
sort(sapply(prefs, function(x)sum(times[grepl(x,names(times))])),decreasing=T)
## with hidden:
## Time difference of 1.51792 mins
## Ivariances     Labels     Imeans Rvariances     Rmeans    Ihidden          q 
##  26.372392  21.895887  21.326742   8.408666   7.174604   4.418221   0.340741 
##     Irates     Rrates     alpha0      Rdata 
##   0.143012   0.053464   0.015036   0.009979 



## Cmodel$plotGraph()

library('nimble')

## Dirichlet-process mixture of product kernels
## Two kinds of variates, two variates per kind:
## real (kernel: gaussian, prior: normal + compound-gamma)
## categorical (kernel: categorical distr, prior: Dirichlet)

nvars <- 10 # switch to include/exclude categorical components
ndata <- 600
ncens <- 5
nclusters <- 64
##
set.seed(229)
rdata <- cbind(
dt2$AVDEL30MIN_neuro[1:ndata],
dt2$CATANIMSC_neuro[1:ndata],
dt2$GDTOTAL_gds[1:ndata],
dt2$ANARTERR_neuro[1:ndata],
dt2$RAVLT_immediate[1:ndata],
matrix(rnorm(ndata*(nvars-ncens), mean=0, sd=2), nrow=ndata)
)
##
breaks0 <- cbind(
(dt2$AVDEL30MIN_neuro[1:ndata] - 0)/(15-0),
(dt2$CATANIMSC_neuro[1:ndata] - 0)/(60-0),
(dt2$GDTOTAL_gds[1:ndata] - 0)/(6-0),
(dt2$ANARTERR_neuro[1:ndata] - 0)/(50-0),
(dt2$RAVLT_immediate[1:ndata] - 0)/(75-0)
)
##
seqs <- list( ((0:15)-0)/(15-0), ((0:60)-0)/(60-0), ((0:6)-0)/(6-0), ((0:50)-0)/(50-0), ((0:75)-0)/(75-0) )
bounds <- sapply(seqs, function(x){(x[-1]+x[-length(x)])/2})
##
## censored <- foreach(datum=rdata[,1:ncens,drop=F], breaks=bounds, .combine=cbind)%do%{
##     sapply(datum, function(x){
##         which(x <= c(breaks,+Inf) & x > c(-Inf,breaks))-1
##     })
## }
## rdata2 <- rdata
## rdata[,1:ncens] <- NA
##
nbounds <- sapply(bounds,length)
maxbounds <- max(nbounds)
mbounds <- t(sapply(bounds,function(x){c(x,rep(+Inf,maxbounds-length(x)))}))


## cbind(rdata2,rdata,censored)
## bounds[!is.na(rdata)] <- signs[!is.na(rdata)]*Inf


datapoints <- c(
    list(Rdata = rdata)
)
##
constants <- list(ndata=ndata, nclusters=nclusters, mbounds=mbounds, nvars=nvars, ncens=ncens, maxbounds=maxbounds, ncens1=ncens+1L)
##
initsFunction <- function(){
    c(
        list(
            palphas = c(1/3,1/3,1/3),
            qalpha0 = rep(0.25/nclusters, nclusters), # cluster probabilities
            alpha0 = 2,
            q=rep(1/nclusters,nclusters),
            Labels = rep(1, ndata), # all data to first cluster
            Rrates = rinvgamma(nvars, shape=0.5, scale=1),
            Rmeans = matrix(rnorm(nvars*nclusters, mean=0, sd=3), nrow=nvars, ncol=nclusters),
            Rvariances = matrix(rinvgamma(nvars*nclusters, shape=0.5, rate=1), nrow=nvars, ncol=nclusters),
            hidden=breaks0
        )
    )
}

infmixture <- nimbleCode({
    alpha0 ~ dcat(prob=palphas[1:3])
    qalpha[1:nclusters] <- 2^alpha0 * qalpha0[1:nclusters]
    q[1:nclusters] ~ ddirch(alpha=qalpha[1:nclusters])
    ##
    for(variate in 1:nvars){
        Rrates[variate] ~ dinvgamma(shape=0.5, scale=1)
    }
    for(cluster in 1:nclusters){
        for(variate in 1:nvars){# real variates
            Rmeans[variate, cluster] ~ dnorm(mean=0, var=3*3)
            Rvariances[variate, cluster] ~ dinvgamma(shape=0.5, rate=Rrates[variate])
        }
    }
    ##
    for(datum in 1:ndata){
        Labels[datum] ~ dcat(prob=q[1:nclusters])
        ##
        for(variate in 1:ncens){# real variates
            Rdata[datum, variate] ~ dinterval(hidden[datum, variate], mbounds[variate,maxbounds])
            hidden[datum, variate] ~ dnorm(mean=Rmeans[variate, Labels[datum]], var=Rvariances[variate, Labels[datum]])
        }
        for(variate in ncens1:nvars){
            Rdata[datum, variate] ~ dnorm(mean=Rmeans[variate, Labels[datum]], var=Rvariances[variate, Labels[datum]])
        }
    }
})

model <- nimbleModel(code=infmixture, name='model',
                     constants=constants,
                     inits=initsFunction(),
                     data=datapoints,
                     dimensions=c(
                         list(hidden=c(ndata,ncens),
                              Labels=ndata,
                              q=nclusters,
                              Rrates=nvars,
                              Rmeans=c(nvars, nclusters),
                              Rvariances=c(nvars, nclusters))
                     )
                     )

Cmodel <- compileNimble(model, showCompilerOutput=FALSE)


confmodel <- configureMCMC(Cmodel,
                           monitors=c('q','alpha0',
                                      ## 'Labels',
                                      ## 'Rdata[2, 1]',
                                      ## 'Rdata[4, 2]',
                                      'Rrates', 'Rmeans', 'Rvariances'
                           ))

mcsampler <- buildMCMC(confmodel)
Cmcsampler <- compileNimble(mcsampler, resetFunctions = TRUE)

##
timestart <- Sys.time()
Cmodel$setInits(initsFunction())
output <- Cmcsampler$run(niter=1100, thin=1, reset=T, resetMV=TRUE, nburnin=100, time=T)
timestop <- Sys.time() - timestart
timestop
## ncens=5, no hidden:
## Time difference of 1.50273 mins
## Rvariances     Rmeans     Labels      Rdata          q     Rrates     alpha0 
##  35.381435  30.354673  20.059967   2.242052   0.291502   0.232479   0.015150 
## ncens=5, hidden:
## Time difference of 1.56104 mins
## Rvariances     Rmeans     Labels     hidden          q     Rrates     alpha0 
##  36.541409  30.998275  20.630430   2.342141   0.323528   0.235517   0.015481 

dpsamples <- as.matrix(Cmcsampler$mvSamples)
times <- Cmcsampler$getTimes()
names(times) <- sapply(confmodel$getSamplers(),function(x)x$target)
## sum(times[grepl('^Rdata',names(times))])
##
prefs <- unique(sub('^([^[]+)(\\[.*\\])', '\\1', names(times)))
sort(sapply(prefs, function(x)sum(times[grepl(x,names(times))])),decreasing=T)



## Cmodel$plotGraph()

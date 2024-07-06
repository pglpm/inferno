library('nimble')

## Dirichlet-process mixture of product kernels
## Two kinds of variates, two variates per kind:
## real (kernel: gaussian, prior: normal + compound-gamma)
## categorical (kernel: categorical distr, prior: Dirichlet)

set.seed(229)
ncopies <- 10
bounds <- rbind(c(-1, 0, 1)+10, c(-1, 0, 1)-10)
rint <- rbind(c(1,2,3,4), c(4,3,2,1))
rcon <- rbind(c(-1.5,-0.5,0.5,1.5)+10, rev(c(-1.5,-0.5,0.5,1.5))-10)
rint <- matrix(rep(rint,ncopies),nrow=nrow(rint))
rcon <- matrix(rep(rcon,ncopies),nrow=nrow(rcon))
##
datapoints <- list(Rint = rint-1)

constants <- list(ndata=ncol(rint), bounds=bounds, nvars=nrow(rint), nbounds=ncol(bounds))
##
initsFunction <- function(){
        list(Rcon=rcon,
            Rmeans = c(-1,1)
        )
}

infmixture <- nimbleCode({
    for(variate in 1:nvars){# real variates
        Rmeans[variate] ~ dnorm(mean=0, var=1)
    }
    ##
    for(datum in 1:ndata){
        for(variate in 1:nvars){
            Rint[variate, datum] ~ dinterval(Rcon[variate, datum], bounds[variate, 1:nbounds])
            Rcon[variate, datum] ~ dnorm(mean=Rmeans[variate], var=1)
        }
    }
})

model <- nimbleModel(code=infmixture, name='test',
                     constants=constants,
                     inits=initsFunction(),
                     data=datapoints,
                     dimensions=c(
                         list(Rint=dim(rint),
                              Rcon=dim(rcon),
                              Rmeans=2
                     )
                     ))
                     

Cmodel <- compileNimble(model, showCompilerOutput=FALSE)


confmodel <- configureMCMC(Cmodel,
                           monitors=c('Rmeans','Rint','Rcon')
                           )

mcsampler <- buildMCMC(confmodel)
Cmcsampler <- compileNimble(mcsampler, resetFunctions = TRUE)

##
Cmodel$setInits(initsFunction())
output <- Cmcsampler$run(niter=1000, thin=1, reset=T, resetMV=TRUE, nburnin=0, time=T)
dpsamples <- as.matrix(Cmcsampler$mvSamples)

pdff('testplotscensor5')
for(var in colnames(dpsamples)){
    tplot(y=dpsamples[,var],main=paste0(var,signif(range(dpsamples[,var]),2),collapse=', '))
}
dev.off()


times <- Cmcsampler$getTimes()
names(times) <- sapply(confmodel$getSamplers(),function(x)x$target)
## sum(times[grepl('^Rdata',names(times))])
##
prefs <- unique(sub('^([^[]+)(\\[.*\\])', '\\1', names(times)))
sort(sapply(prefs, function(x)sum(times[grepl(x,names(times))])),decreasing=T)
##     Rcon   Rmeans 
## 0.050685 0.006235 



## Cmodel$plotGraph()

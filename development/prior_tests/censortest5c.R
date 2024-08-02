library('nimble')

## Dirichlet-process mixture of product kernels
## Two kinds of variates, two variates per kind:
## real (kernel: gaussian, prior: normal + compound-gamma)
## categorical (kernel: categorical distr, prior: Dirichlet)

set.seed(229)
ndata <- 700*8
half1 <- 8
half2 <- 32
bounds <- list(c((-half1):half1)+10, c((-half2):half2)-20)
rint <- list(c(1:(2*half1+2)), c(1:(2*half2+2)))
rcon <- list(((-half1-0.5):(half1+0.5))+10,((-half2-0.5):(half2+0.5))-20)
##
rint <- sapply(1:length(rint),function(x)rep(rint[[x]],ndata)[1:ndata])
rcon <- sapply(1:length(rcon),function(x)rep(rcon[[x]],ndata)[1:ndata])
rleft <- (sapply(1:ncol(rint),function(x){
    c(-Inf,bounds[[x]])[rint[,x]]
}))
rright <- (sapply(1:ncol(rint),function(x){
    c(bounds[[x]],+Inf)[rint[,x]]
}))
maxbounds <- max(sapply(bounds,length))
Rbounds <- sapply(1:length(bounds),function(x){
    c(bounds[[x]],rep(+Inf,maxbounds))[1:maxbounds]
})
##
datapoints <- list(Rint=rint-1)
##
constants <- list(ndata=nrow(rint), nvars=ncol(rint), Rbounds=Rbounds, maxbounds=maxbounds)
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
            Rint[datum, variate] ~ dinterval(Rcon[datum, variate], Rbounds[1:maxbounds, variate])
            Rcon[datum, variate] ~ dnorm(mean=Rmeans[variate], var=1)
        }
    }
})
##
model <- nimbleModel(code=infmixture, name='testinterval',
                     constants=constants,
                     inits=initsFunction(),
                     data=datapoints,
                     dimensions=c(
                         list(Rint=dim(rint),
                              Rcon=dim(rcon),
                              Rbounds=dim(Rbounds),
                              Rmeans=2
                     )
                     ))
Cmodel <- compileNimble(model, showCompilerOutput=FALSE)
##
confmodel <- configureMCMC(Cmodel,
                           monitors=c('Rmeans','Rcon','Rint')
                           )
mcsampler <- buildMCMC(confmodel)
Cmcsampler <- compileNimble(mcsampler, resetFunctions = TRUE)

##
Cmodel$setInits(initsFunction())
output <- Cmcsampler$run(niter=1000, thin=1, reset=T, resetMV=TRUE, nburnin=0, time=T)
dpsamples <- as.matrix(Cmcsampler$mvSamples)

pdff('testplotscensor5ci')
for(var in colnames(dpsamples)){
##    if(!grepl('Rright',var) && !grepl('Rleft',var)){
    if(grepl('\\[1[],]',var) || grepl('\\[2[],]',var)){
        tplot(y=dpsamples[,var],main=paste0(var,signif(range(dpsamples[,var]),4),collapse=', '),ylim=range(dpsamples[,var])+c(-1e-6,+1e-6))
        }
}
dev.off()


times <- Cmcsampler$getTimes()
names(times) <- sapply(confmodel$getSamplers(),function(x)x$target)
## sum(times[grepl('^Rdata',names(times))])
##
prefs <- unique(sub('^([^[]+)(\\[.*\\])', '\\1', names(times)))
sort(sapply(prefs, function(x)sum(times[grepl(x,names(times))])),decreasing=T)
## ndata=700*8
## 8.456406 0.732347
## 8.496815 0.778227 


## ndata=700*4
##     Rcon   Rmeans 
## 3.158189 0.213872 

## ndata=100
##     Rcon   Rmeans 
## 0.056316 0.006519 
##
## ndata=700
##     Rcon   Rmeans 
## 0.465951 0.052905 
##
## ndata=700*8
##     Rcon   Rmeans 
## 6.048791 0.385774 



## cmodel$plotGraph()

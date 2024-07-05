library('nimble')

## Dirichlet-process mixture of product kernels
## Two kinds of variates, two variates per kind:
## real (kernel: gaussian, prior: normal + compound-gamma)
## categorical (kernel: categorical distr, prior: Dirichlet)

set.seed(229)
ndata <- 700
half1 <- 8
half2 <- 32
bounds <- list(c((-half1):half1)+10, c((-half2):half2)-20)
rint <- list(c(1:(2*half1+2)), c(1:(2*half2+2)))
rcon <- list(((-half1-0.5):(half1+0.5))+10,((-half2-0.5):(half2+0.5))-20)
##
rint <- sapply(1:length(rint),function(x)rep(rint[[x]],ndata)[1:ndata])
rcon <- sapply(1:length(rcon),function(x)rep(rcon[[x]],ndata)[1:ndata])
rint[2,1] <- NA
rint[3,2] <- NA
rleft <- (sapply(1:ncol(rint),function(x){
    out <- c(-Inf,bounds[[x]])[rint[,x]]
    out[is.na(out)] <- -Inf
    out
}))
rright <- (sapply(1:ncol(rint),function(x){
    out <- c(bounds[[x]],+Inf)[rint[,x]]
    out[is.na(out)] <- +Inf
    out
}))
##
datapoints <- list(Rleft=rleft, Rright=rright, Raux=1L+0L*rint)
##
constants <- list(ndata=nrow(rleft), nvars=ncol(rleft))
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
            Raux[datum, variate] ~ dconstraint(Rcon[datum, variate] > Rleft[datum, variate] & Rcon[datum, variate] < Rright[datum, variate])
            Rcon[datum, variate] ~ dnorm(mean=Rmeans[variate], var=1)
        }
    }
})
##
model <- nimbleModel(code=infmixture, name='test',
                     constants=constants,
                     inits=initsFunction(),
                     data=datapoints,
                     dimensions=c(
                         list(Rleft=dim(rleft),
                              Rright=dim(rright),
                              Rcon=dim(rcon),
                              Raux=dim(rint),
                              Rmeans=2
                     )
                     ))
Cmodel <- compileNimble(model, showCompilerOutput=FALSE)
##
confmodel <- configureMCMC(Cmodel,
                           monitors=c('Rmeans','Rcon','Rleft', 'Rright', 'Raux')
                           )
mcsampler <- buildMCMC(confmodel)
Cmcsampler <- compileNimble(mcsampler, resetFunctions = TRUE)

##
Cmodel$setInits(initsFunction())
output <- Cmcsampler$run(niter=1000, thin=1, reset=T, resetMV=TRUE, nburnin=0, time=T)
dpsamples <- as.matrix(Cmcsampler$mvSamples)

pdff('testplotscensor6dcon')
for(var in colnames(dpsamples)){
##    if(!grepl('Rright',var) && !grepl('Rleft',var)){
    if((grepl('\\[1[],]',var) || grepl('\\[2[],]',var) || grepl('\\[3[],]',var)) && (grepl('Rmean',var) || grepl('Rcon',var))){
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
##     Rcon   Rmeans 
## 7.186929 0.882888
## 7.039068 0.864236
## 7.112684 0.876044 


## ndata=700*4
##     Rcon   Rmeans 
## 2.370685 0.275082

## ndata=100
##     Rcon   Rmeans 
## 0.050580 0.006521 
##
## ndata=700
##     Rcon   Rmeans 
## 0.433552 0.055872 
##
## ndata=700*8
##     Rcon   Rmeans 
## 6.669616 0.491349 

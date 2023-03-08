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
library('nimble')

datapoints <- list(cont = c(1, 0, NA, NA),
                   lefts = c(0.5, -0.5, 1.5, 3),
                   rights = c(1.5, 0.5, 2.5, +Inf),
                   aux = c(1, 1, 1, 1)
                   )

constants <- list(dvalues=c(1,2,3,4))

initsFunction <- function(){
    list(cont = c(NA, NA, 2, 4))
}

finitemix <- nimbleCode({
    Rmean ~ dnorm(0,3)
    Rrate ~ dinvgamma(shape=1, rate=1)
    Rvar ~ dinvgamma(shape=1, rate=Rrate)
    for(d in 1:4){
        aux[d] ~ dconstraint(cont[d] >= lefts[d] & cont[d] <= rights[d])
        cont[d] ~ dnorm(mean=Rmean, var=Rvar)
    }
})

finitemixnimble <- nimbleModel(code=finitemix, name='finitemixnimble1',
                               constants=constants,
                               inits=initsFunction(),
                               data=datapoints
                               )

Cfinitemixnimble <- compileNimble(finitemixnimble, showCompilerOutput=FALSE)

confnimble <- configureMCMC(Cfinitemixnimble,
                            monitors=c('cont', 'Rmean', 'Rvar', 'Rrate')
                            )

confnimble$getSamplerExecutionOrder()
confnimble$setSamplerExecutionOrder(rev(confnimble$getSamplerExecutionOrder()))

mcsampler <- buildMCMC(confnimble)
Cmcsampler <- compileNimble(mcsampler, resetFunctions = TRUE)

Cfinitemixnimble$setInits(initsFunction())
todelete <- Cmcsampler$run(niter=1000, thin=1, thin2=1, nburnin=0, time=FALSE, reset=TRUE, resetMV=TRUE)
mcsamples <- as.matrix(Cmcsampler$mvSamples)

pdff('testInimble')
for(i in 1:ncol(mcsamples)){
    tplot(y=mcsamples[,i], ylab=colnames(mcsamples)[i])
}
dev.off()


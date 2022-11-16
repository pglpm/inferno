#------------------------------------------------------------------------
#------------------------------------------------------------------------
#------------------------------------------------------------------------
#  user defined functions


dkernelProduct <- nimbleFunction(
  run = function(x = double(1),
                 mean = double(1),
                 var = double(1),
                 p = double(1),
                 q = double(1),
                 log = integer(0, default=0)) 
  {
    returnType(double(0))
    
    ldens <- dnorm(x[1], mean = mean[1], sd = sqrt(var[1]), log=TRUE) +
      dnorm(x[2], mean = mean[2], sd = sqrt(var[2]), log=TRUE) + 
      dcat(x[3], prob=p[1:4], log=TRUE) + 
      dcat(x[4], prob=q[1:4], log=TRUE)
    
    if(log) return(ldens)
    else return(exp(ldens)) 
  }
)
cdkernelProduct <- compileNimble(dkernelProduct)

rkernelProduct <- nimbleFunction(
  run = function(n = integer(0),
                 mean = double(1),
                 var = double(1),
                 p = double(1),
                 q = double(1)) 
  {
    returnType(double(1))
    
    out <- nimNumeric(4)
    
    out[1] <- rnorm(1, mean = mean[1], sd = sqrt(var[1]))
    out[2] <- rnorm(1, mean = mean[2], sd = sqrt(var[2]))
    out[3] <- rcat(1, prob=p[1:4])
    out[4] <- rcat(1, prob=q[1:4])
    
    return(out)
  }
)
crkernelProduct <- compileNimble(rkernelProduct)

ndata <- 5

infmixture <- nimbleCode({
  Labels[1:ndata] ~ dCRP(alpha, size=ndata)
  ##
  ## the nested inverse-gammas are equivalent to a compound-gamma for the variance
  for(variate in 1:2){
    Rrates[variate] ~ dinvgamma(shape=0.5, scale=1)
  }
  ## cluster parameters:
  for(cluster in 1:ndata){
    for(variate in 1:2){# real variates
      Rmeans[cluster, variate] ~ dnorm(mean=0, var=3*3)
      Rvariances[cluster, variate] ~ dinvgamma(shape=0.5, rate=Rrates[variate])
    }
    Pprobs[cluster, 1:4] ~ ddirch(alpha=Calphas[1, 1:4]) # categorical u
    Qprobs[cluster, 1:4] ~ ddirch(alpha=Calphas[2, 1:4]) # categorical v
  }
  ##
  for(datum in 1:ndata){
    data[datum, 1:4] ~ dkernelProduct(mean = Rmeans[Labels[datum], 1:2], 
                                      var = Rvariances[Labels[datum], 1:2],
                                      p = Pprobs[Labels[datum], 1:4],
                                      q = Qprobs[Labels[datum], 1:4])
  }
})


set.seed(1)
datapoints <- list(data = cbind(matrix(rnorm(ndata*2, mean=0, sd=1), nrow=ndata), 
                                matrix(sample(1:4, ndata*2, prob=rep(1, 4), replace=T), nrow=ndata)))
Calphas = matrix(1, nrow=2, ncol=4)

constants <- list(ndata=ndata, Calphas = Calphas)

initsFunction <- function(){
  c(
    list(
      alpha = 1,
      Labels = rep(1, ndata), # all data to first cluster
      Rrates = rinvgamma(2, shape=0.5, scale=1),
      Rmeans = matrix(rnorm(2*ndata, mean=0, sd=3), nrow=ndata),
      Rvariances = matrix(rinvgamma(2*ndata, shape=0.5, rate=1), nrow=ndata),
      Pprobs = rbind(rdirch(1, Calphas[1,]), rdirch(1, Calphas[1,]), rdirch(1, Calphas[1,]), rdirch(1, Calphas[1,]), rdirch(1, Calphas[1,])),
      Qprobs = rbind(rdirch(1, Calphas[2,]), rdirch(1, Calphas[2,]), rdirch(1, Calphas[2,]), rdirch(1, Calphas[2,]), rdirch(1, Calphas[2,]))
    )
  )
}



model <- nimbleModel(code=infmixture, name='model',
                     constants=constants,
                     inits=initsFunction(),
                     data=datapoints,
                     dimensions=c(
                       list(Labels=ndata,
                            Rrates=2,
                            Rmeans=c(ndata, 2),
                            Rvariances=c(ndata, 2),
                            Pprobs = c(ndata, 4),
                            Qprobs = c(ndata, 4))
                     )
)
Cmodel <- compileNimble(model, showCompilerOutput=FALSE)


confmodel <- configureMCMC(Cmodel,
                           monitors=c('Labels',
                                      'Rrates', 'Rmeans', 'Rvariances', 'Pprobs', 'Qprobs'))

mcsampler <- buildMCMC(confmodel)
Cmcsampler <- compileNimble(mcsampler, resetFunctions = TRUE)

Cmodel$setInits(initsFunction())
output <- runMCMC(Cmcsampler, niter=100, thin=1) 

dpsamples <- getSamplesDPmeasure(Cmcsampler)



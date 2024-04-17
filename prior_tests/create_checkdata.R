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
    
sdovermad2 <- 0.5/qnorm(0.75)

Qfunction <- readRDS('Qfunction512.rsd')

set.seed(321)
ndata <- 105
ncomponents <- 3
alpha0 <- 2^((-3):3)
rmean0 <- 0
rvar0 <- (2*sdovermad2)^2
rshapein0 <- 1 # large scales
rshapeout0 <- 1 # small scales
hwidth <- 2 # number of powers of 2 to consider in either direction
rvarscales <- (sdovermad2 * 2^((-hwidth):hwidth))^2
##
dmin <- -3 # left boundary cens
dmax <- 3 # right boundary cens
intn <- 32 # number of integers
ibreaks <- Qfunction(seq(0,1,length.out=intn+1))
catn <- 8 # number of categories
##
compweights <- extraDistr::rdirichlet(n=1, alpha=rep(ncomponents,ncomponents))
rmus <- rnorm(n=ncomponents, mean=0, sd=2*sdovermad2)
rsigmas <- sqrt(extraDistr::rbetapr(n=ncomponents, shape1=1, shape2=1, scale=sample(rvarscales, ncomponents, replace=T)))
dmus <- rnorm(n=ncomponents, mean=0, sd=2*sdovermad2)
dsigmas <- sqrt(extraDistr::rbetapr(n=ncomponents, shape1=1, shape2=1, scale=sample(rvarscales, ncomponents, replace=T)))
imus <- rnorm(n=ncomponents, mean=0, sd=2*sdovermad2)
isigmas <- sqrt(extraDistr::rbetapr(n=ncomponents, shape1=1, shape2=1, scale=sample(rvarscales, ncomponents, replace=T)))
bqs <- rbeta(n=ncomponents, shape1=1, shape2=1)
cqs <- extraDistr::rdirichlet(n=ncomponents, alpha=rep(1,catn))
##
##
comps <- sample(1:ncomponents, ndata, prob=compweights, replace=T)
checkdata <- data.table(
    bin=extraDistr::rbern(n=ndata, prob=bqs[comps]),
    cat=extraDistr::rcat(n=ndata, prob=cqs[comps,]),
    rea=rnorm(n=ndata, mean=rmus[comps], sd=rsigmas[comps]),
    del=pmin(pmax(
        rnorm(n=ndata, mean=dmus[comps], sd=dsigmas[comps])
      , dmin), dmax),
    int=sapply(rnorm(n=ndata, mean=imus[comps], sd=isigmas[comps]),
           function(xxx){sum(xxx > ibreaks)})
)
for(i in 1:5){checkdata[ndata+i-5,i] <- NA}
##
saveRDS(checkdata,'checkdatawithNA.rds')
saveRDS(checkdata[1:(ndata-5)],'checkdatawithoutNA.rds')

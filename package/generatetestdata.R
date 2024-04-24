set.seed(16)
ndata <- 15
## Adding NAs in different places
## to check that missing data are handled correctly

## Continuous variate
Rdata <- rnorm(n=ndata, mean=0, sd=1)
Rdata[2] <- NA

## Censored variate, [-1, 1]
Cdata <- rnorm(n=ndata, mean=0, sd=1.5)
Cdata[Cdata <= -1] <- -1
Cdata[Cdata >= 1] <- 1
Cdata[3] <- NA

## Continuous-rounded variate 0.1
Ddata <- round(rnorm(n=ndata, mean=0, sd=1), 1)
Ddata[4] <- NA

## Ordinal variate, 7 values 1--7
Odata <- sample(1:7, ndata, replace=T, prob=1:7)
Odata[5] <- NA

## Nominal variate, 5 values 'A'-'E'
Ndata <- sample(LETTERS[1:5], ndata, replace=T, prob=LaplacesDemon::rdirichlet(alpha=rep(1,5),n=1))
Ndata[6] <- NA

## Binary variate
Bdata <- sample(c('no','yes'), ndata, replace=T, prob=1:2)
Bdata[7] <- NA

## Continuous, positive
Pdata <- exp(rnorm(n=ndata, mean=0, sd=1))
Pdata[8] <- NA

##
testdata <- data.table(Rvar=Rdata,Cvar=Cdata,Dvar=Ddata,Ovar=Odata,Nvar=Ndata,Bvar=Bdata,Pvar=Pdata)
fwrite(testdata, 'testdata.csv')

metadata <- list(
    name=c('Rvar', 'Cvar', 'Dvar','Ovar','Nvar','Bvar','Pvar'),
    type=c('continuous', 'continuous', 'continuous', 'ordinal', 'nominal', 'binary', 'continuous'),
    Nvalues=c(Inf, Inf, Inf, 7, 5, 2, Inf),
    rounding=c(0, 0, 0.1, NA, NA, NA, 0),
    domainmin=c(-Inf, -1, -Inf, 1, NA, NA, 0),
    domainmax=c(+Inf, +1, +Inf, 7, NA, NA, +Inf),
    minincluded=c(FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE),
    maxincluded=c(FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE),
    centralvalue=c(0, 0, 0, 5, NA, NA, 1),
    lowvalue=c(-0.7, -1, -0.7, 4, NA, NA, 0.5),
    highvalue=c(0.7, 1, 0.7, 7, NA, NA, 2.0),
    plotmin=c(-3, -1, -3, NA, NA, NA, 0),
    plotmax=c(+3, +1, +3, NA, NA, NA, 10),
    V1=c(NA, NA, NA, NA, 'A', 'no', NA),
    V2=c(NA, NA, NA, NA, 'B', 'yes', NA),
    V3=c(NA, NA, NA, NA, 'C', NA, NA),
    V4=c(NA, NA, NA, NA, 'D', NA, NA),
    V5=c(NA, NA, NA, NA, 'E', NA, NA)
)
fwrite(metadata, 'metatestdata.csv')

auxmeta <- buildauxmetadata(data=testdata, metadata=metadata, file=F)


source('bnpi.R')
test <- inferpopulation(data=testdata, metadata=as.data.table(metadata), outputdir='__debug', nsamples=10, nchains=2, cleanup=F, parallel=2)

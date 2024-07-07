set.seed(16)
ndata <- 15# 15
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

## Nominal variate 2, 3 values 'a'-'c'
Ndata2 <- sample(letters[1:3], ndata, replace=T, prob=LaplacesDemon::rdirichlet(alpha=rep(1,3),n=1))
Ndata2[9] <- NA

## Continuous variate, bounded domain
Tdata <- plogis(rnorm(n=ndata, mean=0, sd=1))
Tdata[10] <- NA


##
testdata <- data.frame(Rvrt=Rdata,Cvrt=Cdata,Dvrt=Ddata,Ovrt=Odata,Nvrt=Ndata,Bvrt=Bdata,Pvrt=Pdata,Nvrt2=Ndata2,Tvrt=Tdata)
write.csv(testdata, paste0('testdata_', ndata, '.csv'),
          row.names=FALSE, quote=FALSE, na='')

metadata <- list(
    name=c('Rvrt', 'Cvrt', 'Dvrt','Ovrt','Nvrt','Bvrt','Pvrt','Nvrt2','Tvrt'),
    type=c('continuous', 'continuous', 'continuous', 'ordinal', 'nominal', 'binary', 'continuous','nominal','continuous'),
    Nvalues=c(Inf, Inf, Inf, 7, 5, 2, Inf, 3, Inf),
    rounding=c(0, 0, 0.1, NA, NA, NA, 0, NA, NA),
    domainmin=c(-Inf, -1, -Inf, 1, NA, NA, 0, NA, 0),
    domainmax=c(+Inf, +1, +Inf, 7, NA, NA, +Inf, NA, 1),
    minincluded=c(FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE),
    maxincluded=c(FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE),
    centralvalue=c(0, 0, 0, 5, NA, NA, 1, NA, 0.5),
    lowvalue=c(-0.7, -1, -0.7, 4, NA, NA, 0.5, NA, 0.34),
    highvalue=c(0.7, 1, 0.7, 7, NA, NA, 2.0, NA, 0.66),
    plotmin=c(-3, -1, -3, NA, NA, NA, 0, NA, 0),
    plotmax=c(+3, +1, +3, NA, NA, NA, 10, NA, 1),
    V1=c(NA, NA, NA, NA, 'A', 'no', NA, 'a', NA),
    V2=c(NA, NA, NA, NA, 'B', 'yes', NA, 'b', NA),
    V3=c(NA, NA, NA, NA, 'C', NA, NA, 'c', NA),
    V4=c(NA, NA, NA, NA, 'D', NA, NA, NA, NA),
    V5=c(NA, NA, NA, NA, 'E', NA, NA, NA, NA)
)
write.csv(metadata, 'metatestdata.csv', row.names=FALSE, quote=FALSE, na='')

## source('buildauxmetadata.R')
## auxmeta <- buildauxmetadata(data=testdata, metadata=metadata, file=F)



set.seed(16)
ndata <- 30
## Adding NAs in different places
## to check that missing data are handled correctly

## Continuous variate
Rdata <- rnorm(n=ndata, mean=2, sd=3)

## Closed-domain vvariate, [-1, 1]
Cdata <- rnorm(n=ndata, mean=0, sd=1.5)
Cdata[Cdata <= -1] <- -1
Cdata[Cdata >= 1] <- 1

## Continuous, positive
Pdata <- exp(rnorm(n=ndata, mean=1, sd=2))

## Continuous-rounded variate 0.1
Ddata <- round(rnorm(n=ndata, mean=2, sd=3), 1)

## Ordinal variate, 7 values 1--7
Odata <- sample(1:7, ndata, replace=T, prob=1:7)

## Nominal variate, 5 values 'A'-'E'
Ndata <- sample(LETTERS[1:5], ndata, replace=T, prob=LaplacesDemon::rdirichlet(alpha=rep(1,5),n=1))

## Binary variate
Bdata <- sample(c('no','yes'), ndata, replace=T, prob=1:2)

## Nominal variate 2, 3 values 'a'-'c'
Ndata2 <- sample(letters[1:3], ndata, replace=T, prob=LaplacesDemon::rdirichlet(alpha=rep(1,3),n=1))

## Continuous variate, bounded open domain
Tdata <- plogis(rnorm(n=ndata, mean=1, sd=2))


##
testdata <- data.frame(Rvrt=Rdata,
                       Cvrt=Cdata,
                       Pvrt=Pdata,
                       Tvrt=Tdata,
                       Dvrt=Ddata,
                       Ovrt=Odata,
                       Nvrt=Ndata,
                       Bvrt=Bdata,
                       Nvrt2=Ndata2)
for(i in 1:ncol(testdata)){testdata[nrow(testdata)+1-i,i] <- NA}

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
    V1=c(NA, NA, NA, '1', 'A', 'no', NA, 'a', NA),
    V2=c(NA, NA, NA, '2', 'B', 'yes', NA, 'b', NA),
    V3=c(NA, NA, NA, '3', 'C', NA, NA, 'c', NA),
    V4=c(NA, NA, NA, '4', 'D', NA, NA, NA, NA),
    V5=c(NA, NA, NA, '5', 'E', NA, NA, NA, NA),
    V6=c(NA, NA, NA, '6', NA, NA, NA, NA, NA),
    V7=c(NA, NA, NA, '7', NA, NA, NA, NA, NA)
)
write.csv(metadata, 'metatestdata.csv', row.names=FALSE, quote=FALSE, na='')

## source('buildauxmetadata.R')
## auxmeta <- buildauxmetadata(data=testdata, metadata=metadata, file=F)



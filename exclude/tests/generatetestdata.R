set.seed(16)
pdff('probabilities_testdata')
ndata <- 1e4
## Adding NAs in different places
## to check that missing data are handled correctly

## Continuous variate
means <- c(0, 3)
sds <- c(1, 0.5)
clus <- sample(1:2, ndata, prob=c(0.75, 0.25), replace=TRUE)
Rdata <- rnorm(n=ndata, mean=means[clus], sd=sds[clus])
thist(Rdata, plot=T, border.alpha=1, xlab='Rvrt', ylab=NA)

## Closed-domain vvariate, [-1, 1]
clus <- sample(1:2, ndata, prob=c(0.75, 0.25), replace=TRUE)
Cdata <- rnorm(n=ndata, mean=means[clus], sd=sds[clus])
Cdata[Cdata <= -1] <- -1
Cdata[Cdata >= 4] <- 4
thist(Cdata, plot=T, border.alpha=1, xlab='Cvrt', ylab=NA)

## Continuous, positive
clus <- sample(2:1, ndata, prob=c(0.75, 0.25), replace=TRUE)
Rpdata <- 2^(rnorm(n=ndata, mean=means[clus], sd=sds[clus]))
thist(Rpdata, plot=T, border.alpha=1, xlab='Rpvrt', ylab=NA)

## Continuous variate, bounded open domain
clus <- sample(1:2, ndata, prob=c(0.75, 0.25), replace=TRUE)
Rudata <- plogis(rnorm(n=ndata, mean=means[clus], sd=sds[clus]))
thist(Rudata, plot=T, border.alpha=1, xlab='Ruvrt', ylab=NA)

## Binary variate
Bdata <- sample(c('no','yes'), ndata, replace=T, prob=c(0.25, 0.75))
tplot(x=c('no','yes'), y=c(0.25,0.75), type='b', ylim=c(0,1), xlab='Bvrt', ylab=NA)

## Nominal variate, 5 values 'A'-'E'
## probs <- LaplacesDemon::rdirichlet(alpha=rep(1,5),n=1)
Ndata <- sample(LETTERS[1:5], ndata, replace=T,
                prob=c(0.06, 0.11, 0.34, 0.19, 0.3))
thist(Ndata, plot=T, border.alpha=1, xlab='Nvrt', ylab=NA)

## Ordinal variate, 7 values 1--7
clus <- sample(1:2, ndata, prob=c(0.75, 0.25), replace=TRUE)
probs <- round(((Rdata-(-3))/(4+3))*6 + 1)
probs <- probs[probs >= 1 & probs <= 7]
probs <- as.vector(table(probs))
Odata <- sample(1:7, ndata, replace=T, prob=probs)
tplot(x=1:7, y=probs, type='b', ylim=c(0,1), xlab='Ovrt', ylab=NA)

## Nominal variate, 4 values 'a'-'d'
## the values are completely determined by Rdata and Bdata
## this allows us to test the mutual information
N2data <- letters[1 + (Rdata > 0) + 2 * (Bdata == 'no')]
thist(N2data, plot=T, border.alpha=1, xlab='N2vrt', ylab=NA)
dev.off()

## ## Continuous-rounded variate 0.1
## Ddata <- round(rnorm(n=ndata, mean=2, sd=3), 1)

##
testdata <- data.frame(Rvrt=Rdata,
                       Cvrt=Cdata,
                       Rpvrt=Rpdata,
                       Ruvrt=Rudata,
                       Bvrt=Bdata,
                       Nvrt=Ndata,
                       Ovrt=Odata,
                       N2vrt=N2data)
for(i in 1:ncol(testdata)){testdata[nrow(testdata)+1-i,i] <- NA}

## modelfreeinference::buildmetadata(testdata,file='metadata_test_custom.csv')

for(ndata in c(30,150)){
  write.csv(testdata[1:ndata,,drop=F],
            paste0('data_test_custom_', ndata, '.csv'),
            row.names=FALSE, quote=FALSE, na='')
}



## metadata <- list(
##     name=c('Rvrt', 'Cvrt', 'Dvrt','Ovrt','Nvrt','Bvrt','Pvrt','Nvrt2','Tvrt'),
##     type=c('continuous', 'continuous', 'continuous', 'ordinal', 'nominal', 'binary', 'continuous','nominal','continuous'),
##     Nvalues=c(Inf, Inf, Inf, 7, 5, 2, Inf, 3, Inf),
##     rounding=c(0, 0, 0.1, NA, NA, NA, 0, NA, NA),
##     domainmin=c(-Inf, -1, -Inf, 1, NA, NA, 0, NA, 0),
##     domainmax=c(+Inf, +1, +Inf, 7, NA, NA, +Inf, NA, 1),
##     minincluded=c(FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE),
##     maxincluded=c(FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE),
##     centralvalue=c(0, 0, 0, 5, NA, NA, 1, NA, 0.5),
##     lowvalue=c(-0.7, -1, -0.7, 4, NA, NA, 0.5, NA, 0.34),
##     highvalue=c(0.7, 1, 0.7, 7, NA, NA, 2.0, NA, 0.66),
##     plotmin=c(-3, -1, -3, NA, NA, NA, 0, NA, 0),
##     plotmax=c(+3, +1, +3, NA, NA, NA, 10, NA, 1),
##     V1=c(NA, NA, NA, '1', 'A', 'no', NA, 'a', NA),
##     V2=c(NA, NA, NA, '2', 'B', 'yes', NA, 'b', NA),
##     V3=c(NA, NA, NA, '3', 'C', NA, NA, 'c', NA),
##     V4=c(NA, NA, NA, '4', 'D', NA, NA, NA, NA),
##     V5=c(NA, NA, NA, '5', 'E', NA, NA, NA, NA),
##     V6=c(NA, NA, NA, '6', NA, NA, NA, NA, NA),
##     V7=c(NA, NA, NA, '7', NA, NA, NA, NA, NA)
## )
## write.csv(metadata, 'metatestdata.csv', row.names=FALSE, quote=FALSE, na='')

## source('buildauxmetadata.R')
## auxmeta <- buildauxmetadata(data=testdata, metadata=metadata, file=F)



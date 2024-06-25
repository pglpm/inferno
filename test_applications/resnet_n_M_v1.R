dat <- read.csv(file='../test_applications/_resnet18_1703078471_SPH_9_BH_n2_M10_BH_n4_M8_BH_n4_M10_BH_n4_M12_BH_n6_M10_CUT_15000_values.csv', na.strings = '')
##
dat$Predictions <- NULL
dat$Confidence <- NULL
dat$Correct_Conf <- NULL
dat$Output <- NULL
dat$Truth <- NULL
dat <- dat[dat$Model != "PP13-Sphaleron-THR9-FRZ15-NB0-NSUBPALL",]
##
dat$n <- NA
dat$M <- NA
##
modvalues <- c("BH_n2_M10", "BH_n4_M8", "BH_n4_M10", "BH_n4_M12", "BH_n6_M10")
nvalues <- c(2,4,4,4,6)
Mvalues <- c(10,8,10,12,10)
for(i in seq_along(modvalues)){
  totake <- dat$Model == modvalues[i]
  dat$n[totake] <- nvalues[i]
  dat$M[totake] <- Mvalues[i]
}
dat$Model <- NULL


setwd('../R')
source('bnpi.R')

## buildmetadata(data=dat, file='../test_applications/metadata_n_M', backup_files = F)

metadatafile <- "../test_applications/metadata_n.csv"
##
set.seed(16)
ntrain <- 100
trainpoints <- sort(sample(1:nrow(dat), ntrain))
##
stopparallel()
startparallel(4)
##
outputdir <- file.path('..', 'test_applications', '_resnet_n')
mcoutput <- inferpopulation(data = dat[trainpoints,],
                           metadata = metadatafile,
                           outputdir = outputdir,
                           nsamples = 120, nchains = 120,
                           appendtimestamp = FALSE,
                           appendinfo = FALSE, cleanup = FALSE)

ndomain <- seq(0, 8, by=2)
nprior <- rep(1/length(ndomain), length(ndomain))
Xcols <- paste0('Output_class_', 0:5)

ntest <- 100
testpoints <- sort(sample(setdiff(1:nrow(dat), trainpoints), ntest))

stopparallel()
startparallel(4)
nlikelihood <- samplesFDistribution(Y = dat[testpoints,Xcols],
                                    X = cbind(n=rep(ndomain,
                                                    each=length(testpoints))),
                                    mcoutput = mcoutput)

dim(nlikelihood) <- c(ntest, length(ndomain),
                      prod(dim(nlikelihood))/(ntest*length(ndomain)))

nposterior <- apply(nlikelihood, c(1,3), function(x){
  x <- x*nprior
  x/sum(x)
})

pdff(file.path(outputdir, 'n_posterior'))
for(datum in seq_along(testpoints)){
  tplot(
  x = ndomain, y = nposterior[,datum,],
  xlab = 'n', ylab = 'probability',
  ylim = c(0, 1),
  lty = 1, # solid lines
  lwd = 1, # thin
  col = 7, # grey
  alpha = 0.8
) # quite transparent
##
## plot P(Y|X)
tplot(
  x = ndomain, y = rowMeans(nposterior[,datum,]),
  type = 'b', # line+points
  lty = 1, # solid
  lwd = 3, # thicker
  col = 1, # blue
  add = T
) # add to previous plot
  truevalue <- dat[testpoints[datum], 'n']
  text(
    x = 2.5, y = 1, adj = 0.5,
    labels = paste0(
      'Datapoint ', testpoints[datum],
      ' - True: ', truevalue
    )
  )
  tplot(
    x = dat[testpoints[datum], 'n'], y = 0,
    type = 'p', pch = 2, col=2, add = T
  )
}
dev.off()







Y <- cbind(
  Output_class_0 = +1,
  Output_class_1 = -1,
  Output_class_2 = -1,
  Output_class_3 = -1,
  Output_class_4 = -1,
  Output_class_5 = -1
)
##
stopparallel()
startparallel(4)
nlikelihood <- samplesFDistribution(Y = Y, X = cbind(n=ndomain),
                                    mcoutput = mcoutput)
##
nprior <- rep(1, length(X))/length(X)
nposterior <- nlikelihood * nprior
##
nposterior <- apply(nposterior, 2, function(x)x/sum(x))

## plot samples first
tplot(
  x = ndomain, y = nposterior,
  xlab = 'n', ylab = 'probability',
  ylim = c(0, 1),
  lty = 1, # solid lines
  lwd = 1, # thin
  col = 7, # grey
  alpha = 0.8
) # quite transparent
##
## plot P(Y|X)
tplot(
  x = ndomain, y = rowMeans(nposterior),
  type = 'b', # line+points
  lty = 1, # solid
  lwd = 3, # thicker
  col = 1, # blue
  add = T
) # add to previous plot
dev.off() # close pdf



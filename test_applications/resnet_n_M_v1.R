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

metadatafile <- "../test_applications/metadata_n_M.csv"
##
set.seed(16)
ntrain <- 100
trainpoints <- sort(sample(1:nrow(dat), ntrain))
##
stopparallel()
startparallel(4)
##
outputdir <- file.path('..', 'test_applications', '_resnet_nM')
mcoutput <- inferpopulation(data = dat[trainpoints,],
                           metadata = metadatafile,
                           outputdir = outputdir,
                           nsamples = 6, nchains = 6, parallel=2,
                           appendtimestamp = FALSE,
                           appendinfo = FALSE, cleanup = FALSE)
stopparallel()


mctraces <- readRDS('/home/pglpm/repositories/bayes_nonparametric_inference/test_applications/_resnet_nM/MCtraces.rds')

  traces <- mctraces[apply(mctraces, 1, function(x) { all(is.finite(x)) }), ]



nfmc <- readRDS('/home/pglpm/repositories/bayes_nonparametric_inference/test_applications/_resnet_nM/NONFINITEmcsamples--6_3-2-i0.rds')

allmcsamples <- readRDS('/home/pglpm/repositories/bayes_nonparametric_inference/test_applications/_resnet_nM/test_allmcsamples--06.rds')

tokeep <- readRDS('/home/pglpm/repositories/bayes_nonparametric_inference/test_applications/_resnet_nM/test_tokeep--06.rds')

traces <- readRDS('/home/pglpm/repositories/bayes_nonparametric_inference/test_applications/_resnet_nM/test_traces--6.rds')


test <- mcsubset(allmcsamples,tokeep)


test <- readRDS('/home/pglpm/repositories/bayes_nonparametric_inference/test_applications/_resnet_nM/NONFINITEmcsamples--06_3-2-i0.rds')

pdff('../test_applications/test')
for(nam in names(test)){
  if(grepl('mean', nam, fixed = TRUE)){
    vrt <- test[[nam]]
  }else{
    vrt <- log10(test[[nam]])
  }
  if(length(dim(vrt)) > 2){
    for(j in seq_len(dim(vrt)[1])){
      tplot(x=1:(dim(vrt)[3]),
            y=t(vrt[j,,]),
            type='l', lwd=1, lty=1, col=7,
            main=nam)
    }
  } else{
      tplot(x=1:(dim(vrt)[2]),
            y=t(vrt),
            type='l', lwd=1, lty=1, col=7,
            main=nam)
  }
}
dev.off()




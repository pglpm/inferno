#### parameters and hyperparameters for inferpopulation()
niter0 <- 1024L # initial iterations to try
#### Hyperparameters
nclusters <- 64L # ****
minalpha <- -3L
maxalpha <- 3L
Rshapelo <- 0.5
Rshapehi <- 0.5
Cshapelo <- 0.5
Cshapehi <- 0.5
Dshapelo <- 0.5
Dshapehi <- 0.5
Oshapelo <- 0.5
Oshapehi <- 0.5
Bshapelo <- 1
Bshapehi <- 1
Olambda <- 16^2
##
nalpha <- length(minalpha:maxalpha)
npoints <- nrow(data)

#### other options
Alphatoslice <- TRUE
Ktoslice <- FALSE
RWtoslice <- FALSE
##
showdata <- TRUE # 'histogram' 'scatter' FALSE TRUE
plotmeans <- TRUE # plot frequency averages
totsamples <- 'all' # 'all' number of samples if plotting frequency averages
showsamples <- 100 # number of samples to show.
showquantiles <- c(1,31)/32 # quantiles to show
showhyperparametertraces <- TRUE ##
showsamplertimes <- FALSE ##
family <- 'Palatino'

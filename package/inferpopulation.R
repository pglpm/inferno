#' Description
#'
#' @param data data.table object or filepath ##Can we make it so that this HAS to be data.table?
#' @param metadata
#' @param outputdir String, path to output file folder ## Rename to outputPrefix, also addSuffix?
#' @param nsamples Integer, nr of desired MC samples
#' @param nchains Integer, nr of MC chains
#' @param nsamplesperchain Integer, nr of MC samples per chain
#' @param parallel, Bool
#' @param niterini, Int
#' @param miniter, Int
#' @param maxiter, 
#' @param thinning ??
#' @param plottraces
#' @param showKtraces Bool, when true, it saves the Alpha parameter
#'  during sampling and plots its trace and histogram at the end.
#'  Keeping it to FALSE (default) saves a little computation time.
#' @param showAlphatraces Bool, : when true, it saves the K parameter
#'  more frequently during sampling and plots its trace and histogram
#'  at the end. Keeping it to FALSE (default) saves a little
#'  computation time.
#' @param seed Integer, random number generator seed. If left as default NULL,
#'  a random seed based on the system clock is used in the set.seed() function
#' @param loglikelihood Positive integer or FALSE, default FALSE
#' @param subsampledata ??
#' @param useOquantiles Bool, ??
#' @param output Bool, ??
#' @param cleanup Bool, default TRUE, removes files that can be used for debugging
#' @return ??
#' @import foreach doParallel doRNG data.table LaplacesDemon
inferpopulation <- function(data, metadata, outputdir, nsamples = 0,
                            nchains = 0, nsamplesperchain = 0, parallel = TRUE,
                            niterini = 1024, miniter = 0, maxiter = +Inf,
                            thinning = 0, plottraces = TRUE,
                            showKtraces = FALSE, showAlphatraces = FALSE,
                            seed = NULL, loglikelihood = FALSE, subsampledata,
                            useOquantiles = FALSE, output = FALSE,
                            cleanup = TRUE, timestampdir = TRUE) {

  # Imports
  source('buildauxmetadata.R')
  source('samplesFDistribution.R')
  source('plotFsamples.R')
  source('tplotfunctions.R')
  source('vtransform.R')
  source('proposeburnin.R')
  source('proposethinning.R')
  source('mcsubset.R')
  
  ##################################################
  #### Define functions
  ##################################################
  # Create function printtime to format printing of time
  printtime <- function(tim) {
    paste0(signif(tim, 2), ' ', attr(tim, 'units'))
  }
  # Specific print message function
  printnull <- function(message, outcon) {
    sink(NULL, type = 'message')
    message(message, appendLF = FALSE)
    flush.console()
    sink(outcon, type = 'message')
  }

  cat('\n')

  ##################################################
  #### Setup parallel processing and start timer
  ##################################################

  #### Determine the status of parallel processing
  if (is.logical(parallel) && parallel) {
    if (foreach::getDoParRegistered()) {
      cat(
        'Using already registered', foreach::getDoParName(), 'with',
        foreach::getDoParWorkers(), 'workers\n'
      )
      ncores <- foreach::getDoParWorkers()
    } else {
      cat('No parallel backend registered.\n')
      ncores <- 1
    }
  } else if (is.numeric(parallel) && parallel >= 2) {
    if (foreach::getDoParRegistered()) {
      ncores <- min(foreach::getDoParWorkers(), parallel)
      cat(
        'Using already registered', foreach::getDoParName(), 'with',
        foreach::getDoParWorkers(), 'workers\n'
      )
    } else {
      cl <- parallel::makeCluster(parallel)
      doParallel::registerDoParallel(cl)
      cat(
        'Registered', foreach::getDoParName(), 'with',
        foreach::getDoParWorkers(), 'workers\n'
      )
      ncores <- parallel
    }
  } else {
    cat('No parallel backend registered.\n')
    ncores <- 1
  }

  #### Consistency checks for numbers of samples, chains, cores
  if (nsamples != 0 && nchains != 0 && nsamplesperchain != 0) {
    # Only two out of these three arguments can be given
    stop('Please specify exactly two among "nsamples", "nchains",
         "nsamplesperchain"')
  }
  # nsamples and nchains given
  if (nsamples > 0 && nchains > 0) {
    # Set nchainspercore
    nchainspercore <- ceiling(nchains / ncores)
    if (nchainspercore * ncores > nchains) {
      # Increase nr of chains to have equal nr of chains per core
      nchains <- nchainspercore * ncores
      cat('Increasing number of chains to', nchains, '\n')
    }
    nsamplesperchain <- ceiling(nsamples / nchains)
    if (nsamplesperchain * nchains > nsamples) {
      # Increase nr of samples to have equal nr of samples per chain
      nsamples <- nchains * nsamplesperchain
      cat('Increasing number of samples to', nsamples, '\n')
    }
    # nsamples and nsamplesperchain given
  } else if (nsamples > 0 && nsamplesperchain > 0) {
    # Set nchains
    nchains <- ceiling(nsamples / nsamplesperchain)
    nchainspercore <- ceiling(nchains / ncores)
    if (nchainspercore * ncores > nchains) {
      nchains <- nchainspercore * ncores
    }
    if (nsamplesperchain * nchains > nsamples) {
      nsamples <- nchains * nsamplesperchain
      cat('Increasing number of samples to', nsamples, '\n')
    }
    # nchains and nsamplesperchain given
  } else if (nchains > 0 && nsamplesperchain > 0) {
    nchainspercore <- ceiling(nchains / ncores)
    if (nchainspercore * ncores > nchains) {
      nchains <- nchainspercore * ncores
      cat('Increasing number of chains to', nchains, '\n')
    }
    # Set nsamples
    nsamples <- nchains * nsamplesperchain
  } else {
    stop('Make sure to set two of  "nsamples", "nchains", "nsamplesperchain"
          to positive integer values')
  }
  # We don't need more cores than chains
  if (nchains < ncores) {
    ncores <- nchains
  }

  # Parallellisation if more than one core
  if (ncores < 2) {
    `%dochains%` <- `%do%`
  } else {
    `%dochains%` <- `%dorng%`
  }

  #### Start timer
  timestart0 <- Sys.time()

  ##################################################
  #### Read and process metadata and data
  ##################################################

  #### Build auxiliary metadata object; we'll save it later
  if (is.character(metadata) && file.exists(metadata)) {
    metadata <- data.table::fread(metadata, na.strings = '')
  }
  metadata <- data.table::as.data.table(metadata)
  auxmetadata <- buildauxmetadata(data = data, metadata = metadata)
  cat('Calculating auxiliary metadata\n')

  #### Read dataset
  datafile <- NULL
  # Check if data is missing, or set to some other 'non-value'
  if (missing(data) || is.null(data) || is.logical(data) && data == FALSE) {
    message('Missing data: calculating prior distribution')
    data <- data.table::as.data.table(matrix(NA,
      nrow = 1, ncol = nrow(auxmetadata),
      dimnames = list(NULL, auxmetadata[['name']])
    ))
    loglikelihood <- FALSE
  }
  # Data could also be a file
  if (is.character(data)) {
    datafile <- paste0(sub('.csv$', '', data), '.csv')
    if (file.exists(data)) {
      data <- data.table::fread(datafile, na.strings = '')
    } else {
      stop('Cannot find data file')
    }
  }
  # Set data variable to data.table format
  data <- data.table::as.data.table(data)

  if (!missing(subsampledata) && is.numeric(subsampledata)) {
    ## @@TODO: use faster and memory-saving subsetting
    data <- data[sample(seq_len(nrow(data)),
                   min(subsampledata, nrow(data)),
                   replace = FALSE
                 ), ]
  }

  #### Correct loglikelihood argument if necessary
  # Find which rows have no missing entries
  dataNoNa <- which(apply(data, 1,
    function(xx) {
      !any(is.na(xx))
    }
  )
  )
  nRowsNoNa <- length(dataNoNa)
  if (is.numeric(loglikelihood) && loglikelihood > 1 && nRowsNoNa > 1) {
    # Make sure loglikelihood is not larger than the
    # nr of rows without missing data
    loglikelihood <- min(round(loglikelihood), nRowsNoNa)
  } else if (is.logical(loglikelihood) && loglikelihood && nRowsNoNa > 1) {
    # If loglikelihood is set to TRUE,
    # set length equal to nr of rows without missing data
    loglikelihood <- nRowsNoNa
  } else {
    loglikelihood <- FALSE
  }

  #### Check consistency of variate names
  if (!all(auxmetadata[['name']] %in% colnames(data))) {
    stop('Missing variates in data file. Check the metadata file.')
  }
  # Drop variate columns that are not in the metadata file
  if (!all(colnames(data) %in% auxmetadata[['name']])) {
    subvar <- intersect(colnames(data), auxmetadata[['name']])
    data <- data[, subvar, with = FALSE]
  }

  ##################################################
  #### Various internal parameters
  ##################################################

  #### Hyperparameters (from hyperparameters.R)
  nclusters <- 64L # ****
  minalpha <- -3L
  maxalpha <- 3L
  Rshapelo <- 0.5
  Rshapehi <- 0.5
  Rvarm1 <- 3L^2L
  Cshapelo <- 0.5
  Cshapehi <- 0.5
  Cvarm1 <- 3L^2L
  Dshapelo <- 0.5
  Dshapehi <- 0.5
  Dvarm1 <- 3L^2L
  Oshapelo <- 0.5
  Oshapehi <- 0.5
  Ovarm1 <- 3L^2L
  Bshapelo <- 1L
  Bshapehi <- 1L

  nalpha <- length(minalpha:maxalpha)
  npoints <- nrow(data)

  #### Other options
  Alphatoslice <- TRUE
  Ktoslice <- TRUE
  RWtoslice <- FALSE
  ##
  ## plotmeans <- TRUE # plot frequency averages
  ## totsamples <- 'all' # 'all' number of samples if plotting frequency averages
  showsamples <- 100 # number of samples to show.
  showquantiles <- c(1, 31) / 32 # quantiles to show
  nclustersamples <- 128 ## number of samples of Alpha and K
  showsamplertimes <- FALSE ##
  family <- 'Palatino'

  ##################################################
  #### Folder setup
  ##################################################

  ## append time and sampling info to name of output directory
  if (timestampdir) {
    timestamp <- paste0('-V', nrow(auxmetadata), '-D',
                        (if (npoints == 1 && all(is.na(data))) {
                          0
                        } else {
                          npoints
                        }), '-K', nclusters, '-S', nsamples)
  } else {
    timestamp <- NULL
  }

  if (missing(outputdir) || outputdir == TRUE) {
    outputdir <- paste0('_output_', datafile)
    outputdir <- paste0(sub('.csv$', '', outputdir))
  }

  nameroot <- paste0(outputdir, timestamp)
  ##
  dirname <- paste0(nameroot, '/')
  # Create output directory if it does not exist
  dir.create(dirname)
  # Print information
  cat(
    '\n', paste0(rep('*', max(nchar(dirname), 26)), collapse = ''),
    '\n Saving output in directory\n', dirname, '\n',
    paste0(rep('*', max(nchar(dirname), 26)), collapse = ''), '\n'
  )
  nameroot <- basename(nameroot)
  ## This is in case we need to add some extra specifier to the output files
  ## all 'dashnameroot' can be deleted in a final version
  dashnameroot <- NULL

  ## Save copy of metadata in directory
  fwrite(metadata, file = paste0(dirname, 'metadata.csv'))

 #####################################################
  #### Sort different kinds of variates
  #####################################################

  # R:rounded C:censored continous, D:discretized , O: ordinal,
  # N: nominal, B: binary
  nVars <- list() # How many variates of each type
  varsName <- list() # The names for variates of each type
  for (atype in c('R', 'C', 'D', 'O', 'N', 'B')) {
    nVars[[atype]] <- length(auxmetadata[mcmctype == atype, name])
    varsName[[atype]] <- auxmetadata[mcmctype == atype, name]
  }
  cat('\n ', nVars$R, nVars$C, nVars$D, nVars$O, nVars$N, nVars$B, '\n')

  #### CONSTANTS OF MONTE-CARLO SAMPLER
  # These constants are available in the Nimble environment
  # They don't have to be accessed by constants$varname
  constants <- c(
    list(
      nVars = nVars,
      nclusters = nclusters,
      npoints = npoints,
      nalpha = nalpha,
      probalpha0 = rep(1 / nalpha, nalpha),
      basealphas = rep((2^(minalpha - 1L)) / nclusters, nclusters)
    ),
    if (nVars$R > 0) { # continuous
      list(
        Rn = nVars$R, # This indexing variable is needed internally
        Rmean1 = rep(0, 1),
        Rvarm1 = rep(Rvarm1, 1),
        Rvar1 = rep(1, 1),
        Rshapelo = rep(Rshapelo, 1),
        Rshapehi = rep(Rshapehi, 1)
      )
    },
    if (nVars$C > 0) { # censored
      Cleft <- vtransform(data[, varsName$C, with = FALSE], auxmetadata,
                          Cout = 'left', useOquantiles = useOquantiles)
      Cright <- vtransform(data[, varsName$C, with = FALSE], auxmetadata,
                           Cout = 'right', useOquantiles = useOquantiles)
      list(
        Cn = nVars$C, # This indexing variable is needed internally
        Cmean1 = rep(0, 1),
        Cvarm1 = rep(Cvarm1, 1),
        Cvar1 = rep(1, 1),
        Cshapelo = rep(Cshapelo, 1),
        Cshapehi = rep(Cshapehi, 1),
        Cleft = Cleft,
        Cright = Cright
      )
    },
    if (nVars$D > 0) { # discretized
      Dleft <- vtransform(data[, varsName$D, with = FALSE], auxmetadata,
                          Dout = 'left', useOquantiles = useOquantiles)
      Dright <- vtransform(data[, varsName$D, with = FALSE], auxmetadata,
                            Dout = 'right', useOquantiles = useOquantiles)
      list(
        Dn = nVars$D, # This indexing variable is needed internally
        Dmean1 = rep(0, 1),
        Dvarm1 = rep(Dvarm1, 1),
        Dvar1 = rep(1, 1),
        Dshapelo = rep(Dshapelo, 1),
        Dshapehi = rep(Dshapehi, 1),
        Dleft = Dleft,
        Dright = Dright
      )
    },
    if (nVars$O > 0) { # ordinal
      Oleft <- vtransform(data[, varsName$O, with = FALSE], auxmetadata,
                          Oout = 'left', useOquantiles = useOquantiles)
      Oright <- vtransform(data[, varsName$O, with = FALSE], auxmetadata,
                           Oout = 'right', useOquantiles = useOquantiles)
      list(
        On = nVars$O, # This indexing variable is needed internally
        Omean1 = rep(0, 1),
        Ovarm1 = rep(Ovarm1, 1),
        Ovar1 = rep(1, 1),
        Oshapelo = rep(Oshapelo, 1),
        Oshapehi = rep(Oshapehi, 1),
        Oleft = Oleft,
        Oright = Oright
      )
    },
    if (nVars$N > 0) { # nominal
      Nmaxn <- max(auxmetadata[mcmctype == 'N', Nvalues])
      Nalpha0 <- matrix(1e-100, nrow = nVars$N, ncol = Nmaxn)
      for (avar in seq_along(varsName$N)) {
        nvalues <- auxmetadata[name == varsName$N[avar], Nvalues]
        ## use Hadamard-like prior: 1/nvalues
        ## other choice is flat prior: 1
        Nalpha0[avar, 1:nvalues] <- 1 / nvalues
      }
      list(
        Nn = nVars$N, # This indexing variable is needed internally
        Nmaxn = Nmaxn,
        Nalpha0 = Nalpha0
      )
    },
    if (nVars$B > 0) { # binary
      list(
        Bn = nVars$B, # This indexing variable is needed internally
        Bshapelo = rep(Bshapelo, 1),
        Bshapehi = rep(Bshapehi, 1)
      )
    }
  ) # End constants

 #### DATAPOINTS
  datapoints <- c(
    if (nVars$R > 0) { # continuous
      list(
        Rdata = vtransform(data[, varsName$R, with = FALSE], auxmetadata,
          useOquantiles = useOquantiles
        )
      )
    },
    if (nVars$C > 0) { # censored
      list(
        Caux = vtransform(data[, varsName$C, with = FALSE], auxmetadata,
          Cout = 'aux', useOquantiles = useOquantiles
        ),
        Clat = vtransform(data[, varsName$C, with = FALSE], auxmetadata,
          Cout = 'lat', useOquantiles = useOquantiles
        )
      )
    },
    if (nVars$D > 0) { # discretized
      list(
        Daux = vtransform(data[, varsName$D, with = FALSE], auxmetadata,
          Dout = 'aux', useOquantiles = useOquantiles
        )
      )
    },
    if (nVars$O > 0) { # ordinal
      list(
        Oaux = vtransform(data[, varsName$O, with = FALSE], auxmetadata,
          Oout = 'aux', useOquantiles = useOquantiles
        )
      )
    },
    if (nVars$N > 0) { # nominal
      list(
        Ndata = vtransform(data[, varsName$N, with = FALSE], auxmetadata,
          Nout = 'numeric', useOquantiles = useOquantiles
        )
      )
    },
    if (nVars$B > 0) { # binary
      list(
        Bdata = vtransform(data[, varsName$B, with = FALSE], auxmetadata,
          Bout = 'numeric', useOquantiles = useOquantiles
        )
      )
    }
  ) # End datapoints

  #### Output 
  if (!exists('Nalpha0')) {
    Nalpha0 <- cbind(1)
  }
  cat(
    'Starting Monte Carlo sampling of', nsamples, 'samples by',
    nchains, 'chains'
  )
  cat(
    '\nin a space of',
    (sum(as.numeric(nVars) * c(2, 2, 2, 2, 0, 1)) +
       sum(Nalpha0 > 2e-100) - nrow(Nalpha0) + 1) * nclusters - 1,
    '(effectively',
    paste0(
      (sum(as.numeric(nVars) * c(
        3 + npoints, 3 + npoints, 3 + npoints,
        3 + npoints, 0, 1 + npoints
      )) +
        sum(Nalpha0 > 2e-100) +
        nrow(Nalpha0) * (npoints - 1) + 1) * nclusters - 1 + nalpha - 1,
      ')'
    ), 'dimensions.\n'
  )
  cat(
    'Using', ncores, 'cores:',
    nsamplesperchain, 'samples per chain,', nchainspercore, 'chains per core.\n'
  )
  cat('Core logs are being saved in individual files.\n')
  cat('\nC-compiling samplers appropriate to the variates
      (package Nimble)\nthis can take tens of minutes with
       many data or variates.\n...\r')

  ## cat('Estimating remaining time...\r')
  ## stopCluster(cluster)
  ## stopImplicitCluster()
  ## registerDoSEQ()
  ## ## cl <- makePSOCKcluster(ncores)
  ## if(ncores > 1){
  ##     cl <- makeCluster(ncores)
  ##     registerDoParallel(cl)
  ## }else{
  ##     registerDoSEQ()
  ## }
  ## toexport <- c('constants', 'datapoints', 'nVars', 'varsName', 'nalpha', 'nclusters')
  ## toexport <- c('vtransform','samplesFDistribution','proposeburnin','proposethinning','plotFsamples')

  ## Set the RNG seed if given by user, or if no seed already exists
  if (!missing(seed) || !exists('.Random.seed')) {
    set.seed(seed)
  }
  ## Save current RNG seed in case needed by user
  saveRDS(.Random.seed,
          file = paste0(dirname, 'rng_seed', dashnameroot, '.rds'))

  #####################################################
  #### BEGINNING OF FOREACH LOOP OVER CORES
  #####################################################

  #Iterate over cores, using 'acore' variable as iterator
  chaininfo <- foreach::foreach(
    acore = 1:ncores, .combine = rbind, .inorder = FALSE,
    .packages = c('khroma', 'foreach', 'rngtools')
  ) %dochains% {

    #We have to source scripts again for each chain to be able to access them.
    source('tplotfunctions.R')
    source('vtransform.R')
    source('samplesFDistribution.R')
    source('proposeburnin.R')
    source('proposethinning.R')
    source('plotFsamples.R')
    source('mcsubset.R')

    # Create log file
    outcon <- file(paste0(
      dirname, '_log', dashnameroot,
      '-', acore, '.log'
    ), open = 'w')
    sink(outcon)
    sink(outcon, type = 'message')
    suppressPackageStartupMessages(library('data.table'))
    suppressPackageStartupMessages(library('nimble'))

    ## Function for diagnostics
    funMCSE <- function(x) {
      if (length(x) >= 1000) {
        LaplacesDemon::MCSE(x, method = 'batch.means')$se
      } else {
        LaplacesDemon::MCSE(x)
      }
    }

    ## Parameter and function to test MCMC convergence
    if (thinning <= 0) {
      multcorr <- 2L
      thinning <- 1L
      thresholdfn <- function(diagnIAT, nsamplesperchain, multcorr,
                              diagnBurn2, diagnThin) {
        ceiling(2 * max(diagnBurn2) +
                  ((nsamplesperchain - 1L) * multcorr *
                     ceiling(max(diagnIAT, diagnThin))))
      }
    } else if (thinning > 0) {
      # These two options do exactly the same thing
      multcorr <- (-2L)
      thresholdfn <- function(diagnIAT, nsamplesperchain, multcorr,
                              diagnBurn2, diagnThin) {
        ceiling(2 * max(diagnBurn2) +
                  ((nsamplesperchain - 1L) * (-multcorr) *
                     ceiling(max(diagnIAT, diagnThin))))
      }
      thinning <- ceiling(thinning)
    } else {
      stop('Invalid "thinning" argument.')
    }

    #### CLUSTER REPRESENTATION OF FREQUENCY SPACE

    ## hierarchical probability structure
    finitemix <- nimble::nimbleCode({
      ## Component weights
      Alpha ~ dcat(prob = probalpha0[1:nalpha])
      alphas[1:nclusters] <- basealphas[1:nclusters] * 2^Alpha
      W[1:nclusters] ~ ddirch(alpha = alphas[1:nclusters])

      ## Probability density for the parameters of the components
      # Loop over clusters
      for (k in 1:nclusters) {
        # Check for different types of variates
        if (nVars$R > 0) { # continuous
          for (v in 1:Rn) {
            Rmean[v, k] ~ dnorm(mean = Rmean1, var = Rvarm1)
            Rrate[v, k] ~ dinvgamma(shape = Rshapehi, rate = Rvar1)
            Rvar[v, k] ~ dinvgamma(shape = Rshapelo, rate = Rrate[v, k])
          }
        }
        if (nVars$C > 0) { # censored
          for (v in 1:Cn) {
            Cmean[v, k] ~ dnorm(mean = Cmean1, var = Cvarm1)
            Crate[v, k] ~ dinvgamma(shape = Cshapehi, rate = Cvar1)
            Cvar[v, k] ~ dinvgamma(shape = Cshapelo, rate = Crate[v, k])
          }
        }
        if (nVars$D > 0) { # discretized
          for (v in 1:Dn) {
            Dmean[v, k] ~ dnorm(mean = Dmean1, var = Dvarm1)
            Drate[v, k] ~ dinvgamma(shape = Dshapehi, rate = Dvar1)
            Dvar[v, k] ~ dinvgamma(shape = Dshapelo, rate = Drate[v, k])
          }
        }
        if (nVars$O > 0) { # ordinal
          for (v in 1:On) {
            Omean[v, k] ~ dnorm(mean = Omean1, var = Ovarm1)
            Orate[v, k] ~ dinvgamma(shape = Oshapehi, rate = Ovar1)
            Ovar[v, k] ~ dinvgamma(shape = Oshapelo, rate = Orate[v, k])
          }
        }
        if (nVars$N > 0) { # nominal
          for (v in 1:Nn) {
            Nprob[v, k, 1:Nmaxn] ~ ddirch(alpha = Nalpha0[v, 1:Nmaxn])
          }
        }
        if (nVars$B > 0) { # binary
          for (v in 1:Bn) {
            Bprob[v, k] ~ dbeta(shape1 = Bshapelo, shape2 = Bshapehi)
          }
        }
      }
      ## Probability of data
      for (d in 1:npoints) {
        K[d] ~ dcat(prob = W[1:nclusters])
        ##
        if (nVars$R > 0) { # continuous
          for (v in 1:Rn) {
            Rdata[d, v] ~ dnorm(mean = Rmean[v, K[d]], var = Rvar[v, K[d]])
          }
        }
        if (nVars$C > 0) { # censored
          for (v in 1:Cn) {
            Caux[d, v] ~ dconstraint(Clat[d, v] >= Cleft[d, v] &
                                       Clat[d, v] <= Cright[d, v])
            Clat[d, v] ~ dnorm(mean = Cmean[v, K[d]], var = Cvar[v, K[d]])
          }
        }
        if (nVars$D > 0) { # discretized
          for (v in 1:Dn) {
            Daux[d, v] ~ dconstraint(Dlat[d, v] >= Dleft[d, v] &
                                       Dlat[d, v] < Dright[d, v])
            Dlat[d, v] ~ dnorm(mean = Dmean[v, K[d]], var = Dvar[v, K[d]])
          }
        }
        if (nVars$O > 0) { # ordinal
          for (v in 1:On) {
            Oaux[d, v] ~ dconstraint(Olat[d, v] >= Oleft[d, v] &
                                       Olat[d, v] < Oright[d, v])
            Olat[d, v] ~ dnorm(mean = Omean[v, K[d]], var = Ovar[v, K[d]])
          }
        }
        if (nVars$N > 0) { # nominal
          for (v in 1:Nn) {
            Ndata[d, v] ~ dcat(prob = Nprob[v, K[d], 1:Nmaxn])
          }
        }
        if (nVars$B > 0) { # binary
          for (v in 1:Bn) {
            Bdata[d, v] ~ dbern(prob = Bprob[v, K[d]])
          }
        }
      }
    }) # end finitemix NimbleCode


    #### INITIAL-VALUE FUNCTION
    initsfn <- function() {
      Alpha <- sample(1:nalpha, 1, prob = constants$probalpha0, replace = TRUE)
      W <- nimble::rdirch(n = 1, alpha = constants$basealphas * 2^Alpha)
      outlist <- list(
        Alpha = Alpha,
        W = W,
        K = rep(sample(rep(which(W > 0), 2), 1, replace = TRUE), npoints)
      )
      ##
      if (nVars$R > 0) { # continuous
        Rrate <- matrix(
          nimble::rinvgamma(
            n = nVars$R * nclusters,
            shape = constants$Rshapehi,
            rate = constants$Rvar1
          ),
          nrow = nVars$R, ncol = nclusters
        )
        outlist <- c(
          outlist,
          list(
            Rmean = matrix(
              rnorm(
                n = nVars$R * nclusters,
                mean = constants$Rmean1,
                sd = sqrt(constants$Rvarm1)
              ),
              nrow = nVars$R, ncol = nclusters
            ),
            Rrate = Rrate,
            Rvar = matrix(
              nimble::rinvgamma(
                n = nVars$R * nclusters,
                shape = constants$Rshapelo,
                rate = Rrate
              ),
              nrow = nVars$R, ncol = nclusters
            )
          )
        )
      }
      if (nVars$C > 0) { # censored
        Crate <- matrix(
          nimble::rinvgamma(
            n = nVars$C * nclusters,
            shape = constants$Cshapehi,
            rate = constants$Cvar1
          ),
          nrow = nVars$C, ncol = nclusters
        )
        outlist <- c(
          outlist,
          list(
            Cmean = matrix(
              rnorm(
                n = nVars$C * nclusters,
                mean = constants$Cmean1,
                sd = sqrt(constants$Cvarm1)
              ),
              nrow = nVars$C, ncol = nclusters
            ),
            Crate = Crate,
            Cvar = matrix(
              nimble::rinvgamma(
                n = nVars$C * nclusters,
                shape = constants$Cshapelo,
                rate = Crate
              ),
              nrow = nVars$C, ncol = nclusters
            ),
            ## for data with boundary values
            Clat = vtransform(data[, varsName$C, with = FALSE],
              auxmetadata,
              Cout = 'init',
              useOquantiles = useOquantiles
            )
          )
        )
      }
      if (nVars$D > 0) { # discretized
        Drate <- matrix(
          nimble::rinvgamma(
            n = nVars$D * nclusters,
            shape = constants$Dshapehi,
            rate = constants$Dvar1
          ),
          nrow = nVars$D, ncol = nclusters
        )
        outlist <- c(
          outlist,
          list(
            Dmean = matrix(
              rnorm(
                n = nVars$D * nclusters,
                mean = constants$Dmean1,
                sd = sqrt(constants$Dvarm1)
              ),
              nrow = nVars$D, ncol = nclusters
            ),
            Drate = Drate,
            Dvar = matrix(
              nimble::rinvgamma(
                n = nVars$D * nclusters,
                shape = constants$Dshapelo,
                rate = Drate
              ),
              nrow = nVars$D, ncol = nclusters
            ),
            ## for data with boundary values
            Dlat = vtransform(data[, varsName$D, with = FALSE],
              auxmetadata,
              Dout = 'init',
              useOquantiles = useOquantiles
            )
          )
        )
      }
      if (nVars$O > 0) { # ordinal
        Orate <- matrix(
          nimble::rinvgamma(
            n = nVars$O * nclusters,
            shape = constants$Oshapehi,
            rate = constants$Ovar1
          ),
          nrow = nVars$O, ncol = nclusters
        )
        outlist <- c(
          outlist,
          list(
            Omean = matrix(
              rnorm(
                n = nVars$O * nclusters,
                mean = constants$Omean1,
                sd = sqrt(constants$Ovarm1)
              ),
              nrow = nVars$O, ncol = nclusters
            ),
            Orate = Orate,
            Ovar = matrix(
              nimble::rinvgamma(
                n = nVars$O * nclusters,
                shape = constants$Oshapelo,
                rate = Orate
              ),
              nrow = nVars$O, ncol = nclusters
            ),
            ## for data with boundary values
            Olat = vtransform(data[, varsName$O, with = FALSE],
              auxmetadata,
              Oout = 'init',
              useOquantiles = useOquantiles
            )
          )
        )
      }
      if (nVars$N > 0) { # nominal
        outlist <- c(
          outlist,
          list(
            Nprob = aperm(array(sapply(1:nVars$N, function(avar) {
              sapply(1:nclusters, function(aclus) {
                rdirch(n = 1, alpha = Nalpha0[avar, ])
              })
            }), dim = c(Nmaxn, nclusters, nVars$N)))
          )
        )
      }
      if (nVars$B > 0) { # binary
        outlist <- c(
          outlist,
          list(
            Bprob = matrix(
              rbeta(
                n = nVars$B * nclusters,
                shape1 = constants$Bshapelo,
                shape2 = constants$Bshapehi
              ),
              nrow = nVars$B, ncol = nclusters
            )
          )
        )
      }
      ##
      outlist
    } #End initfns

    # Timer
    timecount <- Sys.time()


    ##################################################
    #### NIMBLE SETUP
    ##################################################
    cat('Nimble setup \n')
    finitemixnimble <- nimbleModel(
      code = finitemix, name = 'finitemixnimble1',
      constants = constants,
      data = datapoints,
      inits = initsfn()
    )

    Cfinitemixnimble <- compileNimble(finitemixnimble,
                                      showCompilerOutput = FALSE)
    gc() #garbage collection

    confnimble <- configureMCMC(
      Cfinitemixnimble, # nodes = NULL
      monitors = c(
        'W',
        if (nVars$R > 0) {
          c('Rmean', 'Rvar')
        },
        if (nVars$C > 0) {
          c('Cmean', 'Cvar')
        },
        if (nVars$D > 0) {
          c('Dmean', 'Dvar')
        },
        if (nVars$O > 0) {
          c('Omean', 'Ovar')
        },
        if (nVars$N > 0) {
          c('Nprob')
        },
        if (nVars$B > 0) {
          c('Bprob')
        }
      ),
      ## It is necessary to monitor K to see if all clusters were used
      ## if 'showAlphatraces' is true then
      ## the Alpha-parameter trace is also recorded and shown
      monitors2 = c(if (showAlphatraces) {
        'Alpha'
      }, 'K')
    )

    targetslist <- sapply(confnimble$getSamplers(), function(xx) xx$target)
    nameslist <- sapply(confnimble$getSamplers(), function(xx) xx$name)
    ## cat('\n******** NAMESLIST', nameslist, '\n')

    ## replace Alpha's cat-sampler with slice
    if (Alphatoslice && !('Alpha' %in% targetslist[nameslist == 'posterior_predictive'])) {
      confnimble$removeSamplers('Alpha')
      confnimble$addSampler(target = 'Alpha', type = 'slice')
    }

    ## replace K's cat-sampler with slice
    if (Ktoslice) {
      for (asampler in grep('^K\\[', targetslist, value = TRUE)) {
        if (!(asampler %in% targetslist[nameslist == 'posterior_predictive'])) {
          confnimble$removeSamplers(asampler)
          confnimble$addSampler(target = asampler, type = 'slice')
        }
      }
    }

    ## replace all RW samplers with slice
    if (RWtoslice) {
      for (asampler in targetslist[nameslist == 'RW']) {
        confnimble$removeSamplers(asampler)
        confnimble$addSampler(target = asampler, type = 'slice')
      }
    }

    ## call this to do a first reordering
    mcsampler <- buildMCMC(confnimble)

    #### change execution order for some variates
    samplerorder <- c(
      'K',
      if (nVars$R > 0) {
        c('Rmean', 'Rrate', 'Rvar')
      },
      if (nVars$C > 0) {
        c('Cmean', 'Crate', 'Cvar')
      },
      if (nVars$D > 0) {
        c('Dmean', 'Drate', 'Dvar')
      },
      if (nVars$O > 0) {
        c('Omean', 'Orate', 'Ovar')
      },
      if (nVars$N > 0) {
        c('Nprob')
      },
      if (nVars$B > 0) {
        c('Bprob')
      },
      'W', 'Alpha'
    )
    ##
    neworder <- foreach::foreach(var = samplerorder, .combine = c) %do% {
      grep(
        paste0('^', var, '(\\[.+\\])*$'),
        sapply(confnimble$getSamplers(), function(x) {
          if (!(x$name == 'posterior_predictive')) {
            x$target
          } else {
            NULL
          }
        })
      )
    }
    ##
    ## cat('\n********NEW ORDER',neworder,'\n')
    confnimble$setSamplerExecutionOrder(c(setdiff(
      confnimble$getSamplerExecutionOrder(),
      neworder
    ), neworder))
    print(confnimble)
    ## print(confnimble$getUnsampledNodes())

    #### Compile Monte Carlo sampler
    mcsampler <- buildMCMC(confnimble)
    Cmcsampler <- compileNimble(mcsampler, resetFunctions = TRUE)

    cat('\nSetup time', printtime(Sys.time() - timecount), '\n')

    if (acore == 1) {
      printnull(paste0('\rCompiled core ', acore, '. Estimating
                       remaining time, please be patient...'), outcon)
    }

    cat('Loop over chains')
    ##################################################
    #### LOOP OVER CHAINS (WITHIN ONE CORE)
    ##################################################
    maxusedclusters <- 0

    ## build test data for assessing stationarity
    testdata <- data.table(sapply(seq_len(nrow(auxmetadata)), function(xx) {
      # consider making this function separate
      xx <- auxmetadata[xx, ]
      toadd <- xx[, paste0('mctest', 1:3), with = FALSE]
      if (xx[['mcmctype']] %in% c('B', 'N')) {
        toadd <- xx[1, paste0('V', toadd), with = FALSE]
      }
      toadd
    }))
    colnames(testdata) <- auxmetadata[['name']]

    # Start timer
    starttime <- Sys.time()

    #### Loop over chains in core
    allflagmc <- FALSE
    for (achain in 1:nchainspercore) {
      showsamplertimes0 <- showsamplertimes && (achain == 1)
      ## showAlphatraces0 <- showAlphatraces && (achain==1)
      niter <- min(niterini, maxiter)
      ## ## Experimental: decrease number of iterations based on previous chain
      ## if(achain == 1){
      ##   niter <- niterini
      ## }else{
      ##  niter <- max(min(niterini,lengthmeasure*2), 128)
      ##  }
      nitertot <- 0
      lengthmeasure <- +Inf
      reset <- TRUE
      traces <- NULL
      allmcsamples <- NULL
      allmcsamplesKA <- list(Alpha = NULL, K = NULL)
      mcsamplesKA <- NULL
      flagll <- FALSE
      flagmc <- FALSE
      if (is.numeric(loglikelihood)) {
        llseq <- sort(sample(dataNoNa, loglikelihood))
      }
      gc() #garbage collection
      chainnumber <- (acore - 1L) * nchainspercore + achain
      padchainnumber <- sprintf(paste0('%0', nchar(nchains), 'i'), chainnumber)
      cat(
        '\nChain #', chainnumber,
        '(chain', achain, 'of', nchainspercore, 'for this core)\n'
      )
      ## cat('Seed:', chainnumber + seed, '\n')
      ## commenting out this changes the randomness of each chain to be the same? AURORA
      ## set.seed(chainnumber + seed) 
      Cfinitemixnimble$setInits(initsfn())

      #### WHILE LOOP CONTINUING UNTIL CONVERGENCE
      subiter <- 1L
      while (nitertot < lengthmeasure) {
        cat('Iterations:', niter, '\n')

        #### MONTE-CARLO CALL
        ## If reporting Alpha or K traces,
        ## then save them more frequently
        ## Otherwise just save the last value
        Cmcsampler$run(
          niter = niter, thin = 1,
          thin2 = (if (showAlphatraces || showKtraces) {
            max(1, round(niter / nclustersamples))
          } else {
            niter
          }),
          nburnin = 0, time = showsamplertimes0,
          reset = reset, resetMV = TRUE
        )

        ## iterationAsLastIndex: See sect 7.7 of Nimble manual
        mcsamples <- as.list(Cmcsampler$mvSamples,
          iterationAsLastIndex = TRUE
        )
        ## ## saveRDS(mcsamples,'__mcsamplestest.rds') # for debug
        mcsamplesKA <- as.list(Cmcsampler$mvSamples2,
          iterationAsLastIndex = FALSE
        )
        ## ## saveRDS(mcsamplesKA,'__mcsamplesKAtest.rds') # for debug
        ## 'mcsamplesKA$K' contains the cluster identity of each training datapoint
        ## but we only want the number of distinct clusters used:
        mcsamplesKA$K <- apply(
          mcsamplesKA$K, 1,
          function(xx) length(unique(xx))
        )
        mcsamplesKA$Alpha <- c(mcsamplesKA$Alpha)
        if (showAlphatraces) {
          dim(mcsamplesKA$Alpha) <- NULL # from matrix to vector
        }

        cat('\nMCMC time', printtime(Sys.time() - starttime), '\n')

        ## #### Remove iterations with non-finite values
        ##                 toRemove <- which(!is.finite(mcsamples), arr.ind=TRUE)
        ##                 ##
        ##                 if(length(toRemove) > 0){
        ##                     printnull('\nWARNING: SOME NON-FINITE OUTPUTS\n', outcon)
        ##                     ##
        ##                     flagmc <- TRUE
        ##                     allflagmc <- TRUE
        ##                     if(length(unique(toRemove[,1])) == nrow(mcsamples)){
        ##                         suppressWarnings(sink())
        ##                         suppressWarnings(sink(NULL,type='message'))
        ##                         registerDoSEQ()
        ##                         if(ncores > 1){ stopCluster(cl) }
        ##                         stop('...TOO MANY NON-FINITE OUTPUTS. ABORTING')
        ##                     }else{
        ##                         mcsamples <- mcsamples[-unique(toRemove[,1]),,drop=FALSE]
        ##                     }
        ##                 }
        #### Remove iterations with non-finite values
        toRemove <- sort(unique(unlist(lapply(mcsamples, function(xx) {
          temp <- which(is.na(xx), arr.ind = TRUE)
          temp[, ncol(temp)]
        }))))
        if (length(toRemove) > 0) {
          cat('\nWARNING:', length(toRemove), 'NON-FINITE SAMPLES\n')
          ##
          flagmc <- TRUE
          allflagmc <- TRUE
          saveRDS(mcsamples, file = paste0(
            dirname, '_NONFINITEmcsamples',
            dashnameroot, '--', padchainnumber,
            '_', achain, '-',
            acore, '-i', nitertot, '.rds'
          ))
          if (length(toRemove) == ncol(mcsamples$W)) {
            cat('\n...TOO MANY NON-FINITE OUTPUTS!\n')
            ## printnull('\n...TOO MANY NON-FINITE OUTPUTS!\n', outcon)
            ## suppressWarnings(sink())
            ## suppressWarnings(sink(NULL,type='message'))
            ## ## registerDoSEQ()
            ## if(exists('cl')){ stopCluster(cl) }
            ## stop('...TOO MANY NON-FINITE OUTPUTS. ABORTING')
            mcsamples <- NULL
          } else {
            mcsamples <- mcsubset(mcsamples, -toRemove)
          }
        }


        ##
        if (showsamplertimes0) {
          samplertimes <- Cmcsampler$getTimes()
          names(samplertimes) <- sapply(
            confnimble$getSamplers(),
            function(x) x$target
          )
          sprefixes <- unique(sub(
            '^([^[]+)(\\[.*\\])', '\\1',
            names(samplertimes)
          ))
          cat('\nSampler times:\n')
          print(sort(sapply(sprefixes, function(x) sum(samplertimes[grepl(x, names(samplertimes))])),
            decreasing = TRUE
          ))
        }

        usedclusters <- mcsamplesKA$K[length(mcsamplesKA$K)]
        cat('\nOCCUPIED CLUSTERS:', usedclusters, 'OF', nclusters, '\n')

        #### Diagnostics
        ## Log-likelihood
        diagntime <- Sys.time()
        ##
        ll <- t(
          log(samplesFDistribution(
            Y = testdata, X = NULL,
            mcoutput = c(mcsamples, list(auxmetadata = auxmetadata)),
            jacobian = FALSE,
            useOquantiles = useOquantiles,
            parallel = FALSE,
            silent = TRUE
          ))
        )

        colnames(ll) <- paste0('log-', c('mid', 'lo', 'hi'))

        if (is.numeric(loglikelihood)) {
          lltime <- Sys.time()
          cat('\nCalculating log-likelihood...')
          ll <- cbind(ll,
            'log-ll' = log(samplesFDistribution(
              Y = data[llseq, ], X = NULL,
              mcoutput = c(mcsamples, list(auxmetadata = auxmetadata)),
              jacobian = FALSE,
              useOquantiles = useOquantiles,
              parallel = FALSE, silent = TRUE,
              combine = '+'
            )) / length(llseq)
          )
          cat('Done,\n', printtime(Sys.time() - lltime), '\n')
        }
        ##
        traces <- rbind(traces, 10 / log(10) * ll)
        ## tplot(y=traces, lwd=1, lty=1)
        ## for(i in 1:ncol(traces)){
        ##     tplot(y=traces[,i], main=colnames(traces)[i])
        ## }
        ## dev.off()
        toRemove <- which(!is.finite(traces), arr.ind = TRUE)
        if (length(toRemove) > 0) {
          flagll <- TRUE
          traces <- traces[-unique(toRemove[, 1]), , drop = FALSE]
        }
        ##
        diagnESS <- LaplacesDemon::ESS(traces)
        cat('\nESSs:', paste0(round(diagnESS), collapse = ', '))
        diagnIAT <- apply(traces, 2, function(x) {
          LaplacesDemon::IAT(x)
        })
        cat('\nIATs:', paste0(round(diagnIAT), collapse = ', '))
        diagnBMK <- LaplacesDemon::BMK.Diagnostic(traces[1:(4 * trunc(nrow(traces) / 4)), ], batches = 4)[, 1]
        cat('\nBMKs:', paste0(round(diagnBMK, 3), collapse = ', '))
        diagnMCSE <- 100 * apply(traces, 2, function(x) {
          funMCSE(x) / sd(x)
        })
        cat('\nMCSEs:', paste0(round(diagnMCSE, 2), collapse = ', '))
        diagnStat <- apply(traces, 2, function(x) {
          LaplacesDemon::is.stationary(as.matrix(x, ncol = 1))
        })
        cat('\nStationary:', paste0(diagnStat, collapse = ', '))
        diagnBurn <- apply(traces, 2, function(x) {
          LaplacesDemon::burnin(matrix(x[1:(10 * trunc(length(x) / 10))], ncol = 1))
        })
        cat('\nBurn-in I:', paste0(diagnBurn, collapse = ', '))
        diagnBurn2 <- proposeburnin(traces, batches = 10)
        cat('\nBurn-in II:', diagnBurn2)
        diagnThin <- proposethinning(traces)
        cat('\nProposed thinning:', paste0(diagnThin, collapse = ', '), '\n')

        cat('\nDiagnostics time', printtime(Sys.time() - diagntime), '\n')

        if (is.null(allmcsamples)) {
          # new chain
          allmcsamples <- mcsamples
        } else {
          # continue chain, concat samples
          allmcsamples <- mapply(
            function(xx, yy) {
              temp <- c(xx, yy)
              dx <- dim(xx)[-length(dim(xx))]
              dim(temp) <- c(dx, length(temp) / prod(dx))
              temp
            },
            allmcsamples, mcsamples,
            SIMPLIFY = FALSE
          )
          if (showAlphatraces || showKtraces) {
            ## Concatenate samples of K and Alpha
            allmcsamplesKA <- mapply(
              function(xx, yy) {
                c(xx, yy)
              },
              allmcsamplesKA, mcsamplesKA,
              SIMPLIFY = FALSE
            )
          }
        }

        if (showAlphatraces || showKtraces) {
          ## Concatenate samples of K and Alpha
          allmcsamplesKA <- mapply(function(xx, yy) {
            c(xx, yy)
          }, allmcsamplesKA, mcsamplesKA, SIMPLIFY = FALSE)
        }

        rm(mcsamples)
        gc() #garbage collection
        nitertot <- ncol(allmcsamples$W)


        #########################################
        #### CHECK IF CHAIN MUST BE CONTINUED ####
        #########################################
        lengthmeasure <- max(miniter, min(
          maxiter,
          thresholdfn(
            diagnIAT, nsamplesperchain, multcorr,
            diagnBurn2, diagnThin
          )
        ))


        cat('\nNumber of iterations', nitertot, ', required', lengthmeasure, '\n')
        ##
        if (nitertot < lengthmeasure) {
          ## limit number of iterations per loop, to save memory
          niter <- min(lengthmeasure - nitertot + 1L, niterini)
          subiter <- subiter + 1L
          cat(
            '\nChain #', chainnumber, '.', subiter,
            '(chain', achain, 'of', nchainspercore,
            'for this core): increasing by', niter, '\n'
          )
        }
        reset <- FALSE
      }


      #########################################
      #### SAVE CHAIN ####
      #########################################
      
      ## tokeep <- seq(to=nrow(allmcsamples), length.out=nsamplesperchain, by=max(thinning,multcorr*ceiling(max(diagnIAT,diagnThin)), na.rm=TRUE))
      ## allmcsamples <- allmcsamples[tokeep,,drop=FALSE]
      ## ##
      ## saveRDS(allmcsamples, file=paste0(dirname,'_mcsamples',dashnameroot,'--', padchainnumber,'.rds'))
      ## ## rm(allmcsamples)

      cat('\nKeeping last', nsamplesperchain, 'samples with thinning',
          max(thinning, multcorr * ceiling(max(diagnIAT, diagnThin)),
              na.rm = TRUE), '\n')

      tokeep <- seq(to = ncol(allmcsamples$W), length.out = nsamplesperchain,
                    by = max(thinning,
                             multcorr * ceiling(max(diagnIAT, diagnThin)),
                             na.rm = TRUE))

      saveRDS(mcsubset(allmcsamples, tokeep),
              file = paste0(dirname, '_mcsamples',
                            dashnameroot, '--',
                            padchainnumber, '.rds'))
      rm(allmcsamples)
      ## nitertot <- ncol(allmcsamples$W)

      gc() #garbage collection

      saveRDS(traces[tokeep, ],
              file = paste0(dirname, '_mctraces',
                            dashnameroot, '--',
                            padchainnumber, '.rds'))


      ###############
      #### PLOTS ####
      ###############

      #### Plot Alpha and cluster occupation, if required
      if (showAlphatraces || showKtraces) {
        cat('Plotting cluster and Alpha information.\n')
        pdff(paste0(dirname, '_hyperparams_traces', dashnameroot, '--',
                    padchainnumber, '_', achain, '-', acore),
             apaper = 4)

        if (showKtraces) {
          cat('\nSTATS OCCUPIED CLUSTERS:\n')
          print(summary(allmcsamplesKA$K))
          ##
          tplot(y = allmcsamplesKA$K, ylab = 'occupied clusters',
                xlab = 'iteration', ylim = c(0, nclusters))
          tplot(x = ((-1):nclusters) + 0.5,
                y = tabulate(allmcsamplesKA$K + 1, nbins = nclusters + 1),
                type = 'h', xlab = 'occupied clusters', ylab = NA,
                ylim = c(0, NA))
        }
        if (showAlphatraces) {
          cat('\nSTATS log2(alpha):\n')
          print(summary(allmcsamplesKA$Alpha + minalpha - 1L, na.rm = TRUE))
          tplot(y = allmcsamplesKA$Alpha + minalpha - 1L,
                ylab = bquote(log2(alpha)), xlab = 'iteration',
                ylim = c(minalpha, maxalpha))
          tplot(x = (minalpha:(maxalpha + 1)) - 0.5,
                y = tabulate(allmcsamplesKA$Alpha, nbin = nalpha),
                type = 'h', xlab = bquote(log2(alpha)), ylab = '',
                ylim = c(0, NA))
        }
        dev.off()
      }

      #### Plot diagnostic traces of current chain
      if (plottraces) {
        cat('\nPlotting traces and samples.\n')

        tracegroups <- as.list(seq_len(ncol(traces)))
        names(tracegroups) <- colnames(traces)
        grouplegends <- foreach::foreach(agroup = seq_along(tracegroups)) %do% {
          c(
            paste0('-- STATS ', names(tracegroups)[agroup], ' --'),
            paste0('min ESS = ',
                   signif(min(diagnESS[tracegroups[[agroup]]]), 6)),
            paste0('max IAT = ',
                   signif(max(diagnIAT[tracegroups[[agroup]]]), 6)),
            paste0('max BMK = ',
                   signif(max(diagnBMK[tracegroups[[agroup]]]), 6)),
            paste0('max MCSE = ',
                   signif(max(diagnMCSE[tracegroups[[agroup]]]), 6)),
            paste0('stationary: ',
                   sum(diagnStat[tracegroups[[agroup]]]), '/',
                   length(diagnStat[tracegroups[[agroup]]])),
            paste0('burn: ', signif(diagnBurn2, 6)),
            paste0('max thin = ',
                   signif(max(diagnThin[tracegroups[[agroup]]]), 6))
          )
        }

        ## Plot various info and traces
        colpalette <- seq_len(ncol(traces))
        names(colpalette) <- colnames(traces)
        cat('\nPlotting MCMC traces')
        graphics.off()
        pdff(paste0(dirname, '_mcpartialtraces', dashnameroot, '--',
                    padchainnumber, '_', achain, '-', acore), apaper = 4)

        ## Summary stats
        matplot(1:2, type = 'l', col = 'white',
                main = paste0('Stats chain ', achain),
                axes = FALSE, ann = FALSE)
        # Legends
        legendpositions <- c('topleft', 'topright', 'bottomleft', 'bottomright')
        for (alegend in seq_along(grouplegends)) {
          legend(
            x = legendpositions[alegend], bty = 'n', cex = 1.5,
            legend = grouplegends[[alegend]]
          )
        }
        legend(
          x = 'center', bty = 'n', cex = 1,
          legend = c(
            paste0('Chain ', chainnumber, '_', achain, '-', acore),
            paste0('Occupied clusters: ', usedclusters, ' of ', nclusters),
            paste0('LL:  ( ', signif(mean(traces[, 1]), 3), ' +- ',
                   signif(sd(traces[, 1]), 3), ' ) dHart'),
            'NOTES:',
            if (flagmc) {
              'some non-finite MC outputs'
            },
            if (usedclusters > nclusters - 5) {
              'too many clusters occupied'
            },
            if (flagll) {
              'non-finite values in diagnostics'
            }
          )
        )
        ## Traces of likelihood and cond. probabilities
        par(mfrow = c(1, 1))
        for (avar in colnames(traces)) {
          tplot(
            y = traces[, avar], type = 'l', lty = 1, col = colpalette[avar],
            main = paste0(
              'ESS = ', signif(diagnESS[avar], 3),
              ' | IAT = ', signif(diagnIAT[avar], 3),
              ' | BMK = ', signif(diagnBMK[avar], 3),
              ' | MCSE = ', signif(diagnMCSE[avar], 3),
              ' | stat: ', diagnStat[avar] * 1L,
              ' | burnI: ', diagnBurn[avar],
              ' | burnII: ', diagnBurn2,
              ' | thin: ', diagnThin[avar]
            ),
            ylab = paste0(avar, '/dHart'),
            xlab = 'Monte Carlo sample',
            family = family, mar = c(NA, 6, NA, NA)
          )
        }
        dev.off()
      }

      ## Can this be deleted? AURORA
      ## #### Plot samples from current chain
      ##             if(plotpartialsamples){
      ##                 ## Samples of marginal frequency distributions
      ##                 ## if(!continue){
      ##                 subsamples <- (if(totsamples=='all'){1:nitertot}else{round(seq(1, nitertot, length.out=totsamples))})
      ##                 ## showsubsample <- round(seq(1, length(subsamples), length.out=showsamples))
      ##                 ##
      ##                 cat('\nPlotting samples of frequency distributions')
      ##                 plotFsamples(file=paste0(dirname,'mcmcdistributions',dashnameroot,'--', padchainnumber,'_',achain,'-',acore),
      ##                              mcsamples=mcsubset(allmcsamples,subsamples),
      ##                              auxmetadata=auxmetadata,
      ##                              data=data,
      ##                              plotuncertainty='samples',
      ##                              uncertainty=showsamples,
      ##                              plotmeans=plotmeans,
      ##                              datahistogram=TRUE, datascatter=TRUE,
      ##                              useOquantiles=useOquantiles,
      ##                              parallel=FALSE, silent=TRUE
      ##                              )
      ##          }

      cat('\nMCMC + diagnostics time', printtime(Sys.time() - starttime), '\n')

      #### Print estimated remaining time
      er_time <- (Sys.time() - starttime) / achain * (nchainspercore - achain + 1)
      if (is.finite(er_time) && er_time > 0) {
        printnull(
          paste0(
            '\rSampling. Core ', acore, ' estimated remaining time: ',
            printtime(er_time),
            '                      '
          ),
          outcon
        )
      }

      maxusedclusters <- max(maxusedclusters, usedclusters)
    } #### END LOOP OVER CHAINS (WITHIN ONE CORE)

    ##
    cat('\nTotal time', printtime(Sys.time() - starttime), '\n')

    cbind(maxusedclusters, allflagmc)
  }
  ############################################################
  #### END OF FOREACH-LOOP OVER CORES
  ############################################################
  suppressWarnings(sink())
  suppressWarnings(sink(NULL, type = "message"))

  maxusedclusters <- max(chaininfo[, 1])
  nonfinitechains <- sum(chaininfo[, 2])

  ############################################################
  #### End of all MCMC
  ############################################################


  ############################################################
  #### Join chains
  ############################################################
  ## Save random seeds used in the parallel processing
  if (!is.null(attr(chaininfo, "rng"))) { # parallel processing
    saveRDS(attr(chaininfo, "rng"), file = paste0(dirname, "rng_parallelseeds", dashnameroot, ".rds"))
  }

  ## Read the samples saved by each chain and concatenate them
  joinmc <- function(mc1, mc2) {
    mapply(
      function(xx, yy) {
        temp <- c(xx, yy)
        dx <- dim(xx)[-length(dim(xx))]
        dim(temp) <- c(dx, length(temp) / prod(dx))
        temp
      },
      mc1, mc2,
      SIMPLIFY = FALSE
    )
  }
  mcsamples <- foreach(chainnumber = 1:(ncores * nchainspercore), .combine = joinmc, .multicombine = F) %do% {
    padchainnumber <- sprintf(paste0("%0", nchar(nchains), "i"), chainnumber)

    readRDS(file = paste0(dirname, "_mcsamples", dashnameroot, "--", padchainnumber, ".rds"))
  }

  #### Save all final parameters together with the aux-metadata in one file
  saveRDS(c(mcsamples, list(auxmetadata = auxmetadata)),
    file = paste0(dirname, "Fdistribution", dashnameroot, ".rds")
  )

  traces <- foreach(chainnumber = 1:(ncores * nchainspercore), .combine = rbind) %do% {
    padchainnumber <- sprintf(paste0("%0", nchar(nchains), "i"), chainnumber)

    readRDS(file = paste0(dirname, "_mctraces", dashnameroot, "--", padchainnumber, ".rds"))
  }
  ## traces <- mcsamples[round(seq(1,nrow(mcsamples),length.out=nsamples)),1:3]
  saveRDS(traces, file = paste0(dirname, "MCtraces", dashnameroot, ".rds"))

  cat("\rFinished Monte Carlo sampling.                                         \n")
  gc()

  ############################################################
  #### Final joint diagnostics
  ############################################################
  ## cat('\nSome diagnostics:\n')
  traces <- traces[apply(traces, 1, function(x) {
    all(is.finite(x))
  }), ]
  ## flagll <- nrow(traces) != nrow(traces)

  ## funMCSE <- function(x){LaplacesDemon::MCSE(x, method='batch.means')$se}
  diagnESS <- LaplacesDemon::ESS(traces)
  cat("\nEffective sample size:", round(max(diagnESS, na.rm = T)), "\n")
  ## diagnIAT <- apply(traces, 2, function(x){LaplacesDemon::IAT(x)})
  ## cat('\nIATs:',paste0(round(diagnIAT), collapse=', '))
  ## diagnBMK <- LaplacesDemon::BMK.Diagnostic(traces[1:(4*trunc(nrow(traces)/4)),], batches=4)[,1]
  ## cat('\nBMKs:',paste0(round(diagnBMK,3), collapse=', '))
  ## diagnMCSE <- 100*apply(traces, 2, function(x){funMCSE(x)/sd(x)})
  ## cat('\nMCSEs:',paste0(round(diagnMCSE,2), collapse=', '))
  ## diagnStat <- apply(traces, 2, function(x){LaplacesDemon::is.stationary(as.matrix(x,ncol=1))})
  ## cat('\nStationary:',paste0(diagnStat, collapse=', '))
  ## diagnBurn <- apply(traces, 2, function(x){LaplacesDemon::burnin(matrix(x[1:(10*trunc(length(x)/10))], ncol=1))})
  ## cat('\nBurn-in I:',paste0(diagnBurn, collapse=', '))
  ## diagnBurn2 <- proposeburnin(traces, batches=10)
  ## cat('\nBurn-in II:',diagnBurn2)
  ## diagnThin <- proposethinning(traces)
  ## ## cat(paste0('\nProposed thinning: ',paste0(diagnThin, collapse=', ')))
  ## ##
  ## cat('\n')
  ## ##       )
  ## ##         }
  ## ##

  ## Plot various info and traces
  cat("\nPlotting final Monte Carlo traces.\n")

  ##
  colpalette <- 1:ncol(traces)
  names(colpalette) <- colnames(traces)
  graphics.off()
  pdff(paste0(dirname, "MCtraces", dashnameroot), apaper = 4)
  ## Traces of likelihood and cond. probabilities
  for (avar in colnames(traces)) {
    tplot(
      y = traces[, avar], type = "l", lty = 1, col = colpalette[avar],
      main = paste0(
        "Effective sample size: ", signif(diagnESS[avar], 3)
        ##             ' | IAT = ', signif(diagnIAT[avar], 3),
        ##             ' | BMK = ', signif(diagnBMK[avar], 3),
        ##             ' | MCSE = ', signif(diagnMCSE[avar], 3),
        ##             ' | stat: ', diagnStat[avar],
        ##             ' | burn I: ', diagnBurn[avar],
        ##             ' | burn II: ', diagnBurn2
      ),
      ylab = paste0(avar, "/dHart"), xlab = "sample", family = family, mar = c(NA, 6, NA, NA)
    )
  }

  cat("\nMax number of occupied clusters:", maxusedclusters, "\n")
  if (maxusedclusters > nclusters - 5) {
    cat("Too many clusters occupied\n")
  }

  if (nonfinitechains > 0) {
    cat(nonfinitechains, "chains with some non-finite outputs\n")
  }

  cat("Plotting marginal samples.\n")
  plotFsamples(
    file = paste0(dirname, "plotsamples_Fdistribution", dashnameroot),
    mcoutput = c(mcsamples, list(auxmetadata = auxmetadata)),
    data = data,
    plotuncertainty = "samples",
    uncertainty = showsamples, plotmeans = TRUE,
    datahistogram = TRUE, datascatter = TRUE,
    useOquantiles = useOquantiles,
    parallel = TRUE, silent = TRUE
  )

  cat("Plotting marginal samples with quantiles.\n")
  plotFsamples(
    file = paste0(dirname, "plotquantiles_Fdistribution", dashnameroot),
    mcoutput = c(mcsamples, list(auxmetadata = auxmetadata)),
    data = data,
    plotuncertainty = "quantiles",
    uncertainty = showquantiles, plotmeans = TRUE,
    datahistogram = TRUE, datascatter = TRUE,
    useOquantiles = useOquantiles,
    parallel = TRUE, silent = TRUE
  )

  cat("\nTotal computation time", printtime(Sys.time() - timestart0), "\n")

  if (exists("cl")) {
    cat("\nClosing connections to cores.\n")
    stopCluster(cl)
  }


  #### remove partial files if required
  if (cleanup) {
    cat("Removing temporary output files.\n")
    file.remove(dir(dirname,
      pattern = paste0("^_mcsamples", dashnameroot, "--.*\\.rds$"),
      full.names = TRUE
    ))
    file.remove(dir(dirname,
      pattern = paste0("^_mctraces", dashnameroot, "--.*\\.rds$"),
      full.names = TRUE
    ))
    ## Should we leave the plots of partial traces?
    ## maybe create an additional function argument?
    file.remove(dir(dirname,
      pattern = paste0("^_mcpartialtraces", dashnameroot, "--.*\\.pdf$"),
      full.names = TRUE
    ))
  }

  #### WHAT MIGHT BE ADDED:
  ## - when 'showKtraces' is true:
  ## a histogram over number of clusters over all chains
  ## (for the moment there's one plot per chain)

  cat("Finished.\n\n")

  if (output) {
    mcsamples
  }
}

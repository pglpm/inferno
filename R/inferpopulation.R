#' Monte Carlo computation of posterior distribution of population frequencies
#'
#' @param data data.table object or filepath: datapoints
#' @param metadata Either the name of the csv file containing metadata
#'   of the current dataset, or a data.frame with the metadata
#' @param outputdir String, path to output file folder ## Rename to
#'   outputPrefix, also addSuffix?
#' @param nsamples Integer, nr of desired MC samples
#' @param nchains Integer, nr of MC chains
#' @param nsamplesperchain Integer, nr of MC samples per chain
#' @param parallel, Bool or numeric: whether to use pre-existing parallel
#'   workers, or how many to create and use
#' @param niterini Number of initial (burn-in) MC iterations
#' @param miniter Minimum number of MC iterations to be done
#' @param maxiter Maximum number of MC iterations
#' @param prior Bool: Calculate the prior distribution of F?
#' @param thinning If NULL, let the diagnostics decide the MC thinning; if positive, use this thinning value
#' @param plottraces Bool: plot MC traces of diagnostic values
#' @param showKtraces Bool, when true, it saves the K parameter during
#'   sampling and plots its trace and histogram at the end. Keeping it to
#'   FALSE (default) saves a little computation time.
#' @param showAlphatraces Bool, : when true, it saves the Alpha parameter
#'   more frequently during sampling and plots its trace and histogram at
#'   the end. Keeping it to FALSE (default) saves a little computation
#'   time.
#' @param seed Integer: seed for random number generator. If left as default
#'   NULL, a random seed based on the system clock is used in the
#'   set.seed() function
#' @param lldata Positive integer or Inf (default 12): number of datapoints
#'   to use for loglikelihood calculations; if Inf, use all.
#' @param subsampledata Numeric: use only a subsample of the datapoints in 'data'
#' @param useLquantiles Bool: internal, use metadata quantiles for ordinal variates
#' @param output if string 'directory', return the output directory name;
#'   if string 'mcoutput', return the 'Fdistribution' object;
#'   anything else, no output
#' @param cleanup Bool, default TRUE, removes files that can be used for
#'   debugging
#' @return name of directory containing output files, or Fdistribution object, or empty
#' @export
#' @import parallel foreach doParallel doRNG nimble
#' @importFrom LaplacesDemon MCSE ESS IAT BMK.Diagnostic is.stationary burnin
inferpopulation <- function(data, metadata, outputdir, nsamples = 1200,
                            nchains = 120, nsamplesperchain = 10, parallel = TRUE,
                            seed = NULL, cleanup = TRUE,
                            appendtimestamp = TRUE, appendinfo = TRUE,
                            subsampledata = NULL, output = 'directory',
                            niterini = 1200, miniter = 3600, maxiter = +Inf,
                            prior = missing(data),
                            thinning = NULL, plottraces = TRUE,
                            showKtraces = FALSE, showAlphatraces = FALSE,
                            lldata = 12, useLquantiles = FALSE) {

   cat('\n') # make sure possible error messages start on new line

  ## Set the RNG seed if given by user, or if no seed already exists
  if (!missing(seed) || !exists('.Random.seed')) {
    set.seed(seed)
  }
  currentseed <- .Random.seed



  ##################################################
  #### Argument-consistency checks
  ##################################################

  #### Determine the status of parallel processing
  if (is.logical(parallel) && parallel) {
    if (foreach::getDoParRegistered()) {
      cat('Using already registered', foreach::getDoParName(),
          'with', foreach::getDoParWorkers(), 'workers\n')
      ncores <- foreach::getDoParWorkers()
    } else {
      cat('No parallel backend registered.\n')
      ncores <- 1
    }
  } else if (is.numeric(parallel) && parallel >= 2) {
    if (foreach::getDoParRegistered()) {
      ncores <- min(foreach::getDoParWorkers(), parallel)
      cat('Using already registered', foreach::getDoParName(),
          'with', foreach::getDoParWorkers(), 'workers\n')
    } else {
      ## ##
      ## ## Alternative way to register cores;
      ## ## might need to be used for portability to Windows?
      ## registerDoSEQ()
      ## cl <- parallel::makePSOCKcluster(ncores)
      ## ##
      cl <- parallel::makeCluster(parallel)
      doParallel::registerDoParallel(cl)
      cat('Registered', foreach::getDoParName(),
          'with', foreach::getDoParWorkers(), 'workers\n')
      ncores <- parallel
    }
  } else {
    cat('No parallel backend registered.\n')
    ncores <- 1
  }

#### Consistency checks for numbers of samples, chains, cores
  ## The defaults are 1200 samples from 120 chains, so 10 samples per chain
  ## If the user changes any of these three,
  ## the user must also take responsibility of changing one of the other two
  ## in an appropriate way

  if (!missing(nsamples) && !missing(nchains) &&
      missing(nsamplesperchain)) {
    ## nsamples and nchains given
    ## Doesn't make sense to have more cores than chains
    if(nchains < ncores) {
      ncores <- nchains
    }
    ## Ineffective is some core has fewer chains than others
    nchainspercore <- ceiling(nchains / ncores)
    if (nchainspercore * ncores > nchains) {
      nchains <- nchainspercore * ncores
      cat('Increasing number of chains to', nchains, '\n')
    }
    ## Ineffective to have chains with different required samples
    nsamplesperchain <- ceiling(nsamples / nchains)
    if(nsamplesperchain * nchains > nsamples) {
      nsamples <- nchains * nsamplesperchain
      cat('Increasing number of samples to', nsamples, '\n')
    }
  } else if (!missing(nsamples) && missing(nchains) &&
             !missing(nsamplesperchain)) {
    ## nsamples and nsamplesperchain given
    nchains <- ceiling(nsamples / nsamplesperchain)
    if(nchains < ncores) {
      ncores <- nchains
    }
    nchainspercore <- ceiling(nchains / ncores)
    if(nchainspercore * ncores > nchains) {
      nchains <- nchainspercore*ncores
    }
    if(nsamplesperchain * nchains > nsamples) {
      nsamples <- nchains * nsamplesperchain
      cat('Increasing number of samples to', nsamples, '\n')
    }
  } else if (missing(nsamples) && !missing(nchains) &&
             !missing(nsamplesperchain)) {
    ## nchains and nsamplesperchain given
    if(nchains < ncores){
      ncores <- nchains
    }
    nchainspercore <- ceiling(nchains / ncores)
    if(nchainspercore * ncores > nchains){
      nchains <- nchainspercore * ncores
      cat('Increasing number of chains to',nchains,'\n')
    }
    nsamples <- nchains * nsamplesperchain
  } else  if ( !(missing(nsamples) && missing(nchains) &&
                 missing(nsamplesperchain)) ) {
    ## The user set all these three arguments or only one
    stop('Please specify exactly two among "nsamples", "nchains",
         "nsamplesperchain"')
  }

  if (is.numeric(thinning) && thinning > 0) {
    thinning <- ceiling(thinning)
  } else if (!is.null(thinning)) {
    stop('Invalid "thinning" argument.')
  }

  ## Parallellisation if more than one core
  ## Done now in case ncores was reduced because of nchains argument
  if (ncores < 2) {
    `%dochains%` <- `%do%`
  } else {
    `%dochains%` <- `%dorng%`
  }

  ## Make sure 'niterini' is at least 2
  niterini <- max(2, niterini)
  #### Start timer
  timestart0 <- Sys.time()

  ##################################################
  #### Read and process data and metadata
  ##################################################


#### Read metadata
  ## Check whether argument 'metadata' is a string for a file name
  ## otherwise we assume it's a data.table or similar object
  if (is.character(metadata) && file.exists(metadata)) {
    metadata <- read.csv(metadata, na.strings = '')
  }
  metadata <- as.data.frame(metadata)

#### Dataset
  ## Check if 'data' is given
  if(!(missing(data) || is.null(data))) {
    ## Check if 'data' argument is an existing file
    ## otherwise we assume is an object
    datafile <- NULL
    if (is.character(data)) {
      ## add '.csv' if missing
      datafile <- paste0(sub('.csv$', '', data), '.csv')
      if (file.exists(data)) {
        data <- read.csv(datafile, na.strings = '')
      } else {
        stop('Cannot find data file')
      }
    }
    data <- as.data.frame(data)

    ## Consistency checks for data
    ## They should be moved to an external function

    ## Are data missing variates?
    if (!all(metadata[['name']] %in% colnames(data))) {
      stop('Missing variates in data. Check data- and metadata-files.')
    }

    ## Drop variates in data that are not in the metadata file
    if (!all(colnames(data) %in% metadata[['name']])) {
      cat('Warning: data have additional variates. Dropping them.\n')
      subvar <- intersect(colnames(data), metadata[['name']])
      data <- data[, subvar, drop = FALSE]
      rm(subvar)
    }

    ## Remove empty datapoints
    tokeep <- which(apply(data, 1, function(xx) { !all(is.na(xx)) }))
    if(length(tokeep) == 0 && !prior) {
      stop('Data are given but empty')
    } else if(length(tokeep) < nrow(data)) {
      cat('Warning: data contain empty datapoints. Dropping them.\n')
      data <- data[tokeep, , drop = FALSE]
    }
    rm(tokeep)

    ## Check if the user wants to use a subset of the dataset
    if (is.numeric(subsampledata)) {
      ## @@TODO: find faster and memory-saving subsetting
      cat('Subsampling data, as requested.\n')
      data <- data[sample(seq_len(nrow(data)),
                          min(subsampledata, nrow(data)),
                          replace = FALSE), ]
    }

    npoints <- nrow(data)

  } else {
    ## data not given: we assume user wants prior calculation
    message('Missing data')
    prior <- TRUE
    npoints <- 0
  }

  ## Build auxiliary metadata object; we'll save it later
  cat('Calculating auxiliary metadata\n')
  auxmetadata <- buildauxmetadata(data = data, metadata = metadata)

  #### Output-folder setup
  if (missing(outputdir) || (is.logical(outputdir) && outputdir)) {
    outputdir <- paste0('_output_', sub('.csv$', '', datafile))
  }

  ## append time and info to output directory, if requested
  suffix <- NULL
  if (appendtimestamp) {
    suffix <- paste0(suffix, '-',
                     strftime(as.POSIXlt(Sys.time()), '%y%m%dT%H%M%S') )
  }
  if (appendinfo) {
    suffix <- paste0(suffix,
                     '-vrt', nrow(auxmetadata),
                     '_dat',
                     (if (npoints == 1 && all(is.na(data))) {
                        0
                      } else {
                        npoints
                      }),
                     ## '-K', nclusters, # unimportant for user
                     '_smp', nsamples)
  }
  dirname <- paste0(outputdir, suffix)
  ##
  # Create output directory if it does not exist
  dir.create(dirname, showWarnings = FALSE)
  # Print information
  cat('\n', paste0(rep('*', max(nchar(dirname), 26)), collapse = ''),
    '\n Saving output in directory\n', dirname, '\n',
    paste0(rep('*', max(nchar(dirname), 26)), collapse = ''), '\n')

  ## This is in case we need to add some extra specifier to the output files
  ## all 'dashnameroot' can be deleted in a final version
  dashnameroot <- NULL

  ## Save copy of metadata in directory
  write.csv(metadata, file = file.path(dirname, 'metadata.csv'),
            row.names = FALSE, quote = FALSE, na = '')

  ## Save initial RNG seed in case needed by user
  saveRDS(currentseed,
          file = file.path(dirname, paste0('rng_seed', dashnameroot, '.rds')))


#### Check if user wants to calculate prior
  if (prior) {
    message('CALCULATING PRIOR DISTRIBUTION')
    ## if data is not empty, we use it to create data for likelihood
    if(!is.null(data)) {
      ## if "testdata" is moved into the for-loop,
      ## then each chain uses a different set of testdata
      for(achain in seq_len(nchains)) {
        testdata <- data[sort(sample(npoints, min(lldata, npoints))), ,
                         drop = FALSE]
        saveRDS(testdata,
                file = file.path(dirname, paste0('_testdata_', achain, '.rds')))
      }
    } else {
      ## no data available: construct one datapoint from the metadata info
      testdata <- data.frame(1:3)[,-1]
      for(xx in seq_len(nrow(auxmetadata))) {
        xx <- as.list(auxmetadata[xx, ])
        toadd <- unlist(xx[paste0('mctest', 1:3)])
        if (xx[['mcmctype']] %in% c('B', 'N', 'O')) {
          toadd <- unlist(xx[paste0('V', toadd)])
        }
        testdata <- cbind(testdata, toadd)
      }
      colnames(testdata) <- auxmetadata[['name']]
      for(achain in seq_len(nchains)) {
        saveRDS(testdata,
                file = file.path(dirname, paste0('_testdata_', achain, '.rds')))
      }
    }
    ## create empty dataset: Monte Carlo sampling is non-Markov
    data <- as.data.frame(
      matrix(NA, nrow = 1, ncol = nrow(metadata),
             dimnames = list(NULL, metadata$name))
    )
  } else {
    ## create data for likelihood
    ## if "testdata" is moved into the for-loop,
    ## then each chain uses a different set of testdata
    for(achain in seq_len(nchains)) {
      testdata <- data[sort(sample(npoints, min(lldata, npoints))), ,
                       drop = FALSE]
      saveRDS(testdata,
              file = file.path(dirname, paste0('_testdata_', achain, '.rds')))
    }
  }
  rm(testdata)

## #### Select and save data for loglikelihood
##   ## Find which datapoints have not all missing entries
##   dataNoNa <- which(apply(data, 1, function(xx) { !all(is.na(xx)) }))
##   ndataNoNa <- length(dataNoNa)
## 
##   if (ndataNoNa > 0) {
##     ## if "testdata" is moved into the for-loop,
##     ## then each chain uses a different set of testdata
##     for(achain in seq_len(nchains)) {
##       testdata <- data[sort(sample(dataNoNa, min(lldata, ndataNoNa))),]
##       saveRDS(testdata,
##               file = file.path(dirname, paste0('_testdata_', achain, '.rds')))
##     }
##   } else {
##       ## no data available: construct one datapoint from the metadata info
##       testdata <- data.table(sapply(seq_len(nrow(auxmetadata)), function(xx) {
##       # consider making this function separate
##       xx <- auxmetadata[xx, ]
##       toadd <- xx[, paste0('mctest', 1:3), with = FALSE]
##       if (xx[['mcmctype']] %in% c('B', 'N', 'O)) {
##         toadd <- xx[1, paste0('V', toadd), with = FALSE]
##       }
##       toadd
##     }))
##     colnames(testdata) <- auxmetadata[['name']]
##     for(achain in seq_len(nchains)) {
##       saveRDS(testdata,
##           file = file.path(dirname, paste0('_testdata_', achain, '.rds')))
##     }
##     rm(testdata)
##     }



  ##################################################
  #### Various internal parameters
  ##################################################

#### Hyperparameters and other internal parameters
  ## source('hyperparameters.R') doesn't seem to work
  nclusters <- 64L
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
  Lshapelo <- 0.5
  Lshapehi <- 0.5
  Lvarm1 <- 3L^2L
  Bshapelo <- 1L
  Bshapehi <- 1L

  nalpha <- length(minalpha:maxalpha)
  npoints <- nrow(data)

  #### Other options
  Alphatoslice <- TRUE # FALSE typically leads to underflow
  Ktoslice <- TRUE # FALSE typically leads to underflow
  RWtoslice <- FALSE
  changeSamplerOrder <- TRUE
  ##
  ## plotmeans <- TRUE # plot frequency averages
  showsamples <- 100 # number of samples to show.
  plotDisplayedQuantiles <- c(1, 31) / 32 # quantiles to show
  nclustersamples <- 128 # number of samples of Alpha and K
  showsamplertimes <- FALSE ##
  family <- 'Palatino' # font family in plots


  ##################################################
  #### Define functions
  ##################################################

  ## Load auxiliary functions from external files
  ## source('samplesFDistribution.R')
  ## source('plotFsamples.R')
  ## source('tplotfunctions.R')
  ## source('vtransform.R')
  ## ## These are used within the foreach core loop
  ## source('proposeburnin.R')
  ## source('proposethinning.R')
  ## source('mcsubset.R')

  # function printtime to format printing of time
  printtime <- function(tim) {
    paste0(signif(tim, 2), ' ', attr(tim, 'units'))
  }
  ## We need to send some messages to the log files, others to the user.
  ## This is done by changing output sink:
  print2user <- function(msg, outcon) {
    sink(NULL, type = 'message')
    message(msg, appendLF = FALSE)
    flush.console()
    sink(outcon, type = 'message')
  }


#####################################################
#### Prepare arguments for Nimble model
#####################################################

  ## R: continuous open domain
  ## C: continuous closed domain (or censored)
  ## D: continuous rounded
  ## L: latent
  ## O: ordinal
  ## N: nominal
  ## B: binary
  ## number and names of variates of each type
  vn <- vnames <- list(R=NULL, C=NULL, D=NULL, L=NULL, O=NULL, N=NULL, B=NULL)

  for (atype in names(vn)) {
    vnames[[atype]] <- auxmetadata[auxmetadata$mcmctype == atype, 'name']
    vn[[atype]] <- length(vnames[[atype]])
  }

  ## ## REMOVE soon: previous version
  ## vn <- list() # How many variates of each type
  ## vnames <- list() # The names for variates of each type
  ## for (atype in c('R', 'C', 'D', 'L', 'O', 'N', 'B')) {
  ##   vnames[[atype]] <- auxmetadata[auxmetadata$mcmctype == atype, 'name']
  ##   vn[[atype]] <- length(vnames[[atype]])
  ## }

  #### CONSTANTS OF NIMBLE MODEL
  ## These constants are available in the Nimble environment
  ## They don't have to be accessed by constants$varname
  ## vn$R + vn$C + vn$D + vn$L
  constants <- c(
    list(
      nclusters = nclusters,
      npoints = npoints,
      nalpha = nalpha,
      probalpha0 = rep(1 / nalpha, nalpha),
      basealphas = rep((2^(minalpha - 1L #- 1L + vn$R + vn$C + vn$D + vn$L
      )) / nclusters, nclusters)
    ),
    if (vn$R > 0) { # continuous open domain
      list(
        Rn = vn$R, # This indexing variable is needed internally
        Rmean1 = rep(0, 1),
        Rvarm1 = rep(Rvarm1, 1),
        Rvar1 = rep(1, 1),
        Rshapelo = rep(Rshapelo, 1),
        Rshapehi = rep(Rshapehi, 1)
      )
    },
    if (vn$C > 0) { # continuous closed domain
      list(
        Cn = vn$C, # This indexing variable is needed internally
        Cmean1 = rep(0, 1),
        Cvarm1 = rep(Cvarm1, 1),
        Cvar1 = rep(1, 1),
        Cshapelo = rep(Cshapelo, 1),
        Cshapehi = rep(Cshapehi, 1),
        Cleft = vtransform(data[, vnames$C, drop = FALSE],
                           auxmetadata = auxmetadata,
                           Cout = 'left',
                           useLquantiles = useLquantiles),
        Cright = vtransform(data[, vnames$C, drop = FALSE],
                            auxmetadata = auxmetadata,
                            Cout = 'right',
                            useLquantiles = useLquantiles),
        Clatinit = vtransform(data[, vnames$C, drop = FALSE],
                              auxmetadata,
                              Cout = 'init',
                          useLquantiles = useLquantiles
                          )
      )
      ## Cleft & Cright are as many as the datapoints
      ## so we do not create copies outside of Nimble
      ## to save RAM
    },
    if (vn$D > 0) { # continuous rounded
      list(
        Dn = vn$D, # This indexing variable is needed internally
        Dmean1 = rep(0, 1),
        Dvarm1 = rep(Dvarm1, 1),
        Dvar1 = rep(1, 1),
        Dshapelo = rep(Dshapelo, 1),
        Dshapehi = rep(Dshapehi, 1),
        Dleft = vtransform(data[, vnames$D, drop = FALSE],
                           auxmetadata = auxmetadata,
                           Dout = 'left',
                           useLquantiles = useLquantiles),
        Dright = vtransform(data[, vnames$D, drop = FALSE],
                            auxmetadata = auxmetadata,
                            Dout = 'right',
                            useLquantiles = useLquantiles),
        Dlatinit = vtransform(data[, vnames$D, drop = FALSE],
                              auxmetadata,
                              Dout = 'init',
                              useLquantiles = useLquantiles)
      )
    },
    if (vn$L > 0) { # latent
      list(
        Ln = vn$L, # This indexing variable is needed internally
        Lmean1 = rep(0, 1),
        Lvarm1 = rep(Lvarm1, 1),
        Lvar1 = rep(1, 1),
        Lshapelo = rep(Lshapelo, 1),
        Lshapehi = rep(Lshapehi, 1),
        Lleft = vtransform(data[, vnames$L, drop = FALSE],
                           auxmetadata = auxmetadata,
                           Lout = 'left',
                           useLquantiles = useLquantiles),
        Lright = vtransform(data[, vnames$L, drop = FALSE],
                            auxmetadata = auxmetadata,
                            Lout = 'right',
                            useLquantiles = useLquantiles),
        Llatinit = vtransform(data[, vnames$L, drop = FALSE],
                              auxmetadata,
                              Lout = 'init',
                              useLquantiles = useLquantiles)
      )
    },
    if (vn$O > 0) { # ordinal
      Omaxn <- max(auxmetadata[auxmetadata$mcmctype == 'O', 'Nvalues'])
      Oalpha0 <- matrix(1e-100, nrow = vn$O, ncol = Omaxn)
      for (avar in seq_along(vnames$O)) {
        nvalues <- auxmetadata[auxmetadata$name == vnames$O[avar], 'Nvalues']
        ## we choose a flatter hyperprior for ordinal variates
          Oalpha0[avar, 1:nvalues] <- 1
      }
      ##
      list(
        On = vn$O, # This indexing variable is needed internally
        Omaxn = Omaxn,
        Oalpha0 = Oalpha0
      )
    },
    if (vn$N > 0) { # nominal
      Nmaxn <- max(auxmetadata[auxmetadata$mcmctype == 'N', 'Nvalues'])
      Nalpha0 <- matrix(1e-100, nrow = vn$N, ncol = Nmaxn)
      for (avar in seq_along(vnames$N)) {
        nvalues <- auxmetadata[auxmetadata$name == vnames$N[avar], 'Nvalues']
        ## we choose a Hadamard-like hyperprior for nominal variates
          Nalpha0[avar, 1:nvalues] <- 1 / nvalues
      }
      ##
      list(
        Nn = vn$N, # This indexing variable is needed internally
        Nmaxn = Nmaxn,
        Nalpha0 = Nalpha0
      )
    },
    if (vn$B > 0) { # binary
      list(
        Bn = vn$B, # This indexing variable is needed internally
        Bshapelo = rep(Bshapelo, 1),
        Bshapehi = rep(Bshapehi, 1)
      )
    }
  ) # End constants

 #### DATAPOINTS
  datapoints <- c(
    if (vn$R > 0) { # continuous open domain
      list(
        Rdata = vtransform(data[, vnames$R, drop = FALSE],
                           auxmetadata = auxmetadata,
                           Rout = 'normalized'
                           useLquantiles = useLquantiles
                           )
      )
    },
    if (vn$C > 0) { # continuous closed domain
      list(
        Caux = vtransform(data[, vnames$C, drop = FALSE],
                          auxmetadata = auxmetadata,
                          Cout = 'aux',
                          useLquantiles = useLquantiles),
        Clat = vtransform(data[, vnames$C, drop = FALSE],
                          auxmetadata = auxmetadata,
                          Cout = 'lat',
                          useLquantiles = useLquantiles)
      )
    },
    if (vn$D > 0) { # continuous rounded
      list(
        Daux = vtransform(data[, vnames$D, drop = FALSE],
                          auxmetadata = auxmetadata,
                          Dout = 'aux',
                          useLquantiles = useLquantiles)
      )
    },
    if (vn$L > 0) { # latent
      list(
        Laux = vtransform(data[, vnames$L, drop = FALSE],
                          auxmetadata = auxmetadata,
                          Lout = 'aux',
                          useLquantiles = useLquantiles)
      )
    },
    if (vn$O > 0) { # nominal
      list(
        Odata = vtransform(data[, vnames$O, drop = FALSE],
                           auxmetadata = auxmetadata,
                           Oout = 'numeric',
                           useLquantiles = useLquantiles)
      )
    },
    if (vn$N > 0) { # nominal
      list(
        Ndata = vtransform(data[, vnames$N, drop = FALSE],
                           auxmetadata = auxmetadata,
                           Nout = 'numeric',
                           useLquantiles = useLquantiles)
      )
    },
    if (vn$B > 0) { # binary
      list(
        Bdata = vtransform(data[, vnames$B, drop = FALSE],
                           auxmetadata = auxmetadata,
                           Bout = 'numeric',
                           useLquantiles = useLquantiles)
      )
    }
  ) # End datapoints

  #### Output information to user
  if (!exists('Nalpha0')) {
    Nalpha0 <- cbind(1)
  }
  if (!exists('Oalpha0')) {
    Oalpha0 <- cbind(1)
  }
  cat(
    'Starting Monte Carlo sampling of', nsamples, 'samples by',
    nchains, 'chains'
  )
  cat(
    '\nin a space of',
    (sum(as.numeric(vn) * c(2, 2, 2, 2, 0, 0, 1)) +
     sum(Nalpha0 > 2e-100) - nrow(Nalpha0) + 1 +
     sum(Oalpha0 > 2e-100) - nrow(Oalpha0) + 1 ) * nclusters - 1,
    '(effectively',
    paste0(
      (sum(as.numeric(vn) * c(
        3 + npoints, 3 + npoints, 3 + npoints,
        3 + npoints, 0, 0, 1 + npoints
      )) +
        sum(Nalpha0 > 2e-100) +
      nrow(Nalpha0) * (npoints - 1) + 1 +
        sum(Oalpha0 > 2e-100) +
      nrow(Oalpha0) * (npoints - 1) + 1 ) * nclusters - 1 + nalpha - 1,
      ')'
    ), 'dimensions.\n'
  )
  cat(
    'Using', ncores, 'cores:',
    nsamplesperchain, 'samples per chain,', nchainspercore, 'chains per core.\n'
  )
  cat('Core logs are being saved in individual files.\n')
  cat('\nC-compiling samplers appropriate to the variates (package Nimble)\n')
  cat('this can take tens of minutes with many data or variates.\n...\r')


  #####################################################
  #### BEGINNING OF FOREACH LOOP OVER CORES
  #####################################################
  #Iterate over cores, using 'acore' variable as iterator
  chaininfo <- foreach(acore = 1:ncores,
                                .combine = rbind,
                                .inorder = FALSE,
                                ##.packages = c('modelfreeinference'),
                                .noexport = c('data')
                                ##.packages = c('khroma', 'foreach', 'rngtools', 'data.table', 'nimble')
    ) %dochains% {
      ## Create log file
      ## Redirect diagnostics and service messages there
      outcon <- file(file.path(dirname,
        paste0('log', dashnameroot,
        '-', acore, '.log')
      ), open = 'w')
      sink(outcon)
      sink(outcon, type = 'message')

      cat('Log core', acore)
      cat(' - Current time:',
          strftime(as.POSIXlt(Sys.time()), '%Y-%m-%d %H:%M:%S'))
      cat('\n')

      ## suppressPackageStartupMessages(library('data.table'))
      suppressPackageStartupMessages(require('nimble'))
      ## requireNamespace("nimble", quietly = TRUE)
      ##library('nimble')

      ##         cat('\n***TEST00***\n')
      ## mcsamples <- readRDS('mcoutput.rds')
      ## mcsamples$auxmetadata <- NULL
      ## chainnumber <- (acore - 1L) * nchainspercore + achain
      ## testdata <- readRDS(file = file.path(dirname,
      ##                 paste0('_testdata_', chainnumber, '.rds')))
      ## cat('\n***TEST0***\n')
      ## str(testdata)
      ## print(testdata)
      ##   ll <- samplesFDistribution(
      ##       Y = testdata, X = NULL,
      ##       mcoutput = c(mcsamples, list(auxmetadata = auxmetadata)),
      ##       jacobian = FALSE,
      ##       useLquantiles = useLquantiles,
      ##       parallel = FALSE,
      ##       combine = '+',
      ##       silent = TRUE
      ##   )
      ## str(ll)
      ## cat('\n***TEST01***\n')
      ## stop('HERE')


      ## We have to source scripts again for each chain to be able to access them.
      ## source('plotFsamples.R')
      ## source('samplesFDistribution.R')
      ## source('tplotfunctions.R')
      ## source('util_vtransform.R')
      ## source('util_proposeburnin.R')
      ## source('util_proposethinning.R')
      ## source('util_mcsubset.R')
      ## source('util_mcmclength.R')


      ## Function for diagnostics
      ## it corrects a little bug in LaplacesDemon::MCSE
      funMCSE <- function(x) {
        if (length(x) >= 1000) {
          LaplacesDemon::MCSE(x, method = 'batch.means')$se
        } else {
          LaplacesDemon::MCSE(x)
        }
      }

      ## ## Not needed?
      ## printtime <- function(tim){paste0(signif(tim,2),' ',attr(tim,'units'))}
      ## print2user <- function(message, outcon){
      ##     sink(NULL,type='message')
      ##     message(message, appendLF=FALSE)
      ##     flush.console()
      ##     sink(outcon,type='message')
      ## }

    #### CLUSTER REPRESENTATION OF FREQUENCY SPACE

    ## hierarchical probability structure
    finitemix <- nimbleCode({
      ## Component weights
      Alpha ~ dcat(prob = probalpha0[1:nalpha])
      alphas[1:nclusters] <- basealphas[1:nclusters] * 2^Alpha
      W[1:nclusters] ~ ddirch(alpha = alphas[1:nclusters])

      ## Probability density for the parameters of the components
      # Loop over clusters
      for (k in 1:nclusters) {
        ## Probability distributions of parameters
        ## of the different variate types
        if (vn$R > 0) { # continuous open domain
          for (v in 1:Rn) {
            Rmean[v, k] ~ dnorm(mean = Rmean1, var = Rvarm1)
            Rrate[v, k] ~ dinvgamma(shape = Rshapehi, rate = Rvar1)
            Rvar[v, k] ~ dinvgamma(shape = Rshapelo, rate = Rrate[v, k])
          }
        }
        if (vn$C > 0) { # continuous closed domain
          for (v in 1:Cn) {
            Cmean[v, k] ~ dnorm(mean = Cmean1, var = Cvarm1)
            Crate[v, k] ~ dinvgamma(shape = Cshapehi, rate = Cvar1)
            Cvar[v, k] ~ dinvgamma(shape = Cshapelo, rate = Crate[v, k])
          }
        }
        if (vn$D > 0) { # continuous rounded
          for (v in 1:Dn) {
            Dmean[v, k] ~ dnorm(mean = Dmean1, var = Dvarm1)
            Drate[v, k] ~ dinvgamma(shape = Dshapehi, rate = Dvar1)
            Dvar[v, k] ~ dinvgamma(shape = Dshapelo, rate = Drate[v, k])
          }
        }
        if (vn$L > 0) { # latent
          for (v in 1:Ln) {
            Lmean[v, k] ~ dnorm(mean = Lmean1, var = Lvarm1)
            Lrate[v, k] ~ dinvgamma(shape = Lshapehi, rate = Lvar1)
            Lvar[v, k] ~ dinvgamma(shape = Lshapelo, rate = Lrate[v, k])
          }
        }
        if (vn$O > 0) { # ordinal
          for (v in 1:On) {
            Oprob[v, k, 1:Omaxn] ~ ddirch(alpha = Oalpha0[v, 1:Omaxn])
          }
        }
        if (vn$N > 0) { # nominal
          for (v in 1:Nn) {
            Nprob[v, k, 1:Nmaxn] ~ ddirch(alpha = Nalpha0[v, 1:Nmaxn])
          }
        }
        if (vn$B > 0) { # binary
          for (v in 1:Bn) {
            Bprob[v, k] ~ dbeta(shape1 = Bshapelo, shape2 = Bshapehi)
          }
        }
      }
      ## Probability of data
      for (d in 1:npoints) {
        K[d] ~ dcat(prob = W[1:nclusters])
        ##
        if (vn$R > 0) { # continuous open domain
          for (v in 1:Rn) {
            Rdata[d, v] ~ dnorm(mean = Rmean[v, K[d]], var = Rvar[v, K[d]])
          }
        }
        if (vn$C > 0) { # continuous closed domain
          for (v in 1:Cn) {
            Caux[d, v] ~ dconstraint(Clat[d, v] >= Cleft[d, v] &
                                       Clat[d, v] <= Cright[d, v])
            Clat[d, v] ~ dnorm(mean = Cmean[v, K[d]], var = Cvar[v, K[d]])
          }
        }
        if (vn$D > 0) { # continuous rounded
          for (v in 1:Dn) {
            Daux[d, v] ~ dconstraint(Dlat[d, v] >= Dleft[d, v] &
                                       Dlat[d, v] < Dright[d, v])
            Dlat[d, v] ~ dnorm(mean = Dmean[v, K[d]], var = Dvar[v, K[d]])
          }
        }
        if (vn$L > 0) { # latent
          for (v in 1:Ln) {
            Laux[d, v] ~ dconstraint(Llat[d, v] >= Lleft[d, v] &
                                       Llat[d, v] < Lright[d, v])
            Llat[d, v] ~ dnorm(mean = Lmean[v, K[d]], var = Lvar[v, K[d]])
          }
        }
        if (vn$O > 0) { # nominal
          for (v in 1:On) {
            Odata[d, v] ~ dcat(prob = Oprob[v, K[d], 1:Omaxn])
          }
        }
        if (vn$N > 0) { # nominal
          for (v in 1:Nn) {
            Ndata[d, v] ~ dcat(prob = Nprob[v, K[d], 1:Nmaxn])
          }
        }
        if (vn$B > 0) { # binary
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
        ## ## assign all points to an unsystematically chosen cluster
        K = rep(sample(rep(which(W > 0), 2), 1, replace = TRUE), npoints)
        ## ## or:
        ## ## assign all points to the most probable cluster
        ## K = rep(which.max(W), npoints)
        ## ## assign all points to the least probable cluster
        ## K = rep(which(W == min(W[W > 0]))[1], npoints)
        ## ## or:
        ## ## distribute points unsystematically among clusters
        ## K = sample(rep(which(W > 0), 2), npoints, replace = TRUE)
        ## ## or:
        ## ## distribute points unsystematically among M=2 clusters
        ## K = sample(sample(rep(which(W > 0), 2), 2, replace = TRUE),
        ##           npoints, replace = TRUE)
      )
      ##
      if (vn$R > 0) { # continuous open domain
        Rrate <- matrix(
          nimble::rinvgamma(
            n = vn$R * nclusters,
            shape = constants$Rshapehi,
            rate = constants$Rvar1
          ),
          nrow = vn$R, ncol = nclusters
        )
        outlist <- c(
          outlist,
          list(
            Rmean = matrix(
              rnorm(
                n = vn$R * nclusters,
                mean = constants$Rmean1,
                sd = sqrt(constants$Rvarm1)
              ),
              nrow = vn$R, ncol = nclusters
            ),
            Rrate = Rrate,
            Rvar = matrix(
              nimble::rinvgamma(
                n = vn$R * nclusters,
                shape = constants$Rshapelo,
                rate = Rrate
              ),
              nrow = vn$R, ncol = nclusters
            )
          )
        )
      }
      if (vn$C > 0) { # ccontinuous closed domain
        Crate <- matrix(
          nimble::rinvgamma(
            n = vn$C * nclusters,
            shape = constants$Cshapehi,
            rate = constants$Cvar1
          ),
          nrow = vn$C, ncol = nclusters
        )
        outlist <- c(
          outlist,
          list(
            Cmean = matrix(
              rnorm(
                n = vn$C * nclusters,
                mean = constants$Cmean1,
                sd = sqrt(constants$Cvarm1)
              ),
              nrow = vn$C, ncol = nclusters
            ),
            Crate = Crate,
            Cvar = matrix(
              nimble::rinvgamma(
                n = vn$C * nclusters,
                shape = constants$Cshapelo,
                rate = Crate
              ),
              nrow = vn$C, ncol = nclusters
            ),
            ## for data with boundary values
            Clat = constants$Clatinit
            ## Clat = vtransform(data[, vnames$C, with = FALSE],
            ##   auxmetadata, Cout = 'init',
            ##   useLquantiles = useLquantiles)
          )
        )
      }
      if (vn$D > 0) { # continuous rounded
        Drate <- matrix(
          nimble::rinvgamma(
            n = vn$D * nclusters,
            shape = constants$Dshapehi,
            rate = constants$Dvar1
          ),
          nrow = vn$D, ncol = nclusters
        )
        outlist <- c(
          outlist,
          list(
            Dmean = matrix(
              rnorm(
                n = vn$D * nclusters,
                mean = constants$Dmean1,
                sd = sqrt(constants$Dvarm1)
              ),
              nrow = vn$D, ncol = nclusters
            ),
            Drate = Drate,
            Dvar = matrix(
              nimble::rinvgamma(
                n = vn$D * nclusters,
                shape = constants$Dshapelo,
                rate = Drate
              ),
              nrow = vn$D, ncol = nclusters
            ),
            ## for data with boundary values
            Dlat = constants$Dlatinit
            ## Dlat = vtransform(data[, vnames$D, with = FALSE],
            ##   auxmetadata, Dout = 'init',
            ##   useLquantiles = useLquantiles)
          )
        )
      }
      if (vn$L > 0) { # latent
        Lrate <- matrix(
          nimble::rinvgamma(
            n = vn$L * nclusters,
            shape = constants$Lshapehi,
            rate = constants$Lvar1
          ),
          nrow = vn$L, ncol = nclusters
        )
        outlist <- c(
          outlist,
          list(
            Lmean = matrix(
              rnorm(
                n = vn$L * nclusters,
                mean = constants$Lmean1,
                sd = sqrt(constants$Lvarm1)
              ),
              nrow = vn$L, ncol = nclusters
            ),
            Lrate = Lrate,
            Lvar = matrix(
              nimble::rinvgamma(
                n = vn$L * nclusters,
                shape = constants$Lshapelo,
                rate = Lrate
              ),
              nrow = vn$L, ncol = nclusters
            ),
            ## for data with boundary values
            Llat = constants$Llatinit
            ## Llat = vtransform(data[, vnames$L, with = FALSE],
            ##   auxmetadata, Lout = 'init',
            ##   useLquantiles = useLquantiles)
          )
        )
      }
      if (vn$O > 0) { # ordinal
        outlist <- c(
          outlist,
          list(
            Oprob = aperm(array(sapply(1:vn$O, function(avar) {
              sapply(1:nclusters, function(aclus) {
                nimble::rdirch(n = 1, alpha = Oalpha0[avar, ])
              })
            }), dim = c(Omaxn, nclusters, vn$O)))
          )
        )
      }
      if (vn$N > 0) { # nominal
        outlist <- c(
          outlist,
          list(
            Nprob = aperm(array(sapply(1:vn$N, function(avar) {
              sapply(1:nclusters, function(aclus) {
                nimble::rdirch(n = 1, alpha = Nalpha0[avar, ])
              })
            }), dim = c(Nmaxn, nclusters, vn$N)))
          )
        )
      }
      if (vn$B > 0) { # binary
        outlist <- c(
          outlist,
          list(
            Bprob = matrix(
              rbeta(
                n = vn$B * nclusters,
                shape1 = constants$Bshapelo,
                shape2 = constants$Bshapehi
              ),
              nrow = vn$B, ncol = nclusters
            )
          )
        )
      }
      ##
      outlist
    } #End initsfns

    # Timer
    timecount <- Sys.time()


    ##################################################
    #### NIMBLE SETUP
    ##################################################
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
        if (vn$R > 0) {
          c('Rmean', 'Rvar')
        },
        if (vn$C > 0) {
          c('Cmean', 'Cvar')
        },
        if (vn$D > 0) {
          c('Dmean', 'Dvar')
        },
        if (vn$L > 0) {
          c('Lmean', 'Lvar')
        },
        if (vn$O > 0) {
          c('Oprob')
        },
        if (vn$N > 0) {
          c('Nprob')
        },
        if (vn$B > 0) {
          c('Bprob')
        }
      ),
      ## It is necessary to monitor K to see if all clusters were used
      ## if 'showAlphatraces' is true then
      ## the Alpha-parameter trace is also recorded and shown
      monitors2 = c(if (showAlphatraces) { 'Alpha' },
                    'K')
    )
      ## ## Uncomment to debug Nimble (in case of Nimble updates)
      ## print(confnimble$getUnsampledNodes())
      ## confnimble$printSamplers(executionOrder=TRUE)


    targetslist <- sapply(confnimble$getSamplers(), function(xx) xx$target)
    nameslist <- sapply(confnimble$getSamplers(), function(xx) xx$name)
    ## cat('\n******** NAMESLIST', nameslist, '\n')

    ## replace Alpha's cat-sampler with slice
      if (Alphatoslice &&
          !('Alpha' %in% targetslist[nameslist == 'posterior_predictive'])) {
        confnimble$replaceSampler(target='Alpha', type='slice')
        ## ## Old replacement method, didn't work in previous Nimble
        ## confnimble$removeSamplers('Alpha')
        ## confnimble$addSampler(target = 'Alpha', type = 'slice')
    }

    ## replace K's cat-sampler with slice
    if (Ktoslice) {
      for (asampler in grep('^K\\[', targetslist, value = TRUE)) {
        if (!(asampler %in% targetslist[nameslist == 'posterior_predictive'])) {
          confnimble$replaceSampler(target=asampler, type='slice')
          ## ## Old replacement method, didn't work in previous Nimble
          ## confnimble$removeSamplers(asampler)
          ## confnimble$addSampler(target = asampler, type = 'slice')
        }
      }
    }

    ## replace all RW samplers with slice
    ## Should find a way to do this faster
    ## testreptime <- Sys.time()
    if (RWtoslice) {
      for (asampler in targetslist[nameslist == 'RW']) {
        confnimble$replaceSampler(target=asampler, type='slice')
        ## ## Old replacement method, didn't work in previous Nimble
        ## confnimble$removeSamplers(asampler)
        ## confnimble$addSampler(target = asampler, type = 'slice')
      }
    }

    ## ## Uncomment when debugging Nimble
    ## print(confnimble$getUnsampledNodes())

    #### change execution order for some variates
    if (changeSamplerOrder) {
      ## call this to do a first reordering of the samplers
      mcsampler <- buildMCMC(confnimble)

      samplerorder <- c(
      'K',
      if (vn$R > 0) {
        c('Rmean', 'Rrate', 'Rvar')
      },
      if (vn$C > 0) {
        c('Cmean', 'Crate', 'Cvar')
      },
      if (vn$D > 0) {
        c('Dmean', 'Drate', 'Dvar')
      },
      if (vn$L > 0) {
        c('Lmean', 'Lrate', 'Lvar')
      },
      if (vn$O > 0) {
        c('Oprob')
      },
      if (vn$N > 0) {
        c('Nprob')
      },
      if (vn$B > 0) {
        c('Bprob')
      },
      'W', 'Alpha'
    )
    ##
    neworder <- foreach(var = samplerorder, .combine = c) %do% {
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
    ## ## Uncomment for debugging
    ## cat('\n********NEW ORDER',neworder,'\n')
    confnimble$setSamplerExecutionOrder(c(setdiff(
      confnimble$getSamplerExecutionOrder(),
      neworder
    ), neworder))

    }

#### Compile Monte Carlo sampler
    print(confnimble)
    mcsampler <- buildMCMC(confnimble)
    ## print(confnimble$getUnsampledNodes())
    Cmcsampler <- compileNimble(mcsampler, resetFunctions = TRUE)

    cat('\nSetup time', printtime(Sys.time() - timecount), '\n')

    ## Inform user that compilation is done, if core 1:
    if (acore == 1) {
      print2user(paste0('\rCompiled core ', acore,
                       '. Estimating remaining time, please be patient...'),
                outcon)
    }

    ## cat('Loop over chains')
    ##################################################
    #### LOOP OVER CHAINS (WITHIN ONE CORE)
    ##################################################
    maxusedclusters <- 0


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
      ##  niter <- max(min(niterini,requirediter*2), 128)
      ##  }
      nitertot <- 0
      requirediter <- +Inf
      reset <- TRUE
      traces <- NULL
      allmcsamples <- NULL
      allmcsamplesKA <- list(Alpha = NULL, K = NULL)
      mcsamplesKA <- NULL
      flagll <- FALSE
      flagmc <- FALSE
      gc() #garbage collection
      chainnumber <- (acore - 1L) * nchainspercore + achain
      padchainnumber <- sprintf(paste0('%0', nchar(nchains), 'i'), chainnumber)
      cat(
        '\nChain #', chainnumber,
        '(chain', achain, 'of', nchainspercore, 'for this core)\n'
      )
      ## Read data to be used in log-likelihood
      testdata <- readRDS(file = file.path(dirname,
                      paste0('_testdata_', chainnumber, '.rds')))

      ## Initial values for this chain
      ## random seed is taken care of by %doRNG%
      Cfinitemixnimble$setInits(initsfn())

      #### WHILE LOOP CONTINUING UNTIL CONVERGENCE
      subiter <- 1L
      while (nitertot < requirediter) {
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
            max(2, floor(niter / 2))
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

        if (showAlphatraces) {
          dim(mcsamplesKA$Alpha) <- NULL # from matrix to vector
        }

        cat('\nCurrent time:',
            strftime(as.POSIXlt(Sys.time()), '%Y-%m-%d %H:%M:%S'))
        cat('\nMCMC time', printtime(Sys.time() - starttime), '\n')

        ## #### Remove iterations with non-finite values
        ## ## old version
        ##                 toRemove <- which(!is.finite(mcsamples), arr.ind=TRUE)
        ##                 ##
        ##                 if(length(toRemove) > 0){
        ##                     print2user('\nWARNING: SOME NON-FINITE OUTPUTS\n', outcon)
        ##                     ##
        ##                     flagmc <- TRUE
        ##                     allflagmc <- TRUE
        ##                     if(length(unique(toRemove[,1])) == nrow(mcsamples)){
        ##                         suppressWarnings(sink())
        ##                         suppressWarnings(sink(NULL,type='message'))
        ##                         registerDoSEQ()
        ##                         if(ncores > 1){ parallel::stopCluster(cl) }
        ##                         stop('...TOO MANY NON-FINITE OUTPUTS. ABORTING')
        ##                     }else{
        ##                         mcsamples <- mcsamples[-unique(toRemove[,1]),,drop=FALSE]
        ##                     }
        ##                 }
#### Remove iterations with non-finite values
        if(any(!is.finite(unlist(mcsamples)))) {
          toRemove <- sort(unique(unlist(lapply(mcsamples, function(xx) {
            temp <- which(!is.finite(xx), arr.ind = TRUE)
            temp[, ncol(temp)]
          }))))

          cat('\nWARNING:', length(toRemove), 'NON-FINITE SAMPLES\n')
          ##
          flagmc <- TRUE
          allflagmc <- TRUE
          saveRDS(mcsamples, file = file.path(dirname,
            paste0('NONFINITEmcsamples',
            dashnameroot, '--', padchainnumber,
            '_', achain, '-',
            acore, '-i', nitertot, '.rds')
          ))
          if (length(toRemove) == ncol(mcsamples$W)) {
            cat('\n...TOO MANY NON-FINITE OUTPUTS!\n')
            ## print2user('\n...TOO MANY NON-FINITE OUTPUTS!\n', outcon)
            ## suppressWarnings(sink())
            ## suppressWarnings(sink(NULL,type='message'))
            ## ## registerDoSEQ()
            ## if(exists('cl')){ parallel::stopCluster(cl) }
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
          print(sort(sapply(sprefixes,
                            function(x) sum(samplertimes[grepl(x, names(samplertimes))])),
            decreasing = TRUE
          ))
        }

        ## Check how many clusters were used at the last step
        ## usedclusters <- mcsamplesKA$K[length(mcsamplesKA$K)]
        usedclusters <- max(mcsamplesKA$K)
        cat('\nUSED CLUSTERS:', usedclusters, 'OF', nclusters, '\n')
        #### Diagnostics
        ## Log-likelihood
        diagntime <- Sys.time()
        ##
        ll <- cbind(
          log(samplesFDistribution(
            Y = testdata, X = NULL,
            mcoutput = c(mcsamples, list(auxmetadata = auxmetadata)),
            jacobian = FALSE,
            useLquantiles = useLquantiles,
            parallel = FALSE,
            combine = '+',
            silent = TRUE
          )) / nrow(testdata)
        )

        ## colnames(ll) <- paste0('log-', c('mid', 'lo', 'hi'))
        ## colnames(ll) <- paste0('log-', seq_len(ncol(ll)))
       colnames(ll) <- paste0('log-F of ', nrow(testdata), ' data')

        ## if (is.numeric(loglikelihood)) {
        ##   lltime <- Sys.time()
        ##   cat('\nCalculating log-likelihood...')
        ##   ll <- cbind(ll,
        ##     'log-ll' = log(samplesFDistribution(
        ##       Y = lldata, X = NULL,
        ##       ## Y = data[llseq, ], X = NULL,
        ##       mcoutput = c(mcsamples, list(auxmetadata = auxmetadata)),
        ##       jacobian = FALSE,
        ##       useLquantiles = useLquantiles,
        ##       parallel = FALSE, silent = TRUE,
        ##       combine = '+'
        ##     )) / nrow(lldata)
        ##   )
        ##   cat('Done,\n', printtime(Sys.time() - lltime), '\n')
        ## }
        ##
        traces <- rbind(traces, ll * 10 / log(10))

        toRemove <- which(!is.finite(traces), arr.ind = TRUE)
        if (length(toRemove) > 0) {
          flagll <- TRUE
          cleantraces <- traces[-unique(toRemove[, 1]), , drop = FALSE]
        } else {
          cleantraces <- traces
        }
        ##
        diagnESS <- LaplacesDemon::ESS(cleantraces)
        cat('\nESSs:', paste0(round(diagnESS), collapse = ', '))
        diagnIAT <- apply(cleantraces, 2, function(x) {
          LaplacesDemon::IAT(x)
        })
        cat('\nIATs:', paste0(round(diagnIAT), collapse = ', '))
        diagnBMK <- LaplacesDemon::BMK.Diagnostic(cleantraces[1:(4 * trunc(nrow(cleantraces) / 4)), , drop = FALSE], batches = 4)[, 1]
        cat('\nBMKs:', paste0(round(diagnBMK, 3), collapse = ', '))
        diagnMCSE <- 100 * apply(cleantraces, 2, function(x) {
          funMCSE(x) / sd(x)
        })
        cat('\nMCSEs:', paste0(round(diagnMCSE, 2), collapse = ', '))
        diagnStat <- apply(cleantraces, 2, function(x) {
          LaplacesDemon::is.stationary(as.matrix(x, ncol = 1))
        })
        cat('\nStationary:', paste0(diagnStat, collapse = ', '))
        diagnBurn <- apply(cleantraces, 2, function(x) {
          LaplacesDemon::burnin(matrix(x[1:(10 * trunc(length(x) / 10))], ncol = 1))
        })
        cat('\nBurn-in I:', paste0(diagnBurn, collapse = ', '))
        diagnBurn2 <- proposeburnin(cleantraces, batches = 10)
        cat('\nBurn-in II:', diagnBurn2)
        diagnThin <- proposethinning(cleantraces)
        cat('\nProposed thinning:', paste0(diagnThin, collapse = ', '), '\n')
        rm(cleantraces)

        cat('\nDiagnostics time', printtime(Sys.time() - diagntime), '\n')

        if (is.null(allmcsamples)) {
          ## new chain
          allmcsamples <- mcsamples
          allmcsamplesKA <- mcsamplesKA
        } else {
          ## continue chain, concat samples
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

        rm(mcsamples)
        gc() #garbage collection

        nitertot <- ncol(allmcsamples$W)


        ##########################################
        #### CHECK IF CHAIN MUST BE CONTINUED ####
        ##########################################

        calcIterThinning <- mcmclength(nsamplesperchain = nsamplesperchain,
                                       thinning=thinning,
                                       diagnESS=diagnESS, diagnIAT=diagnIAT,
                                       diagnBMK=diagnBMK, diagnMCSE=diagnMCSE,
                                       diagnStat=diagnStat, diagnBurn=diagnBurn,
                                       diagnBurn2=diagnBurn2, diagnThin=diagnThin)

        requirediter <- max(miniter,
                             min(maxiter, calcIterThinning$iter))
        currentthinning <- calcIterThinning$thinning

        cat('\nNumber of iterations', nitertot, ', required', requirediter, '\n')
        ##
        if (nitertot < requirediter) {
          ## limit number of iterations per loop, to save memory
          niter <- min(requirediter - nitertot + 1L, niterini)
          subiter <- subiter + 1L
          cat(
            '\nChain #', chainnumber, '- chunk', subiter,
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
          currentthinning, '\n')

      tokeep <- seq(to = ncol(allmcsamples$W),
                    length.out = nsamplesperchain,
                    by = currentthinning)

      saveRDS(mcsubset(allmcsamples, tokeep),
              file = file.path(dirname,
                               paste0('_mcsamples',
                                      dashnameroot, '--',
                                      padchainnumber, '.rds'))
                               )
      rm(allmcsamples)
      ## nitertot <- ncol(allmcsamples$W)

      gc() #garbage collection

      saveRDS(traces[tokeep, , drop = FALSE],
              file = file.path(dirname,
                               paste0('_mctraces',
                            dashnameroot, '--',
                            padchainnumber, '.rds')
                            ))


      ###############
      #### PLOTS ####
      ###############

      #### Plot Alpha and cluster occupation, if required
      if (showAlphatraces || showKtraces) {
        cat('Plotting cluster and Alpha information.\n')
        pdff(file.path(dirname,
                       paste0('hyperparams_traces', dashnameroot, '--',
                    padchainnumber, '_', achain, '-', acore)),
             apaper = 4)

        if (showKtraces) {
          cat('\nSTATS USED CLUSTERS:\n')
          print(summary(allmcsamplesKA$K))
          ##
          tplot(y = allmcsamplesKA$K, ylab = 'used clusters',
                xlab = 'iteration', ylim = c(0, nclusters))
          tplot(x = ((-1):nclusters) + 0.5,
                y = tabulate(allmcsamplesKA$K + 1, nbins = nclusters + 1),
                type = 'h', xlab = 'used clusters', ylab = NA,
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

        tracegroups <- as.list(seq_len(min(4, ncol(traces))))
        names(tracegroups) <- colnames(traces)[1:min(4, ncol(traces))]
        grouplegends <- foreach(agroup = seq_along(tracegroups)) %do% {
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
        ## colpalette <- 1:6 #seq_len(ncol(traces))
        ## names(colpalette) <- colnames(traces)[1:min(6, ncol(traces))]
        cat('\nPlotting MCMC traces')
        graphics.off()
        pdff(file.path(dirname,
                       paste0('_mcpartialtraces', dashnameroot, '--',
                              padchainnumber, '_', achain, '-', acore)
                       ), apaper = 4)

        ## Summary stats
        matplot(1:2, type = 'l', col = 'white',
                main = paste0('Stats chain ', achain),
                axes = FALSE, ann = FALSE)
        # Legends
        legendpositions <- c('topleft', 'topright', 'bottomleft', 'bottomright')
        for (alegend in seq_along(grouplegends)) {
          legend(x = legendpositions[alegend], bty = 'n', cex = 1.5,
            legend = grouplegends[[alegend]]
          )
        }
        legend(x = 'center', bty = 'n', cex = 1,
          legend = c(
            paste0('Chain ', chainnumber, '_', achain, '-', acore),
            paste0('Used clusters: ', usedclusters, ' of ', nclusters),
            paste0('LL:  ( ', signif(mean(traces[, 1]), 3), ' +- ',
                   signif(sd(traces[, 1]), 3), ' ) dHart'),
            'NOTES:',
            if (flagmc) {
              'some non-finite MC outputs'
            },
            if (usedclusters > nclusters - 5) {
              'too many clusters used'
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
            y = traces[is.finite(traces[, avar]), avar],
            type = 'l', lty = 1, col = 1,
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

      ## #### Plot samples from current chain
      ## ## Can be uncommented if debugging is needed
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
      ##                              plotvariability='samples',
      ##                              nsamples=showsamples,
      ##                              plotmeans=plotmeans,
      ##                              datahistogram=TRUE, datascatter=TRUE,
      ##                              useLquantiles=useLquantiles,
      ##                              parallel=FALSE, silent=TRUE
      ##                              )
      ##          }

      cat('\nCurrent time:', strftime(as.POSIXlt(Sys.time()), '%Y-%m-%d %H:%M:%S'))
      cat('\nMCMC + diagnostics time', printtime(Sys.time() - starttime), '\n')

      #### Print estimated remaining time
      remainingTime <- (Sys.time() - starttime) / achain * (nchainspercore - achain + 1)
      if (is.finite(remainingTime) && remainingTime > 0) {
        print2user(
          paste0(
            '\rSampling. Core ', acore, ' estimated remaining time: ',
            printtime(remainingTime),
            '                      '
          ),
          outcon
        )
      }

      maxusedclusters <- max(maxusedclusters, usedclusters)
    } #### END LOOP OVER CHAINS (WITHIN ONE CORE)

    ##
    cat('\nCurrent time:',
        strftime(as.POSIXlt(Sys.time()), '%Y-%m-%d %H:%M:%S'))
    cat('\nTotal time', printtime(Sys.time() - starttime), '\n')

    cbind(maxusedclusters, allflagmc)
    }
############################################################
#### END OF FOREACH-LOOP OVER CORES
############################################################
  ## Close output to log files
  suppressWarnings(sink())
  suppressWarnings(sink(NULL, type = 'message'))

  maxusedclusters <- max(chaininfo[, 1])
  nonfinitechains <- sum(chaininfo[, 2])

  ############################################################
  #### End of all MCMC
  ############################################################


  ############################################################
  #### Join chains
  ############################################################

  ## Save random seeds used in the parallel processing
  if (!is.null(attr(chaininfo, 'rng'))) { # parallel processing
    saveRDS(attr(chaininfo, 'rng'),
            file = file.path(dirname,
                             paste0('rng_parallelseeds', dashnameroot, '.rds')
                             ))
  }

  ## Read the samples saved by each chain and concatenate them
  ## could be moved to a separate file?
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
  mcsamples <- foreach(chainnumber = 1:(ncores * nchainspercore),
                                .combine = joinmc, .multicombine = FALSE) %do% {
    padchainnumber <- sprintf(paste0('%0', nchar(nchains), 'i'), chainnumber)

    readRDS(file = file.path(dirname,
                             paste0('_mcsamples', dashnameroot, '--',
                                    padchainnumber, '.rds')
                             ))
  }

  #### Save all final parameters together with the aux-metadata in one file
  saveRDS(c(mcsamples, list(auxmetadata = auxmetadata, nchains=nchains)),
          file = file.path(dirname,
                           paste0('Fdistribution', dashnameroot, '.rds')
  ))

  traces <- foreach(chainnumber = 1:(ncores * nchainspercore),
                             .combine = rbind) %do% {
    padchainnumber <- sprintf(paste0('%0', nchar(nchains), 'i'), chainnumber)

    readRDS(file = file.path(dirname,
                             paste0('_mctraces', dashnameroot, '--',
                                    padchainnumber, '.rds')
                             ))
  }
  saveRDS(traces, file = file.path(dirname,
                                   paste0('MCtraces', dashnameroot, '.rds')
                                   ))

  cat('\rFinished Monte Carlo sampling.                                 \n')
  gc() #garbage collection

  ############################################################
  #### Final joint diagnostics
  ############################################################

  ## cat('\nSome diagnostics:\n')
  traces <- traces[apply(traces, 1, function(x) { all(is.finite(x)) }), ,
                   drop = FALSE]
  ## flagll <- nrow(traces) != nrow(traces)
  ## funMCSE <- function(x){LaplacesDemon::MCSE(x, method='batch.means')$se}
  diagnESS <- LaplacesDemon::ESS(traces)
  cat('\nEFFECTIVE SAMPLE SIZE > ', nchains,
      '(estim.: ', round(max(diagnESS, na.rm = TRUE)), ')\n')

  ## ## Commented, but we may reconsider displaying this
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
  cat('\nPlotting final Monte Carlo traces.\n')

  ##
  ## colpalette <- seq_len(ncol(traces))
  ## names(colpalette) <- colnames(traces)
  graphics.off()
  pdff(file.path(dirname,
                 paste0('MCtraces', dashnameroot)
                 ), apaper = 4)
  ## Traces of likelihood and cond. probabilities
  for (avar in colnames(traces)) {
    ## Do not join separate chains in the plot
    division <- (if(nrow(traces) > nchains) nchains else 1)
    tplot(x = matrix(seq_len(nsamples), ncol = division),
          y = matrix(traces[, avar], ncol = division),
          type = 'l', lty = 1,
          col = 1:6, # to evidence consecutive chains
          ## col = colpalette[avar], # original, one color per trace
      main = paste0(
        'Effective sample size > ', nchains,
        ' (estimated: ', signif(diagnESS[avar], 3), ')'
        ##             ' | IAT = ', signif(diagnIAT[avar], 3),
        ##             ' | BMK = ', signif(diagnBMK[avar], 3),
        ##             ' | MCSE = ', signif(diagnMCSE[avar], 3),
        ##             ' | stat: ', diagnStat[avar],
        ##             ' | burn I: ', diagnBurn[avar],
        ##             ' | burn II: ', diagnBurn2
      ),
      ylab = paste0(avar, '/dHart'), xlab = 'sample', family = family,
      mar = c(NA, 6, NA, NA)
    )
  }

  cat('\nMax number of used clusters:', maxusedclusters, '\n')
  if (maxusedclusters > nclusters - 5) {
    cat('TOO MANY CLUSTERS USED!\n')
    cat('Consider re-running with increased "nclusters" parameter\n')
  }

  if (nonfinitechains > 0) {
    cat(nonfinitechains, 'chains with some non-finite outputs\n')
  }

  cat('Plotting marginal samples.\n')
  plotFsamples(
    file = file.path(dirname,
                     paste0('plotsamples_Fdistribution', dashnameroot)),
    mcoutput = c(mcsamples, list(auxmetadata = auxmetadata)),
    data = data,
    plotvariability = 'samples',
    nFsamples = showsamples, plotmeans = TRUE,
    datahistogram = TRUE, datascatter = TRUE,
    useLquantiles = useLquantiles,
    parallel = TRUE, silent = TRUE
  )

  cat('Plotting marginal samples with quantiles.\n')
  plotFsamples(
    file = file.path(dirname,
                     paste0('plotquantiles_Fdistribution', dashnameroot)),
    mcoutput = c(mcsamples, list(auxmetadata = auxmetadata)),
    data = data,
    plotvariability = 'quantiles',
    nFsamples = plotDisplayedQuantiles, plotmeans = TRUE,
    datahistogram = TRUE, datascatter = TRUE,
    useLquantiles = useLquantiles,
    parallel = TRUE, silent = TRUE
  )

  cat('\nTotal computation time', printtime(Sys.time() - timestart0), '\n')

  if (exists('cl')) {
    cat('\nClosing connections to cores.\n')
    parallel::stopCluster(cl)
  }


#### Remove partial files if required
  ## Partial files are identified by an initial underscore "_"
    ## Should we leave the plots of partial traces?
    ## maybe create an additional argument to let the user decide?
  if (cleanup) {
    cat('Removing temporary output files.\n')
    file.remove(dir(dirname,
      pattern = paste0('^_.*\\..*$'),
      full.names = TRUE
    ))
  }

  #### WHAT MIGHT BE ADDED:
  ## - when 'showKtraces' is true:
  ## a histogram over number of clusters over all chains
  ## (for the moment there's one plot per chain)

  cat('Finished.\n\n')

    ## What should we output? how about the full name of the output dir?
  if (is.character(output) && output == 'directory') {
    dirname
  } else if (is.character(output) && output == 'mcoutput') {
    readRDS(file.path(dirname,
                      paste0('Fdistribution', dashnameroot, '.rds')
                      ))
  }
}

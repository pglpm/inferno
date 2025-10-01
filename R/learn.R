#' Monte Carlo computation of posterior probability distribution
#'
#' @description
#' Core function to compute the posterior probability distribution of the variates conditional on the given data.
#'
#' @details
#' This function takes as main inputs a set of data and metadata, and computes the probability distribution for new data. Its computation can also be interpreted as an estimation of the frequencies of the variates in the *whole population*, beyond the sample data. The probability distribution is not assumed to be Gaussian or of any other specific shape. The computation is done via Markov-chain Monte Carlo.
#'
#' This function creates an object, contained in a `learnt.rds` file, which is used in all subsequent probabilistic computations. Other information about the computation is provided in logs and plots, saved in a directory specified by the user.
#'
#' See `vignette('inferno_start')` for an introductory example.
#'
#' @param data A dataset, given as a [base::data.frame()] or as a file path to a CSV file.
#' @param metadata A [`metadata`] object, given either as a data.frame object, or as a file pa to a CSV file.
#' @param auxdata A larger dataset, given as a base::data.frame() or as a file path to a CSV file. Such a dataset would be too many to use in the Monte Carlo sampling, but can be used to calculate hyperparameters.
#' @param outputdir Character: path to folder where the output should be saved. If omitted or `NULL` (default), a directory is created that has the same name as the data file but with suffix "`_output_`". If `FALSE`, a directory is created in the temporary-directory space.
#' @param nsamples Integer: number of desired Monte Carlo samples. Default 3600.
#' @param nchains Integer: number of Monte Carlo chains. Default 4.
#' @param nsamplesperchain Integer: number of Monte Carlo samples per chain.
#' @param parallel Logical or positive integer or cluster object. `TRUE`: use roughly half of available cores; `FALSE`: use serial computation; integer: use this many cores. It can also be a cluster object previously created with [parallel::makeCluster()]; in this case the parallel computation will use this object.
#' @param seed Integer: use this seed for the random number generator. If missing or `NULL` (default), do not set the seed.
#' @param cleanup Logical: remove diagnostic files at the end of the computation? Default `TRUE`.
#' @param appendtimestamp Logical: append a timestamp to the name of the output directory `outputdir`? Default `TRUE`.
#' @param appendinfo Logical: append information about dataset and Monte Carlo parameters to the name of the output directory `outputdir`? Default `TRUE`.
#' @param output Character: if `'directory'`, return the output directory name as `VALUE`; if character `'learnt'`, return the `learnt` object containing the parameters obtained from the Monte Carlo computation. Any other value: `VALUE` is `NULL`.
#' @param subsampledata Integer: use only a subset of this many datapoints for the Monte Carlo computation.
#' @param prior Logical: Calculate the prior distribution?
#' @param startupMCiterations Integer: number of initial Monte Carlo iterations. Default 3600.
#' @param minMCiterations Integer: minimum number of Monte Carlo iterations to be doneby a chain. Default 0.
#' @param maxMCiterations Integer: Do at most this many Monte Carlo iterations per chain. Default `Inf`.
#' @param maxhours Numeric: approximate time limit, in hours, for the Monte Carlo computation to last. Default `Inf`.
#' @param ncheckpoints Integer: number of datapoints to use for checking when the Monte Carlo computation should end. If `NULL`, this is equal to number of variates + 2. If Inf, use all datapoints. Default 12.
#' @param maxrelMCSE Numeric positive: desired maximal *relative Monte Carlo Standard Error* of calculated probabilities with respect to their variability with new data. Default `+Inf`, so that `minESS` is used instead. `maxrelMCSE` is related to `minESS` by `maxrelMCSE = 1/sqrt(minESS + initES)`.
#' @param minESS Numeric positive: desired minimal Monte Carlo *Expected Sample Size*. If `NULL`, it is equal to the final `nsamplesperchain`. Default 400. `minESS` is related to `maxrelMCSE` by `minESS = 1/maxrelMCSE^2 - initES`.
#' @param initES Numeric positive: number of initial  *Expected Samples* to discard.
#' @param thinning Integer: thin out the Monte Carlo samples by this value. If `NULL` (default): let the diagnostics decide the thinning value.
#' @param plottraces Logical: save plots of the Monte Carlo traces of diagnostic values? Default `TRUE`.
#' @param showKtraces Logical: save plots of the Monte Carlo traces of the K parameter? Default `FALSE`.
#' @param showAlphatraces Logical: save plots of the Monte Carlo traces of the Alpha parameter? Default `FALSE`.
#' @param hyperparams List: hyperparameters of the prior.
#'
#' @return Name of directory containing output files, or learnt object, or `NULL`, depending on argument `output`.
#'
#' @examples
#'
#' ## Create dataset with 10 points of variate 'V' for demonstration
#' dataset <- data.frame(V = rnorm(n = 10))
#'
#' ## Create metadatafile
#' metadata <- data.frame(
#'     name = 'V',
#'     type = 'continuous'
#' )
#'
#' @import parallel nimble
#'
#' @export
learn <- function(
    data,
    metadata,
    auxdata = NULL,
    outputdir = NULL,
    nsamples = 3600,
    nchains = 8,
    nsamplesperchain = 450,
    parallel = TRUE,
    seed = NULL,
    cleanup = TRUE,
    appendtimestamp = TRUE,
    appendinfo = TRUE,
    output = 'directory',
    subsampledata = NULL,
    prior = missing(data) || is.null(data),
    startupMCiterations = 3600,
    minMCiterations = 0,
    maxMCiterations = +Inf,
    maxhours = +Inf,
    ncheckpoints = 12,
    maxrelMCSE = +Inf, ## Gong-Flegal: 0.038, Z=1000: 0.076, Z=400: 0.12
    minESS = 450, ## Gong-Flegal: 0.038, Z=1000: 0.076, Z=400: 0.12
    initES = 2,
    thinning = NULL,
    plottraces = !cleanup,
    showKtraces = FALSE,
    showAlphatraces = FALSE,
    hyperparams = list(
        ncomponents = 64,
        minalpha = -4,
        maxalpha = 4,
        byalpha = 1,
        Rshapelo = 0.5,
        Rshapehi = 0.5,
        Rvarm1 = 3^2,
        Cshapelo = 0.5,
        Cshapehi = 0.5,
        Cvarm1 = 3^2,
        Dshapelo = 0.5,
        Dshapehi = 0.5,
        Dvarm1 = 3^2,
        Bshapelo = 1,
        Bshapehi = 1,
        Dthreshold = 1,
        tscalefactor = 4.266,
        Oprior = 'Hadamard',
        Nprior = 'Hadamard',
        avoidzeroW = NULL, # NULL: Turek's, TRUE: eps non-conj., FALSE: conj.
        initmethod = 'datacentre',
        Qerror = pnorm(c(-1, 1))
        ## Qerror = c(0.055, 0.945) # pnorm(c(-1, 1))
    )
) {

##################################################
#### Various internal parameters
##################################################

#### Hyperparameters and other internal parameters
    hyperparamsDefaults <- list(
        ncomponents = 64,
        minalpha = -4,
        maxalpha = 4,
        byalpha = 1,
        Rshapelo = 0.5,
        Rshapehi = 0.5,
        Rvarm1 = 3^2,
        Cshapelo = 0.5,
        Cshapehi = 0.5,
        Cvarm1 = 3^2,
        Dshapelo = 0.5,
        Dshapehi = 0.5,
        Dvarm1 = 3^2,
        Bshapelo = 1,
        Bshapehi = 1,
        Dthreshold = 1,
        tscalefactor = 4.266,
        Oprior = 'Hadamard',
        Nprior = 'Hadamard',
        avoidzeroW = NULL,
        initmethod = 'datacentre',
        Qerror = pnorm(c(-1, 1))
        ## Qerror = c(0.055, 0.945) # pnorm(c(-1, 1))
    )
    Qlo <- 0.055 # (100 - 89) / 200 # 0.055
    Qhi <- 0.945 # (100 + 89) / 200 # 0.945

    ## Allow user to specify hyperparameters only partially:
    ## missing ones are given default values above
    ## Check for unknown names first
    if(length(
        setdiff(names(hyperparams), names(hyperparamsDefaults))
    ) > 0){
        stop('Unknown hyperparameters:',
            setdiff(names(hyperparams), names(hyperparamsDefaults))
        )
    }

    for(aname in names(hyperparamsDefaults)){
        if(is.null(hyperparams[[aname]])){
            hyperparams[[aname]] <- hyperparamsDefaults[[aname]]
        }
        assign(aname, hyperparams[[aname]])
    }

#### Other options
    Alphatoslice <- FALSE # FALSE typically leads to underflow
    Ktoslice <- FALSE # FALSE typically leads to underflow
    RWtoslice <- FALSE
    changeSamplerOrder <- TRUE
    ##
    showsamples <- 100 # number of samples to show.
    plotDisplayedQuantiles <- c(5.5, 94.5)/100 # c(1, 31) / 32 # quantiles to show
    ncomponentsamples <- 128 # number of samples of Alpha and K
    showsamplertimes <- FALSE ##
    family <- 'Palatino' # font family in plots


#### Start timer
    timestart0 <- Sys.time()

    cat('\n') # make sure possible error messages start on new line

    ## Set the RNG seed if given by user, or if no seed already exists
    if (!is.null(seed) || !exists('.Random.seed')) {
        set.seed(seed)
    }
    currentseed <- .Random.seed

##################################################
#### Argument-consistency checks
##################################################

#### Consistency checks for numbers of samples, chains, cores
    ## The defaults are 3600 samples from 60 chains, so 60 samples per chain
    ## The user can choose any two
    ## nsamples = nchains * nsamplesperchain

    if(!missing(nchains) && !missing(nsamplesperchain) &&
           missing(nsamples)){
        ## given: nchains, nsamplesperchain
        nsamples <- nchains * nsamplesperchain
    } else if (missing(nchains) && !missing(nsamplesperchain) &&
                   !missing(nsamples)){
        ## given: nsamples, nsamplesperchain
        nchains <- ceiling(nsamples / nsamplesperchain)
        if(nsamples != nchains * nsamplesperchain){
            nsamples <- nchains * nsamplesperchain
            cat('Increasing number of samples to', nsamples,
                'to comply with given "nsamplesperchain"\n')
        }
    } else if (!missing(nchains) && missing(nsamplesperchain) &&
                   !missing(nsamples)){
        ## given: nsamples, nchains
        nsamplesperchain <- ceiling(nsamples / nchains)
        if(nsamples != nchains * nsamplesperchain){
            nsamples <- nchains * nsamplesperchain
            cat('Increasing number of samples to', nsamples,
                'to comply with given "nchains"\n')
        }
    } else if (!(missing(nchains) && missing(nsamplesperchain) &&
                     missing(nsamples))){
        stop('Please specify exactly two among "nsamples", "nchains", "nsamplesperchain"')
    }


#### Requested parallel processing
    ## NB: doesn't make sense to have more cores than chains
    closeexit <- FALSE
    if ('cluster' %in% class(parallel)){
        ## user provides a cluster object
        cl <- parallel
    } else if (isTRUE(parallel)) {
        ## user wants us to register a parallel backend
        ## and to choose number of cores
        ncores <- max(1,
            min(nchains, floor(parallel::detectCores() / 2)))
        cl <- parallel::makeCluster(ncores)
        ## doParallel::registerDoParallel(cl)
        closeexit <- TRUE
        cat('Registered', ncores, 'workers\n\n')
    } else if (isFALSE(parallel)) {
        ## user wants us not to use parallel cores
        ncores <- 1
        cl <- parallel::makeCluster(ncores)
    } else if (is.numeric(parallel) &&
                   is.finite(parallel) && parallel >= 1) {
        ## user wants us to register 'parallel' # of cores
        ncores <- min(nchains, parallel)
        cl <- parallel::makeCluster(ncores)
        closeexit <- TRUE
        cat('Registered', ncores, 'workers\n\n')
    } else {
        stop("Unknown value of argument 'parallel'")
    }

    ## Close parallel connections if any were opened
    if(closeexit) {
        closecoresonexit <- function(){
            cat('\nClosing connections to cores.\n')
            parallel::stopCluster(cl)
            ## parallel::setDefaultCluster(NULL)
        }
        on.exit(closecoresonexit())
    }

    minchainspercore <- nchains %/% ncores
    coreswithextrachain <- nchains %% ncores

    if (is.numeric(thinning) && thinning > 0) {
        thinning <- ceiling(thinning)
    } else if (!is.null(thinning)) {
        stop('Invalid "thinning" argument.')
    }

    ## Parallellisation is done in any case,
    ## so that objects from the Monte Carlo simulation are not left in memory.

    ## Make sure 'startupMCiterations' is at least equal to 2
    startupMCiterations <- max(2, startupMCiterations)

    ## maxrelMCSE and minESS
    maxrelMCSE <- min(maxrelMCSE, 1 / sqrt(minESS + initES))
    minESS <- (1 / maxrelMCSE)^2 - initES

    if(all(!is.finite(c(maxrelMCSE, maxMCiterations)))) {
        stop('Impossible stopping rule')
    }

##################################################
#### Read and process data and metadata
##################################################


#### Read metadata
    ## Check whether argument 'metadata' is a string for a file name
    ## otherwise we assume it's a data.frame or similar object
    if (is.character(metadata) && file.exists(metadata)) {
        metadata <- read.csv(metadata,
            na.strings = '', stringsAsFactors = FALSE,
            colClasses=c(
                name = 'character',
                type = 'character',
                domainmin = 'numeric',
                domainmax = 'numeric',
                datastep = 'numeric',
                minincluded = 'character',
                maxincluded = 'character'
            ))
    }
    metadata <- as.data.frame(metadata)

    ## eliminate possible empty V-columns
    for(i in intersect(paste0('V', 11:3), colnames(metadata))){
        if(all(is.na(metadata[, i]))){
            metadata <- metadata[, -which(colnames(metadata) == i), drop=FALSE]
        }
    }

#### Dataset
    ## Check if 'data' is given
    if(!(missing(data) || is.null(data))) {
        ## Check if 'data' argument is an existing file
        ## otherwise we assume it is an object
        datafile <- NULL
        if (is.character(data)) {
            ## add '.csv' if missing
            datafile <- paste0(sub('.csv$', '', data), '.csv')
            if (file.exists(datafile)) {
                data <- read.csv(datafile,
                    na.strings = '', stringsAsFactors = FALSE, tryLogical = FALSE)
            } else {
                stop('Cannot find data file')
            }
        }
        data <- as.data.frame(data)
        rownames(data) <- NULL

        ## convert factors to strings if necessary
        if(any(sapply(data, is.factor))){
            cat('Converting factors to characters\n')
            . <- sapply(data, is.factor)
            data[, .] <- lapply(data[, ., drop = FALSE], as.character)
        }

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

#### Auxiliary dataset
    ## used to extract information about hyperparameters
    ## Check if 'auxdata' is given
    if(!(missing(auxdata) || is.null(auxdata))) {
        ## Check if 'auxdata' argument is an existing file
        ## otherwise we assume it is an object
        auxdatafile <- NULL
        if (is.character(auxdata)) {
            ## add '.csv' if missing
            auxdatafile <- paste0(sub('.csv$', '', auxdata), '.csv')
            if (file.exists(auxdata)) {
                auxdata <- read.csv(auxdatafile,
                    na.strings = '', stringsAsFactors = FALSE, tryLogical = FALSE)
            } else {
                stop('Cannot find auxdata file')
            }
        }
        auxdata <- as.data.frame(auxdata)

        ## Consistency checks for auxdata
        ## They should be moved to an external function

        ## Are auxdata missing variates?
        if (!all(metadata[['name']] %in% colnames(auxdata))) {
            stop('Missing variates in auxdata. Check auxdata- and metadata-files.')
        }

        ## Drop variates in auxdata that are not in the metadata file
        if (!all(colnames(auxdata) %in% metadata[['name']])) {
            cat('Warning: auxdata have additional variates. Dropping them.\n')
            subvar <- intersect(colnames(auxdata), metadata[['name']])
            auxdata <- auxdata[, subvar, drop = FALSE]
            rm(subvar)
        }

        ## Remove empty datapoints
        tokeep <- which(apply(auxdata, 1, function(xx) { !all(is.na(xx)) }))
        if(length(tokeep) == 0 && !prior) {
            stop('Auxdata are given but empty')
        } else if(length(tokeep) < nrow(auxdata)) {
            cat('Warning: auxdata contain empty datapoints. Dropping them.\n')
            auxdata <- auxdata[tokeep, , drop = FALSE]
        }
        rm(tokeep)

    } else {
        ## auxdata not given
        auxdata <- NULL
    }


    ## Build auxiliary metadata object; we'll save it later
    cat('Calculating auxiliary metadata\n')
    auxmetadata <- buildauxmetadata(
        data = (if (is.null(auxdata)) {data} else {auxdata}),
        metadata = metadata,
        Dthreshold = Dthreshold,
        tscalefactor = tscalefactor
    )
    ## print(auxmetadata) # for debugging

    cat('\nLearning: ', npoints, 'datapoints, ',
        nrow(auxmetadata), 'variates\n')

#### Output-folder setup
    if(isFALSE(outputdir)){
        ## Use a temporary directory
        dirname <- tempdir()
    } else {
        if (is.null(outputdir)) {
            outputdir <- paste0('_output_', sub('.csv$', '', datafile))
        }

            ## append time and info to output directory, if requested
            suffix <- NULL
            if (appendtimestamp) {
                suffix <- paste0(suffix, '-',
                    format(Sys.time(), '%y%m%dT%H%M%S') )
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
                    ## '-K', ncomponents, # unimportant for user
                    '_smp', nsamples)
            }
            dirname <- paste0(outputdir, suffix)
            ##
            ## Create output directory if it does not exist
            dir.create(dirname, showWarnings = FALSE)
    }
    ## Print information
    cat('\n', paste0(rep('*', max(nchar(dirname), 26)), collapse = ''),
        '\n Saving output in directory\n', dirname, '\n',
        paste0(rep('*', max(nchar(dirname), 26)), collapse = ''), '\n')

    ## This is in case we need to add some extra specifier to the output files
    ## all 'dashnameroot' can be deleted in a final version
    dashnameroot <- NULL

    ## Save copy of metadata and auxmetadata in directory
    write.csv(metadata, file = file.path(dirname, 'metadata.csv'),
        row.names = FALSE, quote = FALSE, na = '')
    saveRDS(auxmetadata,
        file = file.path(dirname, paste0('___auxmetadata', dashnameroot, '.rds')))

    ## Save initial RNG seed in case needed by user
    saveRDS(currentseed,
        file = file.path(dirname, paste0('rng_seed', dashnameroot, '.rds')))
    parallel::clusterSetRNGStream(cl = cl, currentseed)
#### number of checkpoints for Monte Carlo stopping rule
    if(is.null(ncheckpoints)) {
        ncheckpoints <- nrow(auxmetadata) + 1
    }
    ncheckpoints <- round(ncheckpoints)
    if (ncheckpoints < 1) {
        stop('"ncheckpoints" must be > 0')
    }

    ## if data is not empty, we use it to create data for likelihood
    if(!is.null(data)) {
        ## testdata is created and saved outside of the parallel processes
        ## so that the dataset does not need to be exported to them
        ## (using extra memory)
        ## Each chain uses a different set of testdata
        for(achain in 0:nchains) {
            pointsid <- sort(sample(seq_len(npoints), min(ncheckpoints, npoints)))
            testdata <- util_prepPcheckpoints(
                x = data[pointsid, , drop = FALSE],
                auxmetadata = auxmetadata,
                pointsid = pointsid
            )
            saveRDS(testdata,
                file = file.path(dirname, paste0('___testdata_', achain, '.rds')))
        }
        rm(testdata)
        rm(pointsid)
    }

#### Check if user wants to calculate prior
    if (prior) {
        message('CALCULATING PRIOR DISTRIBUTION')
        if(is.null(data)) {
            ## no data available: construct one datapoint from the metadata info
            testdata <- as.matrix(
                as.data.frame(
                    lapply(seq_len(nrow(auxmetadata)),
                        function(ii){
                            if(auxmetadata[ii, 'mcmctype'] %in%
                                   c('R', 'C', 'D', 'L')) {
                                rnorm(n = ncheckpoints, mean = 0, sd = 2)
                            } else {
                                sample(seq_len(auxmetadata[ii, 'Nvalues']),
                                    ncheckpoints,
                                    replace = TRUE) +
                                    auxmetadata[ii, 'indexpos']
                            }
                        }
                    ),
                    row.names = seq_len(ncheckpoints),
                    col.names = auxmetadata$name)
            )
            ## ## original alternative way. Delete soon
            ## testdata <- data.frame(1:3)[,-1]
            ## for(xx in seq_len(nrow(auxmetadata))) {
            ##   xx <- as.list(auxmetadata[xx, ])
            ##   toadd <- unlist(xx[paste0('mctest', 1:3)])
            ##   if (xx[['mcmctype']] %in% c('B', 'N', 'O')) {
            ##     toadd <- unlist(xx[paste0('V', toadd)])
            ##   }
            ##   testdata <- cbind(testdata, toadd)
            ## }
            ## colnames(testdata) <- auxmetadata[['name']]
            for(achain in 0:nchains) {
                saveRDS(testdata,
                    file = file.path(dirname, paste0('___testdata_', achain, '.rds')))
            }
            rm(testdata)
        }
        ## create empty dataset: Monte Carlo sampling is non-Markov
        data <- as.data.frame(
            matrix(NA, nrow = 1, ncol = nrow(metadata),
                dimnames = list(NULL, metadata$name))
        )
    }

    ## #### Select and save data for loglikelihood
    ##   ## Find which datapoints have not all missing entries
    ##   dataNoNa <- which(apply(data, 1, function(xx) { !all(is.na(xx)) }))
    ##   ndataNoNa <- length(dataNoNa)
    ##
    ##   if (ndataNoNa > 0) {
    ##     ## if "testdata" is moved into the for-loop,
    ##     ## then each chain uses a different set of testdata
    ##     for(achain in seq_len(nchains)) {
    ##       testdata <- data[sort(sample(dataNoNa, min(ncheckpoints, ndataNoNa))),]
    ##       saveRDS(testdata,
    ##               file = file.path(dirname, paste0('___testdata_', achain, '.rds')))
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
    ##           file = file.path(dirname, paste0('___testdata_', achain, '.rds')))
    ##     }
    ##     rm(testdata)
    ##     }




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
    vn <- vnames <- list(R=NULL, C=NULL, D=NULL, O=NULL, N=NULL, B=NULL)

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
    nalpha <- length(seq(minalpha, maxalpha, by = byalpha))
    probalpha0 <- (1:nalpha)^2.25
    probalpha0 <- probalpha0/sum(probalpha0)
    constants <- c(
        list(
            ncomponents = ncomponents,
            npoints = npoints,
            nalpha = nalpha,
            alphabase = sqrt(2),
            probalpha0 = probalpha0,
            dirchalphas = rep((2^(minalpha - 0.5)) / ncomponents, ncomponents)
        ),
        if (vn$R > 0) { # continuous open domain
            list(
                Rn = vn$R,
                Rmean1 = rep(0, 1),
                Rvarm1 = rep(Rvarm1, 1),
                Rvar1 = rep(1, 1),
                Rshapelo = rep(Rshapelo, 1),
                Rshapehi = rep(Rshapehi, 1)
            )
        },
        if (vn$C > 0) { # continuous closed domain
            list(
                Cn = vn$C,
                Cmean1 = rep(0, 1),
                Cvarm1 = rep(Cvarm1, 1),
                Cvar1 = rep(1, 1),
                Cshapelo = rep(Cshapelo, 1),
                Cshapehi = rep(Cshapehi, 1),
                Cleft = as.matrix(vtransform(data[, vnames$C, drop = FALSE],
                    auxmetadata = auxmetadata,
                    Cout = 'left', logjacobianOr = NULL)),
                Cright = as.matrix(vtransform(data[, vnames$C, drop = FALSE],
                    auxmetadata = auxmetadata,
                    Cout = 'right', logjacobianOr = NULL)),
                Clatinit = as.matrix(vtransform(data[, vnames$C, drop = FALSE],
                    auxmetadata,
                    Cout = 'init', logjacobianOr = NULL))
            )
            ## Cleft & Cright are as many as the datapoints
            ## so we do not create copies outside of Nimble
            ## to save RAM
        },
        if (vn$D > 0) { # continuous rounded
            list(
                Dn = vn$D,
                Dmean1 = rep(0, 1),
                Dvarm1 = rep(Dvarm1, 1),
                Dvar1 = rep(1, 1),
                Dshapelo = rep(Dshapelo, 1),
                Dshapehi = rep(Dshapehi, 1),
                Dleft = as.matrix(vtransform(data[, vnames$D, drop = FALSE],
                    auxmetadata = auxmetadata,
                    Dout = 'left', logjacobianOr = NULL)),
                Dright = as.matrix(vtransform(data[, vnames$D, drop = FALSE],
                    auxmetadata = auxmetadata,
                    Dout = 'right', logjacobianOr = NULL)),
                Dlatinit = as.matrix(vtransform(data[, vnames$D, drop = FALSE],
                    auxmetadata,
                    Dout = 'init', logjacobianOr = NULL))
            )
        },
        ## if (vn$L > 0) { # latent
        ##     list(
        ##         Ln = vn$L,
        ##         Lmean1 = rep(0, 1),
        ##         Lvarm1 = rep(Lvarm1, 1),
        ##         Lvar1 = rep(1, 1),
        ##         Lshapelo = rep(Lshapelo, 1),
        ##         Lshapehi = rep(Lshapehi, 1),
        ##         Lleft = as.matrix(vtransform(data[, vnames$L, drop = FALSE],
        ##             auxmetadata = auxmetadata,
        ##             Lout = 'left')),
        ##         Lright = as.matrix(vtransform(data[, vnames$L, drop = FALSE],
        ##             auxmetadata = auxmetadata,
        ##             Lout = 'right')),
        ##         Llatinit = as.matrix(vtransform(data[, vnames$L, drop = FALSE],
        ##             auxmetadata,
        ##             Lout = 'init'))
        ##     )
        ## },
        if (vn$O > 0) { # ordinal
            Ocards <- auxmetadata[auxmetadata$name %in% vnames$O, 'Nvalues']
            Omaxn <- sum(Ocards)
            Oi <- auxmetadata[auxmetadata$name %in% vnames$O, 'indexpos'] + 1L
            if(vn$O == 1){Oi <- c(Oi, NA_integer_)}
            Of <- Oi + Ocards - 1L
            ## ## we choose a flatter hyperprior for ordinal variates
            ## we choose a Hadamard-like hyperprior for nominal variates
            Oalpha0 <- if(Oprior == 'Hadamard') {
                as.vector(unlist(sapply(Ocards, function(acard){
                    rep(1 / acard, acard)
                })))
            } else {
                as.vector(unlist(sapply(Ocards, function(acard){
                    rep(1, acard)
                })))
            }
            ##
            list(
                On = vn$O,
                ## Ocards = Ocards,
                Oi = Oi,
                Of = Of,
                Oalpha0 = Oalpha0
            )
        },
        if (vn$N > 0) { # nominal
            Ncards <- auxmetadata[auxmetadata$name %in% vnames$N, 'Nvalues']
            Nmaxn <- sum(Ncards)
            Ni <- auxmetadata[auxmetadata$name %in% vnames$N, 'indexpos'] + 1L
            if(vn$N == 1){Ni <- c(Ni, NA_integer_)}
            Nf <- Ni + Ncards - 1L
            ## ## we choose a flatter hyperprior for ordinal variates
            ## we choose a Hadamard-like hyperprior for nominal variates
            Nalpha0 <- if(Nprior == 'Hadamard') {
                as.vector(unlist(sapply(Ncards, function(acard){
                    rep(1 / acard, acard)
                })))
            } else {
                as.vector(unlist(sapply(Ncards, function(acard){
                    rep(1, acard)
                })))
            }
            ##
            list(
                Nn = vn$N,
                ## Ncards = Ncards,
                Ni = Ni,
                Nf = Nf,
                Nalpha0 = Nalpha0
            )
        },
        if (vn$B > 0) { # binary
            list(
                Bn = vn$B,
                Bshapelo = rep(Bshapelo, 1),
                Bshapehi = rep(Bshapehi, 1)
            )
        },
        if(isTRUE(avoidzeroW)){
            list(epsd = .Machine$double.xmin)
        }
    ) # End constants

    saveRDS(constants,
        file = file.path(dirname, paste0('___constants', dashnameroot, '.rds')))

#### DATAPOINTS
    ## for each list element: rows: data; cols: variates of that element
    datapoints <- c(
        if (vn$R > 0) { # continuous open domain
            list(
                Rdata = as.matrix(vtransform(data[, vnames$R, drop = FALSE],
                    auxmetadata = auxmetadata,
                    Rout = 'normalized', logjacobianOr = NULL))
            )
        },
        if (vn$C > 0) { # continuous closed domain
            list(
                Caux = as.matrix(vtransform(data[, vnames$C, drop = FALSE],
                    auxmetadata = auxmetadata,
                    Cout = 'aux', logjacobianOr = NULL)),
                Clat = as.matrix(vtransform(data[, vnames$C, drop = FALSE],
                    auxmetadata = auxmetadata,
                    Cout = 'lat', logjacobianOr = NULL))
            )
        },
        if (vn$D > 0) { # continuous rounded
            list(
                Daux = as.matrix(vtransform(data[, vnames$D, drop = FALSE],
                    auxmetadata = auxmetadata,
                    Dout = 'aux', logjacobianOr = NULL))
            )
        },
        ## if (vn$L > 0) { # latent
        ##     list(
        ##         Laux = as.matrix(vtransform(data[, vnames$L, drop = FALSE],
        ##             auxmetadata = auxmetadata,
        ##             Lout = 'aux'))
        ##     )
        ## },
        if (vn$O > 0) { # nominal
            list(
                Odata = as.matrix(vtransform(data[, vnames$O, drop = FALSE],
                    auxmetadata = auxmetadata,
                    Oout = 'numeric', logjacobianOr = NULL))
            )
        },
        if (vn$N > 0) { # nominal
            list(
                Ndata = as.matrix(vtransform(data[, vnames$N, drop = FALSE],
                    auxmetadata = auxmetadata,
                    Nout = 'numeric', logjacobianOr = NULL))
            )
        },
        if (vn$B > 0) { # binary
            list(
                Bdata = as.matrix(vtransform(data[, vnames$B, drop = FALSE],
                    auxmetadata = auxmetadata,
                    Bout = 'numeric', logjacobianOr = NULL))
            )
        }
    ) # End datapoints

    saveRDS(datapoints,
        file = file.path(dirname, paste0('___datapoints', dashnameroot, '.rds')))

#### Output information to user
    cat(
        'Starting Monte Carlo sampling of', nsamples, 'samples by',
        nchains, 'chains'
    )

    samplespacedims <- (vn$R * 2 + vn$C * 2 + vn$D * 2 + # means, vars
                            (if(vn$O > 0){Omaxn - vn$O}else{0}) +
                            (if(vn$N > 0){Nmaxn - vn$N}else{0}) +
                            vn$B + # independent probs
                            1) * ncomponents - 1 # component weights

    samplespacexdims <- 1 + # Alpha
        (vn$R + vn$C + vn$D) * ncomponents + # rates
        npoints + # K
        (vn$C + vn$D) * npoints * ncomponents + # latents
        sum(is.na(data)) # missing data


    cat('\nin a space of', samplespacedims,
        '(effectively', paste0(samplespacedims + samplespacexdims, ')'),
        'dimensions.\n')

    cat('Using', ncores, 'cores:',
        nsamplesperchain, 'samples per chain, max',
        minchainspercore + (coreswithextrachain > 0), 'chains per core.\n')
    cat('Requested:   ESS', round(minESS),
        '  rel.MCSE', signif(maxrelMCSE, 3), '\n')
    cat('Core logs are being saved in individual files.\n')
    cat('\nC-compiling samplers appropriate to the variates (package Nimble)\n')
    cat('this can take tens of minutes. Please wait...\n')

    ## outconmain <- file(file.path(dirname,
    ##     paste0('log', dashnameroot,
    ##         '-', 0, '.log')
    ## ), open = 'w')
    ## sink(file = outconmain, type = 'output')
    ## sink(file = outconmain, type = 'message')

    ## restoresink <- function(){
    ##     if(sink.number() > 0) {
    ##         ## Close output to log files
    ##         sink(file = NULL, type = 'output')
    ##         sink(file = NULL, type = 'message')
    ##         close(outconmain)
    ##     }
    ## }
    ## on.exit(restoresink())

    ## function  to format printing of time
    printtimediff <- function(tim) {
        paste0(signif(tim, 2), ' ', attr(tim, 'units'))
    }

#####################################################
#### BEGINNING OF FOREACH LOOP OVER CORES
#####################################################
    ## Parallel execution over cores
    ## parallel:clusterExport(cl = cl,
    ##     varlist = c(
    ##         'dirname',
    ##         'dashnameroot',
    ##         'avoidzeroW',
    ##         'constants',
    ##         'initmethod'
    ##     ),
    ##     envir = .GlobalEnv)

## Total computation time: 29 secs 
## Average preparation & finalization time: 25 secs 
## Average Monte Carlo time per chain: 2.3 secs 
## Max total memory used: approx 750 MB
## Max memory used per core: approx 370 MB

    chaininfo <- parallel::parLapply(
        cl = cl,
        X = 1:ncores,
        fun = workerfun,
        ## arguments to workerfun
        dirname = dirname,
        dashnameroot = dashnameroot,
        avoidzeroW = avoidzeroW,
        initmethod = initmethod,
        constants = constants,
        datapoints = datapoints,
        vn = vn,
        showAlphatraces = showAlphatraces,
        Alphatoslice = Alphatoslice,
        Ktoslice = Ktoslice,
        RWtoslice = RWtoslice,
        changeSamplerOrder = changeSamplerOrder,
        minchainspercore = minchainspercore,
        coreswithextrachain = coreswithextrachain,
        nchains = nchains,
        maxhours = maxhours,
        timestart0 = timestart0,
        showsamplertimes = showsamplertimes,
        startupMCiterations = startupMCiterations,
        maxMCiterations = maxMCiterations,
        showKtraces = showKtraces,
        ncomponents = ncomponents,
        plottraces = plottraces,
        Qlo = Qlo,
        Qhi = Qhi,
        Qerror = Qerror,
        minESS = minESS,
        initES = initES,
        nsamplesperchain = nsamplesperchain,
        minMCiterations = minMCiterations,
        printtimediff = printtimediff,
        ##
        chunk.size = NULL
        )

    chaininfo <- do.call(rbind, chaininfo)
    ## restore output to std
    ## flush(outconmain)
    ## sink(file = NULL, type = 'output')
    ## sink(file = NULL, type = 'message')
############################################################
#### END OF PARALLEL FOREACH OVER CORES
############################################################

    maxusedcomponents <- max(chaininfo[, 'maxusedcomponents'])
    maxiterations <- max(chaininfo[, 'maxiterations'])
    nonfinitechains <- sum(chaininfo[, 'nonfinitechains'])
    stoppedchains <- sum(chaininfo[, 'stoppedchains'])
    headertime <- as.difftime(mean(chaininfo[, 'headertime']), units = 'secs')
    MCtime <- as.difftime(sum(chaininfo[, 'MCtime'])/nchains, units = 'secs')
    maxusedmem <- max(chaininfo[, 'usedmem'])
    totusedmem <- sum(chaininfo[, 'usedmem'])
############################################################
#### End of all MCMC
############################################################


############################################################
#### Join chains
############################################################

    headertimestart <- Sys.time()
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
                dx <- dim(yy)[-length(dim(yy))]
                dim(temp) <- c(dx, length(temp) / prod(dx))
                temp
            },
            mc1, mc2,
            SIMPLIFY = FALSE
        )
    }

## mcsamples0 <- foreach(chainnumber = 1:nchains,
##         .combine = joinmc, .multicombine = FALSE) %do% {
##             padchainnumber <- sprintf(paste0('%0', nchar(nchains), 'i'), chainnumber)
## 
##             readRDS(file = file.path(dirname,
##                 paste0('___mcsamples', dashnameroot, '--',
##                     padchainnumber, '.rds')
##             ))
##         }

    for(chainnumber in 1:nchains){
        padchainnumber <- sprintf(paste0('%0', nchar(nchains), 'i'), chainnumber)
        if(chainnumber == 1){
            mcsamples <- readRDS(file = file.path(dirname,
                paste0('___mcsamples', dashnameroot, '--',
                    padchainnumber, '.rds')
            ))
        } else {
            mcsamples <- joinmc(mcsamples,
                readRDS(file = file.path(dirname,
                    paste0('___mcsamples', dashnameroot, '--',
                        padchainnumber, '.rds')
                ))
            )
        }
    }

    ## chainnumber <- 1
    ## padchainnumber <- sprintf(paste0('%0', nchar(nchains), 'i'), chainnumber)
    ## mcsamples <- readRDS(file = file.path(dirname,
    ##             paste0('___mcsamples', dashnameroot, '--',
    ##                 padchainnumber, '.rds')
    ## ))
    ## if(nchains > 1){
    ##     for(chainnumber in 2:nchains){
    ##         padchainnumber <- sprintf(paste0('%0', nchar(nchains), 'i'), chainnumber)
    ##         mcsamples <- joinmc(mcsamples,
    ##             readRDS(file = file.path(dirname,
    ##             paste0('___mcsamples', dashnameroot, '--',
    ##                 padchainnumber, '.rds')
    ##             ))
    ##         )
    ##     }
    ## }

#### Save all final parameters together with the aux-metadata in one file
    saveRDS(c(mcsamples,
        list(auxmetadata = auxmetadata,
            auxinfo = list(
                nchains = nchains,
                npoints = npoints,
                hyperparams = hyperparams
            ))
    ),
    file = file.path(dirname,
        paste0('learnt', dashnameroot, '.rds')
    ))

    cat('\rFinished Monte Carlo sampling.                                 \n')

    cat('Highest number of Monte Carlo iterations across chains:', maxiterations, '\n')
    cat('Highest number of used mixture components:', maxusedcomponents, '\n')
    if (maxusedcomponents > ncomponents - 5) {
        cat('TOO MANY MIXTURE COMPONENTS USED!\nConsider',
            're-running with increased "ncomponents" parameter\n')
    }

    if (nonfinitechains > 0) {
        cat('\nNote:', nonfinitechains, 'chains had some non-finite outputs\n')
    }
    if (stoppedchains > 0) {
        cat('\nNote:', stoppedchains,
            'chains were stopped',
            'before reaching required precision\nin order',
            'to meet the required time constraints\n')
    }



############################################################
#### Final joint diagnostics
############################################################

    testdata <- readRDS(file = file.path(dirname,
        paste0('___testdata_', 0, '.rds')))
    cat('\nChecking test data\n(', paste0('#', testdata$pointsid), ')\n')

    oktraces <- util_Pcheckpoints(
        testdata = testdata,
        learnt = mcsamples
    )

    oktraces <- cbind(
        exp(rowMeans(log(oktraces), na.rm = TRUE)), # geometric mean
        oktraces
    )

    oktraces <- oktraces[apply(oktraces, 1, function(x) { all(is.finite(x)) }), ,
        drop = FALSE]

    colnames(oktraces) <- c('gmean', testdata$pointsid)

    saveRDS(oktraces, file = file.path(dirname,
        paste0('MCtraces', dashnameroot, '.rds')
    ))

    ## ##############################
    ## ## MCSE, ESS, THINNING, ETC ##
    ## ##############################

    N <- nrow(oktraces)

    jointdiagn <- apply(oktraces, 2, function(atrace) {
### same as within cores

        ## quantiles to monitor
        Xlo <- quantile(x = atrace, probs = Qlo,
            na.rm = FALSE, names = FALSE, type = 6)
        Xhi <- quantile(x = atrace, probs = Qhi,
            na.rm = FALSE, names = FALSE, type = 6)
        ## quantile width
        width <- Xhi - Xlo

        ## CIs for lower and upper quantiles
        temp <- funMCEQ(x = atrace, prob = c(Qlo, Qhi), Qpair = Qerror)
        wXlo <- temp[2, 1] - temp[1, 1]
        wXhi <- temp[2, 2] - temp[1, 2]

        ## Transform samples to normalized ranks, as in Vehtari et al. 2021
        essnrmean <- funESS3(qnorm(
        (rank(atrace, na.last = NA, ties.method = 'average') -
             0.5) / N
        ))

        ## We check: relative error of quantiles and ess of norm-rank-mean
        relmcse <- c(1 / sqrt(essnrmean), wXlo / width, wXhi / width)

        autothinning <- N * max(relmcse)^2

### same as within cores

        c(max(relmcse[-1]), essnrmean, autothinning, width)
    })

    ## Output available diagnostics
    toprint <- list(
        'rel. CI error' = jointdiagn[1, ],
        'ESS' = jointdiagn[2, ],
        'needed thinning' = jointdiagn[3, ],
        'average' = colMeans(oktraces),
        'width' = jointdiagn[4, ]
    )

####
    cat('\n')
    for(i in names(toprint)) {
        thisdiagn <- toprint[[i]]
        cat(paste0('\n', i, ':'),
            if(length(thisdiagn) > 1){
                paste(signif(range(thisdiagn), 3), collapse = ' to ')
            } else {
                signif(thisdiagn, 3)
            }
        )
    }


    ## Plot various info and traces
    cat('\nPlotting final Monte Carlo traces and marginal samples...\n')

    ##
    ## colpalette <- seq_len(ncol(oktraces))
    ## names(colpalette) <- colnames(oktraces)
    ## sink(file = outconmain, type = 'output')
    ## sink(file = outconmain, type = 'message')
##    suppressMessages({
    graphics.off()
    pdf(file.path(dirname,
        paste0('MCtraces', dashnameroot, '.pdf')
    ), height = 8.27, width = 11.69)
    ## Traces of likelihood and cond. probabilities
    for (avar in 1:ncol(oktraces)) {
        ## Do not join separate chains in the plot
        division <- (if(N > nchains) nchains else 1)
        flexiplot(
            y = 10*log10(oktraces[is.finite(oktraces[, avar]), avar]),
            ## x = matrix(seq_len(nsamples), ncol = division),
            ## y = matrix(10*log10(oktraces[, avar]), ncol = division),
            type = 'l', lty = 1, col = 1,
            ## col = 1:6, # to evidence consecutive chains
            main = paste0('#', colnames(oktraces)[avar], ': ',
                paste(
                    names(toprint),
                    sapply(toprint, function(xx){
                        signif(xx[avar], 3)
                    }),
                    collapse = ' | ', sep = ': '
                )),
            cex.main = 1.25,
            ylab = paste0('log_F(#',
                colnames(oktraces)[avar],
                ')/dHart'),
            xlab = 'sample', family = family
            ## mar = c(NA, 6, NA, NA)
        )
        abline(v = seq(from = 1, to = nsamples, by = nsamplesperchain)[-1] - 0.5,
            lty = 2, col = 2, lwd = 1)
    }

    ## outcon <- file(file.path(dirname, 'log-1.log'), open = 'a')
    ## sink(file = outcon, type = 'output')
    ## sink(file = outcon, type = 'message')
    ## cat('Plotting marginal samples.\n')
    plotFsamples(
        filename = file.path(dirname,
            paste0('plotsamples_learnt', dashnameroot)),
        learnt = c(mcsamples, list(auxmetadata = auxmetadata)),
        data = data,
        plotvariability = 'samples',
        nFsamples = showsamples, plotprobability = TRUE,
        datahistogram = TRUE, datascatter = TRUE,
        parallel = NULL, silent = TRUE
    )

    ## cat('Plotting marginal samples with quantiles.\n')
    plotFsamples(
        filename = file.path(dirname,
            paste0('plotquantiles_learnt', dashnameroot)),
        learnt = c(mcsamples, list(auxmetadata = auxmetadata)),
        data = data,
        plotvariability = 'quantiles',
        nFsamples = plotDisplayedQuantiles, plotprobability = TRUE,
        datahistogram = TRUE, datascatter = TRUE,
        parallel = NULL, silent = TRUE
    )
##})
    ## restore output to std
    ## flush(outconmain)
    ## sink(file = NULL, type = 'output')
    ## sink(file = NULL, type = 'message')
    ## close(outconmain)

    ## Close connections
    ## invisible(foreach(acore = 1:ncores) %dochains% {
    ##     outcon <- file(file.path(dirname,
    ##         paste0('log', dashnameroot,
    ##             '-', acore, '.log')
    ##     ), open = 'a')
    ##     sink(file = outcon, type = 'output')
    ##     sink(file = outcon, type = 'message')
    ##     flush(outcon)
    ##     sink(file = NULL, type = 'output')
    ##     sink(file = NULL, type = 'message')
    ##     close(outcon)
    ## })


    totalfinaltime <- difftime(Sys.time(), timestart0, units = 'auto')
    cat('\nTotal computation time:', printtimediff(totalfinaltime), '\n')
    cat('Average preparation & finalization time:',
        printtimediff(
            difftime(Sys.time() + headertime, headertimestart, units = 'auto')
        ), '\n')
    cat('Average Monte Carlo time per chain:',
        printtimediff(
            difftime(headertimestart + MCtime, headertimestart, units = 'auto')
        ), '\n')
    cat('Max total memory used: approx', signif(totusedmem, 2), 'MB\n')
    cat('Max memory used per core: approx', signif(maxusedmem, 2), 'MB\n')

    ## if (exists('cl')) {
    ##     cat('\nClosing connections to cores.\n')
    ##     foreach::registerDoSEQ()
    ##     parallel::stopCluster(cl)
    ## }


#### Remove partial files if required
    ## Partial files are identified by at last three initial underscores "___"
    ## Should we leave the plots of partial traces?
    ## maybe create an additional argument to let the user decide?
    if (cleanup) {
        cat('\nRemoving temporary output files.\n')
        file.remove(dir(dirname,
            pattern = paste0('^___.*\\..*$'),
            full.names = TRUE
        ))
    }

#### WHAT MIGHT BE ADDED:
    ## - when 'showKtraces' is true:
    ## a histogram over number of components over all chains
    ## (for the moment there's one plot per chain)

    cat('Finished.\n')

    ## What should we output? how about the full name of the output dir?
    if (is.character(output) && output == 'directory') {
        dirname
    } else if (is.character(output) && output == 'learnt') {
        readRDS(file.path(dirname,
            paste0('learnt', dashnameroot, '.rds')
        ))
    }
}

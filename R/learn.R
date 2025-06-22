#' Monte Carlo computation of posterior probability distribution conditional on given data
#'
#' @param data A dataset, given as a \code{\link[base]{data.frame}}
#' or as a file path to a csv file.
#' @param metadata A metadata object, given either as a data.frame object,
#' or as a file pa to a csv file.
#' @param auxdata A larger dataset, given as a data.frame
#'   or as a file path to a csv file. Such a dataset
#'   would be too many to use in the Monte Carlo sampling,
#'   but can be used to calculate hyperparameters.
#' @param outputdir Character: path to folder where the output should be saved. If omitted, a directory is created that has the same name as the data file but with suffix "`_output_`".
#' @param nsamples Integer: number of desired Monte Carlo samples. Default 3600.
#' @param nchains Integer: number of Monte Carlo chains. Default 60.
#' @param nsamplesperchain Integer: number of Monte Carlo samples per chain.
#' @param parallel Logical or `NULL` or positive integer: `TRUE`: use roughly half of available cores; `FALSE`: use serial computation; `NULL`: don't do anything (use pre-registered condition); integer: use this many cores. Default `NULL`
#' @param seed Integer: use this seed for the random number generator.
#'   If missing or `NULL` (default), do not set the seed.
#' @param cleanup Logical: remove diagnostic files at the end of the computation?
#'   Default `TRUE`.
#' @param appendtimestamp Logical: append a timestamp to the name of
#'   the output directory `outputdir`? Default `TRUE`.
#' @param appendinfo Logical: append information about dataset and Monte Carlo
#'   parameters to the name of the output directory `outputdir`? Default `TRUE`.
#' @param output Character: if `'directory'`, return the output directory name
#'   as `VALUE`; if string `'learnt'`, return the `'learnt'` object
#'   containing the parameters obtained from the Monte Carlo computation.
#'   Any other value: `VALUE` is `NULL`.
#' @param subsampledata Integer: use only a subset of this many datapoints for
#'   the Monte Carlo computation.
#' @param startupMCiterations Integer: number of initial (burn-in)
#'   Monte Carlo iterations. Default 3600.
#' @param minMCiterations Integer: minimum number of Monte Carlo iterations
#'   to be done. Default 0.
#' @param maxMCiterations Integer: Do at most this many Monte Carlo iterations.
#'   Default `Inf`.
#' @param maxhours Numeric: approximate time limit, in hours, for the
#'   Monte Carlo computation to last. Default `Inf`.
#' @param ncheckpoints Integer: number of datapoints to use
#'   for checking when the Monte Carlo computation should end.
#'   If NULL (default), this is equal to number of variates + 2.
#'   If Inf, use all datapoints.
#' @param relerror Numeric: desired maximal relative error of calculated probabilities
#'   with respect to their variability with new data.
#' @param prior Logical: Calculate the prior distribution?
#' @param thinning Integer: thin out the Monte Carlo samples by this value.
#'   If NULL (default): let the diagnostics decide the thinning value.
#' @param plottraces Logical: save plots of the Monte Carlo traces
#'   of diagnostic values? Default `TRUE`.
#' @param showKtraces Logical: save plots of the Monte Carlo traces
#'   of the K parameter? Default `FALSE`.
#' @param showAlphatraces Logical: save plots of the Monte Carlo traces
#'   of the Alpha parameter? Default `FALSE`.
#' @param hyperparams List: hyperparameters of the prior.
#'
#' @return Name of directory containing output files, or learnt object,
#'   or `NULL`, depending on argument `output`.
#'
#' @import parallel foreach doParallel doRNG nimble
#'
#' @export
learn <- function(
    data,
    metadata,
    auxdata = NULL,
    outputdir = NULL,
    nsamples = 3600,
    nchains = 60,
    nsamplesperchain = 60,
    parallel = NULL,
    seed = NULL,
    cleanup = TRUE,
    appendtimestamp = TRUE,
    appendinfo = TRUE,
    output = 'directory',
    subsampledata = NULL,
    startupMCiterations = 3600,
    minMCiterations = 0,
    maxMCiterations = +Inf,
    maxhours = +Inf,
    ncheckpoints = NULL,
    relerror = 0.05,
    prior = missing(data) || is.null(data),
    thinning = NULL,
    plottraces = TRUE,
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
        Lshapelo = 0.5,
        Lshapehi = 0.5,
        Lvarm1 = 3^2,
        Bshapelo = 1,
        Bshapehi = 1,
        Dthreshold = 1,
        tscalefactor = 2,
        initmethod = 'allcentre'
        ## precluster, prior, allcentre
    )
) {

#### Start timer
    timestart0 <- Sys.time()

    cat('\n') # make sure possible error messages start on new line

    ## Set the RNG seed if given by user, or if no seed already exists
    if (!is.null(seed) || !missing(seed) || !exists('.Random.seed')) {
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
        nsamples <- nchains * nsamplesperchain
    } else if (missing(nchains) && !missing(nsamplesperchain) &&
                   !missing(nsamples)){
        nchains <- ceiling(nsamples / nsamplesperchain)
        if(nsamples != nchains * nsamplesperchain){
            nsamples <- nchains * nsamplesperchain
            cat('Increasing number of samples to', nsamples,
                'to comply with given "nsamplesperchain"\n')
        }
    } else if (!missing(nchains) && missing(nsamplesperchain) &&
                   !missing(nsamples)){
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
    if (isTRUE(parallel)) {
        ## user wants us to register a parallel backend
        ## and to choose number of cores
        ncores <- max(1,
            min(nchains, floor(parallel::detectCores() / 2)))
        cl <- parallel::makeCluster(ncores)
        doParallel::registerDoParallel(cl)
        closeexit <- TRUE
        cat('Registered', foreach::getDoParName(),
            'with', foreach::getDoParWorkers(), 'workers\n')
    } else if (isFALSE(parallel)) {
        ## user wants us not to use parallel cores
        ncores <- 1
        foreach::registerDoSEQ()
    } else if (is.null(parallel)) {
        ## user wants us not to do anything
        ncores <- foreach::getDoParWorkers()
    } else if (is.finite(parallel) && parallel >= 1) {
        ## user wants us to register 'parallal' # of cores
        ncores <- min(nchains, parallel)
        cl <- parallel::makeCluster(ncores)
        doParallel::registerDoParallel(cl)
        closeexit <- TRUE
        cat('Registered', foreach::getDoParName(),
            'with', foreach::getDoParWorkers(), 'workers\n')
    } else {
        stop("Unknown value of argument 'parallel'")
    }

    ## Close parallel connections if any were opened
    if(closeexit) {
        closecoresonexit <- function(){
            cat('\nClosing connections to cores.\n')
            foreach::registerDoSEQ()
            parallel::stopCluster(cl)
            ## parallel::setDefaultCluster(NULL)
            env <- foreach:::.foreachGlobals
            rm(list=ls(name=env), pos=env)
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
    ## Done now in case ncores was reduced because of nchains argument
    if (ncores < 1) {
        `%dochains%` <- `%do%`
    } else {
        `%dochains%` <- `%dorng%`
    }

    ## Make sure 'startupMCiterations' is at least 2
    startupMCiterations <- max(2, startupMCiterations)

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
        Dthreshold = hyperparams$Dthreshold,
        tscalefactor = hyperparams$tscalefactor
        )
    ## print(auxmetadata) # for debugging

    cat('\nLearning: ', npoints, 'datapoints, ',
        nrow(auxmetadata), 'variates\n')

#### Output-folder setup
    if (is.null(outputdir)) {
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
            ## '-K', ncomponents, # unimportant for user
            '_smp', nsamples)
    }
    dirname <- paste0(outputdir, suffix)
    ##
                                        ## Create output directory if it does not exist
    dir.create(dirname, showWarnings = FALSE)
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
        file = file.path(dirname, paste0('_auxmetadata', dashnameroot, '.rds')))

    ## Save initial RNG seed in case needed by user
    saveRDS(currentseed,
        file = file.path(dirname, paste0('rng_seed', dashnameroot, '.rds')))

#### number ofcheckpoints for Monte Carlo stopping rule
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
            testdata <- as.matrix(vtransform(
                data[pointsid, , drop = FALSE],
                auxmetadata = auxmetadata,
                Rout = 'normalized',
                Cout = 'boundisinf',
                Dout = 'normalized',
                Oout = 'numeric',
                Nout = 'numeric',
                Bout = 'numeric',
                logjacobianOr = NULL))
            rownames(testdata) <- pointsid
            saveRDS(testdata,
                file = file.path(dirname, paste0('_testdata_', achain, '.rds')))
        }
    }
    rm(testdata)
    rm(pointsid)

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
                                sample(1:auxmetadata[ii, 'Nvalues'], ncheckpoints,
                                    replace = TRUE)
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
                    file = file.path(dirname, paste0('_testdata_', achain, '.rds')))
            }
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
    ## assign the hyperparameter values to corresponding objects
    ## ncomponents <- hyperparams$ncomponents
    ## minalpha <- hyperparams$minalpha
    ## maxalpha <- hyperparams$maxalpha
    ## byalpha <- hyperparams$byalpha
    ## Rshapelo <- hyperparams$Rshapelo
    ## Rshapehi <- hyperparams$Rshapehi
    ## Rvarm1 <- hyperparams$Rvarm1
    ## Cshapelo <- hyperparams$Cshapelo
    ## Cshapehi <- hyperparams$Cshapehi
    ## Cvarm1 <- hyperparams$Cvarm1
    ## Dshapelo <- hyperparams$Dshapelo
    ## Dshapehi <- hyperparams$Dshapehi
    ## Dvarm1 <- hyperparams$Dvarm1
    ## Lshapelo <- hyperparams$Lshapelo
    ## Lshapehi <- hyperparams$Lshapehi
    ## Lvarm1 <- hyperparams$Lvarm1
    ## Bshapelo <- hyperparams$Bshapelo
    ## Bshapehi <- hyperparams$Bshapehi
    for(aname in names(hyperparams)){
        assign(aname, hyperparams[[aname]])
    }

    nalpha <- length(seq(minalpha, maxalpha, by = byalpha))
    npoints <- nrow(data)

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

    ## function printtime to format printing of time
    printtimediff <- function(tim) {
        paste0(signif(tim, 2), ' ', attr(tim, 'units'))
    }
    printtimeend <- function(tim) {
        format(Sys.time() + tim, format='%Y-%m-%d %H:%M')
    }
    ## We need to send some messages to the log files, others to the user.
    ## This is done by changing output sink:
    print2user <- function(msg, outcon) {
        sink(file = NULL, type = 'message')
        message(msg, appendLF = FALSE)
        flush.console()
        sink(file = outcon, type = 'message')
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
                Dn = vn$D, # This indexing variable is needed internally
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
        ##         Ln = vn$L, # This indexing variable is needed internally
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
            Omaxn <- max(auxmetadata[auxmetadata$mcmctype == 'O', 'Nvalues'])
            Oalpha0 <- matrix(1e-100, nrow = vn$O, ncol = Omaxn)
            for (avar in seq_along(vnames$O)) {
                nvalues <- auxmetadata[auxmetadata$name == vnames$O[avar], 'Nvalues']
                ## ## we choose a flatter hyperprior for ordinal variates
                ## we choose a Hadamard-like hyperprior for nominal variates
                Oalpha0[avar, 1:nvalues] <- 1/nvalues
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

    samplespacedims <- vn$R * 2 * ncomponents +
        vn$C * 2 * ncomponents +
        vn$D * 2 * ncomponents +
        sum(apply(Oalpha0, 1, function(x) sum(x > 2e-17) - 1)) * ncomponents +
        sum(apply(Nalpha0, 1, function(x) sum(x > 2e-17) - 1)) * ncomponents +
        sum(apply(Nalpha0, 1, function(x) sum(x > 2e-17) - 1)) +
        vn$B * ncomponents +
        ncomponents - 1
    samplespacexdims <- 1 + # Alpha
        vn$R * ncomponents + # Rrate
        vn$C * ncomponents + # Crate
        vn$D * ncomponents + # Drate
        1 + # K
        vn$C * npoints * ncomponents + # latent C
        vn$D * npoints * ncomponents + # latent D
        sum(is.na(data)) # missing data


    cat('\nin a space of', samplespacedims,
        '(effectively', paste0(samplespacedims + samplespacexdims, ')'),
        'dimensions.\n')

    cat('Using', ncores, 'cores:',
        nsamplesperchain, 'samples per chain, max',
        minchainspercore + (coreswithextrachain > 0), 'chains per core.\n')
    cat('Core logs are being saved in individual files.\n')
    cat('\nC-compiling samplers appropriate to the variates (package Nimble)\n')
    cat('this can take tens of minutes. Please wait...\r')

    ## ## Needed if method F. for K initialization is used
    ## Ksample <- sample(0:1, 1)

#####################################################
#### BEGINNING OF FOREACH LOOP OVER CORES
#####################################################
    ## Parallel execution over cores

    chaininfo <- foreach(acore = 1:ncores,
        .combine = rbind,
        .inorder = FALSE,
        ##.packages = c('predict'),
        .noexport = c('data')
    ) %dochains% {

        ## Create log file
        ## Redirect diagnostics and service messages there
        outcon <- file(file.path(dirname,
            paste0('log', dashnameroot,
                '-', acore, '.log')
        ), open = 'w')
        sink(file = outcon, type = 'output')
        sink(file = outcon, type = 'message')
        if(acore < 0){
            closecons <- function(){
                ## Close output to log files
                sink(file = NULL, type = 'output')
                sink(file = NULL, type = 'message')
                close(outcon)
            }
            on.exit(closecons())
        }
        usedmem <- sum(gc()[,6])

        ## Timer
        headertimestart <- Sys.time()

        cat('Log core', acore)
        cat(' - Current time:',
            strftime(as.POSIXlt(headertimestart), '%Y-%m-%d %H:%M:%S'))
        cat('\n')

        suppressPackageStartupMessages(require('nimble'))
        ## requireNamespace("nimble", quietly = TRUE)
        ##library('nimble')

#### COMPONENT REPRESENTATION OF FREQUENCY SPACE
#### Dirichlet-process mixture of product-kernels

        ## hierarchical probability structure
        finitemix <- nimbleCode({
            ## Component weights
            Alpha ~ dcat(prob = probalpha0[1:nalpha])
            alphas[1:ncomponents] <- dirchalphas[1:ncomponents] * alphabase^Alpha
            W[1:ncomponents] ~ ddirch(alpha = alphas[1:ncomponents])
            W0[1:ncomponents] <- W[1:ncomponents] + 1e-100

            ## Probability density for the parameters of the components
            for (k in 1:ncomponents) {
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
                ## if (vn$L > 0) { # latent
                ##     for (v in 1:Ln) {
                ##         Lmean[v, k] ~ dnorm(mean = Lmean1, var = Lvarm1)
                ##         Lrate[v, k] ~ dinvgamma(shape = Lshapehi, rate = Lvar1)
                ##         Lvar[v, k] ~ dinvgamma(shape = Lshapelo, rate = Lrate[v, k])
                ##     }
                ## }
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
                K[d] ~ dcat(prob = W0[1:ncomponents])
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
                ## if (vn$L > 0) { # latent
                ##     for (v in 1:Ln) {
                ##         Laux[d, v] ~ dconstraint(Llat[d, v] >= Lleft[d, v] &
                ##                                  Llat[d, v] < Lright[d, v])
                ##         Llat[d, v] ~ dnorm(mean = Lmean[v, K[d]], var = Lvar[v, K[d]])
                ##     }
                ## }
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
                    for (v in 1:Bn) { # Bprob is the probability that Bdata=1
                        Bdata[d, v] ~ dbern(prob = Bprob[v, K[d]])
                    }
                }
            }
        }) # end finitemix NimbleCode


#### INITIAL-VALUE FUNCTION
        if(initmethod == 'precluster'){
            ## pre-clustering, k-mean style
            initsfn <- function() {
                ## Create components centres
                ## distance function
                ## NB: all variances will be initialized to 1
                lpnorm <- function(xx){abs(xx)}
                distances <- matrix(0, nrow = npoints, ncol = ncomponents)
                if (vn$R > 0) { # continuous open domain
                    Rmeans <- matrix(rnorm(
                        n = vn$R * ncomponents,
                        mean = constants$Rmean1,
                        sd = sqrt(constants$Rvarm1)
                    ), nrow = vn$R, ncol = ncomponents)
                    ## distances from datapoints
                    distances <- distances + apply(Rmeans, 2, function(ameans){
                        colSums(lpnorm(t(datapoints$Rdata) - ameans), na.rm = TRUE)
                    })
                }
                if (vn$C > 0) { # continuous closed domain
                    Cmeans <- matrix(rnorm(
                        n = vn$C * ncomponents,
                        mean = constants$Cmean1,
                        sd = sqrt(constants$Cvarm1)
                    ), nrow = vn$C, ncol = ncomponents)
                    ## distances from datapoints
                    distances <- distances + apply(Cmeans, 2, function(ameans){
                        colSums(lpnorm(t(datapoints$Clat) - ameans), na.rm = TRUE)
                    })
                }
                if (vn$D > 0) { # discrete
                    Dmeans <- matrix(rnorm(
                        n = vn$D * ncomponents,
                        mean = constants$Dmean1,
                        sd = sqrt(constants$Dvarm1)
                    ), nrow = vn$D, ncol = ncomponents)
                    ## distances from datapoints
                    distances <- distances + apply(Dmeans, 2, function(ameans){
                        colSums(lpnorm(t(constants$Dlatinit) - ameans), na.rm = TRUE)
                    })
                }
                ## if (vn$L > 0) { # 
                ##     Lmeans <- matrix(rnorm(
                ##         n = vn$L * ncomponents,
                ##         mean = constants$Lmean1,
                ##         sd = sqrt(constants$Lvarm1)
                ##     ), nrow = vn$L, ncol = ncomponents)
                ##     ## distances from datapoints
                ##     distances <- distances + apply(Lmeans, 2, function(ameans){
                ##         colSums(lpnorm(t(constants$Llatinit) - ameans), na.rm = TRUE)
                ##     })
                ## }
                ## if (vn$B > 0) {
                ##     Bprobs <- matrix(rbeta(
                ##         n = vn$B * ncomponents,
                ##         shape1 = Bshapelo,
                ##         shape2 = Bshapehi,
                ##         ), nrow = vn$B, ncol = ncomponents)
                ##     ## distances from datapoints
                ##     distances <- distances + apply(Bprobs, 2, function(ameans){
                ##         colSums(lpnorm(t(datapoints$Bdata) - ameans), na.rm = TRUE)
                ##     })
                ## }

                ## assign datapoints to component with closest centre
                K <- apply(distances, 1, which.min)
                occupied <- unique(K)

                ## recalculate components centres according to their points
                if (vn$R > 0) {
                    Rmeans[, occupied] <- sapply(occupied, function(acomponent){
                        colMeans(datapoints$Rdata[which(K == acomponent), , drop = FALSE],
                            na.rm = TRUE)
                    })
                    Rmeans[, -occupied] <- 0
                }
                if (vn$C > 0) {
                    Cmeans[, occupied] <- sapply(occupied, function(acomponent){
                        colMeans(datapoints$Clat[which(K == acomponent), , drop = FALSE],
                            na.rm = TRUE)
                    })
                    Cmeans[, -occupied] <- 0
                }
                if (vn$D > 0) {
                    Dmeans[, occupied] <- sapply(occupied, function(acomponent){
                        colMeans(constants$Dlatinit[which(K == acomponent), , drop = FALSE],
                            na.rm = TRUE)
                    })
                    Dmeans[, -occupied] <- 0
                }
                ## if (vn$L > 0) { # continuous open domain
                ##     Lmeans[, occupied] <- sapply(occupied, function(acomponent){
                ##         colMeans(constants$Llatinit[which(K == acomponent), , drop = FALSE],
                ##             na.rm = TRUE)
                ##     })
                ##     Lmeans[, -occupied] <- 0
                ## }
                ## if (vn$B > 0) {
                ##     Bprobs[, occupied] <- sapply(occupied, function(acomponent){
                ##         colMeans(datapoints$Bdata[which(K == acomponent), , drop = FALSE],
                ##             na.rm = TRUE)
                ##     })
                ##     Bprobs[, -occupied] <- 0.5
                ## }
                ## Alpha <- sample(1:nalpha, 1, prob = constants$probalpha0, replace = TRUE)
                ## W <- c(rep(rempoints, minpoints), rep(1, ncomponents - minpoints))
                ## W <- W/sum(W)

                outlist <- list(
                    Alpha = round(nalpha/2),
                    W = rep(1/ncomponents, ncomponents),
                    ## ## Assign every point to the closest component centre
                    K = K
                    ## ## Other assignment methods:
                    ## ## A. assign all points to an unsystematically chosen component
                    ## K = rep(sample(rep(which(W > 0), 2), 1), nepoints)
                    ## ## B. distribute points unsystematically among components
                    ## K = sample(rep(which(W > 0), 2), npoints, replace = TRUE)
                    ## ## or:
                    ## ## C. assign all points to the most probable component
                    ## K = rep(which.max(W), npoints)
                    ## ## or:
                    ## ## D. assign all points to the least probable component
                    ## K = rep(which(W == min(W[W > 0]))[1], npoints)
                    ## ## or:
                    ## ## E. distribute points unsystematically among M=2 components
                    ## K = sample(sample(rep(which(W > 0), 2), 2, replace = TRUE),
                    ##           npoints, replace = TRUE)
                    ## ## F. mix methods A. and B.
                    ## K = (if(achain %% 2 == Ksample) {
                    ##        ## ## assign all points to an unsystematically chosen component
                    ##        rep(sample(rep(which(W > 0), 2), 1), npoints)
                    ##      } else {
                    ##        ## distribute points unsystematically among components
                    ##        sample(rep(which(W > 0), 2), npoints, replace = TRUE)
                    ##      })
                )
                ##
                if (vn$R > 0) { # continuous open domain
                    outlist <- c(
                        outlist,
                        list(
                            Rmean = Rmeans,
                            Rrate = matrix(
                                nimble::qinvgamma(p = 0.5,
                                    shape = constants$Rshapehi,
                                    rate = constants$Rvar1),
                                nrow = vn$R, ncol = ncomponents
                            ),
                            Rvar = matrix(1,
                                nrow = vn$R, ncol = ncomponents)
                        )
                    )
                }
                if (vn$C > 0) { # continuous closed domain
                    outlist <- c(
                        outlist,
                        list(
                            Cmean = Cmeans,
                            Crate = matrix(
                                nimble::qinvgamma(p = 0.5,
                                    shape = constants$Cshapehi,
                                    rate = constants$Cvar1),
                                nrow = vn$C, ncol = ncomponents
                            ),
                            Cvar = matrix(1,
                                nrow = vn$C, ncol = ncomponents),
                            ## for data with boundary values
                            Clat = constants$Clatinit
                            ## Clat = vtransform(data[, vnames$C, with = FALSE],
                            ##   auxmetadata, Cout = 'init')
                        )
                    )
                }
                if (vn$D > 0) { # continuous rounded
                    outlist <- c(
                        outlist,
                        list(
                            Dmean = Dmeans,
                            Drate = matrix(
                                nimble::qinvgamma(p = 0.5,
                                    shape = constants$Dshapehi,
                                    rate = constants$Dvar1),
                                nrow = vn$D, ncol = ncomponents
                            ),
                            Dvar = matrix(1,
                                nrow = vn$D, ncol = ncomponents),
                            ## for data with boundary values
                            Dlat = constants$Dlatinit
                            ## Dlat = vtransform(data[, vnames$D, with = FALSE],
                            ##   auxmetadata, Dout = 'init')
                        )
                    )
                }
                ## if (vn$L > 0) { # latent
                ##     outlist <- c(
                ##         outlist,
                ##         list(
                ##             Lmean = Lmeans,
                ##             Lrate = matrix(
                ##                 nimble::qinvgamma(p = 0.5,
                ##                     shape = constants$Lshapehi,
                ##                     rate = constants$Lvar1),
                ##                 nrow = vn$L, ncol = ncomponents
                ##             ),
                ##             Lvar = matrix(1,
                ##                 nrow = vn$L, ncol = ncomponents),
                ##             ## for data with boundary values
                ##             Llat = constants$Llatinit
                ##             ## Llat = vtransform(data[, vnames$L, with = FALSE],
                ##             ##   auxmetadata, Lout = 'init')
                ##         )
                ##     )
                ## }
                if (vn$O > 0) { # ordinal
                    outlist <- c(
                        outlist,
                        list(
                            Oprob = aperm(array(sapply(1:vn$O, function(avar) {
                                sapply(1:ncomponents, function(aclus) {
                                    Oalpha0[avar, ]/sum(Oalpha0[avar, ])
                                    ## nimble::rdirch(n = 1, alpha = Oalpha0[avar, ])
                                })
                            }), dim = c(Omaxn, ncomponents, vn$O)))
                        )
                    )
                }
                if (vn$N > 0) { # nominal
                    outlist <- c(
                        outlist,
                        list(
                            Nprob = aperm(array(sapply(1:vn$N, function(avar) {
                                sapply(1:ncomponents, function(aclus) {
                                    Nalpha0[avar, ]/sum(Nalpha0[avar, ])
                                    ## nimble::rdirch(n = 1, alpha = Nalpha0[avar, ])
                                })
                            }), dim = c(Nmaxn, ncomponents, vn$N)))
                        )
                    )
                }
                if (vn$B > 0) { # binary
                    outlist <- c(
                        outlist,
                        list(
                            ## Bprob = Bprobs
                            Bprob = matrix(0.5, nrow = vn$B, ncol = ncomponents)
                        )
                    )
                }
                ##
                outlist
            } #End initsfns
        } else if(initmethod == 'prior'){
            ## values chosen from prior
            initsfn <- function() {
                Alpha <- sample(1:nalpha, 1, prob = probalpha0[1:nalpha])
                W <- nimble::rdirch(n = 1,
                    alpha = constants$dirchalphas[1:ncomponents] *
                        constants$alphabase^Alpha)
                outlist <- list(
                    Alpha = Alpha,
                    W = W,
                    K = sample(rep(which(W > 0), 2), npoints, replace = TRUE)
                )
                ##
                if (vn$R > 0) { # continuous open domain
                    Rrate <- matrix(
                        nimble::rinvgamma(
                            n = vn$R * ncomponents,
                            shape = constants$Rshapehi,
                            rate = constants$Rvar1),
                        nrow = vn$R, ncol = ncomponents
                    )
                    outlist <- c(
                        outlist,
                        list(
                            Rmean = matrix(rnorm(
                                n = vn$R * ncomponents,
                                mean = constants$Rmean1,
                                sd = sqrt(constants$Rvarm1)
                            ), nrow = vn$R, ncol = ncomponents),
                            Rrate = Rrate,
                            Rvar = matrix(
                                nimble::rinvgamma(
                                    n = vn$R * ncomponents,
                                    shape = constants$Rshapelo,
                                    rate = Rrate),
                                nrow = vn$R, ncol = ncomponents
                            )
                        )
                    )
                }
                if (vn$C > 0) { # continuous closed domain
                    Crate <- matrix(
                        nimble::rinvgamma(
                            n = vn$C * ncomponents,
                            shape = constants$Cshapehi,
                            rate = constants$Cvar1),
                        nrow = vn$C, ncol = ncomponents
                    )
                    outlist <- c(
                        outlist,
                        list(
                            Cmean = matrix(rnorm(
                                n = vn$C * ncomponents,
                                mean = constants$Cmean1,
                                sd = sqrt(constants$Cvarm1)
                            ), nrow = vn$C, ncol = ncomponents),
                            Crate = Crate,
                            Cvar = matrix(
                                nimble::rinvgamma(
                                    n = vn$C * ncomponents,
                                    shape = constants$Cshapelo,
                                    rate = Crate),
                                nrow = vn$C, ncol = ncomponents
                            ),
                            Clat = constants$Clatinit
                        )
                    )
                }
                if (vn$D > 0) { # continuous rounded
                    Drate <- matrix(
                        nimble::rinvgamma(
                            n = vn$D * ncomponents,
                            shape = constants$Dshapehi,
                            rate = constants$Dvar1),
                        nrow = vn$D, ncol = ncomponents
                    )
                    outlist <- c(
                        outlist,
                        list(
                            Dmean = matrix(rnorm(
                                n = vn$D * ncomponents,
                                mean = constants$Dmean1,
                                sd = sqrt(constants$Dvarm1)
                            ), nrow = vn$D, ncol = ncomponents),
                            Drate = Drate,
                            Dvar = matrix(
                                nimble::rinvgamma(
                                    n = vn$D * ncomponents,
                                    shape = constants$Dshapelo,
                                    rate = Drate),
                                nrow = vn$D, ncol = ncomponents
                            ),
                            Dlat = constants$Dlatinit
                        )
                    )
                }
                if (vn$O > 0) { # ordinal
                    outlist <- c(
                        outlist,
                        list(
                            Oprob = aperm(array(sapply(1:vn$O, function(avar) {
                                sapply(1:ncomponents, function(aclus) {
                                    nimble::rdirch(n = 1, alpha = Oalpha0[avar, ])
                                })
                            }), dim = c(Omaxn, ncomponents, vn$O)))
                        )
                    )
                }
                if (vn$N > 0) { # nominal
                    outlist <- c(
                        outlist,
                        list(
                            Nprob = aperm(array(sapply(1:vn$N, function(avar) {
                                sapply(1:ncomponents, function(aclus) {
                                    nimble::rdirch(n = 1, alpha = Nalpha0[avar, ])
                                })
                            }), dim = c(Nmaxn, ncomponents, vn$N)))
                        )
                    )
                }
                if (vn$B > 0) { # nominal
                    outlist <- c(
                        outlist,
                        list(
                            Bprob = matrix(
                                rbeta(n = vn$B * ncomponents,
                                    shape1 = Bshapelo, shape2 = Bshapehi),
                                nrow = vn$B, ncol = ncomponents)
                        )
                    )
                }
                ##
                outlist
            } #End initsfns
        } else { # 'allcentre'
            ## all components equal, all points to first component
            initsfn <- function() {
                outlist <- list(
                    Alpha = round(nalpha/2),
                    W = rep(1/ncomponents, ncomponents),
                    ## ## Assign every point to first component
                    K = rep(1, npoints)
                )
                ##
                if (vn$R > 0) { # continuous open domain
                    outlist <- c(
                        outlist,
                        list(
                            Rmean =  matrix(0, nrow = vn$R, ncol = ncomponents),
                            Rrate = matrix(
                                nimble::qinvgamma(p = 0.5,
                                    shape = constants$Rshapehi,
                                    rate = constants$Rvar1),
                                nrow = vn$R, ncol = ncomponents
                            ),
                            Rvar = matrix(1, nrow = vn$R, ncol = ncomponents)
                        )
                    )
                }
                if (vn$C > 0) { # continuous closed domain
                    outlist <- c(
                        outlist,
                        list(
                            Cmean = matrix(0, nrow = vn$C, ncol = ncomponents),
                            Crate = matrix(
                                nimble::qinvgamma(p = 0.5,
                                    shape = constants$Cshapehi,
                                    rate = constants$Cvar1),
                                nrow = vn$C, ncol = ncomponents
                            ),
                            Cvar = matrix(1,
                                nrow = vn$C, ncol = ncomponents),
                            ## for data with boundary values
                            Clat = constants$Clatinit
                        )
                    )
                }
                if (vn$D > 0) { # continuous rounded
                    outlist <- c(
                        outlist,
                        list(
                            Dmean = matrix(0, nrow = vn$D, ncol = ncomponents),
                            Drate = matrix(
                                nimble::qinvgamma(p = 0.5,
                                    shape = constants$Dshapehi,
                                    rate = constants$Dvar1),
                                nrow = vn$D, ncol = ncomponents
                            ),
                            Dvar = matrix(1,
                                nrow = vn$D, ncol = ncomponents),
                            ## for data with boundary values
                            Dlat = constants$Dlatinit
                        )
                    )
                }
                if (vn$O > 0) { # ordinal
                    outlist <- c(
                        outlist,
                        list(
                            Oprob = aperm(array(sapply(1:vn$O, function(avar) {
                                sapply(1:ncomponents, function(aclus) {
                                    Oalpha0[avar, ]/sum(Oalpha0[avar, ])
                                })
                            }), dim = c(Omaxn, ncomponents, vn$O)))
                        )
                    )
                }
                if (vn$N > 0) { # nominal
                    outlist <- c(
                        outlist,
                        list(
                            Nprob = aperm(array(sapply(1:vn$N, function(avar) {
                                sapply(1:ncomponents, function(aclus) {
                                    Nalpha0[avar, ]/sum(Nalpha0[avar, ])
                                })
                            }), dim = c(Nmaxn, ncomponents, vn$N)))
                        )
                    )
                }
                if (vn$B > 0) { # binary
                    outlist <- c(
                        outlist,
                        list(
                            Bprob = matrix(0.5, nrow = vn$B, ncol = ncomponents)
                        )
                    )
                }
                ##
                outlist
            } #End initsfns
        }

#################################################
#### NIMBLE SETUP
##################################################
        finitemixnimble <- nimbleModel(
            code = finitemix,
            name = 'finitemixnimble1',
            constants = constants,
            data = datapoints,
            ## dimensions = c(
            ##     if (vn$R > 0) {
            ##         list(
            ##             Rdata = c(npoints, vn$R)
            ##             )
            ##     },
            ##     if (vn$C > 0) {
            ##         list(
            ##             Caux = c(npoints, vn$C),
            ##             Clat = c(npoints, vn$C),
            ##             Cleft = c(npoints, vn$C),
            ##             Cright = c(npoints, vn$C),
            ##             Clatinit = c(npoints, vn$C)
            ##             )
            ##     },
            ##     if (vn$D > 0) {
            ##         list(
            ##             Daux = c(npoints, vn$D),
            ##             Dleft = c(npoints, vn$D),
            ##             Dright = c(npoints, vn$D),
            ##             Dlatinit = c(npoints, vn$D)
            ##             )
            ##     },
            ##     if (vn$O > 0) {
            ##         list(
            ##             Odata = c(npoints, vn$O)
            ##             )
            ##     },
            ##     if (vn$N > 0) {
            ##         list(
            ##             Ndata = c(npoints, vn$N)
            ##             )
            ##     },
            ##     if (vn$B > 0) {
            ##         list(
            ##             Bdata = c(npoints, vn$B)
            ##             )
            ##     }
            ## ),
            inits = initsfn()
        )

        Cfinitemixnimble <- compileNimble(finitemixnimble,
            showCompilerOutput = FALSE)
        usedmem <- max(usedmem, sum(gc()[,6])) #garbage collection

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
                ## if (vn$L > 0) {
                ##     c('Lmean', 'Lvar')
                ## },
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
            ## It is necessary to monitor K to see if all components were used
            ## if 'showAlphatraces' is true,
            ## then the Alpha-parameter trace is also recorded and shown
            monitors2 = c(if (showAlphatraces) { 'Alpha' },
                'K')
        )

        ## ## Uncomment to debug Nimble (in case of Nimble updates)
        ## print(confnimble$getUnsampledNodes())
        ## cat('\nEX1\n')
        ## confnimble$printSamplers(executionOrder=TRUE)

        targetslist <- sapply(confnimble$getSamplers(), function(xx) xx$target)
        nameslist <- sapply(confnimble$getSamplers(), function(xx) xx$name)
        ## cat('\nNAMESLIST', nameslist, '\n')

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
                ## ## New replacement method, didn't work in previous Nimble
                confnimble$replaceSampler(target=asampler, type='slice')
                ## ## Old replacement method:
                ## confnimble$removeSamplers(asampler)
                ## confnimble$addSampler(target = asampler, type = 'slice')
            }
        }

        ## ## Uncomment when debugging Nimble
        ## print(confnimble$getUnsampledNodes())
        ## cat('\nEX2\n')
        ## confnimble$printSamplers(executionOrder=TRUE)

#### change execution order for some variates
        if (changeSamplerOrder) {
            ## call this to do a first reordering of the samplers
            ## it places posterior-predictive nodes last
            mcsampler <- buildMCMC(confnimble)

            samplerorder <- c(
                'K',
                'W',
                'Alpha',
                if (vn$R > 0) {
                    c('Rmean', 'Rrate', 'Rvar')
                },
                if (vn$C > 0) {
                    c('Cmean', 'Crate', 'Cvar')
                },
                if (vn$D > 0) {
                    c('Dmean', 'Drate', 'Dvar')
                },
                ## if (vn$L > 0) {
                ##     c('Lmean', 'Lrate', 'Lvar')
                ## },
                if (vn$O > 0) {
                    c('Oprob')
                },
                if (vn$N > 0) {
                    c('Nprob')
                },
                if (vn$B > 0) {
                    c('Bprob')
                }
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
            ## cat('\nNEW ORDER',neworder,'\n')
            confnimble$setSamplerExecutionOrder(c(setdiff(
                confnimble$getSamplerExecutionOrder(),
                neworder
            ), neworder))
            ## cat('\nEX3\n')
            ## confnimble$printSamplers(executionOrder=TRUE)

        }

#### Compile Monte Carlo sampler
        print(confnimble)
        mcsampler <- buildMCMC(confnimble)
        ## print(confnimble$getUnsampledNodes())
        Cmcsampler <- compileNimble(mcsampler, resetFunctions = TRUE)

        cat('\nSetup time',
            printtimediff(difftime(Sys.time(), headertimestart, units = 'auto')),
            '\n')
        headertime <- difftime(Sys.time(), headertimestart, units = 'secs')


        ## Inform user that compilation is done, if core 1:
        if (acore == 1) {
            print2user(paste0('\rCompiled core ', acore, '. ',
                'Number of samplers: ',
                length(confnimble$samplerExecutionOrder), '.       \n',
                'Estimating remaining time, please be patient...'),
                outcon)
        }

        ## cat('Loop over chains')
##################################################
#### LOOP OVER CHAINS (WITHIN ONE CORE)
##################################################
        ## Start timer
        MCtimestart <- Sys.time()

        maxusedcomponents <- 0
        maxiterations <- 0
        ## keep count of chains having non-finite outputs
        nonfinitechains <- 0L
        ## keep count of chains stopped before convergence
        stoppedchains <- 0L
        usedmem <- max(usedmem, sum(gc()[,6])) #garbage collection
#### LOOP OVER CHAINS IN CORE
        nchainsperthiscore <- minchainspercore + (acore <= coreswithextrachain)
        ## print2user(paste0('\ncore ',acore,': ',nchainsperthiscore,'\n'), outcon)

        for (achain in 1:nchainsperthiscore) {

            chainnumber <- achain + minchainspercore * (acore - 1) +
                min(coreswithextrachain, acore - 1)
            padchainnumber <- sprintf(paste0('%0', nchar(nchains), 'i'), chainnumber)

            ## calculate how the remaining time is allotted to remaining chains
            timeleft <- (maxhours -
                             as.double(Sys.time() - timestart0, units = 'hours')) /
                (nchainsperthiscore - achain + 1)

            showsamplertimes0 <- showsamplertimes && (achain == 1)
            ## showAlphatraces0 <- showAlphatraces && (achain==1)
            niter <- min(startupMCiterations, maxMCiterations)
            ## ## Experimental: decrease number of iterations based on previous chain
            ## if(achain == 1){
            ##   niter <- startupMCiterations
            ## }else{
            ##  niter <- max(min(startupMCiterations,requirediter*2), 128)
            ##  }
            nitertot <- availiter <- 0L
            requirediter <- +Inf
            reset <- TRUE
            allmcsamples <- NULL
            allmcsamplesKA <- list(Alpha = NULL, K = NULL)
            flagll <- FALSE
            flagnonfinite <- FALSE
            cat(
                '\nChain #', chainnumber,
                '(chain', achain, 'of', nchainsperthiscore, 'for this core)\n'
            )
            ## Read data to be used in log-likelihood
            testdata <- readRDS(file = file.path(dirname,
                paste0('_testdata_', chainnumber, '.rds')))
            cat('\nDatapoints for testing convergence:\n',
                paste0('#', rownames(testdata)), '\n')

            ## will contain the MC traces of the test points
            traces <- matrix(NA, nrow = 0, ncol = nrow(testdata),
                dimnames = list(NULL, rownames(testdata)))


            ## Initial values for this chain
            ## random seed is taken care of by %doRNG%
            Cfinitemixnimble$setInits(initsfn())

            subiter <- 1L
#### WHILE-LOOP CONTINUING UNTIL CONVERGENCE
            while (requirediter > 0) {
                cat('\nIterations:', niter, '\n')

                ## MONTE-CARLO CALL
                ## If reporting Alpha or K traces,
                ## then save them more frequently
                ## Otherwise just save the last value
                Cmcsampler$run(
                    niter = niter,
                    thin = 1,
                    thin2 = (if (showAlphatraces || showKtraces) {
                        max(1, round(niter / ncomponentsamples))
                    } else {
                        max(2, floor(niter / 2))
                    }),
                    nburnin = 0,
                    time = showsamplertimes0,
                    reset = reset,
                    resetMV = TRUE
                )

                ## iterationAsLastIndex: See sect 7.7 of Nimble manual
                mcsamples <- as.list(Cmcsampler$mvSamples,
                    iterationAsLastIndex = TRUE)
                mcsamplesKA <- as.list(Cmcsampler$mvSamples2,
                    iterationAsLastIndex = FALSE)

                ## 'mcsamplesKA$K' contains the component identity
                ## of each training datapoint, but we only want
                ## the number of distinct components used:
                mcsamplesKA$K <- apply(mcsamplesKA$K, 1,
                    function(xx){length(unique(xx))})

                if (showAlphatraces) {
                    dim(mcsamplesKA$Alpha) <- NULL # from matrix to vector
                }

                cat('\nCurrent time:',
                    strftime(as.POSIXlt(Sys.time()), '%Y-%m-%d %H:%M:%S'))
                cat('\nMCMC time',
                    printtimediff(difftime(Sys.time(), MCtimestart, units = 'auto')),
                    '\n')

#### Remove iterations with non-finite values
                if(any(!is.finite(unlist(mcsamples)))) {
                    toRemove <- sort(unique(unlist(lapply(mcsamples, function(xx) {
                        temp <- which(!is.finite(xx), arr.ind = TRUE)
                        temp[, ncol(temp)]
                    }))))

                    cat('\nWARNING:', length(toRemove), 'NON-FINITE SAMPLES\n')
                    ##
                    flagnonfinite <- TRUE
                    nonfinitechains <- TRUE
                    saveRDS(mcsamples, file = file.path(dirname,
                        paste0('NONFINITEmcsamples',
                            dashnameroot, '--', padchainnumber,
                            '_', achain, '-',
                            acore, '-i', nitertot, '.rds')
                    ))
                    if (length(toRemove) == ncol(mcsamples$W)) {
                        cat('\n...TOO MANY NON-FINITE OUTPUTS!\n')
                        ## ## registerDoSEQ()
                        ## if(exists('cl')){ parallel::stopCluster(cl) }
                        ## stop('...TOO MANY NON-FINITE OUTPUTS. ABORTING')
                        mcsamples <- NULL
                        niter <- 0
                    } else {
                        mcsamples <- mcsubset(mcsamples, -toRemove)
                        niter <- niter - length(toRemove)
                    }
                }

                nitertot <- nitertot + niter
                availiter <- availiter + niter

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

                ## Check how many components were used at the last step
                ## usedcomponents <- mcsamplesKA$K[length(mcsamplesKA$K)]
                usedcomponents <- max(mcsamplesKA$K)
                cat('\nUSED COMPONENTS:', usedcomponents, 'OF', ncomponents, '\n')
#### Diagnostics
                ## Log-likelihood
                diagntime <- Sys.time()
                ##
                ll <- cbind(
                    util_Pcheckpoints(
                        Y = testdata,
                        learnt = mcsamples,
                        auxmetadata = auxmetadata
                    )
                )
                traces <- rbind(traces, ll)

                toRemove <- which(!is.finite(traces), arr.ind = TRUE)
                if (length(toRemove) > 0) {
                    flagll <- TRUE
                    cleantraces <- traces[-unique(toRemove[, 1]), , drop = FALSE]
                } else {
                    cleantraces <- traces
                }

                diagn <- mcmcstop(traces = cleantraces,
                    nsamples = nsamplesperchain,
                    availiter = availiter,
                    relerror = relerror,
                    thinning = thinning)

                ## Output available diagnostics
                for(i in names(diagn$toprint)) {
                    thisdiagn <- diagn$toprint[[i]]
                    cat(paste0('\n', i, ':'),
                        if(length(thisdiagn) > 1){
                            paste(signif(range(thisdiagn), 3), collapse = ' to ')
                        } else {
                            signif(thisdiagn, 3)
                        }
                    )
                }

                rm(cleantraces)

                cat('\nDiagnostics time',
                    printtimediff(difftime(Sys.time(), diagntime, units = 'auto')),
                    '\n')

                if (is.null(allmcsamples)) {
                    ## chain just started
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

                ## to save memory, only keep enough last iterations
                currentthinning <- diagn$proposed.thinning

                enoughiter <- currentthinning * (nsamplesperchain + 1L)

                if(availiter > enoughiter) {
                    allmcsamples <- mcsubset(
                        allmcsamples, -seq_len(availiter - enoughiter)
                    )
                    availiter <- enoughiter
                }


                ## ######################################
                ## ## CHECK IF CHAIN MUST BE CONTINUED ##
                ## ######################################

                ## calcIterThinning <- mcmcstop(
                ##     relerror = relerror,
                ##     nsamplesperchain = nsamplesperchain,
                ##     nitertot = nitertot,
                ##     thinning=thinning,
                ##     diagnESS=diagnESS, diagnIAT=diagnIAT,
                ##     diagnBMK=diagnBMK, diagnMCSE=diagnMCSE,
                ##     diagnStat=diagnStat, diagnBurn=diagnBurn,
                ##     diagnBurn2=diagnBurn2, diagnThin=diagnThin)

                requirediter <- max(minMCiterations - nitertot,
                    min(maxMCiterations - nitertot, diagn$reqiter) )

                cat('\nTotal number of iterations', nitertot,
                    '- required further', requirediter, '\n')

                if(requirediter > 0 && as.double(
                    Sys.time() - timestart0, units = 'hours'
                ) >= timeleft) {
                    cat('but stopping chain owing to lack of time\n')
                    stoppedchains <- stoppedchains + 1L
                    requirediter <- 0
                }

                if (requirediter > 0) {
                    ## limit number of iterations per loop, to save memory
                    niter <- min(requirediter + 1L, startupMCiterations)
                    subiter <- subiter + 1L
                    cat(
                        '\nChain #', chainnumber, '- chunk', subiter,
                        '(chain', achain, 'of', nchainsperthiscore,
                        'for this core): increasing by', niter, '\n'
                    )
                }
                reset <- FALSE
            }
#### END WHILE-LOOP OVER CHUNKS OF ONE CHAIN


            ## ################
            ## ## SAVE CHAIN ##
            ## ################

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

            if(any(tokeep <= 0)){
                cat('\nWARNING: have to reduce thinning owing to time constraints\n')
                tokeep <- round(seq(from = 1,
                    to = ncol(allmcsamples$W),
                    length.out = nsamplesperchain))
            }

            saveRDS(mcsubset(allmcsamples, tokeep),
                file = file.path(dirname,
                    paste0('_mcsamples',
                        dashnameroot, '--',
                        padchainnumber, '.rds'))
            )

            saveRDS(allmcsamplesKA,
                file = file.path(dirname,
                    paste0('_hyperparams_traces',
                        dashnameroot, '--',
                        padchainnumber, '.rds'))
            )


            ## put 'tokeep' in first slot to save only the last nsamplesperchain
            saveRDS(traces[ , , drop = FALSE],
                file = file.path(dirname,
                    paste0('_mcpartialtraces',
                        dashnameroot, '--',
                        padchainnumber, '.rds')
                ))

            ## possibly increase the count of chains with non-finite outputs
            nonfinitechains <- nonfinitechains + flagnonfinite

            ## ###########
            ## ## PLOTS ##
            ## ###########

#### Plot diagnostic traces of current chain
            if (plottraces) {
                cat('\nPlotting traces and samples.\n')

                ## tracegroups <- as.list(seq_len(min(4, ncol(traces))))
                ## names(tracegroups) <- colnames(traces)[1:min(4, ncol(traces))]
                ## grouplegends <- foreach(agroup = seq_along(tracegroups)) %do% {
                ##     c(
                ##         paste0('-- STATS ', names(tracegroups)[agroup], ' --'),
                ##         paste0('min ESS = ',
                ##             signif(min(diagnESS[tracegroups[[agroup]]]), 6)),
                ##         paste0('max IAT = ',
                ##             signif(max(diagnIAT[tracegroups[[agroup]]]), 6)),
                ##         paste0('max BMK = ',
                ##             signif(max(diagnBMK[tracegroups[[agroup]]]), 6)),
                ##         paste0('max MCSE = ',
                ##             signif(max(diagnMCSE[tracegroups[[agroup]]]), 6)),
                ##         paste0('stationary: ',
                ##             sum(diagnStat[tracegroups[[agroup]]]), '/',
                ##             length(diagnStat[tracegroups[[agroup]]])),
                ##         paste0('burn: ', signif(diagnBurn2, 6)),
                ##         paste0('max thin = ',
                ##             signif(max(diagnThin[tracegroups[[agroup]]]), 6))
                ##     )
                ## }

                ## Plot various info and traces
                ## colpalette <- 1:6 #seq_len(ncol(traces))
                ## names(colpalette) <- colnames(traces)[1:min(6, ncol(traces))]
                graphics.off()
                pdf(file.path(dirname,
                    paste0('_mcpartialtraces', dashnameroot, '--',
                        padchainnumber, '_', achain, '-', acore, '.pdf')
                ), height = 8.27, width = 11.69)

                ## Summary stats
                matplot(1:2, type = 'l', col = 'white',
                    main = paste0('Stats chain ', achain),
                    axes = FALSE, ann = FALSE)
                ## Legends
                ## legendpositions <- c('topleft', 'topright', 'bottomleft', 'bottomright')
                ## for (alegend in seq_along(grouplegends)) {
                ##     legend(x = legendpositions[alegend], bty = 'n', cex = 1.5,
                ##         legend = grouplegends[[alegend]]
                ##     )
                ## }
                legend(x = 'left', bty = 'n', cex = 1,
                    legend = c(
                        paste0('Chain ', chainnumber, ' - ',
                            achain, ' of core ', acore),
                        paste0('Test points ',
                            paste0('#', rownames(testdata), collapse=' ')
                        ),
                        paste0('Used components: ', usedcomponents, ' of ', ncomponents),
                        ## paste0('LL:  ( ', signif(mean(traces[, 1]), 3), ' +- ',
                        ##     signif(sd(traces[, 1]), 3), ' ) dHart'),
                        'NOTES:',
                        if (flagnonfinite) {
                            'some non-finite MC outputs'
                        },
                        if (usedcomponents > ncomponents - 5) {
                            'too many components used'
                        },
                        if (flagll) {
                            'non-finite values in diagnostics'
                        }
                    )
                )
                ## Traces of likelihood and cond. probabilities
                par(mfrow = c(1, 1))
                for (avar in 1:ncol(traces)) {
                    flexiplot(
                        y = 10*log10(traces[is.finite(traces[, avar]), avar]),
                        type = 'l', lty = 1, col = 1,
                        main = paste0('#', colnames(traces)[avar], ': ',
                            paste(
                                names(diagn$toprint),
                                sapply(diagn$toprint, function(xx){
                                    signif(xx[avar], 3)
                                }),
                                collapse = ' | ', sep = ': '
                            )),
                        cex.main = 1.25,
                        ylab = paste0('log_F(#',
                            colnames(traces)[avar],
                            ')/dHart'),
                        xlab = 'Monte Carlo sample',
                        family = family
                        ## mar = c(NA, 6, NA, NA)
                    )
                }
                dev.off()
            }

#### Plot Alpha and component occupation, if required
            if (showAlphatraces || showKtraces) {
                cat('Plotting component and Alpha information.\n')
                pdf(file = file.path(dirname,
                    paste0('hyperparams_traces', dashnameroot, '--',
                        padchainnumber, '_', achain, '-', acore, '.pdf')),
                    height = 8.27, width = 11.69)

                if (showKtraces) {
                    cat('\nSTATS USED COMPONENTS:\n')
                    print(summary(allmcsamplesKA$K))
                    ##
                    flexiplot(y = allmcsamplesKA$K, ylab = 'used components',
                        xlab = 'iteration', ylim = c(0, ncomponents))
                    flexiplot(x = 0:ncomponents,
                        y = tabulate(allmcsamplesKA$K + 1, nbins = ncomponents + 1),
                        type = 'l', xlab = 'used components', ylab = NA,
                        ylim = c(0, NA))
                }
                if (showAlphatraces) {
                    cat('\nSTATS alpha:\n')
                    print(summary(allmcsamplesKA$Alpha, na.rm = TRUE))
                    flexiplot(y = allmcsamplesKA$Alpha,
                        ylab = bquote(alpha), xlab = 'iteration',
                        ylim = c(1, nalpha))
                    flexiplot(x = seq(minalpha, maxalpha, by = byalpha),
                        y = tabulate(allmcsamplesKA$Alpha, nbins = nalpha),
                        type = 'l', xlab = bquote(alpha), ylab = '',
                        ylim = c(0, NA))
                }
                dev.off()
            }


            cat('\nCurrent time:', strftime(as.POSIXlt(Sys.time()),
                '%Y-%m-%d %H:%M:%S'))
            cat('\nMCMC + diagnostics time',
                printtimediff(difftime(Sys.time(), MCtimestart, units = 'auto')),
                '\n')

#### Print estimated end time
            endTime <- Sys.time() + 180 +
                ( (nchainsperthiscore + (acore > coreswithextrachain) - achain) *
                difftime(Sys.time(), MCtimestart) / achain )
            print2user(
                paste0(
                    '\rSampling. Core ', acore, ' estimated end time: ',
                    format(endTime, format='%Y-%m-%d %H:%M'),
                    '   '
                ),
                outcon
            )

            maxusedcomponents <- max(maxusedcomponents, usedcomponents)
            maxiterations <- max(maxiterations, nitertot)
        }
#### END LOOP OVER CHAINS (WITHIN ONE CORE)

        ##
        cat('\nCurrent time:',
            strftime(as.POSIXlt(Sys.time()), '%Y-%m-%d %H:%M:%S'))

        cat('\nTotal time',
            printtimediff(difftime(Sys.time(), headertimestart, units = 'auto')),
            '\n')

        ## output information from a core,
        ## passed to the originally calling process
        cbind(
            maxusedcomponents = maxusedcomponents,
            maxiterations = maxiterations,
            nonfinitechains = nonfinitechains,
            stoppedchains = stoppedchains,
            headertime = headertime,
            MCtime = difftime(Sys.time(), MCtimestart, units = 'secs'),
            usedmem = max(usedmem, sum(gc()[,6]))
        )
    }
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
                dx <- dim(xx)[-length(dim(xx))]
                dim(temp) <- c(dx, length(temp) / prod(dx))
                temp
            },
            mc1, mc2,
            SIMPLIFY = FALSE
        )
    }
    mcsamples <- foreach(chainnumber = 1:nchains,
        .combine = joinmc, .multicombine = FALSE) %do% {
            padchainnumber <- sprintf(paste0('%0', nchar(nchains), 'i'), chainnumber)

            readRDS(file = file.path(dirname,
                paste0('_mcsamples', dashnameroot, '--',
                    padchainnumber, '.rds')
            ))
        }

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
        paste0('_testdata_', 0, '.rds')))
    cat('\nChecking test data\n(', paste0('#', rownames(testdata)), ')\n')

    traces <- cbind(
        util_Pcheckpoints(
            Y = testdata,
            learnt = mcsamples,
            auxmetadata = auxmetadata
        )
    )
    traces <- traces[apply(traces, 1, function(x) { all(is.finite(x)) }), ,
        drop = FALSE]
    colnames(traces) <- rownames(testdata)

    saveRDS(traces, file = file.path(dirname,
        paste0('MCtraces', dashnameroot, '.rds')
    ))

    diagn <- mcmcstop(traces = traces,
        nsamples = nsamples,
        availiter = 0,
        relerror = relerror,
        thinning = thinning)

    cat('\n')
    for(i in names(diagn$toprint)) {
        thisdiagn <- diagn$toprint[[i]]
        cat(paste0(i, ':'),
            if(length(thisdiagn) > 1){
                paste0(
                    'min: ', signif(min(thisdiagn), 2),
                    '  max: ', signif(max(thisdiagn), 2),
                    '  mean: ', signif(mean(thisdiagn), 2)
                )
                ## paste(signif(range(thisdiagn), 3), collapse = ' to ')
            } else {
                signif(thisdiagn, 3)
            },
            '\n')
    }

    ## Plot various info and traces
    cat('\nPlotting final Monte Carlo traces and marginal samples...\n')

    ##
    ## colpalette <- seq_len(ncol(traces))
    ## names(colpalette) <- colnames(traces)
    graphics.off()
    pdf(file.path(dirname,
        paste0('MCtraces', dashnameroot, '.pdf')
    ), height = 8.27, width = 11.69)
    ## Traces of likelihood and cond. probabilities
    for (avar in 1:ncol(traces)) {
        ## Do not join separate chains in the plot
        division <- (if(nrow(traces) > nchains) nchains else 1)
        flexiplot(
            y = 10*log10(traces[is.finite(traces[, avar]), avar]),
            ## x = matrix(seq_len(nsamples), ncol = division),
            ## y = matrix(10*log10(traces[, avar]), ncol = division),
            type = 'l', lty = 1, col = 1,
            ## col = 1:6, # to evidence consecutive chains
            main = paste0('#', colnames(traces)[avar], ': ',
                paste(
                    names(diagn$toprint),
                    sapply(diagn$toprint, function(xx){
                        signif(xx[avar], 3)
                    }),
                    collapse = ' | ', sep = ': '
                )),
            cex.main = 1.25,
            ylab = paste0('log_F(#',
                colnames(traces)[avar],
                ')/dHart'),
            xlab = 'sample', family = family
            ## mar = c(NA, 6, NA, NA)
        )
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
    ## sink(file = NULL, type = 'output')
    ## sink(file = NULL, type = 'message')
    ## close(outcon)

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
    ## Partial files are identified by an initial underscore "_"
    ## Should we leave the plots of partial traces?
    ## maybe create an additional argument to let the user decide?
    if (cleanup) {
        cat('\nRemoving temporary output files.\n')
        file.remove(dir(dirname,
            pattern = paste0('^_.*\\..*$'),
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

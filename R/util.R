#' Determine the status of parallel processing. If cores are requested,
#' set up the cores. If cores are registered, keep using those.
#' @param parallel logical or integer: whether to use pre-existing parallel
#'   workers, or how many to create and use. Default `TRUE`.
#' @param silent logical: give warnings or updates in the computation?
#'   Default `FALSE`.
#' @keywords internal
#' @return number of cores
setup_parallel <- function(parallel, silent = FALSE) {
    # Set cluster object to false to start
    cluster <- FALSE
    # Check if 'parallel' argument is either logical or numeric
    if (!is.logical(parallel) && !(parallel %% 1 == 0)) {
        stop('Argument parallel must be TRUE, FALSE or an integer.')
    }
    if (parallel == 0) {parallel <- FALSE}
    ## ## Alternative way to register cores;
    ## ## might need to be used for portability to Windows?
    ## registerDoSEQ()
    ## cluster <- makePSOCKcluster(ncores)

    if (isTRUE(parallel)) {
        if (foreach::getDoParRegistered()) {
            print_pretty(c('Using already registered', foreach::getDoParName(),
                           'with', foreach::getDoParWorkers(), 'workers'), silent)
            ncores <- foreach::getDoParWorkers()
        } else {
            print_pretty('No parallel backend registered.', silent)
            ncores <- 1
        }
    } else if (parallel >= 2) {
        if (foreach::getDoParRegistered()) {
            ncores <- min(foreach::getDoParWorkers(), parallel)
            print_pretty(c('Using already registered', foreach::getDoParName(),
                        'with', foreach::getDoParWorkers(), 'workers'), silent)
            if (parallel > ncores) {
                print_pretty(c('NOTE: fewer pre-registered cores',
                            'than requested in the "parallel" argument.'), silent)
            }
        } else {
            cluster <- parallel::makeCluster(parallel)
            doParallel::registerDoParallel(cluster)
            print_pretty(c('Registered', foreach::getDoParName(),
                    'with', foreach::getDoParWorkers(), 'workers'), silent)
            ncores <- parallel
        }
    } else {
        print_pretty('No parallel backend registered.', silent)
        ncores <- 1
    }
    workers <- list('ncores' = ncores, 'cluster' = cluster)
    return(workers)
}

closecoresonexit <- function(cluster, silent = FALSE) {
    print_pretty('Closing connections to cores.', silent)
    foreach::registerDoSEQ()
    parallel::stopCluster(cluster)
    env <- foreach:::.foreachGlobals
    rm(list = ls(name = env), pos = env)
}

print_pretty <- function(message, silent = FALSE) {
    if (!silent) {
        if (is.list(message)) {
            for (item in list) cat(item, '\n')
        } else {cat(message, '\n')}
    }
}

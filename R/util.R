#' Determine the status of parallel processing. If cores are requested,
#' set up the cores. If cores are registered, keep using those.
#' @param parallel logical or integer: whether to use pre-existing parallel
#'   workers, or how many to create and use. Default `TRUE`.
#' @param silent logical: give warnings or updates in the computation?
#'   Default `FALSE`.
#' @keywords internal
#' @return number of cores
checkParallel <- function(parallel, silent) {

    cluster <- FALSE

    if (is.logical(parallel) && parallel) {
        if (foreach::getDoParRegistered()) {
            if (!silent) {
                cat('Using already registered', foreach::getDoParName(),
                    'with', foreach::getDoParWorkers(), 'workers\n')
            }
            ncores <- foreach::getDoParWorkers()
        } else {
            if (!silent) {
                cat('No parallel backend registered.\n')
            }
            ncores <- 1
        }
    } else if (is.numeric(parallel) && parallel >= 2) {
        if (foreach::getDoParRegistered()) {
            ncores <- min(foreach::getDoParWorkers(), parallel)
            if (!silent) {
                cat('Using already registered', foreach::getDoParName(),
                    'with', foreach::getDoParWorkers(), 'workers\n')
                if(parallel > ncores) {
                    cat('NOTE: fewer pre-registered cores',
                        'than requested in the "parallel" argument.\n')
                }
            }
        } else {
            ## ##
            ## ## Alternative way to register cores;
            ## ## might need to be used for portability to Windows?
            ## registerDoSEQ()
            ## cluster <- makePSOCKcluster(ncores)
            ## ##
            cluster <- parallel::makeCluster(parallel)
            doParallel::registerDoParallel(cluster)
            if (!silent) {
                cat('Registered', foreach::getDoParName(),
                    'with', foreach::getDoParWorkers(), 'workers\n')
            }
            ncores <- parallel
        }
    } else {
        if (!silent) {
            cat('No parallel backend registered.\n')
        }
        ncores <- 1
    }
    workers <- list("ncores" = ncores, "cluster" = cluster)
    return(workers)
}

closecoresonexit <- function(cluster, silent) {
    if(!silent) {
        cat('\nClosing connections to cores.\n')
    }
    foreach::registerDoSEQ()
    parallel::stopCluster(cluster)
    env <- foreach:::.foreachGlobals
    rm(list = ls(name = env), pos = env)
}

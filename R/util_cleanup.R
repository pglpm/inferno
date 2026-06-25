#' Cleanup a learn()-output directory
#'
#' For deeper monitoring of the MCMC, the user can require the 'learn()' function not to clean intermediate MCMC-related files generated during the computation.
#'
#' The present function can be used to remove these intermediate files from the output directory created by 'learn()'.
#'
#' @keywords internal
util_cleanup <- function(path){
    if(all(file.exists(file.path(path, c(
        'learnt.rds',
        'MCtraces.rds',
        'MCtraces.pdf',
        'log-1.log',
        'plotquantiles_learnt.pdf',
        'plotsamples_learnt.pdf'
    ))))){
        invisible(file.remove(dir(path,
            pattern = paste0('^___.*\\..*$'),
            full.names = TRUE
        )))
    } else {
        stop("This doesn't look like a 'learn()' output directory")
    }
}


#' Join '____tempPtraces-' files
#'
#' For deeper monitoring of the MCMC, the user can require the 'learn()' function not to clean intermediate MCMC-related files generated during the computation. The files with prefix '____tempPtraces-' contain chunks of MCMC traces.
#'
#' The present function can be used to join them into a single trace.
#' 
#' @keywords internal
util_joinPtraces <- function(path){
    Plist <- list.files(path = path, pattern = '^____tempPtraces-.*\\.rds$')
    chainlist <- unique(sub(
        pattern = '^____tempPtraces-([0-9]+)-.*',
        replacement = '\\1',
        x = Plist
    ))

    for(achain in chainlist){
        chunklist <- grepv(
            pattern = paste0('^____tempPtraces-', achain, '-.*\\.rds$'),
            x = Plist
        )
        chunks <- sort(as.numeric(unique(sub(
            pattern = paste0('^____tempPtraces-', achain, '-([0-9]+)\\.rds$'),
        replacement = '\\1',
        x = chunklist
        ))))

        Ptraces <- NULL
        for(achunk in chunks){
            Ptraces <- rbind(Ptraces,
                readRDS(file.path(path,
                    paste0('____tempPtraces-', achain, '-', achunk, '.rds')
                ))
            )
        }
        saveRDS(Ptraces, file.path(path,
            paste0('_allPtraces-', achain, '.rds')
        ))
    }
    cat('\nDone.\n')
}

#### Cleanup a learn()-output directory
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


#### Join '____tempPtraces-' files
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

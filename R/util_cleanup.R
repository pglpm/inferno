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


#' Join '____tempPtraces-' files
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

#' Convert learnt with R/C/Dvar to learnt with R/C/Dsd
#' 
#' @keywords internal
util_learntvar2sd <- function(file){
    learnt <- readRDS(file = file)
    ## Swap variances and standard deviations
    if(!is.null(learnt$Rvar)){
        learnt$Rvar <- sqrt(learnt$Rvar)
        names(learnt)[which(names(learnt) == 'Rvar')] <- 'Rsd'
    } else if(!is.null(learnt$Rsd)){
        learnt$Rsd <- learnt$Rsd^2
        names(learnt)[which(names(learnt) == 'Rsd')] <- 'Rvar'
    }
    ##
    if(!is.null(learnt$Cvar)){
        learnt$Cvar <- sqrt(learnt$Cvar)
        names(learnt)[which(names(learnt) == 'Cvar')] <- 'Csd'
    } else if(!is.null(learnt$Csd)){
        learnt$Csd <- learnt$Csd^2
        names(learnt)[which(names(learnt) == 'Csd')] <- 'Cvar'
    }
    ##
    if(!is.null(learnt$Dvar)){
        learnt$Dvar <- sqrt(learnt$Dvar)
        names(learnt)[which(names(learnt) == 'Dvar')] <- 'Dsd'
    } else if(!is.null(learnt$Dsd)){
        learnt$Dsd <- learnt$Dsd^2
        names(learnt)[which(names(learnt) == 'Dsd')] <- 'Dvar'
    }
    ##
    saveRDS(object = learnt, file = file)
}

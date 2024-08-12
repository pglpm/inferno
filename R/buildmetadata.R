#' Build a preliminary metadata file from a given dataset
#'
#' @param data A dataset, given as a \code{\link[base]{data.frame}}
#' or as a file path to a csv file.
#' @param file String: name of csv file where the metadata should be saved;
#'   if `NULL`: output metadata as `VALUE`.
#' @param includevrt Character or `NULL`: name of variates in dataset to be included.
#' @param excludevrt Character or `NULL`: name of variates in dataset to be excluded.
#' @param addsummary2metadata Logical: also output some diagnostic statistics
#'    in the metadata? Default `FALSE`
#' @param backupfiles Logical: rename previous metadata file if it exists?
#' Default `TRUE`
#' @param verbose Logical: output heuristics for each variate? Default `TRUE`.
#'
#' @return If `file = NULL`, a preliminary metadata file is created
#'   and `VALUE` is `NULL`;
#'   otherwise `VALUE` is a data.frame containing the metadata.
#'
#' @export
buildmetadata <- function(
    data,
    file = NULL,
    includevrt = NULL,
    excludevrt = NULL,
    addsummary2metadata = FALSE,
    backupfiles = FALSE,
    verbose = TRUE
) {

    ## Greater Common Denominator function
    ## used to calculate multiplicity of values of continuous variates
    ## and guess whether they are rounded
    gcd <- function(...) {
        suppressWarnings(Reduce(function(a, b) {
            if (b == 0) a else Recall(b, a %% b)
        }, c(...)))
    }

    ## keywords to recognize ordinal variates
    ordinalkeywords <- c('low', 'high', 'very', 'medium', 'strong',
        'poor', 'good', 'always')

    ## collect warning messages to be given at the end
    warninglist <- NULL

    ## Read datafile if it exists
    datafile <- NULL
    if (is.character(data) && file.exists(data)) {
        datafile <- data
        data <- read.csv(datafile, na.strings = '')
    }
    data <- as.data.frame(data)

#### Prepare empty metadata frame
    metadata <- as.data.frame(
        c(list(name = NA,
            type = NA,
            Nvalues = NA,
            gapwidth = NA,
            domainmin = NA,
            domainmax = NA,
            minincluded = NA,
            maxincluded = NA,
            ## lowvalue = NA,
            ## centralvalue = NA,
            ## highvalue = NA,
            plotmin = NA,
            plotmax = NA),
            if (addsummary2metadata) {
                list(datamin = NA,
                    datamax = NA,
                    datamaxrep = NA,
                    dataNvalues = NA)
            }
        )
    )[-1,]

    ## build list of variates, including only those in 'includevrt'
    ## and excluding those in 'excludevrt'
    variatelist <- colnames(data)
    if(!is.null(includevrt)) {
        variatelist <- intersect(variatelist, includevrt)
    }
    if(!is.null(excludevrt)) {
        variatelist <- setdiff(variatelist, excludevrt)
    }

    if(verbose){ cat('\nAnalyzing', length(variatelist), 'variates for',
        nrow(data), 'datapoints.\n') }

    ## Loop over variates (columns) in data frame
    ## cat('\n')
    for (name in variatelist) {
        if(verbose){ cat(paste0('\n* "', name, '" variate:\n')) }
        ## remove missing values
        x <- data[[name]]
        x <- x[!is.na(x)]

        ## check if it's a 'factor' object and transform
        if(is.factor(x)){
            x <- as.character(x)
        }

        unx <- sort(unique(x))
        uniquex <- length(unx)

        ## if this variate has only one value, then it bears no information
        if (uniquex <= 1) {
            if(verbose){
                cat('  Only one value present.\n')
            }
            warninglist <- c(warninglist,
                paste0('\n* "', name, '" variate',
                    ' does not have at least two distinct values.',
                    '\nDiscarded because non-informative.',
                    '\nPlease insert its characteristics by hand in the metadata file.'))
            next
        }

        ## If the values are strings,
        ## check whether they can be reinterpreted as numeric
        if(!is.numeric(unx) &&
               !any(is.na(suppressWarnings(as.numeric(unx))))) {
            x <- as.numeric(x)
            unx <- sort(as.numeric(unx))
        }

        ## If the variate is numeric,
        ## calculate several characteristics used later
        if (is.numeric(unx) && uniquex > 2) {
            ## these two are used for diagnostic purposes
            datamin <- min(x)
            datamax <- max(x)
            ## values strictly within the interior of the domain
            ix <- x[!(x %in% range(x))]
            maxrep <- max(table(ix)) # max of repeated inner values
            meanrep <- mean(table(ix)) # average of repeated inner values
            ## differences between consecutive unique values
            jumpsx <- unique(signif(diff(unx), 3))
            rangex <- diff(range(x))
            multi <- 10^(-min(floor(log10(jumpsx))))
            ## jumps between unique values are integer multiples
            ## of 'jumpquantum'
            jumpquantum <- gcd(jumpsx * multi) / multi
            datavalues <- NULL
        } else {
            datavalues <- as.character(unx)
        }

### Check whether the variate has been first rounded
### and then transformed to logarithmic scale.
### This may lead to problems
        if(is.numeric(unx) && datamin > 0 && uniquex > 2) {
            for(base in c(2, exp(1), 10)) {
                ljumpsx <- unique(signif(diff(log(unx, base)), 3))
                lrangex <- diff(log(range(x), base))
                lmulti <- 10^(-min(floor(log10(ljumpsx))))
                ## jumps between unique values are integer multiples
                ## of 'jumpquantum'
                ljumpquantum <- gcd(ljumpsx * lmulti) / lmulti
                if (!(
                    ljumpquantum / lrangex < 1e-5 || # no visible gaps between values
                        (ljumpquantum / lrangex < 1e-3 &&
                             maxrep <= 100)  # visible gaps, but not many value repetitions
                )) {
                    warninglist <- c(warninglist,
                        paste0('\n* "', name, '" variate',
                            ' appears to have been rounded\n',
                            'and then transformed to logarithmic scale.\n',
                            'This may lead to problems in the inference.\n',
                            'Preferably, transform it back to non-logarithmic scale.'))
                    break
                }
            }
        }

#### Now make educated guess for the type of variate

        if (uniquex == 2) {
            ## Binary variate? (has only two unique values)
            type <- 'binary'
            Nvalues <- 2
            gapwidth <- NA
            domainmin <- NA
            domainmax <- NA
            minincluded <- NA
            maxincluded <- NA
            ## lowvalue <- NA
            ## centralvalue <- NA
            ## highvalue <- NA
            plotmin <- NA
            plotmax <- NA
            names(datavalues) <- paste0('V', 1:2)
            ##
            if(verbose){
                cat('  Two different values detected:\n',
                    paste0('"', datavalues, '"', collapse=', '),
                    '\n')
                cat('  Assuming variate to be BINARY.\n')
            }

        } else if (!is.numeric(x) &&
                       !any(sapply(ordinalkeywords, function(keyw){
                           grepl(keyw, datavalues, fixed = TRUE)
                       }))) {
            ## Nominal variate? (non-numeric values)
            type <- 'nominal'
            Nvalues <- uniquex
            gapwidth <- NA
            domainmin <- NA
            domainmax <- NA
            minincluded <- NA
            maxincluded <- NA
            ## lowvalue <- NA
            ## centralvalue <- NA
            ## highvalue <- NA
            plotmin <- NA
            plotmax <- NA
            names(datavalues) <- paste0('V', 1:Nvalues)
            ##
            if(verbose){
                cat('  - ', Nvalues, 'different non-numeric values detected:\n',
                    paste0('"', datavalues, '"', collapse=', '),
                    '\n')
                cat('  which do not seem to refer to an ordered scale.\n')
                cat('  Assuming variate to be NOMINAL.\n')
            }

        } else if (!is.numeric(x)) {
            ## Nominal variate? (non-numeric values)
            ## Ordinal variate with few numeric values?
            type <- 'ordinal'
            Nvalues <- uniquex
            gapwidth <- NA
            domainmin <- NA
            domainmax <- NA
            minincluded <- NA
            maxincluded <- NA
            ## lowvalue <- NA
            ## centralvalue <- NA
            ## highvalue <- NA
            plotmin <- NA
            plotmax <- NA
            names(datavalues) <- paste0('V', 1:Nvalues)
            ##
            if(verbose){
                cat('  - ', Nvalues, 'different non-numeric values detected:\n',
                    paste0('"', datavalues, '"', collapse=', '),
                    '\n')
                cat('  which seem to refer to an ordered scale.\n')
                cat('  Assuming variate to be ORDINAL.\n')
                cat('  Please appropriately reorder its values in metadata file.\n')
            }
            warninglist <- c(warninglist,
                paste0('\n* "', name, '" variate',
                    ' appears to be ordinal.',
                    '\nPlease appropriately reorder its values in metadata file.'))


#### Now we know that the variate is numeric

        } else if (uniquex <= 10 &&
                       (all(x == round(x)) || jumpquantum == round(jumpquantum))) {
            ## Ordinal variate with few numeric values?
            type <- 'ordinal'
            Nvalues <- uniquex
            ## if the values are spaced by an integer,
            ## no need to use the 'V' fields
            if(jumpquantum == round(jumpquantum)) {
                gapwidth <- jumpquantum
                domainmin <- datamin
                domainmax <- datamax
                datavalues <- NULL
            } else {
                gapwidth <- NA
                domainmin <- NA
                domainmax <- NA
                datavalues <- as.character(unx)
                names(datavalues) <- paste0('V', 1:Nvalues)
            }
            minincluded <- NA
            maxincluded <- NA
            ## lowvalue <- NA
            ## centralvalue <- NA
            ## highvalue <- NA
            plotmin <- NA
            plotmax <- NA
            ##
            if(verbose){
                cat('  - Only', Nvalues, 'different numeric values detected:\n')
                if(jumpquantum == round(jumpquantum)) {
                    cat('from', domainmin, 'to', domainmax,
                        'in steps of', jumpquantum, '\n')
                } else {
                    cat(paste0('"', datavalues, '"', collapse=', '), '\n')
                }
                cat('  Assuming variate to be ORDINAL.\n')
            }

            ## } else if (uniquex <= 10 && ) {
            ##     ## Ordinal variate with many numeric values?
            ##     type <- 'ordinal'
            ##     Nvalues <- uniquex
            ##     gapwidth <- jumpquantum
            ##     domainmin <- datamin
            ##     domainmax <- datamax
            ##     minincluded <- NA
            ##     maxincluded <- NA
            ##     ## lowvalue <- NA
            ##     ## centralvalue <- NA
            ##     ## highvalue <- NA
            ##     plotmin <- NA
            ##     plotmax <- NA
            ##     datavalues <- NULL
            ##     ##
            ##     if(verbose){
            ##         cat(' - ', Nvalues, 'different numeric values detected\n')
            ##         cat('  distance between datapoints is a multiple of integer',
            ##             jumpquantum, '\n')
            ##         cat('  Assuming variate to be ORDINAL.\n')
            ##     }
            ##     warninglist <- c(warninglist,
            ##         paste0('* "', name, '" variate',
            ##             ' appears to be ordinal,',
            ##             '\nbut it could also be a rounded continous variate'))


        } else {
            ## The variate seems continuous
            type <- 'continuous'
            Nvalues <- Inf
            ## preliminary values, possibly modified below
            gapwidth <- NA
            domainmin <- -Inf # signif(datamax - 3 * rangex, 1)
            domainmax <- +Inf # signif(datamax + 3 * rangex, 1)

            ## Q1 <- quantile(x, probs = 0.25, type = 6)
            ## centralvalue <- quantile(x, probs = 0.5, type = 6)
            ## Q3 <- quantile(x, probs = 0.75, type = 6)
            ## ## Borderline case if the first and second quartile have the same value
            ## if (lowvalue == highvalue) {
            ##   lowvalue <- if (sum(x < Q1) > 0) {
            ##                 max(x[x < Q1])
            ##               } else {
            ##                 max(x[x <= Q1])
            ##               }
            ##   highvalue <- if (sum(x > Q3) > 0) {
            ##                  min(x[x > Q3])
            ##                } else {
            ##                  min(x[x >= Q3])
            ##                }
            ## }
            plotmin <- datamin - IQR(x, type=6) / 2
            plotmax <- datamax + IQR(x, type=6) / 2
            datavalues <- NULL
            ##
            if(verbose){
                cat('  - Numeric values between',
                    datamin, 'and', datamax, '\n')
                cat('  Assuming variate to be CONTINUOUS.\n')
            }

            roundedflag <- FALSE
#### Consider subcases for gapwidth
            ## for rounded variates it doesn't matter whether the boundaries
            ## are included in the domain or not
            if (!(
                jumpquantum / rangex < 1e-5 || # no visible gaps between values
                    (jumpquantum / rangex < 1e-3 &&
                         maxrep <= 100)  # visible gaps, but not many value repetitions
            )) {
                roundedflag <- TRUE
                ## The variate seems to have quantized differences
                ## between unique values
                ## hence it might be rounded or integer/ordinal

                ## if(rangex / gapwidth > 256) {
                ##   ## The variate seems rounded,
                ##   ## but no latent-variable representation is needed
                ##   message('\nNOTE: variate ', name, ' is probably rounded,\n',
                ##           'but is treated as not rounded ',
                ##           'owing to its large range of values.\n')
                ##   gapwidth <- 0
                ## } else {
                gapwidth <- jumpquantum
                minincluded <- maxincluded <- NA
                ##
                if(verbose){
                    cat('  - Distance between datapoints',
                        'is a multiple of', jumpquantum, '\n')
                    cat('  Assuming variate to be ROUNDED.\n')
                }
                if(jumpquantum >=1) {
                    warninglist <- c(warninglist,
                        paste0('\n* "', name, '" variate',
                            ' appears to be continuous and rounded,',
                            '\nbut it could also be an ordinal variate'))
                }
                ## domainmin <- signif(datamin - 4 * rangex, 1)
                ## domainmax <- signif(datamax + 4 * rangex, 1)

            } # End rounded case


#### Consider subcases for closedness of domain
            if(!roundedflag){
                minincluded <- FALSE
                maxincluded <- FALSE
            }
            if (sum(x == datamin) > 1) {
                ## seems to have left-closed domain
                domainmin <- datamin
                plotmin <- datamin
                if(!roundedflag){
                    minincluded <- TRUE
                }
                ##
                if(verbose){
                    cat('  - Several datapoints have minimum value', datamin, '\n')
                    cat('  Assuming "domainmin" to be this minimum observed value\n')
                    if(!roundedflag){
                        cat('  and to be included in the domain',
                            '(singular probabilities there).\n')
                    }
                }
            } else if (all(unx > 0)) {
                ## seems to be strictly positive
                domainmin <- 0
                plotmin <- max((domainmin + datamin) / 2, plotmin)
                if(!roundedflag){
                    minincluded <- FALSE
                }
                ##
                if(verbose){
                    cat('  - All values are positive\n')
                    cat('  Assuming "domainmin" to be 0\n')
                    if(!roundedflag){
                        cat('  with 0 excluded from domain.\n')
                    }
                }
            } else if (all(unx >= 0)) {
                ## seems to be non-negative
                domainmin <- 0
                plotmin <- max((domainmin + datamin) / 2, plotmin)
                if(!roundedflag){
                    minincluded <- TRUE
                }
                ##
                if(verbose){
                    cat('  - All values are non-negative\n')
                    cat('  Assuming "domainmin" to be 0\n')
                    if(!roundedflag){
                        cat('  with 0 included in the domain.\n')
                    }
                }
            }

            if (sum(x == datamax) > meanrep) {
                ## seems to right-closed domain
                domainmax <- datamax
                plotmax <- datamax
                if(!roundedflag){
                    maxincluded <- TRUE
                }
                ##
                if(verbose){
                    cat('  - Many datapoints have maximum value,', datamax, '\n')
                    cat('  Assuming "domainmax" to be this maximum observed value\n')
                    if(!roundedflag){
                        cat('  and to be included in the domain',
                            '(singular probabilities there).\n')
                    }
                }
            }

        } # End continuous-variate case

#### Create metadata object
        metadata <- merge(metadata,
            c(
                list(name = name,
                    type = type,
                    Nvalues = Nvalues,
                    gapwidth = gapwidth,
                    domainmin = domainmin,
                    domainmax = domainmax,
                    minincluded = minincluded,
                    maxincluded = maxincluded,
                    ## lowvalue = signif(lowvalue,3),
                    ## centralvalue = signif(centralvalue,3),
                    ## highvalue = signif(highvalue,3),
                    plotmin = signif(plotmin,2),
                    plotmax = signif(plotmax,2)),
                if (addsummary2metadata) {
                    list(datamin = datamin,
                        datamax = datamax,
                        datamaxrep = max(table(x)),
                        dataNvalues = uniquex)
                },
                as.list(datavalues)
            ),
            sort = FALSE, all = TRUE)
    } # End loop over variates
    if(verbose){cat('\n')}

    ## Print warnings
    if(!is.null(warninglist)){
        message(
            '=========',
            '\nWARNINGS',
            ' - please make sure to check these variates in the metadata file:')
        for(awarning in warninglist){message(awarning)}
        cat('=========\n')
    }

                                        # Save to file if the file parameter is set
    if (!is.null(file)) {
                                        # If file is a string, name the metadata file that
        if (is.character(file)) {
            file <- paste0(sub('.csv$', '', file), '.csv')
        }
                                        # If file is not a string, and the data is read from a file,
                                        # save the metadata as metadata_<datafile>.csv
        else if (!is.null(datafile)) {
            file <- paste0(sub('.csv$', '', file), '_metadata.csv')
        }
                                        # If the file already exists, rename the old file as a backup file
                                        # and name this metadata file as intended
        if (backupfiles) {
            if (file.exists(file)) {
                file.rename(from = file, to = paste0(sub('.csv$', '', file), '_bak', format(Sys.time(), '%y%m%dT%H%M%S'),'.csv'))
            }
        }
                                        # Save the file
        write.csv(metadata, file, row.names = FALSE, quote = FALSE, na = '')
        cat('\nSaved proposal metadata file as', paste0('"', file, '"'), '\n')

    } else {
                                        # Else just print to console
        metadata
    }
}

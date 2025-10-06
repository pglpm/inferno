unloadNamespace('inferno')
library('inferno')


## metadatatemplate(penguins, file='meta_penguins')

seed <- 16
parallel <- 6

outputdir <- '__learn_penguins'
learntdir <- learn(
    data = penguins,
    prior = FALSE,
    metadata = 'meta_penguins2.csv',
    outputdir = outputdir,
    appendtimestamp = TRUE,
    appendinfo = TRUE,
    outputvalue = 'directory',
    parallel = parallel,
    ## parameters for short test run:
    subsampledata = 10,
    maxhours = 0,
    nsamplesperchain = 60,
    nchains = parallel * 2,
    ##
    seed = seed
)


testf <- function(){
cl <- parallel::makeCluster(3)
doParallel::registerDoParallel(cl)
outcongen <- file('__log0.log', open = 'w')
writeLines('start lines 1', con = outcongen)
writeLines('start lines 2', con = outcongen)
for(i in 1:3){
    writeLines(text = paste0('core ', i, ' text'), con = outcongen)
}
writeLines(text = 'core replace', con = outcongen)
close(outcongen)
outp <- foreach(i = 1:3) %dopar% {
        outcon <- file(paste0('__log', i, '.log'), open = 'w')
        sink(file = outcon, type = 'output')
        sink(file = outcon, type = 'message')
        library('nimble')
        cat('cat to local\n')
        message('message to local')
        Sys.sleep(sample(1:5, 1))
        ## outcongen <- file('__log0.log', open = 'r+')#, blocking = FALSE)
        ## tochange <- readLines(outcongen, -1)
        ## tochange[length(tochange)] <- paste0('new time of core ', i)
        ## writeLines(text = tochange, con = outcongen)
        ## close(outcongen)
        outcongen <- '__log0.log'
        tochange <- readLines(outcongen, -1)
        tochange[length(tochange)] <- paste0('new time of core ', i)
        writeLines(text = tochange, con = outcongen)
        ## writeLines(text = paste0('new time of core ', i),
        ##     con = outcongen, sep='\r')
        ##
        ##
        ## sink(file = outcongen, type = 'output')
        ## sink(file = outcongen, type = 'message')
        ## cat('cat to gen\n')
        ## message('message to gen')
        ## outcon <- file(paste0('__log', i, '.log'), open = 'a')
        ## sink(file = outcon, type = 'output')
        ## sink(file = outcon, type = 'message')
        cat('cat to local 2\n')
        message('message to local 2')
sink(file = NULL, type = 'output')
sink(file = NULL, type = 'message')
        close(outcon)
    i
    }
## sink(file = NULL, type = 'output')
## sink(file = NULL, type = 'message')
foreach::registerDoSEQ()
parallel::stopCluster(cl)
outcongen <- file('__log0.log', open = 'r+')
tochange <- readLines(outcongen, -1)
writeLines(c(tochange, 'final message'), con = outcongen)
close(outcongen)
outp
}
##
testf()

testf <- function(){
    library('doFuture')
future::plan(cl)
outcongen <- file('__log0.log', open = 'w')
writeLines('start lines 1', con = outcongen)
writeLines('start lines 2', con = outcongen)
for(i in 1:3){
    writeLines(text = paste0('core ', i, ' text'), con = outcongen)
}
writeLines(text = 'core replace', con = outcongen)
close(outcongen)
outp <- foreach(i = 1:3) %dopar% {
        outcon <- file(paste0('__log', i, '.log'), open = 'w')
        sink(file = outcon, type = 'output')
        sink(file = outcon, type = 'message')
        library('nimble')
        cat('cat to local\n')
        message('message to local')
        Sys.sleep(sample(1:5, 1))
        ## outcongen <- file('__log0.log', open = 'r+')#, blocking = FALSE)
        ## tochange <- readLines(outcongen, -1)
        ## tochange[length(tochange)] <- paste0('new time of core ', i)
        ## writeLines(text = tochange, con = outcongen)
        ## close(outcongen)
        outcongen <- '__log0.log'
        tochange <- readLines(outcongen, -1)
        tochange[length(tochange)] <- paste0('new time of core ', i)
        writeLines(text = tochange, con = outcongen)
        ## writeLines(text = paste0('new time of core ', i),
        ##     con = outcongen, sep='\r')
        ##
        ##
        ## sink(file = outcongen, type = 'output')
        ## sink(file = outcongen, type = 'message')
        ## cat('cat to gen\n')
        ## message('message to gen')
        ## outcon <- file(paste0('__log', i, '.log'), open = 'a')
        ## sink(file = outcon, type = 'output')
        ## sink(file = outcon, type = 'message')
        cat('cat to local 2\n')
        message('message to local 2')
sink(file = NULL, type = 'output')
sink(file = NULL, type = 'message')
        close(outcon)
    i
    }
## sink(file = NULL, type = 'output')
## sink(file = NULL, type = 'message')
foreach::registerDoSEQ()
parallel::stopCluster(cl)
outcongen <- file('__log0.log', open = 'r+')
tochange <- readLines(outcongen, -1)
writeLines(c(tochange, 'final message'), con = outcongen)
close(outcongen)
outp
}
##
testf()

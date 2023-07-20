mcsubset <- function(mcsamples, subsamples){
    lapply(mcsamples,function(xx){
        do.call('[',c(list(xx),rep(TRUE,length(dim(xx))-1), list(subsamples), list(drop=FALSE)) )
    })
}

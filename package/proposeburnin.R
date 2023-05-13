## Assess stationarity (from LaplacesDemon::LaplacesDemon)
proposeburnin <- function(x, batches=16){
    if(is.null(dim(x))){x <- cbind(x) }
    lx <- nrow(x)
    if(lx%%batches != 0){
        x <- x[1:(batches * trunc(lx/batches)), ]
    }
    lx2 <- nrow(x)
    HD <- LaplacesDemon::BMK.Diagnostic(x, batches = batches)
    Ind <- 1 * (HD > 0.5)
    BurnIn <- lx
    batch.list <- seq(from = 1, to = lx2, by = floor(lx2/batches))
    for (i in 1:(batches-1)) {
        if (sum(Ind[, i:(batches-1)]) == 0) {
            BurnIn <- batch.list[i] - 1
            break
        }
    }
BurnIn
}

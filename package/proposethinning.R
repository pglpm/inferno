## Assess thinning (from LaplacesDemon::LaplacesDemon)
proposethinning <- function(x){
    if(is.null(dim(x))){ x <- cbind(x) }
    lx <- nrow(x)
    LIV <- ncol(x)
    acf.rows <- trunc(10 * log10(lx))
    acf.temp <- matrix(1, acf.rows, LIV)
    Rec.Thin <- rep(1, LIV)
    names(Rec.Thin) <- colnames(x)
    for (j in 1:LIV) {
        temp0 <- acf(x[, j], lag.max = acf.rows, plot = FALSE)
        if (length(temp0$acf[-1, 1, 1]) == acf.rows) 
            acf.temp[, j] <- abs(temp0$acf[-1, 1, 1])
        ##ESS1[j] <- LaplacesDemon::ESS(x[, j])
        Rec.Thin[j] <- which(acf.temp[, j] <= 0.1)[1] * 1
    }
    Rec.Thin[which(is.na(Rec.Thin))] <- nrow(acf.temp)
    Rec.Thin
}

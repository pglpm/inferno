## Calculates the empirical covariance (scatter matrix) with factor 1/N
mycov <- function(x){x <- t(x)-colMeans(x)
    (x%*%t(x))/ncol(x)}

logit2 <- function(x){x}

## > system.time(test1 <- posteriorx(y,dpobj)) # without foreach
##    user  system elapsed
## 398.185   1.646  20.421
## > system.time(test0 <- posteriorx(y,dpobj)) # with foreach
##    user  system elapsed
## 504.901   1.580  26.114

posteriorx <- function(x,dpobj){
    x <- unname(x)
    if(is.null(dim(x))){dim(x) <- c(1,length(x))}
    n <- dpobj$n
    
    G <- dpobj$mixingDistribution$priorParameters
    df <- G$nu - length(G$mu0) + 1
    sigma <- solve(G$Lambda) * (G$kappa0 + 1)/(G$kappa0 * df)

    ac <- dpobj$alphaChain
    nsamples <- length(ac)

    clterm <- numeric(nrow(x))
    
    for(i in seq_along(ac)){
        clterm <- clterm +
       LikelihoodFunction(dpobj=dpobj, ind=i)(x) / (ac[i] + n)}

    c(clterm) * n/nsamples + mean(ac/(ac+n)) * dmvt(x=x, delta=G$mu0, sigma=sigma, df=df, type='shifted', log=FALSE)
}

## posteriorx <- function(x,dpobj){
##     x <- unname(x)
##     if(is.null(dim(x))){dim(x) <- c(1,length(x))}
##     n <- dpobj$n
    
##     G <- dpobj$mixingDistribution$priorParameters
##     df <- G$nu - length(G$mu0) + 1
##     sigma <- solve(G$Lambda) * (G$kappa0 + 1)/(G$kappa0 * df)

##     ac <- dpobj$alphaChain
##     nsamples <- length(ac)
    
##    c(foreach(i=seq_len(nsamples), .combine='+')%do%{
##         LikelihoodFunction(dpobj=dpobj, ind=i)(x)/(ac[i] + n)
##         }) * n/nsamples + mean(ac/(ac+n)) * dmvt(x=x, delta=G$mu0, sigma=sigma, df=df, type='shifted', log=FALSE)
## }

## posteriorxalt2 <- function(x,dpobj){
##     x <- unname(x)
##     if(is.null(dim(x))){dim(x) <- c(1,length(x))}
##     n <- dpobj$n
    
##     G <- dpobj$mixingDistribution$priorParameters
##     df <- G$nu - length(G$mu0) + 1
##     sigma <- solve(G$Lambda) * (G$kappa0 + 1)/(G$kappa0 * df)

##     ac <- dpobj$alphaChain
    
##     rowMeans(vapply(seq_along(ac), function(i){
##         c(LikelihoodFunction(dpobj=dpobj, ind=i)(x)) / (ac[i] + n)
##     },numeric(nrow(x)))) * n +
##         mean(ac/(ac+n)) * dmvt(x=x, delta=G$mu0, sigma=sigma, df=df, type='shifted', log=FALSE)
## }

## ## More precise for the log-likelihood but slower
## myfitfast <- function(dpobj,its,progressBar){
##     x <- dpobj$data
##     n <- dpobj$n
##     dpobj$alphaChain <- NULL
##     dpobj$likelihoodChain <- NULL
##     dpobj$weightsChain <- NULL
##     dpobj$clusterParametersChain <- NULL
##     dpobj$priorParametersChain <- NULL
##     dpobj$labelsChain <- NULL
##     likelihoodChain <- numeric(its)
##     alphaChain <- numeric(its)
##     weightsChain <- vector('list',its)
##     clusterParametersChain <- vector('list',its)

##     for(i in 1:its){
##         dpobj <- UpdateAlpha(ClusterParameterUpdate(ClusterComponentUpdate(dpobj)))
##         params <- dpobj$clusterParameters
##         alphaChain[i] <- dpobj$alpha
##         weightsChain[[i]] <- dpobj$pointsPerCluster / n
##         clusterParametersChain[[i]] <- params
##         labs <- dpobj$clusterLabels
##         likelihoodChain[i]=foreach(l=seq_len(n), .combine='c')%do%{
##             dmvnorm(x[l,], mean=params$mu[,,labs[l]], sig=params$sig[,,labs[l]], log=TRUE,checkSymmetry=FALSE)}
##     }
##     dpobj$alphaChain <- alphaChain
##     dpobj$likelihoodChain <- likelihoodChain
##     dpobj$weightsChain <- weightsChain
##     dpobj$clusterParametersChain <- clusterParametersChain
##     dpobj$labelsChain <- lapply(weightsChain,function(x){seq_len(length(x))})
    
##     dpobj}
## ## > set.seed(2) ; system.time(test0 <- myfitfast(dp,10000))#with foreach-mvnorm
## ##      user    system   elapsed
## ## 17164.580    35.308   914.489
## ## > set.seed(2) ; system.time(test1 <- myfitfast(dp,10000)) #with LikelihoodDP
## ##      user    system   elapsed
## ## 15054.327    47.341   775.507

myfitfastnoalpha <- function(dpobj,its,progressBar){
    x <- dpobj$data
    n <- dpobj$n
    dpobj$alphaChain <- NULL
    dpobj$likelihoodChain <- NULL
    dpobj$weightsChain <- NULL
    dpobj$clusterParametersChain <- NULL
    dpobj$priorParametersChain <- NULL
    dpobj$labelsChain <- NULL
    likelihoodChain <- numeric(its)
    weightsChain <- vector('list',its)
    clusterParametersChain <- vector('list',its)

    for(i in 1:its){
        dpobj <- ClusterParameterUpdate(ClusterComponentUpdate(dpobj))
        
        weightsChain[[i]] <- dpobj$pointsPerCluster / n
        clusterParametersChain[[i]] <- dpobj$clusterParameters
        likelihoodChain[i]=sum(diag(log(LikelihoodDP(dpobj))))
    }
    dpobj$alphaChain <- rep(dpobj$alpha,its)
    dpobj$likelihoodChain <- likelihoodChain
    dpobj$weightsChain <- weightsChain
    dpobj$clusterParametersChain <- clusterParametersChain
    dpobj$labelsChain <- lapply(weightsChain,function(x){seq_len(length(x))})
    
    dpobj}

myfitfastalpha <- function(dpobj,its,alphaits=ncol(dpobj$data),progressBar){
    x <- dpobj$data
    n <- dpobj$n
    dpobj$alphaChain <- NULL
    dpobj$likelihoodChain <- NULL
    dpobj$weightsChain <- NULL
    dpobj$clusterParametersChain <- NULL
    dpobj$priorParametersChain <- NULL
    dpobj$labelsChain <- NULL
    likelihoodChain <- numeric(its)
    alphaChain <- numeric(its)
    weightsChain <- vector('list',its)
    clusterParametersChain <- vector('list',its)

    for(i in 1:its){
        dpobj <- ClusterParameterUpdate(ClusterComponentUpdate(dpobj))
        if(i %% alphaits == 0){dpobj <- UpdateAlpha(dpobj)}
        alphaChain[i] <- dpobj$alpha
        weightsChain[[i]] <- dpobj$pointsPerCluster / n
        clusterParametersChain[[i]] <- dpobj$clusterParameters
        likelihoodChain[i]=sum(diag(log(LikelihoodDP(dpobj))))
    }
    dpobj$alphaChain <- alphaChain
    dpobj$likelihoodChain <- likelihoodChain
    dpobj$weightsChain <- weightsChain
    dpobj$clusterParametersChain <- clusterParametersChain
    dpobj$labelsChain <- lapply(weightsChain,function(x){seq_len(length(x))})
    
    dpobj}

myfitfast <- function(dpobj,its,progressBar){
    x <- dpobj$data
    n <- dpobj$n
    dpobj$alphaChain <- NULL
    dpobj$likelihoodChain <- NULL
    dpobj$weightsChain <- NULL
    dpobj$clusterParametersChain <- NULL
    dpobj$priorParametersChain <- NULL
    dpobj$labelsChain <- NULL
    likelihoodChain <- numeric(its)
    alphaChain <- numeric(its)
    weightsChain <- vector('list',its)
    clusterParametersChain <- vector('list',its)

    for(i in 1:its){
        dpobj <- UpdateAlpha(ClusterParameterUpdate(ClusterComponentUpdate(dpobj)))
        
        alphaChain[i] <- dpobj$alpha
        weightsChain[[i]] <- dpobj$pointsPerCluster / n
        clusterParametersChain[[i]] <- dpobj$clusterParameters
        likelihoodChain[i]=sum(diag(log(LikelihoodDP(dpobj))))
    }
    dpobj$alphaChain <- alphaChain
    dpobj$likelihoodChain <- likelihoodChain
    dpobj$weightsChain <- weightsChain
    dpobj$clusterParametersChain <- clusterParametersChain
    dpobj$labelsChain <- lapply(weightsChain,function(x){seq_len(length(x))})
    
    dpobj}

## Calculates the prior marginal density of DP mixture distribution at x
priorsamplemargdensity <- function(x,dpobj,dims=1:2){
    if(is.null(nrow(x))){dim(x) <- c(1,length(x))}
    pc <- PriorClusters(dpobj)
    w <- pc$weights
    mu <- pc$params[[1]][1,dims,,drop=FALSE]
    sigma <- pc$params[[2]][dims,dims,,drop=FALSE]
    ## foreach(i=seq_along(w),.combine='+')%do%{
    ##     w[i] * dmvnorm(x,mean=mu[,,i],sigma=rbind(sigma[,,i]),log=FALSE,checkSymmetry = FALSE)}
    clusterterm <- numeric(nrow(x))
    for(i in seq_along(w)){
        clusterterm <- clusterterm +
            w[i] * dmvnorm(x,mean=mu[,,i],sigma=rbind(sigma[,,i]),log=FALSE,checkSymmetry = FALSE)}
    clusterterm
}

## Draws samples from the prior marginal density of DP mixture distribution
priorsamplemargscatter <- function(dpobj,nsamples=100,dims=1:2){
    if(is.null(nrow(x))){dim(x) <- c(1,length(x))}
    pc <- PriorClusters(dpobj)
    w <- pc$weights
    mu <- pc$params[[1]][1,dims,,drop=FALSE]
    sigma <- pc$params[[2]][dims,dims,,drop=FALSE]
    wsamples <- sample(seq_along(w),size=nsamples,replace=TRUE,prob=w)
    foreach(i=wsamples,.combine='rbind')%do%{
        rmvnorm(n=1,mean=mu[,,i],sigma=rbind(sigma[,,i]),checkSymmetry = FALSE)
    }}

fold <- function(x){x <- x+1
   (-1)^(x %/% 2) * ((x %% 2) - 1)
}

## Calculates the posterior marginal density of DP mixture distribution at x
posteriorsamplemargdensity <- function(x,dpobj,dims=1:2){
    if(is.null(nrow(x))){dim(x) <- c(1,length(x))}
    nsamples <- length(dpobj$alphaChain)
    s <- sample(x=seq_len(nsamples),size=1)

    pc <- PosteriorClusters(dpobj,s)
    w <- pc$weights
    mu <- pc$params[[1]][1,dims,,drop=FALSE]
    sigma <- pc$params[[2]][dims,dims,,drop=FALSE]
    ## foreach(i=seq_along(w),.combine='+')%do%{
    ##     w[i] * dmvnorm(x,mean=mu[,,i],sigma=rbind(sigma[,,i]),log=FALSE,checkSymmetry = FALSE)}
    clusterterm <- numeric(nrow(x))
    for(i in seq_along(w)){
        clusterterm <- clusterterm +
            w[i] * dmvnorm(x,mean=mu[,,i],sigma=rbind(sigma[,,i]),log=FALSE,checkSymmetry = FALSE)}
    clusterterm
}

## Draws samples from the posterior marginal density of DP mixture distribution



## myfitc <- function(dpobj,its,progressBar=TRUE){
##     x <- dpobj$data
##     n <- dpobj$n
##     dpobj$alphaChain <- NULL
##     dpobj$likelihoodChain <- NULL
##     dpobj$weightsChain <- NULL
##     dpobj$clusterParametersChain <- NULL
##     dpobj$priorParametersChain <- NULL
##     dpobj$labelsChain <- NULL
##     likelihoodChain <- numeric(its)
##     alphaChain <- numeric(its)
##     weightsChain <- vector('list',its)
##     clusterParametersChain <- vector('list',its)
##     if (progressBar){
##         pb <- txtProgressBar(min=0, max=its, width=50, char="-", style=3)
##     }
##     for(i in 1:its){
##         if (progressBar){setTxtProgressBar(pb, i)}
##         dpobj <- UpdateAlpha(ClusterParameterUpdate(ClusterComponentUpdate(dpobj)))
        
##         alphaChain[i] <- dpobj$alpha
##         weightsChain[[i]] <- dpobj$pointsPerCluster / n
##         clusterParametersChain[[i]] <- dpobj$clusterParameters
##         likelihoodChain[i]=sum(log(LikelihoodFunction(dpobj)(x)))
##     }
##     dpobj$alphaChain <- alphaChain
##     dpobj$likelihoodChain <- likelihoodChain
##     dpobj$weightsChain <- weightsChain
##     dpobj$clusterParametersChain <- clusterParametersChain
##     dpobj$labelsChain <- lapply(weightsChain,function(x){seq_len(length(x))})
    
##     dpobj}

## myfital <- function(dpobj,its,thin=ncol(dpobj$data),progressBar=TRUE){
##     x <- dpobj$data
##     n <- dpobj$n
##     dpobj$alphaChain <- NULL
##     dpobj$likelihoodChain <- NULL
##     dpobj$weightsChain <- NULL
##     dpobj$clusterParametersChain <- NULL
##     dpobj$priorParametersChain <- NULL
##     dpobj$labelsChain <- NULL
## likelihoodChain <- numeric(its)
## alphaChain <- numeric(its)
## weightsChain <- vector('list',its)
##     clusterParametersChain <- vector('list',its)
##     if (progressBar){
##     pb <- txtProgressBar(min=0, max=its, width=50, char="-", style=3)
##   }
## for(i in 1:its){
##     if (progressBar){setTxtProgressBar(pb, i)}
##     dpobj <- ClusterParameterUpdate(ClusterComponentUpdate(dpobj))
##     if(i %% thin == 0){dpobj <- UpdateAlpha(dpobj)}
    
##         alphaChain[i] <- dpobj$alpha
##         weightsChain[[i]] <- dpobj$pointsPerCluster / n
##         clusterParametersChain[[i]] <- dpobj$clusterParameters
##     likelihoodChain[i]=sum(log(LikelihoodFunction(dpobj)(x)))
##     }
##     dpobj$alphaChain <- alphaChain
##     dpobj$likelihoodChain <- likelihoodChain
##     dpobj$weightsChain <- weightsChain
##     dpobj$clusterParametersChain <- clusterParametersChain
##     dpobj$labelsChain <- lapply(weightsChain,function(x){seq_len(length(x))})

##     dpobj}

## myfitalfast <- function(dpobj,its,thin=ncol(dpobj$data)){
##     x <- dpobj$data
##     n <- dpobj$n
##     dpobj$alphaChain <- NULL
##     dpobj$likelihoodChain <- NULL
##     dpobj$weightsChain <- NULL
##     dpobj$clusterParametersChain <- NULL
##     dpobj$priorParametersChain <- NULL
##     dpobj$labelsChain <- NULL
## likelihoodChain <- numeric(its)
## alphaChain <- numeric(its)
## weightsChain <- vector('list',its)
##     clusterParametersChain <- vector('list',its)

##     for(i in 1:its){
##     dpobj <- ClusterParameterUpdate(ClusterComponentUpdate(dpobj))
##     if(i %% thin == 0){dpobj <- UpdateAlpha(dpobj)}
    
##         alphaChain[i] <- dpobj$alpha
##         weightsChain[[i]] <- dpobj$pointsPerCluster / n
##         clusterParametersChain[[i]] <- dpobj$clusterParameters
##     likelihoodChain[i]=sum(log(LikelihoodFunction(dpobj)(x)))
##     }
##     dpobj$alphaChain <- alphaChain
##     dpobj$likelihoodChain <- likelihoodChain
##     dpobj$weightsChain <- weightsChain
##     dpobj$clusterParametersChain <- clusterParametersChain
##     dpobj$labelsChain <- lapply(weightsChain,function(x){seq_len(length(x))})

##     dpobj}

## myfitthin <- function(dpobj,its,thin=ncol(dpobj$data),progressBar=TRUE){
##     x <- dpobj$data
##     n <- dpobj$n
##     dpobj$alphaChain <- NULL
##     dpobj$likelihoodChain <- NULL
##     dpobj$weightsChain <- NULL
##     dpobj$clusterParametersChain <- NULL
##     dpobj$priorParametersChain <- NULL
##     dpobj$labelsChain <- NULL
## likelihoodChain <- numeric(its)
## alphaChain <- numeric(its)
## weightsChain <- vector('list',its)
##     clusterParametersChain <- vector('list',its)
##     if (progressBar){
##     pb <- txtProgressBar(min=0, max=its, width=50, char="-", style=3)
##   }
## kk <- 0
## for(i in 1:(its*thin)){
##     dpobj <- ClusterParameterUpdate(ClusterComponentUpdate(dpobj))

##     if(i %% thin == 0){kk <- kk+1
##         if (progressBar){setTxtProgressBar(pb, kk)}
##         dpobj <- UpdateAlpha(dpobj)
##         alphaChain[kk] <- dpobj$alpha
##         weightsChain[[kk]] <- dpobj$pointsPerCluster / n
##         clusterParametersChain[[kk]] <- dpobj$clusterParameters
##         likelihoodChain[kk]=sum(log(LikelihoodFunction(dpobj)(x)))
##     }
## }
##     dpobj$alphaChain <- alphaChain
##     dpobj$likelihoodChain <- likelihoodChain
##     dpobj$weightsChain <- weightsChain
##     dpobj$clusterParametersChain <- clusterParametersChain
##     dpobj$labelsChain <- lapply(weightsChain,function(x){seq_len(length(x))})

##     dpobj}


## > set.seed(1);system.time(myfit(dp,500,4,progressBar=FALSE))
##    user  system elapsed
## 691.169   1.551  35.308
## > set.seed(1);system.time(myfitc(dp,500,progressBar=FALSE))
##    user  system elapsed
## 680.060   1.635  34.811
## > set.seed(1);system.time(myfitfast(dp,500))
##    user  system elapsed
## 667.131   1.525  34.057
## > set.seed(1);system.time(Fit(dp,500,progressBar=FALSE))
##    user  system elapsed
## 854.869   2.599  43.720


  nalpha <- length(minalpha:maxalpha)

sd <- 3
ncsamples <- 2^14
nsamples <- 2^14
  test3 <- t(sapply(1:nsamples, function(i){
    cat('\r',i)
  alpha <- 2^(minalpha -1L +extraDistr::rcat(1, prob=rep(1/nalpha,nalpha)))/nclusters
  W <- extraDistr::rdirichlet(n=1, alpha=rep(alpha,nclusters))
  means <- rnorm(nclusters, mean = mean, sd = sd)
  sds <- sqrt(
    nimble::rinvgamma(nclusters, shape = shapelo, rate = nimble::rinvgamma(nclusters, shape = shapehi, rate = rate)))
  clu <- extraDistr::rcat(n=ncsamples,prob=W)
    quantile(rnorm(n=ncsamples, mean=means[clu], sd=sds[clu]), (1:3)/4,type=6)
  }
  ))
  cat('\n')





  shapehi <- shapelo <- 0.5
  rate <- 1
  mean <- 0
  minalpha <- -3L
  nalpha <- 7L
  sdlist <- c(0.5,0.75,1,2,3,5)
  nclusters <- 64
  ##
  ncsamples <- 2^14
  nsamples <- 2^13
  bqsamples <- foreach(sd=sdlist, .combine=c)%do%{
    cat('\r',sd,'    ')
    ##tfun <- readRDS(paste0('../exclude/proto_package/xQfunction3600_',sd,'.rds'))
    itfun <- readRDS(paste0('../exclude/proto_package/xinvQfunction3600_',sd,'.rds'))
    list(
      apply(sapply(1:nsamples, function(i){
    alpha <- 2^(minalpha -1L +
                extraDistr::rcat(1, prob=rep(1/nalpha,nalpha)))/nclusters
    W <- extraDistr::rdirichlet(n=1, alpha=rep(alpha,nclusters))
    clu <- extraDistr::rcat(n=ncsamples,prob=W)
    means <- rnorm(nclusters, mean = mean, sd = sd)
    sds <- sqrt(nimble::rinvgamma(nclusters, shape = shapelo,
                                  rate = nimble::rinvgamma(nclusters, shape = shapehi, rate = rate)))
    itfun(quantile(rnorm(n=ncsamples, mean=means[clu], sd=sds[clu]),
    (1:3)/4, type=6))
      } ), 1, summary)
    )
  }
  names(bqsamples) <- sdlist
  cat('\n')




  shapehi <- shapelo <- 0.5
  rate <- 1
  mean <- 0
  minalpha <- -3L
  nalpha <- 7L
  sdlist <- c(0.5,0.75,1,2,3,5)
  nclusters <- 64
  ##
  ncsamples <- 2^14
  nsamples <- 2^13
  bqsamples <- foreach(sd=sdlist, .combine=c)%do%{
    cat('\r',sd,'    ')
    ##tfun <- readRDS(paste0('../exclude/proto_package/xQfunction3600_',sd,'.rds'))
    ## itfun <- readRDS(paste0('../exclude/proto_package/xinvQfunction3600_',sd,'.rds'))
    list(
      apply(sapply(1:nsamples, function(i){
    alpha <- 2^(minalpha -1L +
                extraDistr::rcat(1, prob=rep(1/nalpha,nalpha)))/nclusters
    W <- extraDistr::rdirichlet(n=1, alpha=rep(alpha,nclusters))
    clu <- extraDistr::rcat(n=ncsamples,prob=W)
    means <- rnorm(nclusters, mean = mean, sd = sd)
    sds <- sqrt(nimble::rinvgamma(nclusters, shape = shapelo,
                                  rate = nimble::rinvgamma(nclusters, shape = shapehi, rate = rate)))
    out <- quantile(rnorm(n=ncsamples, mean=means[clu], sd=sds[clu]),
    (1:3)/4, type=6)
    out <- c(out, IQR=out[3]-out[1])
    out
      } ), 1, summary)
    )
  }
  names(bqsamples) <- sdlist
  cat('\n')
## > bqsamples
## $`0.5`
##                  25%         50%         75%     IQR.75%
## Min.    -1113.277227 -7.45658137   -1.355479 5.87231e-04
## 1st Qu.    -1.280903 -0.22217465    0.449029 9.25679e-01
## Median     -0.782100  0.00590250    0.799856 1.48607e+00
## Mean       -1.990482  0.00248794    1.997816 3.98830e+00
## 3rd Qu.    -0.433746  0.23233757    1.272219 2.34410e+00
## Max.        1.431023  6.95035595 1118.560151 2.23184e+03
## 
## $`0.75`
##                  25%          50%         75%     IQR.75%
## Min.    -2790.041937 -11.49355125   -2.629746 4.15174e-04
## 1st Qu.    -1.499926  -0.34991758    0.462216 1.03694e+00
## Median     -0.928929   0.00467820    0.922168 1.70295e+00
## Mean       -2.747234  -0.00241487    2.723542 5.47078e+00
## 3rd Qu.    -0.433768   0.34372026    1.471083 2.56961e+00
## Max.        2.211074   6.40368677 2741.613983 5.53166e+03
## 
## $`1`
##                   25%         50%          75%     IQR.75%
## Min.    -12922.441706 -3.91938414    -3.178151 8.26778e-04
## 1st Qu.     -1.656945 -0.45850858     0.452091 1.11633e+00
## Median      -1.014694  0.01042654     1.035171 1.87779e+00
## Mean        -3.278971  0.00924611     3.306188 6.58516e+00
## 3rd Qu.     -0.436721  0.48675828     1.678858 2.79810e+00
## Max.         3.575628  3.61395878 13088.234105 2.60107e+04
## 
## $`2`
##                 25%        50%        75%     IQR.75%
## Min.    -838.417538 -7.5478537  -7.459911 2.08798e-04
## 1st Qu.   -2.565783 -0.8873834   0.465489 1.37890e+00
## Median    -1.481852  0.0311905   1.498367 2.65073e+00
## Mean      -2.451725  0.0332987   2.503813 4.95554e+00
## 3rd Qu.   -0.389754  0.9703781   2.608194 4.01687e+00
## Max.       7.205059  8.4622620 811.995041 1.65041e+03
## 
## $`3`
##                 25%         50%        75%     IQR.75%
## Min.    -614.933556 -10.4501924  -9.905423 9.79869e-05
## 1st Qu.   -3.530139  -1.4296707   0.386981 1.55108e+00
## Median    -1.960042  -0.0257250   1.959510 3.34927e+00
## Mean      -2.541310  -0.0244554   2.534666 5.07598e+00
## 3rd Qu.   -0.378591   1.4021451   3.504946 5.33021e+00
## Max.      10.682938  11.6026875 619.273284 1.19947e+03
## 
## $`5`
##                25%         50%         75%     IQR.75%
## Min.    -799.15169 -16.0838815 -15.5745817 2.95902e-05
## 1st Qu.   -5.35392  -2.4723106   0.0787854 1.80742e+00
## Median    -2.80732  -0.0399473   2.7927887 4.65180e+00
## Mean      -3.37058  -0.0358134   3.2963577 6.66693e+00
## 3rd Qu.   -0.22837   2.4592434   5.3311362 7.89368e+00
## Max.      17.54059  19.6182092 808.7651859 1.60792e+03






  

    cat('\n')
  sd <- 3
ncsamples <- 2^4
nsamples <- 2^4
  test3 <- t(sapply(1:nsamples, function(i){
  alpha <- 2^(minalpha -1L +extraDistr::rcat(1, prob=rep(1/nalpha,nalpha)))/nclusters
cat(i,alpha*nclusters,'\n')
  W <- extraDistr::rdirichlet(n=1, alpha=rep(alpha,nclusters))
  cat(signif(W,2),'\n')
  means <- rnorm(nclusters, mean = mean, sd = sd)
  sds <- sqrt(
    nimble::rinvgamma(nclusters, shape = shapelo, rate = nimble::rinvgamma(nclusters, shape = shapehi, rate = rate)))
  clu <- extraDistr::rcat(n=ncsamples,prob=W)
    quantile(rnorm(n=ncsamples, mean=means[clu], sd=sds[clu]), (1:3)/4,type=6)
  }
  ))
    cat('\n')

    cat('\n')
  sd <- 3
ncsamples <- 2^14
nsamples <- 2^13
  sest3 <- t(sapply(1:nsamples, function(i){
    cat('\r',i)
  alpha <- 2^(minalpha -1L +extraDistr::rcat(1, prob=rep(1/nalpha,nalpha)))/nclusters
  W <- extraDistr::rdirichlet(n=1, alpha=rep(alpha,nclusters))
  means <- rnorm(nclusters, mean = mean, sd = sd)
    sds <- sqrt( nimble::rinvgamma(nclusters, shape = shapelo,
                                   rate = nimble::rinvgamma(nclusters, shape = shapehi, rate = rate)))
  clu <- extraDistr::rcat(n=ncsamples,prob=W)
    rnorm(n=ncsamples, mean=means[clu], sd=sds[clu])
  }
  ))



  shapehi <- shapelo <- 0.5
  rate <- 1
  mean <- 0
  minalpha <- -3L
  nalpha <- 7L
  sdlist <- c(0.5,0.75,1,2,3,5)
  nclusters <- 64
  ##
  nsamples <- 2^4
  xgrid <- seq(0,1,length.out=512)
  for(sd in sdlist){
    tfun <- readRDS(paste0('../exclude/proto_package/xQfunction3600_',sd,'.rds'))
    dtfun <- readRDS(paste0('../exclude/proto_package/xDQfunction3600_',sd,'.rds'))
  pdff(paste0('justatestquantiles_',sd))
  for(i in 1:nsamples){
    alpha <- 2^(minalpha -1L +
                extraDistr::rcat(1, prob=rep(1/nalpha,nalpha)))/nclusters
    W <- extraDistr::rdirichlet(n=1, alpha=rep(alpha,nclusters))
    means <- rnorm(nclusters, mean = mean, sd = sd)
    sds <- sqrt(
      nimble::rinvgamma(nclusters, shape = shapelo, rate = nimble::rinvgamma(nclusters, shape = shapehi, rate = rate)))
    tplot(x=xgrid,y=sapply(tfun(xgrid),function(x)dtfun(x)*sum(W*dnorm(x, mean=means, sd=sds))))
  }
    dev.off()
    }


  
  pdff('justatest')
  xgrid <- seq(-5,5,length.out=128)
     cat('\n')
  sd <- 3
ncsamples <- 2^4
nsamples <- 2^4
  test3 <- t(sapply(1:nsamples, function(i){
  alpha <- 2^(minalpha -1L +extraDistr::rcat(1, prob=rep(1/nalpha,nalpha)))/nclusters
cat(i,alpha*nclusters,'\n')
  W <- extraDistr::rdirichlet(n=1, alpha=rep(alpha,nclusters))
  means <- rnorm(nclusters, mean = mean, sd = sd)
  sds <- sqrt(
    nimble::rinvgamma(nclusters, shape = shapelo, rate = nimble::rinvgamma(nclusters, shape = shapehi, rate = rate)))
  tplot(x=xgrid,y=sapply(xgrid,function(x)sum(W*dnorm(x, mean=means, sd=sds))))
  }
))
  dev.off()
    cat('\n')


}


if (FALSE) {
  nsamples <- 2^24L
  mean <- 0
  sd <- 1
  shapelo <- shapehi <- 0.5
  rate <- 1
  ##
  means <- rnorm(nsamples, mean = mean, sd = sd)
  sds <- sqrt(nimble::rinvgamma(nsamples, shape = shapelo, rate = nimble::rinvgamma(nsamples, shape = shapehi, rate = rate)))

  Qf <- readRDS('Qfunction8192.rds')

  nint <- 128
  seqnint <- (1:(nint - 1)) / nint
  xvals <- Qf(seqnint)
  range(xvals)

  dsamples <- foreach(x = xvals, .combine = c) %dopar% {
    mean(dnorm(x, mean = means, sd = sds))
  }
  dsamples <- (dsamples + rev(dsamples)) / 2
  ##
  ##
  approxq <- approxfun(x = xvals, y = dsamples, yleft = 0, yright = 0)
  file <- paste0('DQfunction', nint)
  saveRDS(approxq, paste0(file, '.rds'))
  ##

  rg <- 2
  tplot(list(xvals, xvals / max(xvals) * rg), dsamples, xlim = c(-rg, rg))


  dq128 <- readRDS('DQfunction128.rds')
  dq256 <- readRDS('DQfunction256.rds')
  dq512 <- readRDS('DQfunction512.rds')
  dq1024 <- readRDS('DQfunction1024.rds')
  dq2048 <- readRDS('DQfunction2048.rds')
  rg <- 0.1
  xgrid <- seq(-rg, rg, length.out = 1024)
  tplot(xgrid, list(dq128(xgrid), dq256(xgrid), dq512(xgrid), dq1024(xgrid), dq2048(xgrid)), lty = 1, alpha = 0.5)

  rg <- 1000
  xgrid <- seq(0, rg, length.out = 1024 * 1000)
  for (fu in list(dq128, dq256, dq512, dq1024, dq2048)) {
    print(any(diff(fu(xgrid)) > 0))
    print(min(fu(xgrid)))
  }

  rg <- 0.01
  xgrid <- 0.05 + seq(-rg, rg, length.out = 1024)
  tplot(xgrid, Qf(xgrid), lty = 1, alpha = 0.5)





  nint2 <- 1024
  seqnint2 <- (1:(nint2 - 1)) / nint2
  xvals2 <- Qf(seqnint2)
  dsamples2 <- diff(seqnint2) / diff(xvals2)
  dsamples2 <- (dsamples2 + rev(dsamples2)) / 2
  xvals2 <- (xvals2[-1] + xvals2[-length(xvals2)]) / 2
  ##
  rg <- 1
  tplot(xvals2, list(dq128(xvals2), dsamples2), xlim = c(-rg, rg), lty = 1, alpha = 0.5)


  tplot(list(xvals, xvals / max(xvals) * 10), xsamples, xlim = c(-10, 10))



  approxq <- approxfun(x = xvals, y = xsamples, yleft = 0, yright = 0)
  file <- paste0('DQfunction', nint)
  saveRDS(approxq, paste0(file, '.rds'))




  fu <- dq2048
  qq2 <- 0.8
  qq2 <- Qf(qq2)
  ygrid <- seq(0, 1, length.out = 1024)
  xgrid <- Qf(ygrid) - qq2
  tplot(ygrid, dnorm(xgrid, mean = 0, sd = 0.5) / fu(xgrid + qq2), ylim = c(0, NA))

  fu <- dnorm
  qq2 <- 0.8
  qq2 <- qnorm(qq2)
  ygrid <- seq(0, 1, length.out = 1024)
  xgrid <- qnorm(ygrid) - qq2
  tplot(ygrid, dnorm(xgrid, mean = 0, sd = 0.5) / fu(xgrid + qq2), ylim = c(0, NA))



  smoothf <- function(x) {
    (2 * x + c(x[-1], 0) + c(0, x[-length(x)])) / 4
  }

  rg <- 1000
  xgrid <- seq(-rg, rg, length.out = 1024 * 1000)
  fu <- dq2048(xgrid)
  print(any(diff(fu) > 0))
  for (i in 1:100000) {
    fu <- smoothf(fu)
    if (!any(diff(fu[xgrid >= 0]) > 0)) {
      print(i)
      break
    }
  }


  rg <- 0.1
  tplot(xgrid, list(fu, dq2048(xgrid), dq128(xgrid)), xlim = c(-rg, rg))
}

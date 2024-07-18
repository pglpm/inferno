#### identity
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



#### Q
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
    out <- itfun(quantile(rnorm(n=ncsamples, mean=means[clu], sd=sds[clu]),
    (1:3)/4, type=6))
    out <- c(out, IQR=out[3]-out[1])
    out
      } ), 1, summary)
    )
  }
  names(bqsamples) <- sdlist
  cat('\n')
## > bqsamples
## $`0.5`
##                                             IQR
## Min.    0.000000 0.0405973 0.136002 0.000176853
## 1st Qu. 0.175535 0.4148987 0.653243 0.288011988
## Median  0.251946 0.5012984 0.747901 0.449082210
## Mean    0.273177 0.5016073 0.728728 0.455551679
## 3rd Qu. 0.349328 0.5893144 0.826274 0.605976345
## Max.    0.869719 1.0000000 1.000000 1.000000000
## 
## $`0.75`
##                                               IQR
## Min.    0.000460718 0.056247 0.115882 2.51061e-06
## 1st Qu. 0.169714976 0.397927 0.636773 2.64932e-01
## Median  0.256528445 0.500739 0.742575 4.27944e-01
## Mean    0.282742199 0.500072 0.718264 4.35522e-01
## 3rd Qu. 0.370882978 0.603771 0.828201 5.83455e-01
## Max.    0.898664483 0.954968 0.999550 9.99089e-01
## 
## $`1`
##                                               IQR
## Min.    0.000000 0.0489413 0.0741271 0.0000463546
## 1st Qu. 0.172277 0.3897240 0.6161444 0.2348679718
## Median  0.268156 0.5049670 0.7365724 0.4013776619
## Mean    0.296413 0.5047609 0.7079359 0.4115230060
## 3rd Qu. 0.387073 0.6201730 0.8313841 0.5591574437
## Max.    0.928300 0.9567687 1.0000000 1.0000000000
## 
## $`2`
##                                               IQR
## Min.    0.000000 0.0386417 0.0501078 0.0000629624
## 1st Qu. 0.176393 0.3613992 0.5695115 0.1661148499
## Median  0.285760 0.5036462 0.7168730 0.3457604982
## Mean    0.323069 0.5032388 0.6810376 0.3579686801
## 3rd Qu. 0.438164 0.6460534 0.8284001 0.5088187322
## Max.    0.946475 0.9552894 1.0000000 1.0000000000
## 
## $`3`
##                                                IQR
## Min.    0.000000 0.00621539 0.0299442 0.0000118767
## 1st Qu. 0.173026 0.34816704 0.5411501 0.1248861572
## Median  0.297474 0.49962554 0.7008272 0.3093579023
## Mean    0.337101 0.50060337 0.6641995 0.3270985807
## 3rd Qu. 0.466007 0.65228184 0.8255625 0.4785464645
## Max.    0.969978 0.98727757 1.0000000 1.0000000000
## 
## $`5`
##                                               IQR
## Min.    0.000000 0.0117332 0.0148066 0.0000155034
## 1st Qu. 0.171989 0.3283963 0.5059223 0.0932568258
## Median  0.300288 0.4968905 0.6898193 0.2781653174
## Mean    0.342523 0.4961157 0.6479655 0.3054423968
## 3rd Qu. 0.478499 0.6610900 0.8228768 0.4645112428
## Max.    0.977670 0.9829611 1.0000000 1.0000000000


#### log
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
    itfun <- exp
    list(
      apply(sapply(1:nsamples, function(i){
    alpha <- 2^(minalpha -1L +
                extraDistr::rcat(1, prob=rep(1/nalpha,nalpha)))/nclusters
    W <- extraDistr::rdirichlet(n=1, alpha=rep(alpha,nclusters))
    clu <- extraDistr::rcat(n=ncsamples,prob=W)
    means <- rnorm(nclusters, mean = mean, sd = sd)
    sds <- sqrt(nimble::rinvgamma(nclusters, shape = shapelo,
                                  rate = nimble::rinvgamma(nclusters, shape = shapehi, rate = rate)))
    out <- itfun(quantile(rnorm(n=ncsamples, mean=means[clu], sd=sds[clu]),
    (1:3)/4, type=6))
    out <- c(out, IQR=out[3]-out[1])
    out
      } ), 1, summary)
    )
  }
  names(bqsamples) <- sdlist
  cat('\n')
## > bqsamples
## $`0.5`
##              25%         50%      75%     IQR.75%
## Min.    0.000000 5.27103e-59 0.262448 0.000400319
## 1st Qu. 0.278012 7.94012e-01 1.522336 0.878929920
## Median  0.457317 9.96770e-01 2.157162 1.582225134
## Mean    0.503885 5.28199e+64      Inf         Inf
## 3rd Qu. 0.653372 1.25111e+00 3.565110 3.084133859
## Max.    4.800065 4.32701e+68      Inf         Inf
## 
## $`0.75`
##              25%         50%      75%     IQR.75%
## Min.    0.000000 2.26109e-10 0.121063 0.000459467
## 1st Qu. 0.224582 7.08354e-01 1.563078 0.983350920
## Median  0.401778 1.00569e+00 2.477196 1.878080813
## Mean    0.515478 1.18518e+00      Inf         Inf
## 3rd Qu. 0.633600 1.40455e+00 4.358566 3.746640723
## Max.    9.396164 7.31841e+01      Inf         Inf
## 
## $`1`
##               25%         50%       75%      IQR.75%
## Min.     0.000000  0.00133436 0.0670076 0.0000294633
## 1st Qu.  0.184907  0.63018482 1.5693953 1.0322576853
## Median   0.361934  1.01072161 2.8237343 2.2177731232
## Mean     0.563184  1.36641482       Inf          Inf
## 3rd Qu.  0.646973  1.59251764 5.4710128 4.7585420676
## Max.    20.195736 21.95618169       Inf          Inf
## 
## $`2`
##                  25%          50%          75%      IQR.75%
## Min.    3.06370e-214   0.00145239  2.38029e-03  5.00691e-04
## 1st Qu.  7.66598e-02   0.39875266  1.57020e+00  1.08969e+00
## Median   2.22238e-01   1.01293503  4.52260e+00  3.74404e+00
## Mean     1.42264e+00   3.99144296 2.95709e+219 2.95709e+219
## 3rd Qu.  6.53622e-01   2.54813590  1.32468e+01  1.17528e+01
## Max.     2.29886e+02 676.51900516 2.42245e+223 2.42245e+223
## 
## $`3`
##                 25%         50%          75%     IQR.75%
## Min.    0.00000e+00 5.79093e-05  0.000074232 2.89858e-05
## 1st Qu. 3.09851e-02 2.54295e-01  1.487824200 1.09420e+00
## Median  1.45111e-01 1.02832e+00  7.286997731 6.11139e+00
## Mean    1.28556e+01 5.27053e+01          Inf         Inf
## 3rd Qu. 6.94809e-01 4.41779e+00 34.891196470 3.03108e+01
## Max.    1.33236e+04 8.57140e+04          Inf         Inf
## 
## $`5`
##                 25%         50%         75%     IQR.75%
## Min.    0.00000e+00 6.68029e-09 6.78711e-09 2.08577e-10
## 1st Qu. 4.95433e-03 8.30395e-02 1.04581e+00 7.83668e-01
## Median  5.81255e-02 9.31836e-01 1.53210e+01 1.28771e+01
## Mean    5.40848e+03 1.51973e+04         Inf         Inf
## 3rd Qu. 8.31693e-01 1.05061e+01 1.92464e+02 1.67218e+02
## Max.    2.90775e+07 2.94194e+07         Inf         Inf


























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

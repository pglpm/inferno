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
    out <- rnorm(n=ncsamples, mean=means[clu], sd=sds[clu])
    c(quantile(out, (1:3)/4, type=6), IQR=IQR(out), MAD=mad(out))
      } ), 1, summary)
    )
  }
  names(bqsamples) <- sdlist
  cat('\n')
## Id
## $`0.5`
##                  25%          50%         75%         IQR         MAD
## Min.    -2455.281906 -16.45779758   -1.522157 5.16965e-04 3.83178e-04
## 1st Qu.    -1.293719  -0.21954952    0.425965 9.12182e-01 6.49637e-01
## Median     -0.778540   0.00440033    0.803288 1.46699e+00 1.07053e+00
## Mean       -2.722123   0.00693539    2.753016 5.47399e+00 4.04379e+00
## 3rd Qu.    -0.422418   0.23549166    1.315688 2.39788e+00 1.76707e+00
## Max.        1.100349  25.46162442 2428.086709 4.88279e+03 3.61999e+03
## 
## $`0.75`
##                 25%         50%        75%         IQR         MAD
## Min.    -888.742473 -4.88716057  -1.910988 7.61462e-04 5.64744e-04
## 1st Qu.   -1.472718 -0.34297028   0.445659 1.02680e+00 7.09972e-01
## Median    -0.907200  0.00739862   0.894602 1.67631e+00 1.20146e+00
## Mean      -2.313046  0.00638922   2.318621 4.63064e+00 3.40412e+00
## 3rd Qu.   -0.441822  0.34874156   1.471164 2.54571e+00 1.85910e+00
## Max.       2.217239  6.67591747 891.189039 1.77915e+03 1.31827e+03
## 
## $`1`
##                  25%         50%         75%         IQR         MAD
## Min.    -2531.134508 -3.27826562   -2.966249 6.70205e-05 4.96438e-05
## 1st Qu.    -1.656276 -0.45695319    0.435532 1.08462e+00 7.33766e-01
## Median     -1.023775  0.00699411    1.041844 1.86562e+00 1.32472e+00
## Mean       -2.395823  0.00793756    2.418870 4.81388e+00 3.52074e+00
## 3rd Qu.    -0.448847  0.46627767    1.700507 2.80022e+00 2.03974e+00
## Max.        3.005248  3.75440255 2451.168681 4.98170e+03 3.69296e+03
## 
## $`2`
##                  25%         50%         75%         IQR         MAD
## Min.    -2228.716903 -7.11410269   -6.864767 1.43381e-04 1.06304e-04
## 1st Qu.    -2.590929 -0.92053984    0.423323 1.39644e+00 8.78831e-01
## Median     -1.511884 -0.00778223    1.516936 2.65890e+00 1.75596e+00
## Mean       -2.834768 -0.00818035    2.851680 5.68546e+00 4.06740e+00
## 3rd Qu.    -0.463697  0.94014301    2.573440 4.03026e+00 2.81677e+00
## Max.        6.508074 10.50531107 2307.416423 4.52000e+03 3.35442e+03
## 
## Id
## $`3`
##                  25%         50%         75%         IQR         MAD
## Min.    -1860.112859 -11.0169576  -10.679016 8.30752e-04 6.15527e-04
## 1st Qu.    -3.555009  -1.4851725    0.371013 1.56323e+00 9.37258e-01
## Median     -1.979507  -0.0109402    1.955058 3.37023e+00 2.10868e+00
## Mean       -3.308923  -0.0368181    3.236807 6.54454e+00 4.56643e+00
## 3rd Qu.    -0.354466   1.4387305    3.500002 5.32415e+00 3.55291e+00
## Max.        8.965375  14.9074137 1747.839160 3.60726e+03 2.67938e+03
## 
## $`5`
##                 25%         50%        75%         IQR         MAD
## Min.    -598.583898 -19.5330165 -18.688365 1.27167e-04 9.43895e-05
## 1st Qu.   -5.402424  -2.4645901   0.119228 1.72744e+00 9.71448e-01
## Median    -2.909748  -0.0350620   2.806298 4.67615e+00 2.67145e+00
## Mean      -3.556294  -0.0505318   3.416222 6.97133e+00 4.61721e+00
## 3rd Qu.   -0.235264   2.4401380   5.403990 7.79460e+00 5.01956e+00
## Max.      16.528143  17.9977636 472.759445 1.07098e+03 8.09052e+02



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
    itfun <- readRDS(paste0('xinvQfunction3600_',sd,'.rds'))
    list(
      apply(sapply(1:nsamples, function(i){
    alpha <- 2^(minalpha -1L +
                extraDistr::rcat(1, prob=rep(1/nalpha,nalpha)))/nclusters
    W <- extraDistr::rdirichlet(n=1, alpha=rep(alpha,nclusters))
    clu <- extraDistr::rcat(n=ncsamples,prob=W)
    means <- rnorm(nclusters, mean = mean, sd = sd)
    sds <- sqrt(nimble::rinvgamma(nclusters, shape = shapelo,
                                  rate = nimble::rinvgamma(nclusters, shape = shapehi, rate = rate)))
    out <- itfun(rnorm(n=ncsamples, mean=means[clu], sd=sds[clu]))
    c(quantile(out, (1:3)/4, type=6), IQR=IQR(out), MAD=mad(out))
      } ), 1, summary)
    )
  }
  names(bqsamples) <- sdlist
cat('\n')
## Q
## $`0.5`
##              25%        50%      75%         IQR         MAD
## Min.    0.000000 0.00608062 0.150194 0.000186814 0.000138509
## 1st Qu. 0.172971 0.41161611 0.655634 0.294951768 0.192916697
## Median  0.248706 0.49917737 0.751867 0.459100861 0.306713512
## Mean    0.270186 0.49958806 0.730574 0.460341181 0.304385445
## 3rd Qu. 0.345728 0.58647375 0.827703 0.613273957 0.410733775
## Max.    0.839770 0.92887774 1.000000 1.000000000 0.724931767
## 
## $`0.75`
##                 25%      50%      75%         IQR         MAD
## Min.    0.000305386 0.110495 0.129007 0.000358519 0.000265973
## 1st Qu. 0.172306167 0.395612 0.631818 0.263272693 0.163585991
## Median  0.260143628 0.501273 0.739157 0.419786539 0.278599779
## Mean    0.285047180 0.501362 0.715938 0.430845648 0.278743656
## 3rd Qu. 0.368953532 0.604742 0.825069 0.578142895 0.382263304
## Max.    0.881607153 0.957186 0.999678 0.999372456 0.735049190
## 
## $`1`
##              25%       50%      75%          IQR          MAD
## Min.    0.000000 0.0838557 0.100599 0.0000921514 0.0000683179
## 1st Qu. 0.168445 0.3827362 0.603145 0.2250336263 0.1405717113
## Median  0.262757 0.4974913 0.732263 0.3979373465 0.2560869034
## Mean    0.292924 0.4982269 0.702899 0.4099297951 0.2621473407
## 3rd Qu. 0.385255 0.6139892 0.827434 0.5614411662 0.3666670841
## Max.    0.921052 0.9287755 1.000000 1.0000000000 0.7378825357
## 
## $`2`
##              25%       50%       75%        IQR          MAD
## Min.    0.000000 0.0402894 0.0489904 0.00019374 0.0000155096
## 1st Qu. 0.172572 0.3548342 0.5597833 0.16075118 0.0946237322
## Median  0.285517 0.4952934 0.7109782 0.34120642 0.2049664192
## Mean    0.321664 0.4982858 0.6763986 0.35469083 0.2196215338
## 3rd Qu. 0.436445 0.6421251 0.8270260 0.50525526 0.3206806232
## Max.    0.952641 0.9577334 1.0000000 1.00000000 0.7344296937
## 
## Q
## $`3`
##              25%        50%       75%        IQR          MAD
## Min.    0.000000 0.00499336 0.0260572 0.00017074 0.0000640864
## 1st Qu. 0.173673 0.34722078 0.5308651 0.12171627 0.0675326383
## Median  0.298940 0.49726945 0.7015743 0.30586743 0.1732775174
## Mean    0.336147 0.49819724 0.6613702 0.32518154 0.1957611123
## 3rd Qu. 0.459672 0.65108102 0.8227428 0.48216344 0.2988670003
## Max.    0.971607 0.97241515 1.0000000 1.00000000 0.7388856264
## 
## $`5`
##              25%        50%       75%          IQR          MAD
## Min.    0.000000 0.00409817 0.0149372 0.0000599431 0.0000382326
## 1st Qu. 0.178760 0.33739152 0.5171224 0.0917887075 0.0501137789
## Median  0.310841 0.49954438 0.6894968 0.2699893285 0.1380206971
## Mean    0.350638 0.49958762 0.6497577 0.2990835921 0.1738736241
## 3rd Qu. 0.492464 0.66464342 0.8216897 0.4570610620 0.2724802101
## Max.    0.983087 0.98385615 1.0000000 1.0000000000 0.7243791045


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
    out <- itfun(rnorm(n=ncsamples, mean=means[clu], sd=sds[clu]))
    c(qq <- quantile(out, (1:3)/4, type=6), IQR=IQR(out), MAD=mad(out), sqrt((qq[3]-qq[2])*(qq[2]-qq[1])))
      } ), 1, summary)
    )
  }
  names(bqsamples) <- sdlist
  cat('\n')
## Log
## $`0.5`
##              25%         50%      75%         IQR         MAD
## Min.    0.000000 3.15421e-14 0.281885 0.000467121 4.67643e-14
## 1st Qu. 0.273681 7.93567e-01 1.531835 0.900565271 5.69052e-01
## Median  0.454044 9.98846e-01 2.184936 1.597215158 9.18430e-01
## Mean    0.505527 1.08479e+00      Inf         Inf 9.96471e-01
## 3rd Qu. 0.654784 1.26887e+00 3.600595 3.121903537 1.30287e+00
## Max.    8.711117 1.82051e+01      Inf         Inf 2.69909e+01
## 
## $`0.75`
##              25%         50%       75%          IQR         MAD
## Min.    0.000000 6.53324e-02 0.0828285 0.0000583159 4.29604e-05
## 1st Qu. 0.232448 6.99368e-01 1.5826374 1.0070066138 5.59160e-01
## Median  0.400231 1.00344e+00 2.5245315 1.9472536895 9.64238e-01
## Mean    0.515824 2.97009e+07       Inf          Inf 4.40346e+07
## 3rd Qu. 0.635812 1.41393e+00 4.3792231 3.7565377249 1.47154e+00
## Max.    7.651265 2.43310e+11       Inf          Inf 3.60731e+11
## 
## $`1`
##               25%        50%       75%        IQR          MAD
## Min.     0.000000  0.0184723 0.0200941 0.00130665  0.000951792
## 1st Qu.  0.186537  0.6256158 1.5511706 1.01932248  0.523126814
## Median   0.350297  0.9791389 2.7977406 2.23242290  0.966893006
## Mean     0.567280  1.3593252       Inf        Inf  1.329931841
## 3rd Qu.  0.625702  1.5850092 5.3642381 4.70179412  1.664814094
## Max.    34.466069 34.7956124       Inf        Inf 29.377214311
## 
## $`2`
##                  25%         50%         75%          IQR         MAD
## Min.       0.0000000 4.38489e-04  0.00167189  0.000559134 1.67578e-04
## 1st Qu.    0.0718221 3.83064e-01  1.54524792  1.111254234 3.42353e-01
## Median     0.2192162 1.00261e+00  4.61477071  3.870736835 1.05908e+00
## Mean       1.8147808 1.84911e+02         Inf          Inf 2.71756e+02
## 3rd Qu.    0.6385548 2.65251e+00 13.71913282 12.253952219 2.83433e+00
## Max.    2507.5984353 1.47647e+06         Inf          Inf 2.18902e+06
## 
## Log
## $`3`
##                 25%         50%         75%         IQR         MAD
## Min.    0.00000e+00 2.47971e-06 2.60627e-06 2.42912e-07 1.79935e-07
## 1st Qu. 3.15319e-02 2.36556e-01 1.44954e+00 1.03847e+00 2.26630e-01
## Median  1.47985e-01 1.00596e+00 7.04187e+00 5.88993e+00 1.08692e+00
## Mean    2.47405e+01 6.43355e+01         Inf         Inf 6.68471e+01
## 3rd Qu. 7.14161e-01 4.25040e+00 3.36054e+01 2.92318e+01 4.50420e+00
## Max.    1.30007e+05 3.30255e+05         Inf         Inf 3.66165e+05
## 
## $`5`
##                 25%         50%         75%         IQR         MAD
## Min.    0.00000e+00 2.08208e-09 2.72504e-09 1.12488e-09 7.99765e-10
## 1st Qu. 4.83938e-03 9.61178e-02 1.16826e+00 9.00553e-01 9.37142e-02
## Median  5.96916e-02 9.77205e-01 1.64344e+01 1.40723e+01 1.11014e+00
## Mean    5.32045e+03 1.13091e+04         Inf         Inf 9.04111e+03
## 3rd Qu. 7.62024e-01 1.13797e+01 1.96645e+02 1.75829e+02 1.19130e+01
## Max.    2.21538e+07 2.65484e+07         Inf         Inf 1.20896e+07



## Id
## $`3`
##                  25%         50%         75%         IQR         MAD
## Min.    -1860.112859 -11.0169576  -10.679016 8.30752e-04 6.15527e-04
## 1st Qu.    -3.555009  -1.4851725    0.371013 1.56323e+00 9.37258e-01
## Median     -1.979507  -0.0109402    1.955058 3.37023e+00 2.10868e+00
## Mean       -3.308923  -0.0368181    3.236807 6.54454e+00 4.56643e+00
## 3rd Qu.    -0.354466   1.4387305    3.500002 5.32415e+00 3.55291e+00
## Max.        8.965375  14.9074137 1747.839160 3.60726e+03 2.67938e+03

## Log
## $`3`
##                 25%         50%         75%         IQR         MAD
## Min.    0.00000e+00 2.47971e-06 2.60627e-06 2.42912e-07 1.79935e-07
## 1st Qu. 3.15319e-02 2.36556e-01 1.44954e+00 1.03847e+00 2.26630e-01
## Median  1.47985e-01 1.00596e+00 7.04187e+00 5.88993e+00 1.08692e+00
## Mean    2.47405e+01 6.43355e+01         Inf         Inf 6.68471e+01
## 3rd Qu. 7.14161e-01 4.25040e+00 3.36054e+01 2.92318e+01 4.50420e+00
## Max.    1.30007e+05 3.30255e+05         Inf         Inf 3.66165e+05

## Q
## $`3`
##              25%        50%       75%        IQR          MAD
## Min.    0.000000 0.00499336 0.0260572 0.00017074 0.0000640864
## 1st Qu. 0.173673 0.34722078 0.5308651 0.12171627 0.0675326383
## Median  0.298940 0.49726945 0.7015743 0.30586743 0.1732775174
## Mean    0.336147 0.49819724 0.6613702 0.32518154 0.1957611123
## 3rd Qu. 0.459672 0.65108102 0.8227428 0.48216344 0.2988670003
## Max.    0.971607 0.97241515 1.0000000 1.00000000 0.7388856264





















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

expo <- 12
expo2 <- expo/2
mm <- 0
MM <- 100
aa <- (mm*2^expo2-MM*2^-expo2)/(2^expo2-2^-expo2)
bb <- (MM*2^expo2-mm*2^-expo2)/(2^expo2-2^-expo2)
dd <- exp(log(MM-mm)+log1p(2^-expo)-log1p(-2^-expo))
##
## y <- qlogis((x-aa)/dd)
## x <- plogis(y)*dd+aa
##
## p(x) <- dnorm(y)/dlogis(y)/dd
##
xgrid <- seq(mm,MM,length.out=256)
ygrid <- qlogis((xgrid-aa)/dd)
tplot(x=xgrid,
      y=ygrid)

tplot(x=xgrid,
      y=dnorm(ygrid)/dlogis(ygrid)/dd,
      ylim=c(0,NA))

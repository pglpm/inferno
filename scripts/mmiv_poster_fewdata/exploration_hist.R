
myticksh <- function(x, n=10){
    rg <- range(x, na.rm=T)
    ext <- diff(rg)
    deltas <- c(0.1, 0.2, 0.5, 1) * 10^floor(log10(ext))
    delta <- deltas[which.min(abs(sapply(deltas,function(d){diff(rg%/%d)})-n))]
    (((rg[1]%/%delta)-1):(rg[2]%/%delta + 1))*delta
}



data <- rbeta(10000, 3, 2)

delta <- diff(quantile(data,c(1,3)/4))/8
## med <- quantile(data,0.5)
## breaks <- seq(med - delta*ceiling((med-min(data))/delta),
##               med + delta*ceiling((max(data)-med)/delta),
##               by=delta)
breaks <- myticksh(data, n=range(data)/delta)
freq <- FALSE
hi <- hist(data, breaks=breaks, plot=FALSE)
##
x <- hi$mids
##
myplot(xlim=range(hi$breaks), ylim=range(c(0,hi$density)))
for(i in 1:length(hi$mids)){
    polygon(x=rep(hi$breaks[i:(i+1)], each=2),
            y=c(0,rep(hi$density[i],2),0),
            col=paste0(palette()[1],'88'), border=1#'#555555ff'
            )
}
fivenumaxis(1, data)

matplot(x=NA, y=NA, xlim=0:1, ylim=0:1)
matpoints(x=NA, y=1:2)

testx <- NA
xx <- testx[!is.na(testx) && is.finite(testx)]
length(xx)==0


xdata <- c(0.18,2.9)
myticks(xdata)
length(myticks(xdata))
myplot(xdata,xdata)


xdata <- c(-0.1,3)
n <- 10
    rg <- range(xdata, na.rm=T)
    ext <- diff(rg)
deltas <- c(0.1, 0.2, 0.5, 1) * 10^floor(log10(ext))
    delta <- deltas[which.min(abs(sapply(deltas,function(d){diff(rg%/%d) + (rg[2]%%d > 0)} - n)))]
testt <- ((rg[1]%/%delta):((rg[2]%/%delta) + (rg[2]%%delta > 0)))*delta
rbind(da=range(xdata),ti=range(testt))
testt
sapply(deltas, function(delta){
    c(length(((rg[1]%/%delta):((rg[2]%/%delta) + (rg[2]%%delta > 0)))*delta),
      diff(rg%/%delta) + (rg[2]%%delta > 0))})

delta <- deltas[which.min(abs(sapply(deltas,function(d){diff(rg%/%d)+2})-n))]
(((rg[1]%/%delta)+(rg[1]%%delta > 0)):(rg[2]%/%delta))*delta



xdata <- c(-0.1,3)
    rg <- range(xdata, na.rm=T)
    ext <- diff(rg)
    deltas <- c(0.1, 0.2, 0.5, 1) * 10^floor(log10(ext))
## delta <- deltas[which.min(abs(sapply(deltas,function(d){diff(rg%/%d)+(rg[1]%%d == 0)})-n))]
## (((rg[1]%/%delta)+(rg[1]%%delta > 0)):(rg[2]%/%delta))*delta
    deltas <- c(0.1, 0.2, 0.5, 1) * 10^floor(log10(ext))
sapply(deltas, function(d){
c(diff(rg%/%d)+2, (rg[2]%/%d +1)-((rg[1]%/%d))+1)})

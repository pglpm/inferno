
##################################################################
##################################################################
#### Calibration plots
calibrationplots <- function(condfreqs, priorP, divs=2, distfunction=NULL, title=''){
## distfunction <- function(freq1,freq0){
##     sum(freq1*log(freq1/freq0),na.rm=TRUE)
## }
    if(is.null(distfunction)){
        distfunction<- function(freq0,freq1){sum(freq1*log(freq1/freq0),na.rm=TRUE)}
        }
##
#divs <- 1
##
nodes <- seq(0,1,1/divs)
## centres1 <- nodes[-(divs+1)]+1/3/divs
## centres2 <- nodes[-(divs+(0:1))]+2/3/divs
##
pbins <- rbind(
    foreach(i1=1:(divs+1), .combine=rbind)%:%foreach(i3=1:(divs-i1+2), .combine=rbind)%do%{
    p1 <- nodes[i1]
    p3 <- nodes[i3]
    c(p1, 1-p1-p3, p3)
    })
rownames(pbins) <- paste0('bin',1:nrow(pbins))
#pbins <- normalizem(pbins + 1e-6)
##
## pbins <- rbind(
##     foreach(i1=1:divs, .combine=rbind)%:%foreach(i3=1:(divs-i1+1), .combine=rbind)%do%{
##     p1 <- centres1[i1]
##     p3 <- centres1[i3]
##     c(p1, 1-p1-p3, p3)
## },
##     foreach(i3=1:(divs-1), .combine=rbind)%:%foreach(i1=1:(divs-i3), .combine=rbind)%do%{
##     p1 <- centres2[i1]
##     p3 <- centres2[i3]
##     c(p1, 1-p1-p3, p3)
##     })
## rownames(pbins) <- paste0('bin',1:nrow(pbins))
##
##
freqbins <- 0 * pbins
psample <- normalizem(t(t(condfreqs[,,1]) * priorP))
for(sample in 1:nrow(condfreqs)){
    distances <- apply(pbins,1,function(x){distfunction(psample[sample,],x)})
    bin <- which.min(distances)
    outcome <- testdata[sample,bin_RMSD]
    freqbins[bin,outcome] <- freqbins[bin,outcome] + 1
}
##
wfreqbins <- rowSums(freqbins)
rfreqbins <- freqbins/wfreqbins
##
##
pdf(file=paste0(title,'_calibration_plots_',nrow(pbins),'bins.pdf'),height=11.7,width=11.7)
matplot(x=rbind(nodes[1],nodes[divs+1],nodes[1],nodes[1]), y=rbind(nodes[1],nodes[1],nodes[divs+1],nodes[1]), type='l', lty=1, xlim=c(0,1), ylim=c(0,1), col=mygrey, xlab='p(1)', ylab='p(3)', cex.axis=2, cex.lab=2, cex.main=2, main='distribution of posterior probabilities')
psample <- normalizem(t(t(condfreqs[,,1]) * priorP))
matpoints(x=psample[,1],y=psample[,3],type='p',pch=1,col=mygrey)
##
matplot(x=rbind(nodes[1],nodes[divs+1],nodes[1],nodes[1]), y=rbind(nodes[1],nodes[1],nodes[divs+1],nodes[1]), type='l', lty=1, xlim=c(0,1), ylim=c(0,1), col=mygrey, xlab='p(1)', ylab='p(3)', cex.axis=2, cex.lab=2, cex.main=2, main='binning of probability simplex')
psampled <- rdirichlet(n=10000, alpha=rep(1,3))
for(i in 1:10000){
    dists <- apply(pbins,1,function(x){distfunction(psampled[i,],x)})
    cent <- which.min(dists)
    matpoints(x=psampled[i,1],y=psampled[i,3],type='p',pch=20,col=mypalette[(cent%%7)+1])
}
matpoints(x=pbins[,1],y=pbins[,3],type='p',pch=18,cex=3,col='black')
##
for(bin in order(wfreqbins, decreasing=TRUE)){
    matplot(x=rbind(nodes[1],nodes[divs+1],nodes[1],nodes[1]), y=rbind(nodes[1],nodes[1],nodes[divs+1],nodes[1]), type='l', lty=1, xlim=c(0,1), ylim=c(0,1), col=mygrey, xlab='p(1)', ylab='p(3)', cex.axis=2, cex.lab=2, cex.main=2, main=paste0('W = ',wfreqbins[bin]))
    matpoints(x=pbins[,1], y=pbins[,3], type='p', pch=15, col='black')
    for(i in 1:10000){
        dists <- apply(pbins,1,function(x){distfunction(psampled[i,],x)})
        if(bin==which.min(dists)){
            matpoints(x=psampled[i,1],y=psampled[i,3],type='p',pch=20,col=mygrey)
        }
    }
#    matlines(x=rbind(pbins[bin,1],rfreqbins[bin,1]),y=rbind(pbins[bin,3],rfreqbins[bin,3]),type='l',lty=1,col=myblue)
    matpoints(x=pbins[bin,1],y=pbins[bin,3],type='p',pch=18,cex=3,col='black')
    matpoints(x=rfreqbins[bin,1],y=rfreqbins[bin,3],type='p',pch=20,cex=4,col=myred)
}
dev.off()
}

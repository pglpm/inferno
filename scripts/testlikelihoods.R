nfsamples <- nrow(mcsamples)
subsample <- round(seq(1,nfsamples, length.out=64))
##
ad0 <- cbind(0)
ad1 <- cbind(1)
colnames(ad0) <- colnames(ad1) <- predictands

graphics.off()
pdff(paste0(dirname,'prelim_likelihoods-R',basename,'--',mcmcseed,'-',stage),'a4')
for(v in setdiff(unlist(variate),predictands)){#cat(avar)
    contvar <- varinfo[v,'type'] %in% variatetypes[c('R','L','T','S')]
    rg <- varinfo[v,c('plotmin','plotmax')]
    if(contvar){
        Xgrid <- cbind(seq(rg[1], rg[2], length.out=256))
    }else{
        Xgrid <- seq(varinfo[v,'min'], varinfo[v,'max'], length.out=varinfo[v,'n'])
        Xgrid <- cbind(Xgrid[Xgrid >= rg[1] & Xgrid <= rg[2]])
    }
    colnames(Xgrid) <- v
    ##
    plotsamples0 <- samplesFDistribution(Y=Xgrid, X=ad0, mcsamples=mcsamples, varinfo=varinfo, subsamples=round(seq(1,nrow(mcsamples),length.out=nfsamples)), jacobian=TRUE)
    plotsamples1 <- samplesFDistribution(Y=Xgrid, X=ad1, mcsamples=mcsamples, varinfo=varinfo, subsamples=round(seq(1,nrow(mcsamples),length.out=nfsamples)), jacobian=TRUE)
    ##
    datum0 <- data0[Subgroup_num_==0][[v]]
    datum0 <- datum0[!is.na(datum0)]
    datum1 <- data0[Subgroup_num_==1][[v]]
    datum1 <- datum1[!is.na(datum1)]
##
    par(mfrow=c(1,1))
    if(varinfo[v,'type'] == variatetypes['S']){
        interior <- which(Xgrid < varinfo[v,'max'])
        interiordata0 <- which(datum0 < varinfo[v,'max'])
        interiordata1 <- which(datum1 < varinfo[v,'max'])
    }else{
        interior <- 1:length(Xgrid)
        interiordata0 <- 1:length(datum0)
        interiordata1 <- 1:length(datum1)
    }
    ymax <- max(tquant(apply(plotsamples0[interior,subsample],2,function(x){tquant(x,31/32)}),31/32, na.rm=T), tquant(apply(plotsamples1[interior,subsample],2,function(x){tquant(x,31/32)}),31/32, na.rm=T))
    ##
    ##
    tplot(x=Xgrid[interior], y=plotsamples0[interior,subsample], type='l', col=5, alpha=7/8, lty=1, lwd=2, xlab=paste0(v), ylab=paste0('frequency'), ylim=c(0, ymax), xlim=range(Xgrid), family=family)
    tplot(x=Xgrid[interior], y=plotsamples1[interior,subsample], type='l', col=2, alpha=7/8, lty=1, lwd=2, xlab=paste0(v), ylab=paste0('frequency'), ylim=c(0, ymax), xlim=range(Xgrid), family=family,add=T)
        ##
    tplot(x=(Xgrid[interior]), y=list(rowMeans(plotsamples0, na.rm=T)[interior],rowMeans(plotsamples1, na.rm=T)[interior]), type='l', col=c(1,6), alpha=0.25, lty=c(1,2), lwd=4, add=T)
    legend
    ##
    ##
    if(FALSE){
    histo <- thist(datum0[interiordata0], n=round(length(interiordata0)/64))
    histomax <- max(rowMeans(plotsamples)[interior])/max(histo$density)
    tplot(x=histo$breaks, y=histo$density*histomax, col=6, alpha=15/16, border=NA, lty=1, lwd=1, family=family, ylim=c(0,NA), add=TRUE)
    ##
    histo <- thist(datum1[interiordata1], n=round(length(interiordata1)/64))
    histomax <- max(rowMeans(plotsamples)[interior])/max(histo$density)
    tplot(x=histo$breaks, y=histo$density*histomax, col=5, alpha=15/16, border=NA, lty=1, lwd=1, family=family, ylim=c(0,NA), add=TRUE)
    ##
    }
}
dev.off()

###############################################
## patient-dependent example
###############################################


v <- 'RAVLT_immediate'
Vgrid <- cbind(seq(varinfo[v,'min'],varinfo[v,'max'],length.out=varinfo[v,'n']))
colnames(Vgrid) <- v
##
ad0 <- cbind(0)
ad1 <- cbind(1)
colnames(ad0) <- colnames(ad1) <- predictands
##
plotsamples0 <- samplesFDistribution(Y=Vgrid, X=ad0, mcsamples=mcsamples, varinfo=varinfo, subsamples=1:nrow(mcsamples), jacobian=TRUE)
plotsamples1 <- samplesFDistribution(Y=Vgrid, X=ad1, mcsamples=mcsamples, varinfo=varinfo, subsamples=1:nrow(mcsamples), jacobian=TRUE)

dataad <- data0[[predictands]]
dataad <- dataad[!is.na(dataad)]
base1 <- sum(dataad==1)/length(dataad)

subsample <- round(seq(10,nrow(mcsamples), length.out=64))
graphics.off()
pdff('example_personalized_prognosis')
p1 <- base1
toplot <- plotsamples1*p1/(plotsamples1*p1+plotsamples0*(1-p1))
tplot(x=Vgrid, y=rowMeans(toplot), lty=1, col='#000000', alpha=0.25, lwd=2,
      ylim=c(0,1), ylab='probability of conversion to AD',
      xlab=paste0('measured "',v,'"'))
## tplot(x=Vgrid, y=toplot[,subsample], lty=1, col=7, alpha=0.5, lwd=1,
##       ylim=c(0,1), ylab='probability of conversion to AD',
##       xlab=paste0('measured "',v,'"'))
plotquantiles(x=Vgrid, y=apply(toplot,1,function(x)tquant(x,c(1,15)/16)),
              col='#000000', alpha=7/8)
##
p1 <- 0.15
toplot <- plotsamples1*p1/(plotsamples1*p1+plotsamples0*(1-p1))
tplot(x=Vgrid, y=rowMeans(toplot), lty=1, col=3, alpha=0.25, lwd=2,
      ylim=c(0,1), ylab='probability of conversion to AD',
      xlab=paste0('measured "',v,'"'),add=T)
## tplot(x=Vgrid, y=toplot[,subsample], lty=1, col=7, alpha=0.5, lwd=1,
##       ylim=c(0,1), ylab='probability of conversion to AD',
##       xlab=paste0('measured "',v,'"'))
plotquantiles(x=Vgrid, y=apply(toplot,1,function(x)tquant(x,c(1,15)/16)),
              col=3, alpha=7/8)
##
p1 <- 0.75
toplot <- plotsamples1*p1/(plotsamples1*p1+plotsamples0*(1-p1))
tplot(x=Vgrid, y=rowMeans(toplot), lty=1, col=4, alpha=0.25, lwd=2,
      ylim=c(0,1), ylab='probability of conversion to AD',
      xlab=paste0('measured "',v,'"'),add=T)
## tplot(x=Vgrid, y=toplot[,subsample], lty=1, col=7, alpha=0.5, lwd=1,
##       ylim=c(0,1), ylab='probability of conversion to AD',
##       xlab=paste0('measured "',v,'"'))
plotquantiles(x=Vgrid, y=apply(toplot,1,function(x)tquant(x,c(1,15)/16)),
              col=4, alpha=7/8)
##
legend(x=35, y=1.05,
       legend=paste0('patient with ',c(signif(base1,2),0.15,0.75),' base prob.'),
       bty='n',
       lty=1,lwd=4,col=c(7,3,4))
dev.off()


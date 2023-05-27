plotFsamples <- function(file, mcsamples, auxmetadata, dataset, plotmeans=TRUE, nsubsamples=100, showdata='scatter', parallel=TRUE){

    family <- 'Palatino'
    source('pglpm_plotfunctions.R')
    source('vtransform.R')
    source('samplesFDistribution.R')

    if(nsubsamples=='all'){nsubsamples <- nrow(mcsamples)}
    if(plotmeans){
        mcsubsamples <- round(seq(1,nrow(mcsamples),length.out=nsubsamples))
    }else{
        mcsubsamples <- 1:nrow(mcsamples)
    }
    subsamples <- round(seq(1,length(mcsubsamples),length.out=nsubsamples))

    graphics.off()
    pdff(file, apaper=4)
    par(mfrow=c(1,1))

    for(v in auxmetadata[['name']]){
        varinfo <- as.list(auxmetadata[name==v])
        vtype <- varinfo[['mcmctype']]

        if(vtype %in% c('R','D','C','O')){ #
            if(vtype == 'O'){
                Xgrid <- seq(varinfo[['domainmin']], varinfo[['domainmax']],
                             length.out=varinfo[['Nvalues']])
                Xgrid <- cbind(Xgrid[Xgrid >= varinfo[['plotmin']] & Xgrid <= varinfo[['plotmax']]])
            }else{
                Xgrid <- cbind(seq(
                    varinfo[['plotmin']], varinfo[['plotmax']],
                    length.out=256
                ))
            }
            colnames(Xgrid) <- v
            
            xleft <- Xgrid > varinfo[['censormin']]
            xright <- Xgrid < varinfo[['censormax']]

            plotsamples <- samplesFDistribution(Y=Xgrid, X=NULL, mcsamples=mcsamples, auxmetadata=auxmetadata, subsamples=mcsubsamples, jacobian=TRUE, parallel=parallel)
            
            ymax <- tquant(apply(plotsamples[xleft & xright, subsamples, drop=F],
                                 2,function(x){tquant(x,31/32)}),31/32, na.rm=T)

            dataplot <- FALSE
            ## data plots if required
            if(showdata=='histogram' && !(missing(dataset) || is.null(dataset)) && !all(is.na(dataset[[v]]))){
                datum <- dataset[[v]]
                datum <- datum[!is.na(datum)]
                dleft <- datum > varinfo[['censormin']]
                dright <- datum < varinfo[['censormax']]
                if(vtype == 'O'){
                    nh <- (varinfo[['domainmax']]-varinfo[['domainmin']])/(varinfo[['Nvalues']]-1L)/2
                    nh <- seq(varinfo[['domainmin']]-nh, varinfo[['domainmax']]+nh,
                              length.out=varinfo[['Nvalues']]+1L)
                }else{
                    nh <- max(16,round(length(datum[dleft & dright])/64))
                }
                histo <- thist(datum[dleft & dright],
                               n=nh)
                hleft <- sum(!dleft)/length(datum)
                hright <- sum(!dright)/length(datum)
                
                ymax <- max(ymax, histo$density)

                dataplot <- TRUE
            }

            ## Plot F samples
            tplot(x=Xgrid[xleft & xright], y=plotsamples[xleft & xright,subsamples,drop=F],
                  xlim=range(Xgrid), ylim=c(0,ymax),
                  type='l', lty=1, lwd=2,
                  col=5, alpha=7/8,
                  xlab=v,
                  ylab=paste0('frequency',(if(vtype=='O'){''}else{' density'})),
                  family=family)
            if(any(!(xleft & xright))){
                tplot(x=Xgrid[!(xleft & xright)],
                      y=plotsamples[!(xleft & xright),subsamples,drop=F]*ymax,
                      type='p', pch=2, cex=2,
                      col=5, alpha=7/8,
                      family=family, add=TRUE)
            }
            
            ## Plot F means if required
            if(plotmeans){
                tplot(x=Xgrid[xleft & xright], y=rowMeans(plotsamples[xleft & xright,,drop=F], na.rm=T),
                      type='l', lty=1, lwd=4,
                      col=1, alpha=0.25,
                      add=TRUE)
                if(any(!(xleft & xright))){
                    tplot(x=Xgrid[!(xleft & xright)],
                          y=rowMeans(plotsamples, na.rm=T)[!(xleft & xright)]*ymax,
                          type='p', pch=2, cex=2,
                          col=1, alpha=0.25,
                          lty=1, lwd=3,
                          add=TRUE)
                }
            }
            
            if(dataplot){
                histomax <- 1 # max(rowMeans(plotsamples))/max(histo$density)
                tplot(x=histo$mids, y=histo$density*histomax,
                      xlim=range(Xgrid), ylim=c(0,ymax),
                      type='l', lty=1, lwd=4,
                      col=yellow, alpha=2/4, border=darkgrey, border.alpha=3/4,
                      xlab=v,
                      ylab=paste0('frequency',(if(vtype=='O'){''}else{' density'})),
                      family=family, add=TRUE)

                if(any(!(dleft & dright))){
                    tplot(x=c(varinfo[['censormin']],varinfo[['censormax']]), y=c(hleft,hright)*ymax,
                          type='p', pch=0, cex=2,
                          col=7, alpha=0,
                          lty=1, lwd=5,
                          family=family, add=TRUE)
                }
                fiven <- fivenum(datum)
                abline(v=fiven,col=paste0(palette()[c(2,4,5,4,2)], '44'),lwd=4)
            }
            
            ## nominal or binary variate
        }else{ 

            Xgrid <- cbind(unlist(varinfo[paste0('V',1:varinfo[['Nvalues']])]))
            colnames(Xgrid) <- v
            Ngrid <- vtransform(x=Xgrid, auxmetadata=auxmetadata,
                                Nout='numeric', Bout='numeric')

            plotsamples <- samplesFDistribution(Y=Xgrid, X=NULL, mcsamples=mcsamples, auxmetadata=auxmetadata, subsamples=mcsubsamples, jacobian=TRUE, parallel=parallel)
            
            ymax <- tquant(apply(plotsamples[, subsamples, drop=F],
                                 2,function(x){tquant(x,31/32)}),31/32, na.rm=T)

            dataplot <- FALSE
            ## data plots if required
            if(showdata=='histogram' && !(missing(dataset) || is.null(dataset)) && !all(is.na(dataset[[v]]))){
                datum <- dataset[[v]]
                datum <- datum[!is.na(datum)]
                histo <- as.vector(table(factor(datum, levels=Xgrid)))/length(datum)
                
                ymax <- max(ymax, histo)

                dataplot <- TRUE
            }

            ## Plot F samples
            tplot(x=Ngrid, y=plotsamples[,subsamples,drop=F],
                  xlim=range(Ngrid), ylim=c(0,ymax),
                  xticks=Ngrid, xlabels=Xgrid,
                  type='l', lty=1, lwd=2,
                  col=5, alpha=7/8,
                  xlab=v,
                  ylab='frequency',
                  family=family)
            
            ## Plot F means if required
            if(plotmeans){
                tplot(x=Ngrid, y=rowMeans(plotsamples, na.rm=T),
                      type='l', lty=1, lwd=4,
                      col=1, alpha=0.25,
                      add=TRUE)
            }

            if(dataplot){
                histomax <- 1 # max(rowMeans(plotsamples))/max(histo$density)
                tplot(x=Ngrid, y=histo*histomax,
                      xlim=range(Ngrid), ylim=c(0,ymax),
                      xticks=Ngrid, xlabels=Xgrid,
                      type='l', lty=1, lwd=4,
                      col=yellow, alpha=2/4, border=darkgrey, border.alpha=3/4,
                      xlab=v,
                      ylab='frequency',
                      family=family, add=TRUE)
            }
            
        }

        if(showdata=='scatter' && !(missing(dataset) || is.null(dataset)) && !all(is.na(dataset[[v]]))){
            datum <- dataset[[v]]
            datum <- datum[!is.na(datum)]
            if(!(vtype %in% c('R','D','C','O'))){
            datum <- vtransform(x=matrix(datum,ncol=1,nrow=length(datum),dimnames=list(NULL,v)), auxmetadata=auxmetadata,
                       Nout='numeric', Bout='numeric')
            }
            scatteraxis(side=1, n=NA, alpha=0.5, ext=5,
                        x=datum+runif(length(datum),
                                      min=-min(diff(sort(unique(datum))))/4,
                                      max=min(diff(sort(unique(datum))))/4),
                        col=yellow)
            if(vtype %in% c('R','D','C','O')){
                fiven <- fivenum(datum)
                abline(v=fiven,col=paste0(palette()[c(2,4,5,4,2)], '44'),lwd=4)
            }
        }

    }
    dev.off()
}

plotFsamples <- function(mcsamples, varinfoaux, dataset, file, subsamples, showdata){

    source('pglpm_plotfunctions.R')
    source('vtransform.R')
    source('samplesFDistribution.R')

    graphics.off()
    pdff(file, 'a4')
    par(mfrow=c(1,1))

    for(v in varinfoaux[['name']]){

        varinfo <- as.list(varinfoaux[name==v])
        vtype <- varinfo[['mcmctype']]

        if(vtype %in% c('R','D')){ #
            Xgrid <- cbind(seq(
                varinfo[['plotmin']], varinfo[['plotmax']],
                length.out=256
            ))
            colnames(Xgrid) <- v

            plotsamples <- samplesFDistribution(Y=Xgrid, X=NULL, mcsamples=mcsamples, varinfoaux=varinfoaux, subsamples=subsamples, jacobian=TRUE)
            
            ymax <- tquant(apply(plotsamples[,showsubsample],2,function(x){tquant(x,31/32)}),31/32, na.rm=T)

            addplot <- FALSE
            ## data plots if required
            if(showdata=='histogram'){
                datum <- dataset[[v]]
                datum <- datum[!is.na(datum)]
                histo <- thist(datum, n=32)

                ymax <- max(ymax, histo$density)

                histomax <- 1 # max(rowMeans(plotsamples))/max(histo$density)
                tplot(x=histo$mids, y=histo$density*histomax,
                      xlim=range(Xgrid), ylim=c(0,ymax),
                      lty=1, lwd=4,
                      col=yellow, alpha=2/4, border=darkgrey, border.alpha=3/4,
                      xlab=v, ylab='freq. density',
                      family=family)
                addplot <- TRUE
            }
            tplot(x=Xgrid, y=plotsamples[,subsample],
                  xlim=range(Xgrid), ylim=c(0,ymax),
                  type='l', lty=1, lwd=2,
                  col=5, alpha=7/8,
                  xlab=v, ylab='freq. density',
                  family=family, add=addplot)

            tplot(x=Xgrid, y=rowMeans(plotsamples, na.rm=T),
                  type='l', lty=1, lwd=4,
                  col=1, alpha=0.25,
                  add=TRUE)

        }else if(vtype %in% c('C')){ #
            Xgrid <- cbind(seq(
                varinfo[['plotmin']], varinfo[['plotmax']],
                length.out=256
            ))
            colnames(Xgrid) <- v

            plotsamples <- samplesFDistribution(Y=Xgrid, X=NULL, mcsamples=mcsamples, varinfoaux=varinfoaux, subsamples=subsamples, jacobian=TRUE)
            
            ymax <- tquant(apply(plotsamples[,showsubsample],2,function(x){tquant(x,31/32)}),31/32, na.rm=T)

            addplot <- FALSE
            ## data plots if required
            if(showdata=='histogram'){
                datum <- dataset[[v]]
                datum <- datum[!is.na(datum)]
                histo <- thist(datum, n=32)

                ymax <- max(ymax, histo$density)

                histomax <- 1 # max(rowMeans(plotsamples))/max(histo$density)
                tplot(x=histo$mids, y=histo$density*histomax,
                      xlim=range(Xgrid), ylim=c(0,ymax),
                      lty=1, lwd=4,
                      col=yellow, alpha=2/4, border=darkgrey, border.alpha=3/4,
                      xlab=v, ylab='freq. density',
                      family=family)
                addplot <- TRUE
            }
            tplot(x=Xgrid, y=plotsamples[,subsample],
                  xlim=range(Xgrid), ylim=c(0,ymax),
                  type='l', lty=1, lwd=2,
                  col=5, alpha=7/8,
                  xlab=v, ylab='freq. density',
                  family=family, add=addplot)

            tplot(x=Xgrid, y=rowMeans(plotsamples, na.rm=T),
                  type='l', lty=1, lwd=4,
                  col=1, alpha=0.25,
                  add=TRUE)


        }







            
        }else if(vtype=='O'){
            rg <- seq(
                varinfo[['domainmin']], varinfo[['domainmax']],
                length.out=varinfo[['Nvalues']]
            )
            Xgrid <- rg <- cbind(rg[rg >= varinfo[['plotmin']] & rg <= varinfo[['plotmax']]])
        else{
            Xgrid <- rg <- cbind(unlist(varinfo[paste0('V',1:varinfo[['Nvalues']])]))
        }
        
        plotsamples <- samplesFDistribution(Y=Xgrid, X=NULL, mcsamples=mcsamples, varinfoaux=varinfoaux, subsamples=subsamples, jacobian=TRUE)

        ymax <- tquant(apply(plotsamples[,showsubsample],2,function(x){tquant(x,31/32)}),31/32, na.rm=T)

        ## data plots if required
        if(is.character(showdata) || (is.logical(showdata) && showdata)){

            ## as histogram
            if(showdata=='histogram' || (showdata==TRUE && vtype %in% c('O','N','B'))){
                ## better via tabulation
                if(vtype %in% c('O','N','B')){
                    datum <- vtransform(rg, varinfoaux=varinfoaux, Oout='original', Nout='numeric', Bout='numeric')
                    datum <- datum[!is.na(datum)]
                    histo <- thist(datum, n=seq(min(datum)-0.5,max(datum)+0.5,by=1))
                    ymax <- max(ymax, histo$counts)
                    tplot(x=histo$mids,
                          y=histo$counts,
                          ylim=c(0,ylim),
                          xlim=range(rg))
                }else{
                    datum <- dataset[[v]]
                    datum <- datum[!is.na(datum)]
                    pleft <- pright <- 0
                    if(varinfo[['censormin']] >= varinfo[['plotmin']]){
                        ## check left-censored values
                        vleft <- datum <= varinfo[['censormin']]
                        pleft <- sum(vleft)/length(datum)
                        datum <- datum[!vleft]
                    }
                    if(varinfo[['censormax']] <= varinfo[['plotmax']]){
                        ## check right-censored values
                        vright <- datum >= varinfo[['censormax']]
                        pright <- sum(vright)/length(datum)
                        datum <- datum[!vright]
                    }
                    histox <- thist(datum, n=32)
                    histoy <- histox$density
                    histox <- histox$mids
                }
            }

            ymax <- max(ymax, histoy, pleft, pright)

            tplot(x=)

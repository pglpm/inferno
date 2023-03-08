library('data.table')
library('png')
library('foreach')

Tprob <- 2^-6
params <- c('type','min','max','n','location','scale',
                                  'Q1','Q2','Q3','datamin','datamax','plotmin','plotmax',
                                  'hmean','hsd','hshapeout','hshapein','hvarscale')
##
dt <- fread('~/repositories/ADBayes/worldbrain/scripts/ingrid_data_nogds6.csv')
varnames <- colnames(dt)
##
variate <- list(S=c("TRAASCOR_neuro","TRABSCOR_neuro"),
                L=c("AGE","LRHHC_n_long"),
                B=c("Apoe4_","Gender_num_","Subgroup_num_"),
                I=c("ANARTERR_neuro","AVDEL30MIN_neuro","AVDELTOT_neuro","CATANIMSC_neuro","GDTOTAL_gds","RAVLT_immediate")
                )
if(any(sort(unlist(variate))!=sort(varnames))){warning('variate mismatch')}
varinfo <- matrix(NA, nrow=length(unlist(variate)), ncol=length(params),
                  dimnames=list(unlist(variate),params) )

## integer
varinfo[variate$I,'type'] <- 1L
varinfo[variate$I,'n'] <- c(51L,16L,16L,64L,7L,76L)
varinfo[variate$I,'min'] <- 0L
varinfo[variate$I,'max'] <- varinfo[variate$I,'n']-1L
##
varinfo[variate$I,'location'] <- varinfo[variate$I,'min']
varinfo[variate$I,'scale'] <- (varinfo[variate$I,'max']-varinfo[variate$I,'min'])/(varinfo[variate$I,'n']-1L)
varinfo[variate$I,'plotmin'] <- sapply(variate$I,function(v){
    dat <- dt[[v]]
    max(varinfo[v,'min'], min(dat, na.rm=T) - diff(tquant(dat, c(0.25,0.5))))
    })
varinfo[variate$I,'plotmax'] <- sapply(variate$I,function(v){
    dat <- dt[[v]]
    min(varinfo[v,'max'], max(dat, na.rm=T) + diff(tquant(dat, c(0.5,0.75))))
    })

## binary
varinfo[variate$B,'type'] <- 2L
varinfo[variate$B,'n'] <- 2L
varinfo[variate$B,'min'] <- 0L
varinfo[variate$B,'max'] <- varinfo[variate$B,'n']-1L
##
varinfo[variate$B,'location'] <- varinfo[variate$B,'min']
varinfo[variate$B,'scale'] <- varinfo[variate$B,'max']-varinfo[variate$B,'min']
varinfo[variate$B,'plotmin'] <- varinfo[variate$B,'min']
varinfo[variate$B,'plotmax'] <- varinfo[variate$B,'max']

## real
varinfo[variate$R,'type'] <- 0L
varinfo[variate$R,'n'] <- 0
varinfo[variate$R,'min'] <- -Inf
varinfo[variate$R,'max'] <- +Inf
##
varinfo[variate$R,'location'] <- apply(dt[,variate$R,with=F],2,median,na.rm=T)
varinfo[variate$R,'scale'] <- apply(dt[,variate$R,with=F],2,IQR,na.rm=T)*sd2iqr
varinfo[variate$R,'plotmin'] <- apply(dt[,variate$R,with=F],2,function(x){
    min(x, na.rm=T) - diff(tquant(x, c(0.25,0.5)))
    })
varinfo[variate$R,'plotmax'] <- apply(dt[,variate$R,with=F],2,function(x){
    max(x, na.rm=T) + diff(tquant(x, c(0.75,0.5)))
    })

## logarithmic
varinfo[variate$L,'type'] <- -1L
varinfo[variate$L,'n'] <- 0
varinfo[variate$L,'min'] <- 0
varinfo[variate$L,'max'] <- +Inf
##
varinfo[variate$L,'location'] <- apply(log(dt[,variate$L,with=F]),2,median,na.rm=T)
varinfo[variate$L,'scale'] <- apply(log(dt[,variate$L,with=F]),2,IQR,na.rm=T)*sd2iqr
varinfo[variate$L,'plotmin'] <- apply(dt[,variate$L,with=F],2,function(x){
    max(0, min(x, na.rm=T) - diff(tquant(x, c(0.25,0.5))))
    })
varinfo[variate$L,'plotmax'] <- apply(dt[,variate$L,with=F],2,function(x){
    max(x, na.rm=T) + diff(tquant(x, c(0.5,0.75)))
})

## logarithmic censored
varinfo[variate$S,'type'] <- -3L
varinfo[variate$S,'n'] <- 0
varinfo[variate$S,'min'] <- 0
varinfo[variate$S,'max'] <- c(150L,300L)
##
varinfo[variate$S,'location'] <- apply(log(dt[,variate$S,with=F]),2,median,na.rm=T)
varinfo[variate$S,'scale'] <- apply(log(dt[,variate$S,with=F]),2,IQR,na.rm=T)*sd2iqr
varinfo[variate$S,'plotmin'] <- apply(dt[,variate$S,with=F],2,function(x){
    max(0, min(x, na.rm=T) - diff(tquant(x, c(0.25,0.5))))
    })
varinfo[variate$S,'plotmax'] <- sapply(variate$S,function(v){
    dat <- dt[[v]]
    min(varinfo[v,'max'], max(dat, na.rm=T) + diff(tquant(dat, c(0.5,0.75))))
    })

## doubly bounded
varinfo[variate$T,'type'] <- -2L
varinfo[variate$T,'n'] <- Tprob
varinfo[variate$T,'min'] <- 0L
varinfo[variate$T,'max'] <- c(150L,300L)
varinfo[variate$T,'scale'] <- (varinfo[variate$T,'max']-varinfo[variate$T,'min'])/(1-2*Tprob)
varinfo[variate$T,'location'] <- varinfo[variate$T,'min']-varinfo[variate$T,'n']*varinfo[variate$T,'scale']
varinfo[variate$T,'plotmin'] <- sapply(variate$T,function(v){
    dat <- dt[[v]]
    max(varinfo[v,'min'], min(dat, na.rm=T) - diff(tquant(dat, c(0.25,0.5))))
    })
varinfo[variate$T,'plotmax'] <- sapply(variate$T,function(v){
    dat <- dt[[v]]
    min(varinfo[v,'max'], max(dat, na.rm=T) + diff(tquant(dat, c(0.5,0.75))))
    })


##
varinfo[varnames,'Q1'] <- apply(dt[,..varnames], 2, tquant, 0.25)
varinfo[varnames,'Q2'] <- apply(dt[,..varnames], 2, tquant, 0.5)
varinfo[varnames,'Q3'] <- apply(dt[,..varnames], 2, tquant, 0.75)
varinfo[varnames,'datamin'] <- apply(dt[,..varnames], 2, min, na.rm=T)
varinfo[varnames,'datamax'] <- apply(dt[,..varnames], 2, max, na.rm=T)


varinfo[variate$R,'hmean'] <- 0L
varinfo[variate$L,'hmean'] <- 0L
varinfo[variate$S,'hmean'] <- 0L
varinfo[variate$T,'hmean'] <- 0L
varinfo[variate$I,'hmean'] <- 0L
varinfo[variate$B,'hmean'] <- NA
varinfo[variate$C,'hmean'] <- NA
##
varinfo[variate$R,'hsd'] <- 3 #***
varinfo[variate$L,'hsd'] <- 2
varinfo[variate$S,'hsd'] <- 2
varinfo[variate$T,'hsd'] <- 1
varinfo[variate$I,'hsd'] <- 7/8
varinfo[variate$B,'hsd'] <- NA
varinfo[variate$C,'hsd'] <- NA
##
varinfo[variate$R,'hshapeout'] <- 1
varinfo[variate$L,'hshapeout'] <- 1
varinfo[variate$S,'hshapeout'] <- 1
varinfo[variate$T,'hshapeout'] <- 1
varinfo[variate$I,'hshapeout'] <- 1
varinfo[variate$B,'hshapeout'] <- 1
varinfo[variate$C,'hshapeout'] <- 1
##
varinfo[variate$R,'hshapein'] <- 1
varinfo[variate$L,'hshapein'] <- 1
varinfo[variate$S,'hshapein'] <- 1
varinfo[variate$T,'hshapein'] <- 1
varinfo[variate$I,'hshapein'] <- 1
varinfo[variate$B,'hshapein'] <- 1
varinfo[variate$C,'hshapein'] <- 1
##
varinfo[variate$R,'hvarscale'] <- 1
varinfo[variate$L,'hvarscale'] <- 1
varinfo[variate$S,'hvarscale'] <- 1
varinfo[variate$T,'hvarscale'] <- 1/4
varinfo[variate$I,'hvarscale'] <- 1/4
varinfo[variate$B,'hvarscale'] <- NA
varinfo[variate$C,'hvarscale'] <- NA


predictands <- 'Subgroup_num_'
predictors <- setdiff(unlist(variate),'Subgroup_num_')

write.csv(varinfo, 'varinfo.csv')
saveRDS(varinfo, 'varinfo.rds')
write.table(predictors, 'predictors.csv', row.names=F, col.names=F)
# varinfo <- data.matrix(read.csv('varinfo.csv', row.names=1))








summary(dt)

summary(sapply(colnames(dt),function(v){transfdir(data.matrix(dt[,..v]), varinfo)}))

summary(sapply(colnames(dt),function(v){recjacobian(data.matrix(dt[,..v]), varinfo)}))


sapply(colnames(dt),function(v){
    x <- data.matrix(dt[,..v])
    y <- transfinv(transfdir(x,varinfo),varinfo)
    rdiff <- y/x-1
    rdiff[x==0 & y==0] <- 0
    rdiff[x==0 & y!=0] <- Inf
    max(abs(rdiff),na.rm=T)
})

v <- 'GDTOTAL_gds'
varinfo[v,]

transfdir(0:6,v)
transfinv(0:6,v)

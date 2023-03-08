library('data.table')
library('png')
library('foreach')

sd2iqr <- 0.5/qnorm(0.75)
Tprob <- 2^-6
params <- c('name','type','transform','min','max','n','tmin','tmax',
            'location','scale',
            'Q1','Q2','Q3','datamin','datamax','plotmin','plotmax')
##
dt <- fread('~/repositories/ADBayes/worldbrain/scripts/ingrid_data_nogds6.csv')
varnames <- colnames(dt)
##
variate <- list(##S=c("TRAASCOR_neuro","TRABSCOR_neuro"),
                L=c("AGE","LRHHC_n_long"),
                B=c("Apoe4_","Gender_num_","Subgroup_num_"),
                I=c("ANARTERR_neuro","AVDEL30MIN_neuro","AVDELTOT_neuro","CATANIMSC_neuro","GDTOTAL_gds","RAVLT_immediate","TRAASCOR_neuro","TRABSCOR_neuro")
                )
if(any(sort(unlist(variate))!=sort(varnames))){warning('variate mismatch')}
varinfo <- lapply(1:length(params),function(x){
    x <- rep(NA,length(unlist(variate)))
    names(x) <- unlist(variate)
    x
})
names(varinfo) <- params

for(i in variate){varinfo[['name']][i] <- i}


## integer
varinfo[['type']][variate$I] <- 'I'
varinfo[['transform']][variate$I] <- 'probit'
varinfo[['n']][variate$I] <- c(51L,16L,16L,64L,7L,76L,150L,300L)
varinfo[['min']][variate$I] <- c(0L,0L,0L,0L,0L,0L,1L,1L)
varinfo[['max']][variate$I] <- varinfo[['n']][variate$I]-1L + varinfo[['min']][variate$I]
varinfo[['tmin']][variate$I] <- -Inf
varinfo[['tmax']][variate$I] <- +Inf
varinfo[['location']][variate$I] <- varinfo[['min']][variate$I]
varinfo[['scale']][variate$I] <- (varinfo[['max']][variate$I]-varinfo[['min']][variate$I])/(varinfo[['n']][variate$I]-1L)
##
varinfo[['plotmin']][variate$I] <- transf(x=rbind(sapply(variate$I,function(v){
    dat <- dt[[v]]
    max(varinfo[['min']][v], min(dat, na.rm=T) - diff(tquant(dat, c(0.25,0.5))))
})), varinfo=varinfo, Iout='correct')
varinfo[['plotmax']][variate$I] <- transf(x=rbind(sapply(variate$I,function(v){
    dat <- dt[[v]]
    min(varinfo[['max']][v], max(dat, na.rm=T) + diff(tquant(dat, c(0.5,0.75))))
})), varinfo=varinfo, Iout='correct')
varinfo[['Q1']][variate$I] <- transf(x=rbind(apply(dt[,variate$I,with=F], 2, tquant, 0.25)), varinfo=varinfo, Iout='correct')
varinfo[['Q2']][variate$I] <- transf(x=rbind(apply(dt[,variate$I,with=F], 2, tquant, 0.5)), varinfo=varinfo, Iout='correct')
varinfo[['Q3']][variate$I] <- transf(x=rbind(apply(dt[,variate$I,with=F], 2, tquant, 0.75)), varinfo=varinfo, Iout='correct')
varinfo[['datamin']][variate$I] <- transf(x=rbind(apply(dt[,variate$I,with=F], 2, min, na.rm=T)), varinfo=varinfo, Iout='correct')
varinfo[['datamax']][variate$I] <- transf(x=rbind(apply(dt[,variate$I,with=F], 2, max, na.rm=T)), varinfo=varinfo, Iout='correct')


## binary
varinfo[['type']][variate$B] <- 'B'
varinfo[['transform']][variate$B] <- 'identity'
varinfo[['n']][variate$B] <- 2L
varinfo[['min']][variate$B] <- 0L
varinfo[['max']][variate$B] <- varinfo[['n']][variate$B]-1L
varinfo[['tmin']][variate$B] <- -Inf
varinfo[['tmax']][variate$B] <- +Inf
varinfo[['location']][variate$B] <- varinfo[['min']][variate$B]
varinfo[['scale']][variate$B] <- varinfo[['max']][variate$B]-varinfo[['min']][variate$B]
##
varinfo[['plotmin']][variate$B] <- transf(x=rbind(sapply(variate$B,function(v){
    dat <- dt[[v]]
    max(varinfo[['min']][v], min(dat, na.rm=T) - diff(tquant(dat, c(0.25,0.5))))
})), varinfo=varinfo, Bout='correct')
varinfo[['plotmax']][variate$B] <- transf(x=rbind(sapply(variate$B,function(v){
    dat <- dt[[v]]
    min(varinfo[['max']][v], max(dat, na.rm=T) + diff(tquant(dat, c(0.5,0.75))))
})), varinfo=varinfo, Bout='correct')
varinfo[['Q1']][variate$B] <- transf(x=rbind(apply(dt[,variate$B,with=F], 2, tquant, 0.25)), varinfo=varinfo, Bout='correct')
varinfo[['Q2']][variate$B] <- transf(x=rbind(apply(dt[,variate$B,with=F], 2, tquant, 0.5)), varinfo=varinfo, Bout='correct')
varinfo[['Q3']][variate$B] <- transf(x=rbind(apply(dt[,variate$B,with=F], 2, tquant, 0.75)), varinfo=varinfo, Bout='correct')
varinfo[['datamin']][variate$B] <- transf(x=rbind(apply(dt[,variate$B,with=F], 2, min, na.rm=T)), varinfo=varinfo, Bout='correct')
varinfo[['datamax']][variate$B] <- transf(x=rbind(apply(dt[,variate$B,with=F], 2, max, na.rm=T)), varinfo=varinfo, Bout='correct')

## logarithmic
varinfo[['type']][variate$L] <- 'R'
varinfo[['transform']][variate$L] <- 'log'
varinfo[['n']][variate$L] <- 0
varinfo[['min']][variate$L] <- 0
varinfo[['max']][variate$L] <- +Inf
varinfo[['tmin']][variate$L] <- -Inf
varinfo[['tmax']][variate$L] <- +Inf
varinfo[['location']][variate$L] <- apply(log(dt[,variate$L,with=F]),2,median,na.rm=T)
varinfo[['scale']][variate$L] <- apply(log(dt[,variate$L,with=F]),2,IQR,na.rm=T)*sd2iqr
##
varinfo[['plotmin']][variate$L] <- apply(dt[,variate$L,with=F],2,function(x){
    max(0, min(x, na.rm=T) - diff(tquant(x, c(0.25,0.5))))
    })
varinfo[['plotmax']][variate$L] <- apply(dt[,variate$L,with=F],2,function(x){
    max(x, na.rm=T) + diff(tquant(x, c(0.5,0.75)))
})
##
varinfo[['Q1']][variate$L] <- apply(dt[,variate$L,with=F], 2, tquant, 0.25)
varinfo[['Q2']][variate$L] <- apply(dt[,variate$L,with=F], 2, tquant, 0.5)
varinfo[['Q3']][variate$L] <- apply(dt[,variate$L,with=F], 2, tquant, 0.75)
varinfo[['datamin']][variate$L] <- apply(dt[,variate$L,with=F], 2, min, na.rm=T)
varinfo[['datamax']][variate$L] <- apply(dt[,variate$L,with=F], 2, max, na.rm=T)

## logarithmic censored
## varinfo[['type']][variate$S] <- 'O'
## varinfo[['transform']][variate$S] <- 'log'
## varinfo[['n']][variate$S] <- 0
## varinfo[['min']][variate$S] <- 0
## varinfo[['max']][variate$S] <- +Inf
## varinfo[['tmin']][variate$S] <- -Inf
## varinfo[['tmax']][variate$S] <- c(150L,300L)
## varinfo[['location']][variate$S] <- apply(log(dt[,variate$S,with=F]),2,median,na.rm=T)
## varinfo[['scale']][variate$S] <- apply(log(dt[,variate$S,with=F]),2,IQR,na.rm=T)*sd2iqr
## ##
## varinfo[['plotmin']][variate$S] <- apply(dt[,variate$S,with=F],2,function(x){
##     max(0, min(x, na.rm=T) - diff(tquant(x, c(0.25,0.5))))
##     })
## varinfo[['plotmax']][variate$S] <- sapply(variate$S,function(v){
##     dat <- dt[[v]]
##     min(varinfo[['tmax']][v], max(dat, na.rm=T) + diff(tquant(dat, c(0.5,0.75))))
## })

## varinfo[['t']][variate$S] <- lapply(variate$S,function(v){function(x){(log(x)-varinfo[['location']][v])/varinfo[['scale']][v]}})
## varinfo[['ij']][variate$S] <- lapply(variate$S,function(v){function(x){
##     out <- x*varinfo[['scale']][v]
##     out[x >= varinfo[['max']][v]] <- 1L
##     out
## }})



## varinfo[['hmean']][variate$R] <- 0L
## varinfo[['hmean']][variate$L] <- 0L
## varinfo[['hmean']][variate$S] <- 0L
## varinfo[['hmean']][variate$T] <- 0L
## varinfo[['hmean']][variate$I] <- 0L
## varinfo[['hmean']][variate$B] <- NA
## varinfo[['hmean']][variate$C] <- NA
## ##
## varinfo[['hsd']][variate$R] <- 2
## varinfo[['hsd']][variate$L] <- 2
## varinfo[['hsd']][variate$S] <- 2
## varinfo[['hsd']][variate$T] <- 1
## varinfo[['hsd']][variate$I] <- 7/8
## varinfo[['hsd']][variate$B] <- NA
## varinfo[['hsd']][variate$C] <- NA
## ##
## varinfo[['hshapeout']][variate$R] <- 1
## varinfo[['hshapeout']][variate$L] <- 1
## varinfo[['hshapeout']][variate$S] <- 1
## varinfo[['hshapeout']][variate$T] <- 1
## varinfo[['hshapeout']][variate$I] <- 1
## varinfo[['hshapeout']][variate$B] <- 1
## varinfo[['hshapeout']][variate$C] <- 1
## ##
## varinfo[['hshapein']][variate$R] <- 1
## varinfo[['hshapein']][variate$L] <- 1
## varinfo[['hshapein']][variate$S] <- 1
## varinfo[['hshapein']][variate$T] <- 1
## varinfo[['hshapein']][variate$I] <- 1
## varinfo[['hshapein']][variate$B] <- 1
## varinfo[['hshapein']][variate$C] <- 1
## ##
## varinfo[['hvarscale']][variate$R] <- 1
## varinfo[['hvarscale']][variate$L] <- 1
## varinfo[['hvarscale']][variate$S] <- 1
## varinfo[['hvarscale']][variate$T] <- 1/4
## varinfo[['hvarscale']][variate$I] <- 1/4
## varinfo[['hvarscale']][variate$B] <- NA
## varinfo[['hvarscale']][variate$C] <- NA


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

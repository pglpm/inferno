testd <- data.frame(testvrt=c(-7,NA,
    -6,-6,NA,(-5):(-2),NA, 1:4, NA, 5, 5,
    NA, 6)*2.54)
quantile(testd, 0.5, type=6, na.rm=T)
##   50%
## -1.27
summary(testd)
 ##    testvrt
 ## Min.   :-17.78
 ## 1st Qu.:-12.06
 ## Median : -1.27
 ## Mean   : -1.27
 ## 3rd Qu.:  9.53
 ## Max.   : 15.24
 ## NA's   :5

#### discretized
metadata <- data.frame(
    name='testvrt',
    type='continuous',
    Nvalues=Inf,
    gapwidth=2.54,
    domainmin= -6 * 2.54,
    domainmax= 5 * 2.54,
    minincluded=F,
    maxincluded=F,
    plotmin=NA,
    plotmax=NA
)
##
source('util_buildauxmetadata.R', local=T)
auxmetadata <- buildauxmetadata(testd, metadata)
auxmetadata

##
source('util_vtransform.R', local=T)
arglist <- c('init', 'left', 'right', 'aux', 'normalized')
as.data.frame(c(testd,sapply(arglist, function(x)
    vtransform(x=testd, auxmetadata=auxmetadata,
        Dout=x)
)))
##
testout <- vtransform(x=testd, auxmetadata=auxmetadata, Dout='normalized')
source('util_vtransform.R', local=T)
arglist <- c('mi', 'original')
as.data.frame(c(testd,testout,sapply(arglist, function(x)
    vtransform(x=testout, auxmetadata=auxmetadata,
        Dout=x)
)))



#### ordinal
metadata <- data.frame(
    name='testvrt',
    type='ordinal',
    Nvalues=NA,
    gapwidth=2.54,
    domainmin= -7 * 2.54,
    domainmax= 6 * 2.54,
    minincluded=T,
    maxincluded=T,
    plotmin=NA,
    plotmax=NA
)
##
source('util_buildauxmetadata.R', local=T)
auxmetadata <- buildauxmetadata(testd, metadata)
auxmetadata

##
source('util_vtransform.R', local=T)
arglist <- c('init', 'left', 'right', 'aux', 'normalized')
as.data.frame(c(testd,sapply(arglist, function(x)
    vtransform(x=testd, auxmetadata=auxmetadata,
        Dout=x)
)))
##
testout <- vtransform(x=testd, auxmetadata=auxmetadata, Dout='normalized')
source('util_vtransform.R', local=T)
arglist <- c('mi', 'original')
as.data.frame(c(testd,testout,sapply(arglist, function(x)
    vtransform(x=testout, auxmetadata=auxmetadata,
        Dout=x)
)))

#### closed domain
metadata <- data.frame(
    name='testvrt',
    type='continuous',
    Nvalues=Inf,
    gapwidth=0,
    domainmin= -7 * 2.54,
    domainmax= 6 * 2.54,
    minincluded=T,
    maxincluded=F,
    plotmin=NA,
    plotmax=NA
)
##
source('util_buildauxmetadata.R', local=T)
auxmetadata <- buildauxmetadata(testd, metadata)
auxmetadata

##
source('util_vtransform.R', local=T)
arglist <- c('init', 'left', 'right', 'lat', 'aux', 'boundisinf')
as.data.frame(c(testd,sapply(arglist, function(x)
    vtransform(x=testd, auxmetadata=auxmetadata,
        Cout=x)
)))
##
testout <- vtransform(x=testd, auxmetadata=auxmetadata,
    Cout='lat')
source('util_vtransform.R', local=T)
arglist <- c('mi', 'original')
as.data.frame(c(testd,testout,sapply(arglist, function(x)
    vtransform(x=testout, auxmetadata=auxmetadata,
        Cout=x)
)))


#### open domain
metadata <- data.frame(
    name='testvrt',
    type='continuous',
    Nvalues=Inf,
    gapwidth=0,
    domainmin= -7 * 2.54,
    domainmax= 6 * 2.54,
    minincluded=F,
    maxincluded=F,
    plotmin=NA,
    plotmax=NA
)
##
load('sysdata.rda')
source('util_buildauxmetadata.R', local=T)
auxmetadata <- buildauxmetadata(testd, metadata)
auxmetadata

##
source('util_vtransform.R', local=T)
arglist <- c('normalized')
as.data.frame(c(testd,sapply(arglist, function(x)
    vtransform(x=testd, auxmetadata=auxmetadata,
        Rout=x)
)))
##
testout <- vtransform(x=testd, auxmetadata=auxmetadata,
    Rout='normalized')
source('util_vtransform.R', local=T)
arglist <- c('mi', 'original')
as.data.frame(c(testd,testout,sapply(arglist, function(x)
    vtransform(x=testout, auxmetadata=auxmetadata,
        Rout=x)
)))



####
####
testd <- data.frame(testvrt=c(letters[c(1:2,4:5)], NA))

#### nominal
metadata <- data.frame(
    name='testvrt',
    type='nominal',
    Nvalues=Inf,
    rounding=2.54,
    domainmin=-10*2.54,
    domainmax=9*2.54,
    minincluded=F,
    maxincluded=T,
    plotmin=NA,
    plotmax=NA,
    V1='a',V2='b',V3='c',V4='d',V5='e'
)
##
source('util_buildauxmetadata.R', local=T)
auxmetadata <- buildauxmetadata(testd, metadata)
auxmetadata
##
source('util_vtransform.R', local=T)
arglist <- c('numeric')
as.data.frame(c(testd,sapply(arglist, function(x)
    vtransform(x=testd, auxmetadata=auxmetadata,
        Nout=x)
)))
##
testout <- vtransform(x=testd, auxmetadata=auxmetadata,
    Nout='numeric')
source('util_vtransform.R', local=T)
arglist <- c('mi', 'original')
as.data.frame(c(testd,testout,sapply(arglist, function(x)
    vtransform(x=testout, auxmetadata=auxmetadata,
        Nout=x)
)))


#### ordinal
metadata <- data.frame(
    name='testvrt',
    type='ordinal',
    Nvalues=Inf,
    rounding=2.54,
    domainmin=-10*2.54,
    domainmax=9*2.54,
    minincluded=F,
    maxincluded=T,
    plotmin=NA,
    plotmax=NA,
    V1='a',V2='b',V3='c',V4='d',V5='e'
)
##
source('util_buildauxmetadata.R', local=T)
auxmetadata <- buildauxmetadata(testd, metadata)
auxmetadata
##
source('util_vtransform.R', local=T)
arglist <- c('numeric')
as.data.frame(c(testd,sapply(arglist, function(x)
    vtransform(x=testd, auxmetadata=auxmetadata,
        Oout=x)
)))
##
testout <- vtransform(x=testd, auxmetadata=auxmetadata,
    Oout='numeric')
source('util_vtransform.R', local=T)
arglist <- c('mi', 'original')
as.data.frame(c(testd,testout,sapply(arglist, function(x)
    vtransform(x=testout, auxmetadata=auxmetadata,
        Oout=x)
)))



####
####
testd <- data.frame(testvrt=c(letters[1:2], NA))

#### nominal
metadata <- data.frame(
    name='testvrt',
    type='binary',
    Nvalues=Inf,
    rounding=2.54,
    domainmin=-10*2.54,
    domainmax=9*2.54,
    minincluded=F,
    maxincluded=T,
    plotmin=NA,
    plotmax=NA,
    V1='a',V2='b'
)
##
source('util_buildauxmetadata.R', local=T)
auxmetadata <- buildauxmetadata(testd, metadata)
auxmetadata
##
source('util_vtransform.R', local=T)
arglist <- c('numeric')
as.data.frame(c(testd,sapply(arglist, function(x)
    vtransform(x=testd, auxmetadata=auxmetadata,
        Bout=x)
)))
##
testout <- vtransform(x=testd, auxmetadata=auxmetadata,
    Bout='numeric')
source('util_vtransform.R', local=T)
arglist <- c('mi', 'original')
as.data.frame(c(testd,testout,sapply(arglist, function(x)
    vtransform(x=testout, auxmetadata=auxmetadata,
        Bout=x)
)))



#### Tests ordinal & rounded

testd <- data.frame(testvrt=c((-10):(-2), 1:9, NA)*2.54)
diff(unique(testd$testvrt))
quantile(testd$testvrt, 0.5, type=6, na.rm=T)
##  50%
## 1.27

source('buildmetadata.R', local=T)
metadata <- buildmetadata(testd)
##
metadata$domainmin <- min(testd$testvrt, na.rm=T)
metadata$domainmax <- max(testd$testvrt, na.rm=T)
metadata$minincluded <- T
metadata$maxincluded <- F
metadata









testl <- list(a=1,b=2,c=3)
testl2 <- list(a=100,b=2,c=3)
##
testfu <- function(x, testl){
    with(testl, {
        d <- a+x
        c(a,b,d,x)
    })
}
testfu(20, testl)
testfu(20, testl2)


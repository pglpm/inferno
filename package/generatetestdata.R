set.seed(789)
ndata <- 15
##
Rdata <- rnorm(n=ndata, mean=0, sd=1)
Rdata[2] <- NA
##
Cdata <- rnorm(n=ndata, mean=0, sd=3)
Cdata[Cdata <= -1] <- -1
Cdata[Cdata >= 1] <- 1
Cdata[3] <- NA
##
Ddata <- round(rnorm(n=ndata, mean=0, sd=1), 1)
Ddata[4] <- NA
##
Odata <- round(rnorm(n=ndata, mean=0, sd=3))
Odata[5] <- NA
##
Ndata <- round(rnorm(n=ndata, mean=0, sd=3))
Ndata[6] <- NA
nval <- letters[1:length(unique(Ndata[!is.na(Ndata)]))]
names(nval) <- unique(Ndata[!is.na(Ndata)])
Ndata <- nval[as.character(Ndata)]
##
Bdata <- rbinom(n=ndata, prob=0.5, size=1)
Bdata[7] <- NA
##
Rdatal <- exp(rnorm(n=ndata, mean=0, sd=1))
Rdatal[8] <- NA
##
##
testdatafile <- data.table(Rvar=Rdata,Cvar=Cdata,Dvar=Ddata,Ovar=Odata,Nvar=Ndata,Bvar=Bdata,Rvarl=Rdatal)
fwrite(testdatafile, 'testdata.csv')


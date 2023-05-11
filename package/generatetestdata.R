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


## more than one variate of each kind
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
Rdata2 <- rnorm(n=ndata, mean=0, sd=1)
Rdata2[3] <- NA
##
Cdata2 <- rnorm(n=ndata, mean=0, sd=1)
Cdata2[Cdata2 <= -1] <- -2
Cdata2[Cdata2 >= 1] <- 2
Cdata2[4] <- NA
##
Ddata2 <- round(rnorm(n=ndata, mean=0, sd=1), 1)
Ddata2[5] <- NA
##
Odata2 <- round(rnorm(n=ndata, mean=0, sd=2))
Odata2[6] <- NA
##
Ndata2 <- round(rnorm(n=ndata, mean=0, sd=1))
Ndata2[7] <- NA
nval <- letters[1:length(unique(Ndata2[!is.na(Ndata2)]))]
names(nval) <- unique(Ndata2[!is.na(Ndata2)])
Ndata2 <- nval[as.character(Ndata2)]
##
Bdata2 <- rbinom(n=ndata, prob=0.5, size=1)
Bdata2[8] <- NA
##
Rdatal2 <- exp(rnorm(n=ndata, mean=0, sd=1))
Rdatal2[9] <- NA
##
##
testdatafile2 <- data.table(Rvvv=Rdata,Cvvv=Cdata,Dvvv=Ddata,Ovvv=Odata,Nvvv=Ndata,Bvvv=Bdata,Rvvvl=Rdatal,Rvvv2=Rdata2,Cvvv2=Cdata2,Dvvv2=Ddata2,Ovvv2=Odata2,Nvvv2=Ndata2,Bvvv2=Bdata2,Rvvv2l=Rdatal2)
fwrite(testdatafile2, 'testdata2.csv')


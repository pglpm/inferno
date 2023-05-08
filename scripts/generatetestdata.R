set.seed(789)
ndata <- 10
##
Rdata <- rnorm(n=ndata, mean=0, sd=1)
Rdata2 <- rnorm(n=ndata, mean=1, sd=1)
##
Cdata <- rnorm(n=ndata, mean=0, sd=3)
Cdata[Cdata <= -1] <- -1
Cdata[Cdata >= 1] <- 1
##
Ddata <- round(rnorm(n=ndata, mean=0, sd=1),1)
##
Odata <- rpois(n=ndata, lambda=4)
##
Ndata <- letters[rpois(n=ndata, lambda=4)+1]
##
Bdata <- rbinom(n=ndata, prob=0.5, size=1)
##
##
Rdata[2] <- Cdata[3] <- Ddata[4] <- Odata[5] <- Ndata[6] <- Bdata[7] <- Rdata2[8] <- NA
##
testdatafile <- data.table(Rvar=Rdata,Cvar=Cdata,Dvar=Ddata,Ovar=Odata,Nvar=Ndata,Bvar=Bdata,Rvar2=Rdata2)
fwrite(testdatafile, 'testdata.csv')


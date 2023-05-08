set.seed(789)
ndata <- 10
##
Rdata <- rnorm(n=ndata, mean=0, sd=1)
##
Ldata <- rnorm(n=ndata, mean=0, sd=2)
Ldata[Ldata <= -1] <- -1
Ldata[Ldata >= 1] <- 1
##
Odata <- rpois(n=ndata, lambda=4)
##
Ndata <- letters[rpois(n=ndata, lambda=4)+1]
##
Bdata <- rbinom(n=ndata, prob=0.5, size=1)
##
##
Rdata[2] <- Ldata[3] <- Odata[4] <- Ndata[5] <- Bdata[6] <- NA
##
testdatafile <- data.table(Rvar=Rdata,Lvar=Ldata,Ovar=Odata,Nvar=Ndata,Bvar=Bdata)
fwrite(testdatafile, 'testdata.csv')


## Function to generate probability of number of components given alpha,
## rossi2014 eqn (2.5.2)

mmexp <- function(minalpha, maxalpha=-minalpha, byalpha=1, nclusters=64){
    sapply(2^seq(minalpha, maxalpha, by=byalpha),
        function(a){
            exp(as.numeric(log(abs(c(gmp::Stirling1(nclusters,0),
                gmp::Stirling1.all(nclusters))))) +
                (0:nclusters)*log(a) +
                lgamma(a)-lgamma(64+a))})
}

aseq <- function(minalpha, maxalpha=-minalpha, byalpha=1){
    (1:length(seq(minalpha,maxalpha,by=byalpha)))
}

vnorm <- function(x){x/sum(x)}

ma <- -4; Ma <- 4; ba <- 0.5
tplot(y=mmexp(ma,Ma,ba) %*% vnorm(aseq(ma,Ma,ba)^2))

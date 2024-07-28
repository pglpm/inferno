## Function for calculating total number of needed MCMC iterations
mcmclength <- function(nsamplesperchain, nitertot, thinning,
                        diagnESS, diagnIAT, diagnBMK, diagnMCSE,
                        diagnStat, diagnBurn, diagnBurn2, diagnThin) {
  ## This function uses (or can potentially use) various diagnostics from
  ## the LaplacesDemon package in order to determine
  ## the total number of needed MCMC iterations as well as
  ## the needed thinning
  ##
  ## Current method: require:
  ## 1. MCSE/sd > 6.27%
  ## 2. number of iterations larger than the required samples * ESS
  autothinning <- 2 * ceiling(max(diagnIAT, diagnThin))
  list(reqiter = nitertot +
         (if(
         max(diagnMCSE) < 6.2 &&
         nitertot > (autothinning * (nsamplesperchain - 1L))
       ) {
            0
          } else {
            3 * autothinning
          }),
       ##
       thinning = (if (is.null(thinning)) {
                     autothinning
                   } else {
                     thinning
                   })
       )
}

## oldmcmclength <- function(nsamplesperchain, nitertot, thinning,
##                         diagnESS, diagnIAT, diagnBMK, diagnMCSE,
##                         diagnStat, diagnBurn, diagnBurn2, diagnThin) {
##   ## This function uses (or can potentially use) various diagnostics from
##   ## the LaplacesDemon package in order to determine
##   ## the total number of needed MCMC iterations as well as
##   ## the needed thinning
##   ## Current method: no. of iterations is given by:
##   ## 2 * diagnosed burn-in +
##   ## 2 * (maximum between needed thinning and integrated autocorrelation time) *
##   ## number of independent samples needed (minus last one)
##   list(reqiter = ceiling(3 * max(diagnBurn2) +
##           ( 2 * ceiling(max(diagnIAT, diagnThin)) *
##            (nsamplesperchain - 1L) ) ),
##        thinning = (if (is.null(thinning)) {
##          2 * ceiling(max(diagnIAT, diagnThin))
##        } else {
##          thinning
##        })
##        )
## }

mixq <- function(x){sum(parmList$q[1,] * pnorm(x, mean=parmList$meanC[1,1,], sd=1/sqrt(parmList$tauC[1,1,])))}

fn <- function(x){(log(sum(parmList$q[1,] * pnorm(x, mean=parmList$meanC[1,1,], sd=1/sqrt(parmList$tauC[1,1,])))) - log(0.95))^2}
xq1 <- optim(0, fn, gr = function(x) pracma::grad(fn, x), 
method = "L-BFGS-B",
lower = -Inf, upper = Inf,
control = list(factr = 1e-10,
maxit = 100))


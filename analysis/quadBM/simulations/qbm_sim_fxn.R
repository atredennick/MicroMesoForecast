##  R Script for the Population Percent Cover Change Model



####
####  Define Population Growth Function
####
growFunc <- function(N, int, slope, clims, climcovs, tau){
  mu <- int+slope*log(N)+sum(clims*climcovs)
  newN <- rlnormTrunc(1, meanlog = mu, sdlog = tau, min = 0, max = 1)
  return(newN)
}

growFunc_noClim <- function(N, int, slope, tau){
  mu <- int+slope*log(N)
  newN <- rlnormTrunc(1, meanlog = mu, sdlog = tau, min = 0, max = 1)
  return(newN)
}
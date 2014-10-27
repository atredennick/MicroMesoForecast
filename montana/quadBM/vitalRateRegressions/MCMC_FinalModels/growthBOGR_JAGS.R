model{
  #growth process and likelihood
  for(i in 1:nObs){
    logit(gMu[i]) <- intYr[yrs[i]] + intG[grp[i]] + size[i]*beta[yrs[i]] 
    p[i] <- gMu[i] * tau
    q[i] <- (1 - gMu[i]) * tau
    Y[i] ~ dbeta(p[i], q[i])
  }
  
  #priors
  for(g in 1:nGrp){
    intG[g] ~ dnorm(0, intVarG)
  }
  
  for(y in 1:nYrs){
    beta[y] ~ dnorm(betaMu, betaVar)
    intYr[y] ~ dnorm(interceptMu, intVarY)
  }
  
  interceptMu ~ dnorm(0, .001)
  betaMu ~ dnorm(0, .001)
  betaVar ~ dgamma(0.001, 0.001)
  intVarY ~ dgamma(0.001, 0.001)
  intVarG ~ dgamma(2, 0.5)
  t0 ~ dnorm(0, .01) 
  tau <- exp(t0)
}
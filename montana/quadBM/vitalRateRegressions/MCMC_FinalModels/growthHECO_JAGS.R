model{
  #growth process and likelihood
  for(i in 1:nObs){
    logit(gMu[i]) <- intYr[yrs[i]] + intG[grp[i]] + size[i]*beta[yrs[i]] + temp1*TmeanSpr1[i] + temp2*TmeanSpr2[i] + rain1*ppt1[i] + rain2*ppt2[i] 
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
  
#   temp1 ~ dnorm(0,.001)
#   temp2 ~ dnorm(0,.001)
#   rain1 ~ dnorm(0,.001)
#   rain2 ~ dnorm(0,.001)
  temp1 ~ dunif(-10,10)
  temp2 ~ dunif(-10,10)
  rain1 ~ dunif(-10,10)
  rain2 ~ dunif(-10,10)
  interceptMu ~ dnorm(0, .001)
  betaMu ~ dnorm(0, .001)
  betaVar ~ dgamma(0.001, 0.001)
  intVarY ~ dgamma(0.001, 0.001)
  intVarG ~ dgamma(2, 0.5)
  t0 ~ dnorm(0, .01) 
  tau <- exp(t0)
}
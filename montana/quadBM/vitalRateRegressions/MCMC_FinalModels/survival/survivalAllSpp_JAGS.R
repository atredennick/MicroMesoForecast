model{
  #likelihood and process model
  for(i in 1:nObs){
    logit(p[i]) <- intYr[yrs[i]] + intG[grp[i]] + size[i]*beta[yrs[i]] + temp1*TmeanSpr1[i] + temp2*TmeanSpr2[i] + rain1*ppt1[i] + rain2*ppt2[i] 
    S[i] ~ dbern(p[i])
  }
  
  #priors
  for(g in 1:nGrp){
    intG[g] ~ dnorm(0, intVarG)
  }
  
  for(y in 1:nYrs){
    beta[y] ~ dnorm(betaMu, betaVar)
    intYr[y] ~ dnorm(interceptMu, intVarY)
  }
  
  temp1 ~ dnorm(0,1e-6)
  temp2 ~ dnorm(0,1e-6)
  rain1 ~ dnorm(0,1e-6)
  rain2 ~ dnorm(0,1e-6)
  interceptMu ~ dnorm(0, 1e-6)
  betaMu ~ dnorm(0, 1e-6)
  betaVar ~ dgamma(0.001, 0.001)
  intVarY ~ dgamma(0.001, 0.001)
  intVarG ~ dgamma(2, 0.5) 
}
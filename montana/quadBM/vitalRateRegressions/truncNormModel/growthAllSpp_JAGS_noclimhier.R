model{
  #likelihood and process model
  for(i in 1:nObs){
    gMu[i] <- intYr[spp[i],yrs[i]] + intG[spp[i],grp[i]] + beta[spp[i],yrs[i]]*X[i] + temp1[spp[i]]*TmeanSpr1[i] + temp2[spp[i]]*TmeanSpr2[i] + rain1[spp[i]]*ppt1[i] + rain2[spp[i]]*ppt2[i] + rainLag[spp[i]]*pptLag[i] +
      TmeanSpr1[i]*ppt1[i]*climInt1[spp[i]] + TmeanSpr2[i]*ppt2[i]*climInt2[spp[i]]
    C[i] ~ dlnorm(gMu[i], tau[spp[i]]) T(0,1)
  }
  
  #priors
  for(s in 1:nSpp){
    betaSpp[s] ~ dnorm(0, 1e-6)
    temp1[s] ~ dnorm(0, 1e-6)
    temp2[s] ~ dnorm(0, 1e-6)
    rain1[s] ~ dnorm(0, 1e-6)
    rain2[s] ~ dnorm(0, 1e-6)
    rainLag[s] ~ dnorm(0, 1e-6)
    climInt1[s] ~ dnorm(0, 1e-6)
    climInt2[s] ~ dnorm(0, 1e-6)
    intercept[s] ~ dnorm(0, 1e-6)
    intVaryY[s] ~ dgamma(0.001, 0.001)
    betaVar[s] ~ dgamma(0.001, 0.001)
    tau[s] ~ dgamma(0.001, 0.001)
#     t0[s] ~ dnorm(0, .01) 
#     tau[s] <- exp(t0[s])
    intVarG[s] ~ dgamma(2, 0.5) 
    for(y in 1:nYrs){
      intYr[s,y] ~ dnorm(intercept[s], intVaryY[s])
      beta[s,y] ~ dnorm(betaSpp[s], betaVar[s])
    }
    for(g in 1:nGrp){
      intG[s,g] ~ dnorm(0, intVarG[s])
    }
  }
  
#   temp1Mu ~ dnorm(0,1e-6)
#   temp2Mu ~ dnorm(0,1e-6)
#   rain1Mu ~ dnorm(0,1e-6)
#   rain2Mu ~ dnorm(0,1e-6)
#   rainLagMu ~ dnorm(0,1e-6)
#   climInt1Mu ~ dnorm(0,1e-6)
#   climInt2Mu ~ dnorm(0,1e-6)
#   temp1Var ~ dgamma(0.001, 0.001)
#   temp2Var ~ dgamma(0.001, 0.001)
#   rain1Var ~ dgamma(0.001, 0.001)
#   rain2Var ~ dgamma(0.001, 0.001)
#   rainLagVar ~ dgamma(0.001, 0.001)
#   climInt1Var ~ dgamma(0.001, 0.001)
#   climInt2Var ~ dgamma(0.001, 0.001)
}
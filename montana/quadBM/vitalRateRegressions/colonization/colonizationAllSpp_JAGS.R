model{
  #likelihood and process model
  for(i in 1:nObs){
    logit(p[i]) <- intercept[spp[i]] + intG[spp[i],grp[i]] + temp1[spp[i]]*TmeanSpr1[i] + temp2[spp[i]]*TmeanSpr2[i] + rain1[spp[i]]*ppt1[i] + rain2[spp[i]]*ppt2[i]  
    C[i] ~ dbern(p[i])
  }
  
  #priors
  for(s in 1:nSpp){
    temp1[s] ~ dnorm(temp1Mu, temp1Var)
    temp2[s] ~ dnorm(temp2Mu, temp2Var)
    rain1[s] ~ dnorm(rain1Mu, rain1Var)
    rain2[s] ~ dnorm(rain2Mu, rain2Var)
    intercept[s] ~ dnorm(0, 1e-6)
    intVarG[s] ~ dgamma(2, 0.5)
    for(g in 1:nGrp){
      intG[s,g] ~ dnorm(0, intVarG[s])
    }
#     intVaryY[s] ~ dgamma(0.001, 0.001)
    
#     for(y in 1:nYrs){
#       intYr[s,y] ~ dnorm(intercept[s], intVaryY[s])
#     }
  }
  
  temp1Mu ~ dnorm(0,1e-6)
  temp2Mu ~ dnorm(0,1e-6)
  rain1Mu ~ dnorm(0,1e-6)
  rain2Mu ~ dnorm(0,1e-6)
#   interceptMu ~ dnorm(0, 1e-6)
#   intVarY ~ dgamma(0.001, 0.001)
#   intVarSpp ~ dgamma(0.001, 0.001)
  temp1Var ~ dgamma(0.001, 0.001)
  temp2Var ~ dgamma(0.001, 0.001)
  rain1Var ~ dgamma(0.001, 0.001)
  rain2Var ~ dgamma(0.001, 0.001) 
}
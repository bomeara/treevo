GetBMRatePrior<-function(phy, traits, timeStep) {
  print("BM rate prior is an exponential distribution with a mean value approximately equal to the likelihood estimation")
  GetBrownianSDRate<-function(phy, traits, timeStep) { #conversion from continuous rate to discrete
    continuous.time.sigma.squared<-fitContinuous.hacked(phy, traits)$Trait1$beta[1]
    numSteps <- getSimulationSplits(phy)[1, 1] / timeStep
    discrete.time.sigma.squared<-continuous.time.sigma.squared * max(branching.times(phy)) / numSteps
    discrete.time.sd<-sqrt(discrete.time.sigma.squared)
    return(discrete.time.sd)
  }
  LikelihoodRateEst<-GetBrownianSDRate(phy, traits, timeStep)
  intrinsicPriorsValues<-LikelihoodRateEst^-1
  print(paste("LikelihoodRateEst =", LikelihoodRateEst, "; Mean of exp distribution for prior is ~", mean(rexp(10000, LikelihoodRateEst^-1)), "; exp rate =",  LikelihoodRateEst^-1))
  #set prior distribution mean to BM estimate
  return(c(intrinsicPriorsValues))
}

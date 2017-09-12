#' Get Brownian Motion Rate Prior
#' 
#' This function automatically calculates prior distributions for
#' the rate of trait evolution under the Brownian Motion (BM) model.
#' 
#' Returns a matrix of prior values that can be used in the \code{doRun} functions.
#' Builds on functions in \code{geiger} to estimate distribution.
#' 

#' @inheritParams doSimulation
#' @inheritParams doRun_prc

#' @param timeStep time in a single iteration of the discrete-time simulation

#' @return Returns a matrix of prior values

#' @author Brian O'Meara and Barb Banbury

# @references O'Meara and Banbury, unpublished
# @keywords GetBMRatePrior

#' @examples
#' 
#' data(simRun)
#' 
#' GetBMRatePrior(phy=simPhy, traits=simChar, timeStep=1)
#' 
#' 

#' @name GetBMRatePrior
#' @rdname GetBMRatePrior
#' @export
GetBMRatePrior<-function(phy, traits, timeStep) {
  print("BM rate prior is an exponential distribution with a mean value approximately equal to the likelihood estimation")
  GetBrownianSDRate<-function(phy, traits, timeStep) { #conversion from continuous rate to discrete
    if(is.null(names(traits))){
		names(traits) <- colnames(traits)
		}
	continuous.time.sigma.squared <- fitContinuous(phy, traits)$opt$sigsq
    numSteps <- getSimulationSplits(phy)[1, 1] / timeStep
    discrete.time.sigma.squared<-continuous.time.sigma.squared * max(branching.times(phy)) / numSteps
    discrete.time.sd<-sqrt(discrete.time.sigma.squared)
    return(discrete.time.sd)
  }
  LikelihoodRateEst<-GetBrownianSDRate(phy, traits, timeStep)
  intrinsicPriorsValues<-LikelihoodRateEst^-1
  print(paste("LikelihoodRateEst =", LikelihoodRateEst, "; Mean of exp distribution for prior is ~", 
	mean(rexp(10000, LikelihoodRateEst^-1)), "; exp rate =",  LikelihoodRateEst^-1))
  #set prior distribution mean to BM estimate
  return(c(intrinsicPriorsValues))
}

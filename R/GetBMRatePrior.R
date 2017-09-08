#' Get BM Rate Prior
#' 
#' This function automatically calculates prior distributions for BM model of
#' evolution
#' 
#' Returns a matrix of prior values that can be used in the doRun functions.
#' Builds on functions in geiger to estimate distribution.
#' 
#' @param phy Tree (Phylogenetic tree in phylo format)
#' @param traits data matrix with rownames equal to phy
#' @param timeStep time in a single iteration of the discrete-time simulation
#' @return Returns a matrix of prior values
#' @author Brian O'Meara and Barb Banbury

# @references O'Meara and Banbury, unpublished
# @keywords GetBMRatePrior

#' @examples
#' 
#' example(simRun)
#' GetBMRatePrior(phy, traits, timeStep)
#' 
#' 

#' @name GetBMRatePrior
#' @rdname GetBMRatePrior
#' @export
GetBMRatePrior<-function(phy, traits, timeStep) {
  print("BM rate prior is an exponential distribution with a mean value approximately equal to the likelihood estimation")
  GetBrownianSDRate<-function(phy, traits, timeStep) { #conversion from continuous rate to discrete
    if(is.null(names(traits)))
      names(traits) <- colnames(traits)
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

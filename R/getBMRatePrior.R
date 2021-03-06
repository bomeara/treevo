#' Get Brownian Motion Rate Prior
#' 
#' This function automatically calculates prior distributions for
#' the rate of trait evolution under the Brownian Motion (BM) model on a
#' discrete time-scale, at a given \code{timeStep}, in the sense that that
#' variable is used with other \code{TreEvo} functions like \code{doRun_prc}.
#' 
#' Returns a matrix of prior values that can be used in the \code{doRun} functions.
#' Builds on functions in \code{phylolm} to estimate distribution.
#' 

#' @inheritParams doSimulation
#' @inheritParams simulateWithPriors
#' @inheritParams doRun

#' @param timeStep time in a single iteration of the discrete-time simulation

#' @return Returns a matrix of prior values

#' @author Brian O'Meara and Barb Banbury

# @references O'Meara and Banbury, unpublished
# @keywords getBMRatePrior

#' @examples
#' 
#' data(simRunExample)
#' 
#' #timeStep = 0.1 -> effectively ~100 steps over the tree length
#' 
#' getBMRatePrior(phy = simPhyExample, traits = simCharExample, 
#'    timeStep = 0.01, verbose = TRUE)
#' 
#' 

#' @name getBMRatePrior
#' @rdname getBMRatePrior
#' @export
getBMRatePrior <- function(phy, traits, timeStep, verbose = TRUE){
  if(verbose){
    message("BM rate prior is an exponential distribution with a mean value approximately equal to the likelihood estimation")
    }
  LikelihoodRateEst <- GetBrownianSDRate(phy=phy, traits=traits, timeStep=timeStep)
  expectedRateEst <- LikelihoodRateEst^-1
  if(verbose){
    message(paste0("LikelihoodRateEst = ",
		LikelihoodRateEst,
		"; Mean of exp distribution for prior is ~", 
        mean(rexp(10000, LikelihoodRateEst^-1)),
		"; expected rate estimate= ",
		LikelihoodRateEst^-1))
		}
  #set prior distribution mean to BM estimate
  return(expectedRateEst)
  }


GetBrownianSDRate <- function(phy, traits, timeStep) { 
	#conversion from continuous rate to discrete
    if(is.null(names(traits))){
        names(traits) <- colnames(traits)
        }
    #
    continuous.time.sigma.squared <- numeric()
    for(i in 1:ncol(traits)){
        trait <- traits[, i]
        names(trait) <- rownames(traits)
        continuous.time.sigma.squared[i] <- getBM(phy = phy, dat = trait)$beta
        }
    #
    numSteps <- getSimulationSplits(phy)[1, 1] / timeStep
    discrete.time.sigma.squared <- continuous.time.sigma.squared * max(node.depth.edgelength(phy)) / numSteps
    discrete.time.sd <- sqrt(discrete.time.sigma.squared)
    return(discrete.time.sd)
	}

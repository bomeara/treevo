#extrinsic models
#note that these work for univariate, but need to be generalized for multivariate
#otherstates has one row per taxon, one column per state
#states is a vector for each taxon, with length = nchar

#' Extrinsic Character Evolution Models
#' 
#' Functions describing various models of 'extrinsic' evolution (i.e. evolutionary processes
#' dependent on factors extrinsic to the evolving lineage, such as environmental change, or
#' other evolving lineages that interact with the lineage in question (competitors, predators, etc).
#' 
#' The following extrinsic models are:
#' 
#' \code{nullExtrinsic} describes a model of no extrinsic character change.
#' It has no parameters, really.
#' 
#' \code{nearestNeighborDisplacementExtrinsic} describes a model of extrinsic trait evolution where character
#' values of a focal taxon depend on the values of closest relatives on the tree (e.g. competitive exclusion).
#' The input parameters for this model are:
#' \code{nearestNeighborDisplacementExtrinsic} with parameters \code{params = sd, springK, maximumForce}
#' 
#' \code{everyoneDisplacementExtrinsic} describes a model of extrinsic trait evolution where the character
#' values of a focal taxon depend on the values of all co-extant relatives on the simulated tree.
#' The input parameters for this model are:
#' \code{everyoneDisplacementExtrinsic} with parameters \code{params = sd, springK, maximumForce}
#' 
#' \code{ExponentiallyDecayingPushExtrinsic} describes a model of extrinsic trait evolution where the character
#' values of a focal taxon is 'pushed' away from other taxa with similar values, but the force of that 'push'
#' exponentially decays as lineages diverge and their character values become less similar.
#' The input parameters for this model are:
#' \code{ExponentiallyDecayingPushExtrinsic} with parameters \code{params = sd, maximumForce, halfDistance}
#' 

#' @inheritParams intrinsicModels

#' @param selfstates Vector of current trait values for the population of
#' interest. May be multiple for some models, but generally expected to be
#' only a single value. Multivariate \code{TreEvo} is not yet supported.

#' @param otherstates Matrix of current trait values for all concurrent taxa/populations other
#' than the one of interest, with one row for each taxon, and a column for each trait.
#' May be multiple states per taxa/populations for some models, but generally expected to be
#' only a single value. Multivariate \code{TreEvo} is not yet supported.

#' @aliases abcmodels.extrinsic

#' @seealso Intrinsic models are described at \code{\link{intrinsicModels}}.

#' @return
#' A vector of values representing character displacement of that lineage over a single time step.


#' @author Brian O'Meara and Barb Banbury

#' @examples
#
#' \donttest{
#
#' set.seed(1)
#' # Examples of simulations with various extrinsic models (and null intrinsic model)
#' tree <- rcoal(20)
#' # get realistic edge lengths
#' tree$edge.length <- tree$edge.length*20
#' 
#' #No trait evolution except due to
#'        # character displacement due to nearest neighbor taxon
#' char <- doSimulation(
#'     phy = tree, 
#'     intrinsicFn = nullIntrinsic, 
#'     extrinsicFn = nearestNeighborDisplacementExtrinsic, 
#'     startingValues = c(10), #root state
#'     intrinsicValues = c(0), 
#'     extrinsicValues = c(0.1, 0.1, 0.1), 
#'  generation.time = 100000)
#' 
#' #Similarly, no trait evolution except due to
#'        # character displacement from all other taxa in the clade
#' char <- doSimulation(
#'     phy = tree, 
#'     intrinsicFn = nullIntrinsic, 
#'     extrinsicFn = everyoneDisplacementExtrinsic, 
#'     startingValues = c(10), #root state
#'     intrinsicValues = c(0), 
#'     extrinsicValues = c(0.1, 0.1, 0.1), 
#'     generation.time = 100000)
#' 
#' # A variant where force of character displacement decays exponentially
#'         # as lineages become more different
#' char <- doSimulation(
#'     phy = tree, 
#'     intrinsicFn = nullIntrinsic, 
#'     extrinsicFn = ExponentiallyDecayingPushExtrinsic, 
#'     startingValues = c(10), #root state
#'     intrinsicValues = c(0), 
#'     extrinsicValues = c(0.1, 0.1, 2), 
#'     generation.time = 100000)
#
#' }
#

#' @name extrinsicModels
#' @rdname extrinsicModels
#' @export
nullExtrinsic <- function(params, selfstates, otherstates, timefrompresent) {
    newdisplacement <- 0*selfstates
    return(newdisplacement)
}

#' @rdname extrinsicModels
#' @export
nearestNeighborDisplacementExtrinsic <- function(params, selfstates, otherstates, timefrompresent) {
    #params[1] is sd, params[2] is springK, params[3] is maxforce
    repulsorTaxon <- which.min(abs(otherstates-selfstates))
    repulsorValue <- otherstates[repulsorTaxon]
    sd <- params[1]
    springK <- params[2]
    maxforce <- params[3]
    localsign <- sign(selfstates[1]- repulsorValue)
    #message(abs((selfstates[1]-repulsorValue)))
    if(localsign == 0) {
        localsign = sign(rpgm::rpgm.rnorm(n = 1))    
    }
    newdisplacement <- rpgm::rpgm.rnorm(n = 1, mean = localsign*min(c(abs(springK/((selfstates[1]-repulsorValue)*(selfstates[1]-repulsorValue))), maxforce), na.rm = TRUE), sd = sd)
    return(newdisplacement)
}


#' @rdname extrinsicModels
#' @export
everyoneDisplacementExtrinsic <- function(params, selfstates, otherstates, timefrompresent) {
    #this is set up for one character only right now
    #params[1] is sd, params[2] is springK, params[3] is maxforce
    sd <- params[1]
    springK  <- params[2]
    maxforce <- params[3]
    netforce <- 0
    for (i in 1:length(otherstates)) {
            localsign <- sign(selfstates[1]-otherstates[i])
            if(localsign == 0) {
                localsign = sign(rpgm::rpgm.rnorm(n = 1))    
            }
            netforce <- netforce+localsign*min(c(abs(springK/((selfstates[1]-otherstates[i])*(selfstates[1]-otherstates[i]))), maxforce), na.rm = TRUE)
    }
    newdisplacement <- rpgm::rpgm.rnorm(n = 1, mean = netforce, sd = sd)
    return(newdisplacement)
}


#' @rdname extrinsicModels
#' @export
ExponentiallyDecayingPushExtrinsic <- function(params, selfstates, otherstates, timefrompresent) {
    #params[1] is sd, params[2] is maxForce when character difference = 0, params[3] is half
        # distance (the phenotypic distance at which repulsion is half maxForce)
    repulsorTaxon <- which.min(abs(otherstates-selfstates))
    repulsorValue <- otherstates[repulsorTaxon]
    sd <- params[1]
    maxForce <- params[2]
    halfDistance <- params[3] #is like half life
    rate <- log(2, base = exp(1))/ halfDistance
    localsign <- sign(selfstates[1]- repulsorValue)
    if(localsign == 0) {  #to deal with case of identical values
        localsign = sign(rpgm::rpgm.rnorm(n = 1))    
    }
    newdisplacement <- rpgm::rpgm.rnorm(n = 1, mean = maxForce*localsign*exp(-1*rate*abs((selfstates[1]-repulsorValue))), sd = sd)
    return(newdisplacement)
}

#Extra functions for calculating Exponential Decay Push priors
qexpOptimality <- function(x, q, diff) {
  return(abs(qexp(q, x)-diff))
}

EstimateRate <- function(diff, quantile = 0.95) {
  result <- optimize(f = qexpOptimality, interval = c(0, 100000), q = quantile, diff = diff)
  return(result$minimum)
}

GetExpPushPriors <- function(numSteps, phy, data) {
  #returns a matrix with exponential rates for the three Exponential push priors
  timeStep <- 1/numSteps  #out of doRun_rej code
  sd <- getBMRatePrior(phy, data, timeStep)  #new TreEvo function
  data.sort <- sort(data[, 1])
  data.min.diff <- min(abs(data.sort[1:(length(data.sort)-1)]-data.sort[2:(length(data.sort))]))
  data.max.diff <- abs(max(data.sort)-min(data.sort))
  half.distance.prior <- EstimateRate(data.min.diff)
  max.push.prior <- EstimateMaxForce(data.max.diff, halfDistance = qexp(0.5, half.distance.prior), numSteps)
  return(matrix(c(sd, sd, half.distance.prior, half.distance.prior, max.push.prior, max.push.prior), nrow = 2))  #each one is doubled so they will exp rate param will fit into intrinsicPriorValues matrix in doRun_rej
}

OneStepPush <- function(maxForce, halfDistance, state) {
  rate <- log(2, base = exp(1))/ halfDistance
  push <- maxForce*2*exp(-1*rate*abs(state))
  return(push)
}

MultiStepPush <- function(maxForce, halfDistance, numSteps) {
  state <- 0
  for (i in sequence(numSteps)) {
    state <- OneStepPush(maxForce, halfDistance, state)+state
  }
  return(state)
}

ForceOptimality <- function(x, halfDistance, numSteps, diff) {
  return(abs(diff - MultiStepPush(x, halfDistance, numSteps)))
}

EstimateMaxForce <- function(diff, halfDistance, numSteps) {
  result <- optimize(f = ForceOptimality, interval = c(0, 100000), halfDistance = halfDistance, numSteps = numSteps, diff = diff)
  return(result$minimum)
}

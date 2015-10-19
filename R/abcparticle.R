#create abcparticle class object


#' abcparticle
#' 
#' An internal TreEvo function that creates a list of objects in class
#' abcparticle
#' 
#' 
#' @param id Simulation ID
#' @param generation Simulation generation
#' @param weight Simulation weight
#' @param distance Simulation distance
#' @param startingValues Vector of parameter estimates for startingStates
#' @param intrinsicValues Vector of parameter estimates for intrinsicValues
#' @param extrinsicValues Vector of parameter estimates for extrinsicValues
#' @author Brian O'Meara and Barb Banbury
#' @references O'Meara and Banbury, unpublished
abcparticle <- function( id=NA, generation=NA, weight=NA, distance=NA, startingValues=NA, intrinsicValues=NA, extrinsicValues=NA ) {
	particle <- list(id=id, generation=generation, weight=weight, distance=distance, 
		startingValues=startingValues, intrinsicValues=intrinsicValues, 
		extrinsicValues=extrinsicValues)
	class(particle) <- "abcparticle"
	return(particle)
}

initializeStatesFromMatrices <- function(particle, startingPriorsValues, startingPriorsFns, intrinsicPriorsValues, intrinsicPriorsFns, extrinsicPriorsValues, extrinsicPriorsFns) {
	particle$startingValues  <- rep(NA,length=dim(startingPriorsValues)[2])
	particle$intrinsicValues <- rep(NA,length=dim(intrinsicPriorsValues)[2])
	particle$extrinsicValues <- rep(NA,length=dim(extrinsicPriorsValues)[2])
	for (j in sequence(dim(startingPriorsValues)[2])) {
		particle$startingValues[j] <- pullFromPrior(priorValues=startingPriorsValues[, j], priorFn=startingPriorsFns[j])
	}
	for (j in sequence(dim(intrinsicPriorsValues)[2])) {
		particle$intrinsicValues[j] <- pullFromPrior(priorValues=intrinsicPriorsValues[, j], priorFn=intrinsicPriorsFns[j])
	}
	for (j in sequence(dim(extrinsicPriorsValues)[2])) {
		particle$extrinsicValues[j] <- pullFromPrior(priorValues=extrinsicPriorsValues[, j], priorFn=extrinsicPriorsFns[j])
	}
	return(particle)
}

mutateStates <- function(particle, startingPriorsValues, startingPriorsFns, intrinsicPriorsValues, intrinsicPriorsFns, extrinsicPriorsValues, extrinsicPriorsFns, standardDevFactor) {
	for (j in sequence(dim(startingPriorsValues)[2])) {
		particle$startingValues[j] <- mutateState(startingState=particle$startingValues[j], standardDevFactor=standardDevFactor, priorValues=startingPriorsValues[, j], priorFn=startingPriorsFns[j])
	}
	for (j in sequence(dim(intrinsicPriorsValues)[2])) {
		particle$intrinsicValues[j] <- mutateState(startingState=particle$intrinsicValues[j], standardDevFactor=standardDevFactor, priorValues=intrinsicPriorsValues[, j], priorFn=intrinsicPriorsFns[j])
	}
	for (j in sequence(dim(extrinsicPriorsValues)[2])) {
		particle$extrinsicValues[j] <- mutateState(startingState=particle$extrinsicValues[j], standardDevFactor=standardDevFactor, priorValues=extrinsicPriorsValues[, j], priorFn=extrinsicPriorsFns[j])
	}
	return(particle)
}

simulateTips <- function(particle, splits, phy, intrinsicFn, extrinsicFn, timeStep) {
	newtips<-convertTaxonFrameToGeigerData(doSimulation(splits, intrinsicFn, extrinsicFn, particle$startingValues, particle$intrinsicValues, particle$extrinsicValues, timeStep), phy)
	return(newtips)
}

#create abcparticleS3 class object
abcparticleS3 <- function( id=NA, generation=NA, weight=NA, distance=NA, startingStates=NA, intrinsicValues=NA, extrinsicValues=NA ) {
	particle <- list(id=id, generation=generation, weight=weight, distance=distance, 
		startingStates=startingStates, intrinsticValues=intrinsicValues, 
		extrinsicValues=extrinsicValues)
	class(particle) <- "abcparticleS3"
	return(particle)
}

initializeStatesFromMatrices <- function(particle, startingPriorsValues, startingPriorsFns, intrinsicPriorsValues, intrinsicPriorsFns, extrinsicPriorsValues, extrinsicPriorsFns) {
	particle$startingStates  <- rep(NA,length=dim(startingPriorsValues)[2])
	particle$intrinsicStates <- rep(NA,length=dim(intrinsicPriorsValues)[2])
	particle$extrinsicStates <- rep(NA,length=dim(extrinsicPriorsValues)[2])
	for (j in sequence(dim(startingPriorsValues)[2])) {
		particle$startingStates[j] <- pullFromPrior(priorValues=startingPriorsValues[, j], priorFn=startingPriorsFns[j])
	}
	for (j in sequence(dim(intrinsicPriorsValues)[2])) {
		particle$intrinsicStates[j] <- pullFromPrior(priorValues=intrinsicPriorsValues[, j], priorFn=intrinsicPriorsFns[j])
	}
	for (j in sequence(dim(extrinsicPriorsValues)[2])) {
		particle$extrinsicStates[j] <- pullFromPrior(priorValues=extrinsicPriorsValues[, j], priorFn=extrinsicPriorsFns[j])
	}
	return(particle)
}

mutateStates <- function(particle, startingPriorsValues, startingPriorsFns, intrinsicPriorsValues, intrinsicPriorsFns, extrinsicPriorsValues, extrinsicPriorsFns, standardDevFactor) {
	for (j in sequence(dim(startingPriorsValues)[2])) {
		particle$startingStates[j] <- mutateState(startingState=particle$startingStates[j], standardDevFactor=standardDevFactor, priorValues=startingPriorsValues[, j], priorFn=startingPriorsFns[j])
	}
	for (j in sequence(dim(intrinsicPriorsValues)[2])) {
		particle$intrinsicStates[j] <- mutateState(startingState=particle$intrinsicStates[j], standardDevFactor=standardDevFactor, priorValues=intrinsicPriorsValues[, j], priorFn=intrinsicPriorsFns[j])
	}
	for (j in sequence(dim(extrinsicPriorsValues)[2])) {
		particle$extrinsicStates[j] <- mutateState(startingState=particle$extrinsicStates[j], standardDevFactor=standardDevFactor, priorValues=extrinsicPriorsValues[, j], priorFn=extrinsicPriorsFns[j])
	}
	return(particle)
}

simulateTips <- function(particle, splits, phy, intrinsicFn, extrinsicFn, timeStep) {
	newtips<-convertTaxonFrameToGeigerData(doSimulation(splits, intrinsicFn, extrinsicFn, particle$startingStates, particle$intrinsicValues, particle$extrinsicValues, timeStep), phy)
	return(newtips)
}
#create abcparticle class object


#  abcparticle
# 
#  An internal TreEvo function that creates a list of objects in class
#  abcparticle
# 
# 

#  @param id Simulation ID

#  @param generation Simulation generation

#  @param weight Simulation weight

#  @param distance Simulation distance

#  @param startingValues Vector of parameter estimates for startingStates

#  @param intrinsicValues Vector of parameter estimates for intrinsicValues

#  @param extrinsicValues Vector of parameter estimates for extrinsicValues

#  @author Brian O'Meara and Barb Banbury

# @references O'Meara and Banbury, unpublished

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

# simulateTips <- function(particle, splits, phy, intrinsicFn, extrinsicFn, timeStep) {
# 	newtips<-convertTaxonFrameToGeigerData(doSimulation(splits, intrinsicFn, extrinsicFn, particle$startingValues, particle$intrinsicValues, particle$extrinsicValues, timeStep), phy)
# 	return(newtips)
# }

simulateTips <- function(particle, taxon.df, phy, intrinsicFn, extrinsicFn, timeStep) {
	newtips<-doSimulationWithPossibleExtinction(taxon.df, intrinsicFn, extrinsicFn, particle$startingValues, particle$intrinsicValues, particle$extrinsicValues, timeStep)
	return(newtips)
}

#  Mutate Character State
#  
#  This function mutates the character state of a given taxon by one discrete
#  time step.
#  
#  

#  @param startingState Character state prior to mutating.

#  @param standardDevFactor Standard deviation.

#  @param priorFn Shape of the prior distribution. Must be one of
#  "fixed", "uniform", "normal", "lognormal", "gamma", or "exponential".

#  @param priorValues Vector of parameter Values for the prior function.


#  @author Brian O'Meara and Barb Banbury

# @references O'Meara and Banbury, unpublished


#  @examples
#  
#  data(simRun)
#  
#  mutateState(startingState, standardDevFactor, priorValues, priorFn)
#  

#  @name mutateState
#  @rdname mutateState
#  @export
mutateState<-function(startingState, standardDevFactor, priorFn, priorValues) {
	newState<-NA
	minBound=-Inf
	maxBound=Inf
	validNewState<-FALSE  #was lowercase, but not recognised 
	priorFn<-match.arg(arg=priorFn,
		choices=c("fixed", "uniform", "normal", "lognormal", "gamma", "exponential"),several.ok=FALSE);
	if (priorFn=="fixed" || priorFn=="uniform") {
		minBound<-min(priorValues)
		maxBound<-max(priorValues)
	}
	else if (priorFn=="lognormal" || priorFn=="gamma" || priorFn=="exponential") {
		minBound<-0
	}
	
	sdToUse<-standardDevFactor
	if (priorFn=="fixed") {
		sdToUse<-0
	}
	else if (priorFn=="uniform") {
		sdToUse<-standardDevFactor*(abs(max(priorValues)-min(priorValues)))
			if (sdToUse<0){
				print(paste("priorFN=", priorFn, "standardDevFactor=", standardDevFactor, "range(priorValues)=", range(priorValues)))
			}
	}
	else if (priorFn=="normal") {
		sdToUse<-standardDevFactor*priorValues[2]
			if (sdToUse<0){
				print(paste("priorFN=", priorFn, "standardDevFactor=", standardDevFactor, "range(priorValues)=", range(priorValues)))
			}
	}
	else if (priorFn=="lognormal") {
		sdToUse<-standardDevFactor*priorValues[2]
			if (sdToUse<0){
				print(paste("priorFN=", priorFn, "standardDevFactor=", standardDevFactor, "range(priorValues)=", range(priorValues)))
			}
	}
	else if (priorFn=="gamma") {
		sdToUse<-standardDevFactor*sqrt(priorValues[1]*priorValues[2]*priorValues[2])
			if (sdToUse<0){
				print(paste("priorFN=", priorFn, "standardDevFactor=", standardDevFactor, "range(priorValues)=", range(priorValues)))
			}
	}
	else if (priorFn=="exponential") {
		sdToUse<-standardDevFactor/priorValues[1]
			if (sdToUse<0){
				print(paste("priorFN=", priorFn, "standardDevFactor=", standardDevFactor, "range(priorValues)=", range(priorValues)))
			}
	}
	else {
		stop(priorFn," was not a recognized prior function")
	}

	while(!validNewState) {
		newState<-rnorm(n=1, mean=startingState, sd=sdToUse)
		validNewState<-TRUE
		if(is.na(newState)) {
			print(paste("MUTATESTATE_ERROR: newState = ",newState," sdToUse=",sdToUse," startingState=",startingState," priorFn=",priorFn," startingState=",startingState," priorValues=\n",sep=""))
			print(priorValues)
		}
		if (newState<minBound){
			validNewState<-FALSE
			}
		if (newState>maxBound){
			validNewState<-FALSE
			}	
#		if (!validNewState)	{
			#cat("newState ",newState," does not fit into one of the bounds (", minBound, "--", maxBound, ")\n")
#		}	
	}
	
newState
}

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

abcparticle <- function( id = NA, generation = NA, weight = NA, distance = NA, parentid = NA, startingValues = NA, intrinsicValues = NA, extrinsicValues = NA ) {
    particle <- list(
        id = id, 
        generation = generation, 
        weight = weight, 
        distance = distance, 
        parentid = parentid, 
        startingValues = startingValues, 
        intrinsicValues = intrinsicValues, 
        extrinsicValues = extrinsicValues)
    #
    #
    # note that the Priors Values are matrices, where each column is a different model parameter, 
        # and the rows represent different parameters of the prior for that parameter
    #
    class(particle) <- "abcparticle"
    return(particle)
    }

initializeStatesFromMatrices <- function(particle, 
    startingPriorsValues, startingPriorsFns, 
    intrinsicPriorsValues, intrinsicPriorsFns, 
    extrinsicPriorsValues, extrinsicPriorsFns) {
    #
    # note that the Priors Values are matrices, where each column is a different model parameter, 
        # and the rows represent different parameters of the prior for that parameter
    #
    particle$startingValues  <- rep(NA, length = length(startingPriorsValues))
    particle$intrinsicValues <- rep(NA, length = length(intrinsicPriorsValues))
    particle$extrinsicValues <- rep(NA, length = length(extrinsicPriorsValues))
    for (j in 1:length(startingPriorsValues)) {
        particle$startingValues[j] <- pullFromPrior(priorValues = startingPriorsValues[[j]], priorFn = startingPriorsFns[j])
        }
    for (j in 1:length(intrinsicPriorsValues)) {
        particle$intrinsicValues[j] <- pullFromPrior(priorValues = intrinsicPriorsValues[[j]], priorFn = intrinsicPriorsFns[j])
        }
    for (j in 1:length(extrinsicPriorsValues)) {
        particle$extrinsicValues[j] <- pullFromPrior(priorValues = extrinsicPriorsValues[[j]], priorFn = extrinsicPriorsFns[j])
        }
    return(particle)
    }
 
mutateStates <- function(particle, 
    startingPriorsValues, startingPriorsFns, 
    intrinsicPriorsValues, intrinsicPriorsFns, 
    extrinsicPriorsValues, extrinsicPriorsFns, 
    standardDevFactor) {
    #
    # note that the Priors Values are matrices, where each column is a different model parameter, 
        # and the rows represent different parameters of the prior for that parameter
    #
    for (j in 1:length(startingPriorsValues)) {
        particle$startingValues[j] <- mutateState(startingState = particle$startingValues[j], standardDevFactor = standardDevFactor, 
            priorValues = startingPriorsValues[[j]], priorFn = startingPriorsFns[j])
        }
    for (j in 1:length(intrinsicPriorsValues)) {
        particle$intrinsicValues[j] <- mutateState(startingState = particle$intrinsicValues[j], standardDevFactor = standardDevFactor, 
            priorValues = intrinsicPriorsValues[[j]], priorFn = intrinsicPriorsFns[j])
        }
    for (j in 1:length(extrinsicPriorsValues)) {
        particle$extrinsicValues[j] <- mutateState(startingState = particle$extrinsicValues[j], standardDevFactor = standardDevFactor, 
            priorValues = extrinsicPriorsValues[[j]], priorFn = extrinsicPriorsFns[j])
        }
    return(particle)
    }
 
# simulateTips <- function(particle, splits, phy, intrinsicFn, extrinsicFn, timeStep) {
#     newtips <- convertTaxonFrameToGeigerData(doSimulation(splits, intrinsicFn, extrinsicFn, particle$startingValues, particle$intrinsicValues, particle$extrinsicValues, timeStep), phy)
#     return(newtips)
# } 

#simulateTips <- function(particle, taxonDF, phy, intrinsicFn, extrinsicFn, timeStep) {
#    newtips <- doSimulation(taxonDF = taxonDF, 
#        intrinsicFn = intrinsicFn, extrinsicFn = extrinsicFn, startingValues = particle$startingValues, 
#        intrinsicValues = particle$intrinsicValues, extrinsicValues = particle$extrinsicValues, timeStep = timeStep)
#    return(newtips)
#    }

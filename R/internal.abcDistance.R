#  Calculate ABC Distances
#
#  This function uses results from a PLS model to calculate the distance from
#  each sim to the original data.
#
#  This function runs a PLS regression on each free parameter in the model, unlike
#  a true multivariate PLS regression.  For ABC, this seems to result in much
#  better results, without one parameter dominating the combined variance.
#

#  @param summaryValuesMatrix Matrix of summary statistics from simulations

#  @param originalSummaryValues Summary statistics for original data.

#  @param pls.model.list indexing of a list of pls models for different params

#  @return Returns euclidean distance of each simulation's summary values to
#  the original summary stats.

#  @author Brian O'Meara and Barb Banbury

# @references O'Meara and Banbury, unpublished

# @keywords abcDistance

#  @examples
#
#
#  \donttest{
#
#  data(simRunExample)
#
#  # example simulation
#
#  simDataParallel <- parallelSimulateWithPriors(
#    nrepSim = 10, multicore = FALSE, coreLimit = 1,
#    phy = simPhy,
#    intrinsicFn = brownianIntrinsic,
#    extrinsicFn = nullExtrinsic,
#    startingPriorsFns = "normal",
#    startingPriorsValues = matrix(c(mean(simChar[,1]), sd(simChar[,1]))),
#    intrinsicPriorsFns = c("exponential"),
#    intrinsicPriorsValues = matrix(c(10, 10), nrow = 2, byrow = FALSE),
#    extrinsicPriorsFns = c("fixed"),
#    extrinsicPriorsValues = matrix(c(0, 0), nrow = 2, byrow = FALSE),
#    timeStep = 10^-6,
#    checkpointFile = NULL, checkpointFreq = 24,
#    verbose = FALSE,
#    freevector = NULL, taxonDF = NULL,
#    niter.brown = 25, niter.lambda = 25, niter.delta = 25, niter.OU = 25, niter.white = 25)
#
#  nParFree <- sum(attr(simDataParallel,"freevector"))
#
#  #separate the simulation results: 'true' generating parameter values from the summary values
#  trueFreeValuesMat <- simDataParallel[,1:nParFree]
#  summaryValuesMat <- simDataParallel[,-1:-nParFree]
#
#  PLSmodel <- returnPLSModel(trueFreeValuesMatrix = trueFreeValuesMat,
#  	  	summaryValuesMatrix = summaryValuesMat,
#     		validation = "CV", scale = TRUE, variance.cutoff = 95)
#
#  abcDistance(summaryValuesMatrix = summaryValuesMat, originalSummaryValues, pls.model.list = PLSmodel)
#
#  }
#

# # respective lines from doRun_prc:
#	  pls.model.list <- apply(trueFreeValuesMatrix, 2, returnPLSModel, summaryValuesMatrix = summaryValuesMatrix, validation = validation,
#		scale = scale, variance.cutoff = variance.cutoff)
#	  originalSummaryValues <- summaryStatsLong(phy, traits, niter.brown = 200, niter.lambda = 200, niter.delta = 200, niter.OU = 200, niter.white = 200)
#	  distanceVector <- abcDistance(summaryValuesMatrix, originalSummaryValues, pls.model.list)

#  @name abcDistance
#  @rdname abcDistance
#  @export
abcDistance <- function(summaryValuesMatrix, originalSummaryValues, pls.model.list) {
  abcDistancesRaw <- sapply(sequence(length(pls.model.list)),
	SingleParameterPLSDistanceSquaredFixedPLS,
		pls.model.list = pls.model.list, summaryValuesMatrix = summaryValuesMatrix,
		originalSummaryValues = originalSummaryValues, scale = scale
		)
  if (class(abcDistancesRaw) != "matrix") { #it must be a vector, but apply likes matrices
  	abcDistancesRaw <- matrix(abcDistancesRaw, nrow = 1)
  }
  abcDistancesRawTotal <- apply(abcDistancesRaw, 1, sum)
  abcDistances <- sqrt(abcDistancesRawTotal) #Euclid rules.
  return(abcDistances)
}

SingleParameterPLSDistanceSquaredFixedPLS <- function(index, pls.model.list, summaryValuesMatrix, originalSummaryValues, scale = TRUE) {
  pls.model <- pls.model.list[[index]]
  summaryValues.transformed <- PLSTransform(summaryValuesMatrix, pls.model)
  originalSummaryValues.transformed <- PLSTransform(originalSummaryValues, pls.model)
  distanceByRow <- function(x,originalSummaryValues.transformed) {
    return(dist(matrix(c(x,originalSummaryValues.transformed),byrow = TRUE,nrow = 2))[1])
  }
  if (class(summaryValues.transformed) != "matrix") {
  	summaryValues.transformed <- matrix(summaryValues.transformed,nrow = 1)
  }
  raw.distances <- apply(summaryValues.transformed, 1, distanceByRow, originalSummaryValues.transformed = originalSummaryValues.transformed)
  return(raw.distances^2)
}

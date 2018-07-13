#  Partial Least Squares Rejection
#
#  This function performs the ABC-rejection analysis using an input simulation
#  data. Particles are accepted is they fall sufficiently close to the target
#  data (within the tolerance). Distances are calculated using \code{abcDistance}.
#

#  @param summaryValuesMatrix Matrix of summary statistics from simulations.

#  @param trueFreeValuesMatrix Matrix of true free values from simulations.

#  @inheritParams doSimulation
#  @inheritParams doRun

#  @param abcTolerance Proportion of accepted simulations.

#  @param verbose option to message progress to screen.

#  @param validation Cross Validation procedure for ABC.

#  @param scale Scale for \code{pls.model.list}.

#  @param variance.cutoff Variance cutoff for \code{pls.model.list}.

#  @return Returns a list of the particle data frame and ABC distances.

#  @author Brian O'Meara and Barb Banbury

# @references O'Meara and Banbury, unpublished
# @keywords PLSRejection doRun doRun_rej abc

# @examples
#
# PLSRejection(summaryValuesMatrix, 
#   trueFreeValuesMatrix, 
# 	 phy = simPhy, traits = simChar, abcTolerance = results[[1]]$abcTolerance)
#

# This is a difficult function to write an example for, as mainly an algorithm for doRun_rej; unclear if even needs to be exported
# Will internalize and see if it makes any difference


#  @name PLSRejection
#  @rdname PLSRejection
#  @export
PLSRejection <- function(summaryValuesMatrix, trueFreeValuesMatrix, 
      phy, traits, abcTolerance,
      verbose = TRUE, validation = "CV", scale = TRUE, variance.cutoff = 95) {
  ##############################################################################
  originalSummaryValues <- summaryStatsLong(phy = phy, traits = traits)
  if (verbose) {
    message("Done getting originalSummaryValues")
  }
  abcDistancesRaw <- sapply(sequence(dim(trueFreeValuesMatrix)[2]), SingleParameterPLSDistanceSquared,
         summaryValuesMatrix = summaryValuesMatrix, 
         trueFreeValuesMatrix = trueFreeValuesMatrix, originalSummaryValues = originalSummaryValues, 
		 validation = validation, scale = scale, variance.cutoff = variance.cutoff)
  abcDistancesRawTotal <- apply(abcDistancesRaw, 1, sum)
  abcDistances <- sqrt(abcDistancesRawTotal) #Euclid rules.
  #
  # here's where we do abc - calculate quantile on distances
      # drop what distances are greater than that cutoff
  abcQuantile <- quantile(abcDistances, prob = abcTolerance)
  acceptedParticles <- trueFreeValuesMatrix[
    which(abcDistances <= abcQuantile)
	, , drop = FALSE] 
  #
  acceptedDistances <- abcDistances[which(abcDistances <= quantile(abcDistances, prob = abcTolerance))]
  #
  particleDataFrame <- data.frame(cbind(
	rep(1, dim(acceptedParticles)[1]), 
	as.vector(which(abcDistances <= quantile(abcDistances, prob = abcTolerance))), 
	seq(1:dim(acceptedParticles)[1]), 
	rep(0, dim(acceptedParticles)[1]), 
	acceptedDistances, rep(1, 
	dim(acceptedParticles)[1]), 
	acceptedParticles
	))
  colnames(particleDataFrame) <- c("generation", "attempt", "id", "parentid", "distance", 
	"weight", paste("param", seq(dim(trueFreeValuesMatrix)[2]), sep = ""))
  return(list(particleDataFrame = particleDataFrame, abcDistances = abcDistances))
  }

SingleParameterPLSDistanceSquared <- function(index, summaryValuesMatrix, trueFreeValuesMatrix, 
		originalSummaryValues, validation = "CV", scale = TRUE, variance.cutoff = 95) {
  trueFreeValuesMatrix <- trueFreeValuesMatrix[, index]
  pls.model <- returnPLSModel(trueFreeValuesMatrix, summaryValuesMatrix,
    validation = validation, scale = scale, variance.cutoff = variance.cutoff)
  summaryValues.transformed <- PLSTransform(summaryValuesMatrix, pls.model)
  originalSummaryValues.transformed <- PLSTransform(originalSummaryValues, pls.model)
  #
  distanceByRow <- function(x, originalSummaryValues.transformed) {
    return(dist(matrix(c(x, originalSummaryValues.transformed), byrow = TRUE, nrow = 2))[1])
	}
  #
  raw.distances <- apply(summaryValues.transformed, 1, distanceByRow,
    originalSummaryValues.transformed = originalSummaryValues.transformed
	)
  return(raw.distances^2)
  }

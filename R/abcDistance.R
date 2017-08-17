#' Calculate ABC Distances
#' 
#' This function uses results from a PLS model to calculate the distance from
#' each sim to the original data.
#' 
#' This function runs a PLS regression on each free parameter in the model, unlike
#' a true multivariate PLS regression.  For ABC, this seems to result in much
#' better results, without one parameter dominating the combined variance.
#' 

#' @param summaryValuesMatrix Matrix of summary statistics from simulations

#' @param originalSummaryValues Summary statistics for original data.

#' @param pls.model.list indexing of a list of pls models for different params

#' @return Returns euclidean distance of each simulation's summary values to
#' the original summary stats.

#' @author Brian O'Meara and Barb Banbury

# @references O'Meara and Banbury, unpublished

# @keywords abcDistance

#' @examples
#' 
#' #abcDistance(summaryValuesMatrix, originalSummaryValues, pls.model.list)
#' 
abcDistance<-function(summaryValuesMatrix, originalSummaryValues, pls.model.list) {
  abcDistancesRaw<-sapply(sequence(length(pls.model.list)), SingleParameterPLSDistanceSquaredFixedPLS, pls.model.list=pls.model.list, summaryValuesMatrix=summaryValuesMatrix, originalSummaryValues=originalSummaryValues, scale=scale)
  if (class(abcDistancesRaw)!="matrix") { #it must be a vector, but apply likes matrices
  	abcDistancesRaw<-matrix(abcDistancesRaw, nrow=1)
  }
  abcDistancesRawTotal<-apply(abcDistancesRaw, 1, sum)
  abcDistances<-sqrt(abcDistancesRawTotal) #Euclid rules.
  return(abcDistances)
}

SingleParameterPLSDistanceSquaredFixedPLS<-function(index, pls.model.list, summaryValuesMatrix, originalSummaryValues, scale=TRUE) {
  pls.model <- pls.model.list[[index]]
  summaryValues.transformed<-PLSTransform(summaryValuesMatrix, pls.model)
  originalSummaryValues.transformed<-PLSTransform(originalSummaryValues, pls.model)
  distanceByRow<-function(x,originalSummaryValues.transformed) {
    return(dist(matrix(c(x,originalSummaryValues.transformed),byrow=TRUE,nrow=2))[1])
  }
  if (class(summaryValues.transformed)!="matrix") {
  	summaryValues.transformed<-matrix(summaryValues.transformed,nrow=1)
  }
  raw.distances<-apply(summaryValues.transformed, 1, distanceByRow, originalSummaryValues.transformed=originalSummaryValues.transformed)
  return(raw.distances^2)
}

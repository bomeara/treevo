abcDistance<-function(summaryValuesMatrix, originalSummaryValues, pls.model.list) {
  abcDistancesRaw<-sapply(sequence(length(pls.model.list)), SingleParameterPLSDistanceSquaredFixedPLS, summaryValuesMatrix=summaryValuesMatrix, 
                          originalSummaryValues=originalSummaryValues, scale=scale)
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
  raw.distances<-apply(summaryValues.transformed, 1, distanceByRow, originalSummaryValues.transformed=originalSummaryValues.transformed)
  return(raw.distances^2)
}
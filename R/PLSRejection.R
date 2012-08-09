PLSRejection<-function(summaryValuesMatrix, trueFreeValuesMatrix, phy, traits, abcTolerance, verbose=TRUE) {
  
  pls.model<-returnPLSModel(trueFreeValuesMatrix[,1],summaryValuesMatrix)
  originalSummaryValues<-summaryStatsLong(phy=phy, data=traits)
  abcDistancesRaw<-sapply(sequence(dim(trueFreeValuesMatrix)[2]), SingleParameterPLSDistanceSquared, summaryValuesMatrix=summaryValuesMatrix, 
         trueFreeValuesMatrix=trueFreeValuesMatrix, originalSummaryValues=originalSummaryValues)
  abcDistancesRawTotal<-apply(abcDistancesRaw, 1, sum)
  abcDistances<-sqrt(abcDistancesRawTotal) #Euclid rules.
}

SingleParameterPLSDistanceSquared<-function(index, summaryValuesMatrix, trueFreeValuesMatrix, originalSummaryValues) {
  trueFreeValuesMatrix<-trueFreeValuesMatrix[,index]
  pls.model<-returnPLSModel(trueFreeValuesMatrix,summaryValuesMatrix)
  summaryValues.transformed<-PLSTransform(summaryValuesMatrix, pls.model)
  originalSummaryValues.transformed<-PLSTransform(originalSummaryValues, pls.model)
  distanceByRow<-function(x,originalSummaryValues.transformed) {
    return(dist(matrix(c(x,originalSummaryValues.transformed),byrow=TRUE,nrow=2))[1])
  }
  raw.distances<-apply(summaryValues.transformed, 1, distanceByRow, originalSummaryValues.transformed=originalSummaryValues.transformed)
  return(raw.distances^2)
}
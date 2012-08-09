PLSRejection<-function(summaryValuesMatrix, trueFreeValuesMatrix, phy, traits, abcTolerance, verbose=TRUE) {
  originalSummaryValues<-summaryStatsLong(phy=phy, data=traits)
  if (verbose) {
    print("Done getting originalSummaryValues") 
  }
  abcDistancesRaw<-sapply(sequence(dim(trueFreeValuesMatrix)[2]), SingleParameterPLSDistanceSquared, summaryValuesMatrix=summaryValuesMatrix, 
         trueFreeValuesMatrix=trueFreeValuesMatrix, originalSummaryValues=originalSummaryValues)
  abcDistancesRawTotal<-apply(abcDistancesRaw, 1, sum)
  abcDistances<-sqrt(abcDistancesRawTotal) #Euclid rules.
  
  acceptedParticles<-trueFreeValuesMatrix[which(abcDistances<=quantile(abcDistances, prob=abcTolerance)), ] #here's where we diy abc
  acceptedDistances<-abcDistances[which(abcDistances<=quantile(abcDistances, prob=abcTolerance))]
  
  particleDataFrame<-data.frame(cbind(rep(1, dim(acceptedParticles)[1]), as.vector(which(abcDistances<=quantile(abcDistances, prob=abcTolerance))), seq(1:dim(acceptedParticles)[1]), rep(0, dim(acceptedParticles)[1]), acceptedDistances, rep(1, dim(acceptedParticles)[1]), acceptedParticles))
  colnames(particleDataFrame)<-c("generation", "attempt", "id", "parentid", "distance", "weight",  paste("param", seq(dim(trueFreeValuesMatrix)[2]),sep=""))
  
  return(list(particleDataFrame=particleDataFrame, abcDistances=abcDistances))
}

SingleParameterPLSDistanceSquared<-function(index, summaryValuesMatrix, trueFreeValuesMatrix, originalSummaryValues) {
  print(paste("my index is ",index))
  trueFreeValuesMatrix<-trueFreeValuesMatrix[,index]
  pls.model<-returnPLSModel(trueFreeValuesMatrix,summaryValuesMatrix)
  print("got pls.model")
  summaryValues.transformed<-PLSTransform(summaryValuesMatrix, pls.model)
  print("did PLS transform on SV")
  originalSummaryValues.transformed<-PLSTransform(originalSummaryValues, pls.model)
  print("did PLS transform on original SV")
  distanceByRow<-function(x,originalSummaryValues.transformed) {
    return(dist(matrix(c(x,originalSummaryValues.transformed),byrow=TRUE,nrow=2))[1])
  }
  raw.distances<-apply(summaryValues.transformed, 1, distanceByRow, originalSummaryValues.transformed=originalSummaryValues.transformed)
  return(raw.distances^2)
}
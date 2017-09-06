#' Partial Least Squares Rejection
#' 
#' This function automatically calculates prior distributions for the Brownian Motion model of
#' evolution.
#' 
#' This function performs the ABC-rejection analysis using an input simulation
#' data. Particles are accepted is they fall sufficiently close to the target
#' data (within the tolerance). Distances are calculated using \code{abcDistance}.
#' 

#' @param summaryValuesMatrix Matrix of summary statistics from simulations

#' @param trueFreeValuesMatrix Matrix of true free values from simulations

#' @param phy Tree (Phylogenetic tree in phylo format)

#' @param traits data matrix with rownames equal to phy

#' @param abcTolerance Proportion of accepted simulations

#' @param verbose option to print progress to screen

#' @param validation Cross Validation procedure for abc

#' @param scale scale for pls.model.list

#' @param variance.cutoff variance cutoff for pls.model.list

#' @return Returns a list of the particle data frame and abc distances.

#' @author Brian O'Meara and Barb Banbury

# @references O'Meara and Banbury, unpublished
# @keywords PLSRejection doRun doRun_rej abc

#' @examples
#' 
#' #PLSRejection(summaryValuesMatrix, trueFreeValuesMatrix, phy, traits, abcTolerance)
#' 


#' @name intrinsicModels
#' @rdname intrinsicModels
#' @export
PLSRejection<-function(summaryValuesMatrix, trueFreeValuesMatrix, phy, traits, abcTolerance, verbose=TRUE, validation="CV", scale=TRUE, variance.cutoff=95) {
  originalSummaryValues<-summaryStatsLong(phy=phy, traits=traits)
  if (verbose) {
    print("Done getting originalSummaryValues") 
  }
  abcDistancesRaw<-sapply(sequence(dim(trueFreeValuesMatrix)[2]), SingleParameterPLSDistanceSquared, summaryValuesMatrix=summaryValuesMatrix, 
         trueFreeValuesMatrix=trueFreeValuesMatrix, originalSummaryValues=originalSummaryValues, validation=validation, scale=scale, variance.cutoff=variance.cutoff)
  abcDistancesRawTotal<-apply(abcDistancesRaw, 1, sum)
  abcDistances<-sqrt(abcDistancesRawTotal) #Euclid rules.
  
  acceptedParticles<-trueFreeValuesMatrix[which(abcDistances<=quantile(abcDistances, prob=abcTolerance)), ] #here's where we diy abc
  acceptedDistances<-abcDistances[which(abcDistances<=quantile(abcDistances, prob=abcTolerance))]
  
  particleDataFrame<-data.frame(cbind(rep(1, dim(acceptedParticles)[1]), as.vector(which(abcDistances<=quantile(abcDistances, prob=abcTolerance))), seq(1:dim(acceptedParticles)[1]), rep(0, dim(acceptedParticles)[1]), acceptedDistances, rep(1, dim(acceptedParticles)[1]), acceptedParticles))
  colnames(particleDataFrame)<-c("generation", "attempt", "id", "parentid", "distance", "weight",  paste("param", seq(dim(trueFreeValuesMatrix)[2]),sep=""))
  
  return(list(particleDataFrame=particleDataFrame, abcDistances=abcDistances))
}

SingleParameterPLSDistanceSquared<-function(index, summaryValuesMatrix, trueFreeValuesMatrix, originalSummaryValues, validation="CV", scale=TRUE, variance.cutoff=95) {
  trueFreeValuesMatrix<-trueFreeValuesMatrix[,index]
  pls.model<-returnPLSModel(trueFreeValuesMatrix,summaryValuesMatrix, validation=validation, scale=scale, variance.cutoff=variance.cutoff)
  summaryValues.transformed<-PLSTransform(summaryValuesMatrix, pls.model)
  originalSummaryValues.transformed<-PLSTransform(originalSummaryValues, pls.model)
  distanceByRow<-function(x,originalSummaryValues.transformed) {
    return(dist(matrix(c(x,originalSummaryValues.transformed),byrow=TRUE,nrow=2))[1])
  }
  raw.distances<-apply(summaryValues.transformed, 1, distanceByRow, originalSummaryValues.transformed=originalSummaryValues.transformed)
  return(raw.distances^2)
}

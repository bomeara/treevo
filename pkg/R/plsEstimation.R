plsEstimation<-function(trueFreeValues, boxcoxSummaryValues, vipthresh) {
#----------------- Find best set of summary stats to use for this problem. (Start) -----------------
			#Use mixOmics to to find the optimal set of summary stats. Store this info in the todo vector. Note that this uses a different package (mixOmics rather than pls than that used by Weggman et al. because this package can calculate variable importance in projection and deals fine with NAs)
	library("mixOmics")
	plsResult<-pls(Y=trueFreeValues, X=boxcoxSummaryValues)
	vipResult<-vip(plsResult)
	todo<-rep(1, dim(boxcoxSummaryValues)[2]) #initialize the vector that indicates which summary stats to include
			
	summaryIndexOffset=0 #since R excludes invariant columns from regression, this offests so we don't try to extract from these columns
	nearZeroVarVector<-mixOmics:::nearZeroVar(boxcoxSummaryValues)
	nearZeroVarVector<-nearZeroVarVector$Position
	#print(nearZeroVarVector)
	for (summaryIndex in 1:dim(boxcoxSummaryValues)[2]) {
		if (summaryIndex %in% nearZeroVarVector) {
			summaryIndexOffset=summaryIndexOffset+1
			todo[summaryIndex]<-0 #exclude this summary stat because it lacks variation
		}	
		else if (max(vipResult[summaryIndex-summaryIndexOffset, ]) < vipthresh) {
			todo[summaryIndex]<-0 #exclude this summary stat, because it is too unimportant
		}	
	}
	prunedSummaryValues<-boxcoxSummaryValues[, which(todo>0)]
	prunedPlsResult<-pls(Y=trueFreeValues, X=prunedSummaryValues)
	return(list(prunedSummaryValues=prunedSummaryValues,prunedPlsResult=prunedPlsResult,todo=todo))			
}

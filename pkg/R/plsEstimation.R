plsEstimation<-function(trueFreeValues, boxcoxSummaryValues, ncomp) {
	library("mixOmics", quietly=T)
	plsResult<-pls(Y=trueFreeValues, X=boxcoxSummaryValues, ncomp=ncomp)
	vipResult<-vip(plsResult)
			
	summaryIndexOffset=0 #since R excludes invariant columns from regression, this offests so we don't try to extract from these columns
	nearZeroVarVector<-mixOmics:::nearZeroVar(boxcoxSummaryValues)
	nearZeroVarVector<-nearZeroVarVector$Position
	#print(nearZeroVarVector)
	prunedSummaryValues<-boxcoxSummaryValues
	#prunedPlsResult<-pls(Y=trueFreeValues, X=prunedSummaryValues)  #no need to do another pls(), so here for the prunedPlsResult, we can just use plsResult
	return(list(prunedSummaryValues=prunedSummaryValues,plsResult=plsResult, vipResult=vipResult))			
}

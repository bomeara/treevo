summarizeRejection<-function(summaryValues, trueFreeValues, vipthresh, traits, todo, abcMethod, abcTolerance, jobName){
	#get boxcoxLambda and boxcoxAddition
	library(abc, quietly=T)
	boxcoxEstimates<-boxcoxEstimation(summaryValues)
	boxcoxAddition<-boxcoxEstimates$boxcoxAddition
	boxcoxLambda<-boxcoxEstimates$boxcoxLambda
	boxcoxSummaryValues<-boxcoxEstimates$boxcoxSummaryValues
	#get pls parameters
	save(trueFreeValues, file="tFV")
	save(boxcoxSummaryValues, file="bcSV")
	plsEstimates<-plsEstimation(trueFreeValues, boxcoxSummaryValues, vipthresh)
	prunedPlsResult<-plsEstimates$prunedPlsResult
	prunedSummaryValues<-plsEstimates$prunedSummaryValues
	todo<-plsEstimates$todo
	#do abc() using transformed summary stat for observed data and set of transformed summary stats and true values for simulated data
	originalSummaryStats<-summaryStatsLong(phy=phy, data=traits, todo=todo, jobName=jobName)[which(todo==1)]
	abcResults<-abc(target=originalSummaryStats, param=trueFreeValues, sumstat=prunedSummaryValues, tol=abcTolerance, method=abcMethod)
	save(originalSummaryStats,trueFreeValues, prunedSummaryValues,file=paste("prunedSimulations",jobName,".Rdata",sep=""))

	return(abcResults)
}

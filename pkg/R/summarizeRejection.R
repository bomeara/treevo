summarizeRejection<-function(summaryValuesMatrix, trueFreeValuesMatrix, phy, traits, vipthresh, abcTolerance){
	#get boxcoxLambda and boxcoxAddition
	library(abc, quietly=T)
	boxcoxEstimates<-boxcoxEstimation(summaryValuesMatrix)
	boxcoxAddition<-boxcoxEstimates$boxcoxAddition
	boxcoxLambda<-boxcoxEstimates$boxcoxLambda
	boxcoxSummaryValuesMatrix<-boxcoxEstimates$boxcoxSummaryValuesMatrix
	boxcoxOriginalSummaryStats<-boxcoxTransformation(summaryStatsLong(phy=phy, data=traits),boxcoxAddition, boxcoxLambda)
	plsEstimates<-vector("list", length=dim(trueFreeValuesMatrix)[2])
	plsResult<-vector("list", length=dim(trueFreeValuesMatrix)[2])
	vipResult<-matrix(ncol=dim(trueFreeValuesMatrix)[2], nrow=length(boxcoxOriginalSummaryStats))
  abcDistances<-matrix(ncol=dim(trueFreeValuesMatrix)[2], nrow=dim(trueFreeValuesMatrix)[1])
  for(plsEst in 1:dim(trueFreeValuesMatrix)[2]){
    plsEstimates[[plsEst]]<-plsEstimation(trueFreeValuesMatrix[,plsEst], boxcoxSummaryValuesMatrix, ncomp=1)
    plsResult[[plsEst]]<-plsEstimates[[plsEst]]$plsResult
    vipResult[,plsEst]<-plsEstimates[[plsEst]]$vipResult  
    abcDistances[,plsEst]<-abc(target=boxcoxOriginalSummaryStats[which(vipResult[,plsEst]>vipthresh)],param=trueFreeValuesMatrix[,plsEst],sumstat=boxcoxSummaryValuesMatrix[,which(vipResult[,plsEst]>vipthresh)],tol=1.0,method="rejection")$dist  #why is the last one == 0?
  }
  sumStandardizedDistances<-apply(abcDistances, 1, sum)/max(apply(abcDistances, 1, sum))
	abcResults<-trueFreeValuesMatrix[which(sumStandardizedDistances<abcTolerance),]  #unadj.values
	
  
  #abcResults<-abc(target=boxcoxOriginalSummaryStats[which(vipResult[,plsEst]>vipthresh)],param=trueFreeValuesMatrix[,plsEst],sumstat=boxcoxSummaryValuesMatrix[,which(vipResult[,plsEst]>vipthresh)],tol=1.0,method="rejection") #rerun abc() with kept values and tolerance=1.0?
  
  return(abcResults)
}

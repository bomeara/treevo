
boxcoxplsSummary<-function(summaryValues, plsResult, boxcoxLambda, boxcoxAddition, whichVip) {
	for (summaryValueIndex in 1:length(summaryValues)) {
		summaryValues[summaryValueIndex]<-(summaryValues[summaryValueIndex] + boxcoxAddition[summaryValueIndex])^boxcoxLambda[summaryValueIndex]
	}
	predictResult<-c()
	for(plsIndex in 1:length(plsResult)){
		prunedSummaryValues<-summaryValues[whichVip[,plsIndex]]
		predictSingleResult<-(predict(plsResult[[plsIndex]], matrix(prunedSummaryValues, nrow=1))$predict[, , 1])
		predictResult<-append(predictResult, predictSingleResult)
	}
	
	predictResult
}


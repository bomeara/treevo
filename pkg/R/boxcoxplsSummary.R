
boxcoxplsSummary<-function(summaryValues, plsResult, boxcoxLambda, boxcoxAddition) {
	for (summaryValueIndex in 1:length(summaryValues)) {
		summaryValues[summaryValueIndex]<-(summaryValues[summaryValueIndex] + boxcoxAddition[summaryValueIndex])^boxcoxLambda[summaryValueIndex]
	}
	prunedSummaryValues<-summaryValues
	predictResult<-(predict(plsResult, matrix(prunedSummaryValues, nrow=1))$predict[, , 1])
	predictResult
}


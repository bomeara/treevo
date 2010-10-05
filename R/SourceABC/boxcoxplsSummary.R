
boxcoxplsSummary<-function(todo, summaryValues, prunedPlsResult, boxcoxLambda, boxcoxAddition) {
	for (summaryValueIndex in 1:length(summaryValues)) {
		summaryValues[summaryValueIndex]<-(summaryValues[summaryValueIndex] + boxcoxAddition[summaryValueIndex])^boxcoxLambda[summaryValueIndex]
	}
	prunedSummaryValues<-summaryValues[which(todo>0)]
	predictResult<-(predict(prunedPlsResult, matrix(prunedSummaryValues, nrow=1))$predict[, , 1])
	predictResult
}


boxcoxTransformation<-function(summaryValues, boxcoxAddition, boxcoxLambda) { #yes, a row of summary values
	for (summaryValueIndex in 1:length(summaryValues)) {
        summaryValues[summaryValueIndex] <- (summaryValues[summaryValueIndex] + 
            boxcoxAddition[summaryValueIndex])^boxcoxLambda[summaryValueIndex]
    }
    return(summaryValues)
}

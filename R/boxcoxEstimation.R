boxcoxEstimation<-function(summaryValues){
	library("car", quietly=T)
	boxcoxLambda<-rep(NA, dim(summaryValues)[2])
	boxcoxAddition<-rep(NA, dim(summaryValues)[2])
	for (summaryValueIndex in 1:dim(summaryValues)[2]) {
		boxcoxAddition[summaryValueIndex]<-0
		lowValue<-min(summaryValues[, summaryValueIndex])-4*sd(summaryValues[, summaryValueIndex])
		if (lowValue<=0) {
			boxcoxAddition[summaryValueIndex]<-4*abs(lowValue) #just for some protection against low values, since box.cox needs non-negative values
		}
		summary<-summaryValues[, summaryValueIndex]+boxcoxAddition[summaryValueIndex]
		boxcoxLambda[summaryValueIndex]<-1
		if(sd(summaryValues[, summaryValueIndex])>0) { #box.cox fails if all values are identical
			newLambda<-as.numeric(try(powerTransform(summary,method="Nelder-Mead")$lambda)) #new car uses powerTransform instead of box.cox.powers
			if (!is.na(newLambda)) {
				boxcoxLambda[summaryValueIndex]<-newLambda
			}
		}
		summaryValues[, summaryValueIndex]<-summary^boxcoxLambda[summaryValueIndex]
	}
	return(list(boxcoxAddition=boxcoxAddition,boxcoxLambda=boxcoxLambda,boxcoxSummaryValues=summaryValues))
}

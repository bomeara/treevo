#' Box-Cox Transformation
#' 
#' This function Box-Cox transforms summary values from a single simulation
#' 
#' 
#' @param summaryValues Vector of summary statistics
#' @param boxcoxAddition Vector of boxcox additions from boxcoxEstimation
#' @param boxcoxLambda Vector of boxcox lambda values from boxcoxEstimation
#' @return Returns a vector of Box-Cox transformed summary statistics from a
#' single observation.
#' @author Brian O'Meara and Barb Banbury
# @references O'Meara and Banbury, unpublished; Bates et al. 2009
# @keywords boxcoxTransformation Box-Cox
boxcoxTransformation<-function(summaryValues, boxcoxAddition, boxcoxLambda) { #yes, a row of summary values
	for (summaryValueIndex in 1:length(summaryValues)) {
        summaryValues[summaryValueIndex] <- (summaryValues[summaryValueIndex] + 
            boxcoxAddition[summaryValueIndex])^boxcoxLambda[summaryValueIndex]
    }
    return(summaryValues)
}

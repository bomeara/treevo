#' Box-Cox Transformation for Simulations
#' 
#' \code{boxcoxTransformation} Box-Cox transforms summary values from a single simulation, while
#' \code{boxcoxTransformationMatrix} Box-Cox transforms summary values from multiple simulations.
#' 
#' 

#' @param summaryValuesMatrix Matrix of summary statistics from simulations

#' @param summaryValuesVector Vector of summary statistics from a single simulation

#' @param boxcoxAddition Vector of boxcox additions from boxcoxEstimation

#' @param boxcoxLambda Vector of boxcox lambda values from boxcoxEstimation

#' @return \code{boxcoxTransformation} returns a vector of Box-Cox transformed summary statistics from a
#' single observation. \code{boxcoxTransformationMatrix} returns a matrix of Box-Cox transformed summary statistics with the
#' same dimensions as \code{summaryValuesMatrix}.

#' @author Brian O'Meara and Barb Banbury

# @references O'Meara and Banbury, unpublished; Bates et al. 2009
# @keywords boxcoxEstimation Box-Cox

#' @examples
#' 
#' data(simRun)
#' 
#' boxcoxTransformation(summaryValuesMatrix)
#' 
#' # boxcoxTransformationMatrix needs examples
#'

#' @name boxcoxTransformation
#' @rdname boxcoxTransformation
#' @export
boxcoxTransformation<-function(summaryValues, boxcoxAddition, boxcoxLambda) { #yes, a row of summary values
	for (summaryValueIndex in 1:length(summaryValues)) {
        summaryValues[summaryValueIndex] <- (summaryValues[summaryValueIndex] + 
            boxcoxAddition[summaryValueIndex])^boxcoxLambda[summaryValueIndex]
    }
    return(summaryValues)
}

#' @rdname boxcoxTransformation
#' @export
boxcoxTransformationMatrix<-function (summaryValuesMatrix) {
  #library("car", quietly = T)
  boxcoxLambda <- rep(NA, dim(summaryValuesMatrix)[2])
  boxcoxAddition <- rep(NA, dim(summaryValuesMatrix)[2])
  for (summaryValueIndex in 1:dim(summaryValuesMatrix)[2]) {
    boxcoxAddition[summaryValueIndex] <- 0
    lowValue <- min(summaryValuesMatrix[, summaryValueIndex]) - 4 * sd(summaryValuesMatrix[, summaryValueIndex])
    if (lowValue <= 0) {
      boxcoxAddition[summaryValueIndex] <- 4 * abs(lowValue)
    }
    summary <- summaryValuesMatrix[, summaryValueIndex] + boxcoxAddition[summaryValueIndex]
    boxcoxLambda[summaryValueIndex] <- 1
    if (sd(summaryValuesMatrix[, summaryValueIndex]) > 0) {
      newLambda <- as.numeric(try(powerTransform(summary, method = "Nelder-Mead")$lambda))
      if (!is.na(newLambda)) {
        boxcoxLambda[summaryValueIndex] <- newLambda
      }
    }
    summaryValuesMatrix[, summaryValueIndex] <- summary^boxcoxLambda[summaryValueIndex]
  }
  return(list(boxcoxAddition = boxcoxAddition, boxcoxLambda = boxcoxLambda, 
              boxcoxSummaryValuesMatrix = summaryValuesMatrix))
}

boxcoxEstimation<-function (summaryValuesMatrix) {
  library("car", quietly = T)
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

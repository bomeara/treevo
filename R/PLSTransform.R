#' PLS Transform
#' 
#' Uses results from the PLS model to transform summary values
#' 
#' This function uses the scores from the pls model to transform the summary
#' statistics.
#' 
#' @param summaryValuesMatrix Matrix of summary statistics from simulations
#' @param pls.model results from returnPLSModel
#' @return Returns transformed summary statistics
#' @author Brian O'Meara and Barb Banbury
# @references O'Meara and Banbury, unpublished
# @keywords PLSTransform PLS
#' @examples
#' 
#' #PLSTransform(summaryValuesMatrix, pls.model)
#' 
PLSTransform<-function(summaryValuesMatrix, pls.model) {
  if (class(summaryValuesMatrix)!="matrix") {
    summaryValuesMatrix<-matrix(summaryValuesMatrix,ncol=max(c(1,length(summaryValuesMatrix)),na.rm=TRUE)) 
  }
  predict(pls.model,summaryValuesMatrix,type="scores")
}

#' Partial Least Squares Transform
#' 
#' This function uses results from a Partial Least Squares (PLS) model to transform summary values.
#' 

#' @param summaryValuesMatrix Matrix of summary statistics from simulations

#' @param pls.model Results output by from \code{\link{returnPLSModel}}

#' @return Returns transformed summary statistics

#' @author Brian O'Meara and Barb Banbury

# @references O'Meara and Banbury, unpublished
# @keywords PLSTransform PLS

#' @examples
#' 
#' #PLSTransform(summaryValuesMatrix, pls.model)
#' 

#' @name PLSTransform
#' @rdname PLSTransform
#' @export
PLSTransform<-function(summaryValuesMatrix, pls.model) {
  if (class(summaryValuesMatrix)!="matrix") {
    summaryValuesMatrix<-matrix(summaryValuesMatrix,ncol=max(c(1,length(summaryValuesMatrix)),na.rm=TRUE)) 
  }
  predict(pls.model,summaryValuesMatrix,type="scores")
}

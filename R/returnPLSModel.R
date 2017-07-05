#' Return PLS Model
#' 
#' Returns a PLS Model
#' 
#' This function runs a PLS regression on each free param in the model, unlike
#' a true multivariate PLS regression.  For ABC, this seems to result in much
#' better results, without one parameter dominating the combined variance.
#' 

#' @param trueFreeValuesMatrix Matrix of true free values from simulations

#' @param summaryValuesMatrix Matrix of summary statistics from simulations

#' @param validation Cross Validation procedure for abc

#' @param scale scale for pls.model.list

#' @param variance.cutoff variance cutoff for pls.model.list

#' @return Returns a pls model

#' @author Brian O'Meara and Barb Banbury

#' @references O'Meara and Banbury, unpublished

#' @keywords returnPLSModel PLSTransform PLS

#' @examples
#' 
#' # returnPLSModel(trueFreeValuesMatrix,summaryValuesMatrix,
#' #   validation="CV", scale=scale, variance.cutoff=variance.cutoff)
#' 
returnPLSModel<-function(trueFreeValuesMatrix, summaryValuesMatrix, validation="CV", scale=TRUE, variance.cutoff=95) {
  #note that this assumes that trueFreeValues is for a single param at a time, which works MUCH better
  trueFreeValuesMatrix<-trueFreeValuesMatrix
  if (dim(summaryValuesMatrix)[2]>1) {
    warning("in practice, doing PLS works best if you do each free parameter separately, so one does not dominate")
  }
  if (class(trueFreeValuesMatrix)!="matrix") {
    trueFreeValuesMatrix<-matrix(trueFreeValuesMatrix,nrow=max(c(1,length(trueFreeValuesMatrix)),na.rm=TRUE)) 
  }
  pls.model <- plsr(trueFreeValuesMatrix~summaryValuesMatrix,validation=validation,scale=scale) #scaling is important
  explained.variance <-cumsum(sort(attr(scores(pls.model),"explvar"),decreasing=TRUE))
  ncomp.final<-min(c(as.numeric(which(explained.variance>=variance.cutoff)[1]),length(explained.variance)),na.rm=TRUE) #min is to deal with case of never explaining >95%
  pls.model.final <- plsr(trueFreeValuesMatrix~summaryValuesMatrix,ncomp=ncomp.final, validation="none",scale=scale) #now rerun with the ideal number of components
  return(pls.model.final)
}

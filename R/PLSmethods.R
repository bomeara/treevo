#' Fitting Univariate Partial Least Squares Models to Free Parameters in ABC
#' 
#' Function \code{returnPLSModel} fits a PLS regression (using \link{\code{plsr}}) individually to each freely varying parameter of a model, unlike
#' a true multivariate PLS regression. A secondary step than limits the number
#' of components to those that explain some minimum cumulative percentage
#' of variance (see argument \code{variance.cutoff}).  For ABC, this seems to result in much
#' better results, without one parameter dominating the combined variance.
#' 
#' Function \code{PLSTransform} uses results from a Partial Least Squares (PLS) model fit with \code{returnPLSModel} to transform summary values.
#' 

#' @param trueFreeValuesMatrix Matrix of true free values from simulations.

#' @param summaryValuesMatrix Matrix of summary statistics from simulations.

#' @param validation Character argument controlling what validation procedure is used by \link{\code{plsr}}.
#' Default is \code{"CV"} for cross-validation.

#' @param scale This argument is passed to \link{\code{plsr}}.  It may be a numeric vector, or logical. If numeric vector, 
#' the input is scaled by dividing each variable with the corresponding element of scale. 
#' If \code{scale = TRUE}, the inpus is scaled by dividing each variable by its sample standard deviation. 
#' If cross-validation is selected (the default for \code{returnPLSModel}),
#' scaling by the standard deviation is done for every segment.
 
#' @param variance.cutoff Minimum threshold percentage of variance explained for the
#' number of components included in the final PLS model fit. This value is a 
#' percentage and must be between 0 and 100. Default is 95 percent.

# for \code{pls.model.list} # Used in doc for both scale and variance.cutoff
	# I HAVE NO IDEA WHAT THIS MEANS? WHAT IS pls.model.list???


#' @param pls.model Output from \link{\code{returnPLSModel}}.

#' @seealso
#' Function \code{returnPLSModel} effectively wraps function \link{\code{plsr}} from package \code{mvr}.

#' @return 
#' Function \code{returnPLSModel} returns a PLS model, and function \code{PLSTransform} returns transformed summary statistics.


#' @author Brian O'Meara and Barb Banbury

#' @examples
#' 
#' data(simRun)
#' 
#' PLSmodel<-returnPLSModel(trueFreeValuesMatrix,
#' 	  summaryValuesMatrix=summaryValuesMatrix,
#'    validation="CV", scale=TRUE, variance.cutoff=95)
#' 
#' 
#' PLSTransform(summaryValuesMatrix=summaryValuesMatrix, pls.model=PLSmodel)
#' 

#' @name PLSmethods
#' @rdname PLSmethods
#' @export
returnPLSModel<-function(trueFreeValuesMatrix, summaryValuesMatrix, validation="CV", scale=TRUE, variance.cutoff=95) {
  #note that this assumes that trueFreeValues is for a single param at a time, which works MUCH better
  
  trueFreeValuesMatrix<-trueFreeValuesMatrix
  if (dim(summaryValuesMatrix)[2]>1) {
    warning("in practice, doing PLS works best if you do each free parameter separately, so one does not dominate")
  }
  if (class(trueFreeValuesMatrix)!="matrix") {
    trueFreeValuesMatrix<-matrix(trueFreeValuesMatrix,nrow=max(c(1,length(trueFreeValuesMatrix)),na.rm=TRUE)) 
  }
  #
  #scaling is important
  pls.model <- plsr(trueFreeValuesMatrix~summaryValuesMatrix,validation=validation,scale=scale) 
  explained.variance <-cumsum(sort(attr(scores(pls.model),"explvar"),decreasing=TRUE))
  ncomp.final<-min(c(as.numeric(which(explained.variance>=variance.cutoff)[1]),
	length(explained.variance)),na.rm=TRUE) #min is to deal with case of never explaining >95%
  #
  # now rerun with the ideal number of components
  pls.model.final <- plsr(trueFreeValuesMatrix~summaryValuesMatrix,
	ncomp=ncomp.final, validation="none",scale=scale) 
  return(pls.model.final)
}


#' @rdname PLSmethods
#' @export
PLSTransform<-function(summaryValuesMatrix, pls.model) {
  if (class(summaryValuesMatrix)!="matrix") {
    summaryValuesMatrix<-matrix(summaryValuesMatrix,ncol=max(c(1,length(summaryValuesMatrix)),na.rm=TRUE)) 
  }
  predict(pls.model,summaryValuesMatrix,type="scores")
}
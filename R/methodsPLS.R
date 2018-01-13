#' Fitting Univariate Partial Least Squares Models to Free Parameters in ABC
#' 
#' Function \code{returnPLSModel} fits a PLS regression (using 
#' \code{\link{plsr}}) individually to each freely varying parameter of a model, unlike
#' a true multivariate PLS regression. A secondary step than limits the number
#' of components to those that explain some minimum cumulative percentage
#' of variance (see argument \code{variance.cutoff}).  For ABC, this seems to result in much
#' better results, without one parameter dominating the combined variance.
#' 
#' Function \code{PLSTransform} uses results from a Partial Least Squares (PLS)
#' model fit with \code{returnPLSModel} to transform summary values.
#' 

#' @param trueFreeValuesMatrix Matrix of true free values from simulations.

#' @param summaryValuesMatrix Matrix of summary statistics from simulations.

#' @param validation Character argument controlling what validation procedure is used by \code{\link{plsr}}.
#' Default is \code{"CV"} for cross-validation.

#' @param scale This argument is passed to \code{\link{plsr}}.  It may be a numeric vector, or logical. If numeric vector, 
#' the input is scaled by dividing each variable with the corresponding element of scale. 
#' If \code{scale = TRUE}, the inpus is scaled by dividing each variable by its sample standard deviation. 
#' If cross-validation is selected (the default for \code{returnPLSModel}),
#' scaling by the standard deviation is done for every segment.
 
#' @param variance.cutoff Minimum threshold percentage of variance explained for the
#' number of components included in the final PLS model fit. This value is a 
#' percentage and must be between 0 and 100. Default is 95 percent.

#' @param ... Additional arguments, passed to \code{\link{plsr}}. In particular, if the number of observations is less than
#' 10 (such as for the example below), it may be necessary to pass a revised \code{segments} argument, which 
#' \code{plsr} passes to \code{\link{mvrCv}}, as the default \code{segments} value is 10 (and there cannot be more segments
#' than observations). Note that this should be an atypical situation outside of examples, as the number of observations
#' should often be much greater than 10. 

# for \code{pls.model.list} # Used in doc for both scale and variance.cutoff
	# I HAVE NO IDEA WHAT THIS MEANS? WHAT IS pls.model.list???


#' @param pls.model Output from \code{\link{returnPLSModel}}.

#' @seealso
#' Function \code{returnPLSModel} effectively wraps function \code{\link{plsr}} from package \code{mvr}.

#' @return 
#' Function \code{returnPLSModel} returns a PLS model, and function \code{PLSTransform} returns transformed summary statistics.


#' @author Brian O'Meara and Barb Banbury

#' @examples
#' \donttest{
#' set.seed(1)
#' data(simRunExample)
#' 
#' # example simulation
#' 
#' nSimulations<-6
#' 
#' simDataParallel<-parallelSimulateWithPriors( 
#'   nrepSim=nSimulations, multicore=FALSE, coreLimit=1, 
#'   phy=simPhy,
#'   intrinsicFn=brownianIntrinsic,
#'   extrinsicFn=nullExtrinsic,
#'   startingPriorsFns="normal",
#'   startingPriorsValues=matrix(c(mean(simChar[,1]), sd(simChar[,1]))),
#'   intrinsicPriorsFns=c("exponential"),
#'   intrinsicPriorsValues=matrix(c(10, 10), nrow=2, byrow=FALSE),
#'   extrinsicPriorsFns=c("fixed"),
#'   extrinsicPriorsValues=matrix(c(0, 0), nrow=2, byrow=FALSE), 
#'   generation.time=100000,
#'   checkpointFile=NULL, checkpointFreq=24,
#'   verbose=FALSE,
#'   freevector=NULL, taxon.df=NULL)
#' 
#' nParFree<-sum(attr(simDataParallel,"freevector"))
#' 
#' #separate the simulation results: 'true' generating parameter values from the summary values
#' trueFreeValuesMat<-simDataParallel[,1:nParFree]
#' summaryValuesMat<-simDataParallel[,-1:-nParFree]
#' 
#' PLSmodel<-returnPLSModel(trueFreeValuesMatrix=trueFreeValuesMat,
#' 	  	summaryValuesMatrix=summaryValuesMat,
#'    	validation="CV", scale=TRUE, variance.cutoff=95 , segments = nSimulations)
#' 
#' PLSmodel
#' 
#' PLSTransform(summaryValuesMatrix=summaryValuesMat, pls.model=PLSmodel)
#' } 
#' 



#' @name methodsPLS
#' @rdname methodsPLS
#' @export
returnPLSModel<-function(trueFreeValuesMatrix, summaryValuesMatrix, validation="CV", scale=TRUE, variance.cutoff=95, ...) {
  #note that this assumes that trueFreeValues is for a single param at a time, which works MUCH better
  
  trueFreeValuesMatrix<-trueFreeValuesMatrix
  if (dim(summaryValuesMatrix)[2]>1) {
    warning("in practice, doing PLS works best if you do each free parameter separately, so one parameter does not dominate")
  }
  if (class(trueFreeValuesMatrix)!="matrix") {
    trueFreeValuesMatrix<-matrix(trueFreeValuesMatrix,nrow=max(c(1,length(trueFreeValuesMatrix)),na.rm=TRUE)) 
  }
  #
  #scaling is important
  pls.model <- makeQuiet(plsr(trueFreeValuesMatrix~summaryValuesMatrix,validation=validation,scale=scale,...)) 
  explained.variance <-cumsum(sort(attr(scores(pls.model),"explvar"),decreasing=TRUE))
  ncomp.final<-min(c(as.numeric(which(explained.variance>=variance.cutoff)[1]),
	length(explained.variance)),na.rm=TRUE) #min is to deal with case of never explaining >95%
  #
  # now rerun with the ideal number of components
  pls.model.final <- makeQuiet(plsr(trueFreeValuesMatrix~summaryValuesMatrix,
	ncomp=ncomp.final, validation="none",scale=scale,...)) 
  return(pls.model.final)
}


#' @rdname methodsPLS
#' @export
PLSTransform<-function(summaryValuesMatrix, pls.model) {
  if (class(summaryValuesMatrix)!="matrix") {
    summaryValuesMatrix<-matrix(summaryValuesMatrix,ncol=max(c(1,length(summaryValuesMatrix)),na.rm=TRUE)) 
  }
  predict(pls.model,summaryValuesMatrix,type="scores")
}
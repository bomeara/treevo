#' Box-Cox Transformation for Simulation Summary Values
#' 
#' \code{boxcoxTransformation} Box-Cox transforms summary values from a single simulation, while
#' \code{boxcoxTransformationMatrix} Box-Cox transforms summary values from multiple simulations. The
#' output of \code{boxcoxTransformationMatrix} needs to be provided as input for \code{boxcoxTransformation}.
#' 
#' 

#' @param summaryValuesMatrix Matrix of summary statistics from simulations

#' @param summaryValuesVector Vector of summary statistics from a single simulation

#' @param boxcoxAddition Vector of boxcox additions from \code{boxcoxTransformationMatrix}

#' @param boxcoxLambda Vector of boxcox lambda values from \code{boxcoxTransformationMatrix}

#' @return \code{boxcoxTransformation} returns a vector of Box-Cox transformed summary statistics from a
#' single observation. \code{boxcoxTransformationMatrix} returns a matrix of Box-Cox transformed summary statistics with the
#' same dimensions as \code{summaryValuesMatrix}.

#' @author Brian O'Meara and Barb Banbury

# @references O'Meara and Banbury, unpublished; Bates et al. 2009
# @keywords boxcoxEstimation Box-Cox

#' @examples
#' \donttest{
#' set.seed(1)
#' data(simRunExample)
#' 
#' # example simulation
#' 
#' simDataParallel<-parallelSimulateWithPriors(
#'   nrepSim=5, multicore=FALSE, coreLimit=1,
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
#'   freevector=NULL, taxonDF=NULL)
#' 
#' nParFree<-sum(attr(simDataParallel,"freevector"))
#' 
#' #separate the simulation results: 'true' generating parameter values from the summary values
#' summaryValuesMat<-simDataParallel[,-1:-nParFree]
#' 
#' boxTranMat<-boxcoxTransformationMatrix(summaryValuesMatrix = summaryValuesMat)
#' boxTranMat
#' 
#' boxcoxTransformation(summaryValuesVector=summaryValuesMat[,1],
#'  boxcoxAddition=boxTranMat$boxcoxAddition, boxcoxLambda=boxTranMat$boxcoxLambda)
#' }


#' @name boxcoxTransformation
#' @rdname boxcoxTransformation
#' @export
boxcoxTransformation<-function(summaryValuesVector, boxcoxAddition, boxcoxLambda) {
	#yes, a row of summary values
	for (summaryValueIndex in 1:length(summaryValuesVector)) {
        summaryValuesVector[summaryValueIndex] <- (summaryValuesVector[summaryValueIndex] +
            boxcoxAddition[summaryValueIndex])^boxcoxLambda[summaryValueIndex]
    }
    return(summaryValuesVector)
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
    summaryVM <- summaryValuesMatrix[, summaryValueIndex] + boxcoxAddition[summaryValueIndex]
    boxcoxLambda[summaryValueIndex] <- 1
    if (sd(summaryValuesMatrix[, summaryValueIndex]) > 0) {
		# this alternative is thanks to https://stackoverflow.com/questions/33999512/how-to-use-the-box-cox-power-transformation-in-r/34002020
      #newLambda <- makeQuiet(as.numeric(try(car::powerTransform(summaryVM, method = "Nelder-Mead")$lambda)))
      bc <- MASS::boxcox(variable ~ 1, data=data.frame(variable=summaryVM),
		 lambda=seq(-20, 20, 1/10000), plotit=FALSE)
	  newLambda<-bc$x[which(max(bc$y)==bc$y)[1]]
      if (!is.na(newLambda)) {
        boxcoxLambda[summaryValueIndex] <- newLambda
      }
    }
    summaryValuesMatrix[, summaryValueIndex] <- summaryVM^boxcoxLambda[summaryValueIndex]
  }
  res<-list(boxcoxAddition = boxcoxAddition, boxcoxLambda = boxcoxLambda,
       boxcoxSummaryValuesMatrix = summaryValuesMatrix)
  return(res)
}

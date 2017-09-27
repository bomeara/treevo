#' TreEvo--abc for comparative methods
#'
#' A package for applying Approximate Bayesian Computation to estimating parameters of trait evolution in comparative analyses.
#'
#' \tabular{ll}{ Package: \tab TreEvo\cr Type: \tab Package\cr Version: \tab
#' 0.3.3\cr Date: \tab 2012-07-02\cr License: \tab GPL\cr }
#'

#' @name TreEvo-package

#' @aliases TreEvo-package treevo

#' @docType package

#' @author Brian O'Meara, Barb L. Banbury, David W. Bapst
#'

#' Maintainer: David Bapst <dwbapst@gmail.com>

# @keywords treevo abc

#' @examples
#'
#' \donttest{
#'
#' data(simRun)
#'
#' doRun_prc(
#'   phy = simPhy,
#'   traits = simChar,
#'   intrinsicFn=brownianIntrinsic,
#'   extrinsicFn=nullExtrinsic,
#'   startingPriorsFns="normal",
#'   startingPriorsValues=matrix(c(mean(simChar[,1]), sd(simChar[,1]))),
#'   intrinsicPriorsFns=c("exponential"),
#'   intrinsicPriorsValues=matrix(c(10, 10), nrow=2, byrow=FALSE),
#'   extrinsicPriorsFns=c("fixed"),
#'   extrinsicPriorsValues=matrix(c(0, 0), nrow=2, byrow=FALSE),
#'   TreeYears=1000,
#'   standardDevFactor=0.2,
#'   plot=FALSE,
#'   StartSims=300,
#'   epsilonProportion=0.7,
#'   epsilonMultiplier=0.7,
#'   nStepsPRC=5,
#'   numParticles=100,
#'   jobName="examplerun_prc",
#'   stopRule=FALSE,
#'   multicore=FALSE,
#'   coreLimit=1
#' )
#' 
#' }
#' 

# NAMESPACE IMPORTING

#' @import ape
#' @import stats
#' @import geiger
#' @import phytools

#' @importFrom pls plsr scores 
#' @importFrom corpcor pseudoinverse 
#' @importFrom coda effectiveSize HPDinterval as.mcmc 
#' @importFrom foreach foreach getDoParWorkers '%dopar%' 
#' @importFrom rgl plot3d title3d rgl.viewpoint open3d rgl.material triangles3d lines3d spheres3d 
#' @importFrom partitions blockparts 
#' @importFrom car powerTransform 
#' @importFrom mvtnorm dmvnorm 
#' @importFrom graphics curve layout layout.show legend lines plot plot.new points polygon rect segments symbols text title 
#' @importFrom grDevices dev.new dev.off gray pdf rgb 
#' @importFrom methods as setAs



#'
NULL
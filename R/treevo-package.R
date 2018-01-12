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
#' # example analysis, using data simulated with TreEvo
#' 
#' \donttest{
#' 
#' set.seed(1)
#' # let's simulate some data, and then try to infer the parameters using ABC
#' # get a 30-taxon coalescent tree
#' tree<-rcoal(30)
#' # get realistic edge lengths
#' tree$edge.length<-tree$edge.length*20
#' 
#' genRate<-c(0.01)
#' ancState<-c(10)
#' 
#' #Simple Brownian motion
#' simCharOut<-doSimulation(
#' 	phy=tree,
#' 	intrinsicFn=brownianIntrinsic,
#' 	extrinsicFn=nullExtrinsic,
#' 	startingValues=ancState, #root state
#' 	intrinsicValues=genRate,
#' 	extrinsicValues=c(0),
#' 	generation.time=100000,
#' 	saveHistory=FALSE)
#' 
#' # clean for use with doRun
#' simChar<-simCharOut[,"statesmatrix",drop=FALSE]
#' rownames(simChar)<-tree$tip.label[simCharOut$taxonid]
#' 
#' # NOTE: the example analyses below sample too few particles, 
#' 	# over too few steps, with too few starting simulations
#' 	# - all for the sake of examples that reasonably test the functions
#' 
#' # Please set these values to more realistic levels for your analyses!
#' 
#' results<-doRun_prc(
#'   phy = tree,
#'   traits = simChar,
#'   intrinsicFn=brownianIntrinsic,
#'   extrinsicFn=nullExtrinsic,
#'   startingPriorsFns="normal",
#'   startingPriorsValues=matrix(c(mean(simChar[,1]), sd(simChar[,1]))),
#'   intrinsicPriorsFns=c("exponential"),
#'   intrinsicPriorsValues=matrix(c(10, 10), nrow=2, byrow=FALSE),
#'   extrinsicPriorsFns=c("fixed"),
#'   extrinsicPriorsValues=matrix(c(0, 0), nrow=2, byrow=FALSE),
#'   generation.time=100000,
#'   standardDevFactor=0.2,
#'   plot=FALSE,
#'   StartSims=10,
#'   epsilonProportion=0.7,
#'   epsilonMultiplier=0.7,
#'   nStepsPRC=3,
#'   numParticles=20,
#'   jobName="examplerun_prc",
#'   stopRule=FALSE,
#'   multicore=FALSE,
#'   coreLimit=1
#'   )
#' }
#' 
#' 

# NAMESPACE IMPORTING

#' @import ape
#' @import stats
#' @import phytools

#' @importFrom phylolm phylolm
#' @importFrom pls plsr scores 
#' @importFrom corpcor pseudoinverse 
#' @importFrom coda effectiveSize HPDinterval as.mcmc 
#' @importFrom foreach foreach getDoParWorkers '%dopar%' 
#' @importFrom partitions blockparts 
#' @importFrom MASS boxcox 
#' @importFrom mvtnorm dmvnorm 
#' @importFrom graphics curve layout legend lines plot plot.new points polygon rect segments symbols text title 
#' @importFrom grDevices dev.off gray pdf rgb 
#' @importFrom methods as setAs
#' @importFrom utils capture.output

# @importFrom rgl plot3d title3d rgl.viewpoint open3d rgl.material triangles3d lines3d spheres3d 


#' 
NULL


#' TreEvo--abc for comparative methods
#' 
#' A package for applying Approximate Bayesian Computation to estimating parameters of trait evolution in comparative analyses.
#' 
#' \tabular{ll}{ Package: \tab TreEvo\cr Type: \tab Package\cr Version: \tab
#' 0.3.3\cr Date: \tab 2012-07-02\cr License: \tab GPL\cr }
#' 

#' @name TreEvo-package

#' @aliases TreEvo-package treevo TreEvo Treevo

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
#' # get a 20-taxon coalescent tree
#' tree <- rcoal(20)
#' # get realistic edge lengths
#' tree$edge.length <- tree$edge.length*20
#' 
#' genRateExample <- c(0.01)
#' ancStateExample <- c(10)
#' 
#' #Simple Brownian motion
#' simCharExample <- doSimulation(
#'     phy = tree, 
#'     intrinsicFn = brownianIntrinsic, 
#'     extrinsicFn = nullExtrinsic, 
#'     startingValues = ancStateExample, #root state
#'     intrinsicValues = genRateExample, 
#'     extrinsicValues = c(0), 
#'     generation.time = 100000
#'     )

# 
# # clean for use with doRun
# simCharExample <- simCharOut[, "statesmatrix", drop = FALSE]
# rownames(simCharExample) <- tree$tip.label[simCharOut$taxonid]
#

#' 
#' # NOTE: the example analyses below sample too few particles, 
#'     # over too few steps, with too few starting simulations
#'     # - all for the sake of examples that reasonably test the functions
#' 
#' # Please set these values to more realistic levels for your analyses!
#' 
#' resultsBMExample <- doRun_prc(
#'   phy = tree, 
#'   traits = simCharExample, 
#'   intrinsicFn = brownianIntrinsic, 
#'   extrinsicFn = nullExtrinsic, 
#'   startingPriorsFns = "normal", 
#'   startingPriorsValues = list(c(mean(simCharExample[, 1]), sd(simCharExample[, 1]))), 
#'   intrinsicPriorsFns = c("exponential"), 
#'   intrinsicPriorsValues = list(10), 
#'   extrinsicPriorsFns = c("fixed"), 
#'   extrinsicPriorsValues = list(0), 
#'   generation.time = 100000, 
#'   nRuns = 2, 
#'   nStepsPRC = 3, 
#'   numParticles = 20, 
#'      nInitialSimsPerParam = 10, 
#'   jobName = "examplerun_prc", 
#'   stopRule = FALSE, 
#'   multicore = FALSE, 
#'   coreLimit = 1
#'   )
#' }
#' 
#' 

# NAMESPACE IMPORTING

#' @import ape
#' @import stats
#' @import phytools

#' @importFrom parallel detectCores 
#stopCluster

#' @importFrom phylolm phylolm
#' @importFrom pls plsr scores
#' @importFrom coda effectiveSize HPDinterval as.mcmc
#' @importFrom foreach foreach '%dopar%' registerDoSEQ
#' @importFrom fastmatch '%fin%'
#' @importFrom partitions blockparts
#' @importFrom MASS boxcox
#' @importFrom graphics curve layout legend lines plot plot.new points polygon rect segments symbols text title
#' @importFrom grDevices dev.off gray pdf rgb
#' @importFrom methods as setAs
#' @importFrom utils capture.output
#' @importFrom RcppZiggurat zrnorm

# package rpgm deprecated on CRAN as of Jan 2019 (04-25-19)
# @importFrom rpgm rpgm.rnorm

# package mvr doesn't exist, was it renamed pls?
# @importFrom mvr plsr


# @importFrom plotrix listDepth

# @importFrom corpcor pseudoinverse
# @importFrom mvtnorm dmvnorm
# @importFrom rgl plot3d title3d rgl.viewpoint open3d rgl.material triangles3d lines3d spheres3d


#' 
NULL


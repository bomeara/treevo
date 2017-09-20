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
#' breakThisExample

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

# need to import class gpc.poly from gpclib, and setAs from methods
# @importClassesFrom gpclib gpc.poly
# and triangulate
# @importFrom gpclib triangulate 

# interval functions - DO NOT EXPORT
# @export abcparticle
# @export abctaxon
# @export fitContinuous.hacked
# @export interparticleDistance
# @export pairings
# @export pullFromPrior
# @export summarizeTaxonStates
# @export summaryStatsLong
#  @export sumStatNames
#  @export getSimulationSplits
#  @export getTaxonDFWithPossibleExtinction

# exporting, but would be better to do this in each functions ro2 there, not here...
# the following should all be exported without these statements now (09-06-17)

#  @export compareListIPD
#  @export convertTaxonFrameToGeigerData
#  @export credibleInt
#  @export doSimulation
#  @export doSimulationForPlotting
#  @export mutateState
#  @export PairwiseESS
#  @export PairwiseKS
#  @export parallelSimulation
#  @export parentOffspringPlots
#  @export plotPosteriors
#  @export plotPrior
#  @export simulateData
#  @export ThreeD.ABCplots
#  @export getUnivariatePriorCurve
#  @export getUnivariatePosteriorCurve
#  @export plotUnivariatePosteriorVsPrior
#  @export PLSRejection
#  @export returnPLSModel
#  @export PLSTransform
#  @export GetBMRatePrior
#  @export doSimulationWithPossibleExtinction
#  @export createAbcTaxonFromHeightsRow

#'
NULL
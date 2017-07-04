

#' Simulated data
#'
#' Simulated 30-taxon coalescent tree and simulated character data from BM
#' intrinsic model. Character data was generated using doSimulation.
#'
#'
#' @name simData
#' @aliases simData simPhy simChar
#' @docType data
#' @format simPhy in phylo format and simChar in single column matrix
#' @keywords datasets
NULL



#' TreEvo--abc for comparative methods
#'
#' A package for using approximate Bayesian computation for comparative methods
#' in trait evolution
#'
#' \tabular{ll}{ Package: \tab TreEvo\cr Type: \tab Package\cr Version: \tab
#' 0.3.3\cr Date: \tab 2012-07-02\cr License: \tab GPL\cr }
#'
#' @name TreEvo-package
#' @aliases TreEvo-package treevo
#' @docType package
#' @author Brian O'Meara, Barb L. Banbury
#'
#' Maintainer: Barb Banbury <darwinthesun@@gmail.com>
#' @keywords treevo abc
#' @examples
#'
#' #Continuous character simulation under Brownian motion
#' library(ape)
#' phy<-rcoal(20)
#' char<-doSimulationWithPossibleExtinction(
#' 	taxon.df=getTaxonDFWithPossibleExtinction(phy),
#' 	intrinsicFn=brownianIntrinsic,
#' 	extrinsicFn=nullExtrinsic,
#' 	startingValues=c(30),
#' 	intrinsicValues=c(.01),
#' 	extrinsicValues=c(0),
#' 	timeStep=0.001)
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
#' @importFrom rgl plot3d title3d rgl.viewpoint open3d rgl.material triangles3d lines3d spheres3d triangulate 
#' @importFrom partitions blockparts 
#' @importFrom car powerTransform 
#' @importFrom mvtnorm dmvnorm 
#' @importFrom graphics curve layout layout.show legend lines plot plot.new points polygon rect segments symbols text title 
#' @importFrom grDevices dev.new dev.off gray pdf rgb 
#' @importFrom methods as 

# exporting, but would be better to do this in each functions ro2 there, not here...

#' @export nullIntrinsic
#' @export brownianIntrinsic
#' @export boundaryIntrinsic
#' @export boundaryMinIntrinsic
#' @export autoregressiveIntrinsic
#' @export MinBoundaryAutoregressiveIntrinsic
#' @export autoregressiveIntrinsicTimeSlices
#' @export autoregressiveIntrinsicTimeSlicesConstantMean
#' @export autoregressiveIntrinsicTimeSlicesConstantSigma
#' @export nullExtrinsic
#' @export nearestNeighborDisplacementExtrinsic
#' @export ExponentiallyDecayingPush
#' @export everyoneDisplacementExtrinsic
#' @export abcDistance
#' @export abctaxon
#' @export abcparticle
#' @export BCP
#' @export boxcoxEstimation
#' @export boxcoxTransformation
#' @export compareListIPD
#' @export convertTaxonFrameToGeigerData
#' @export CredInt
#' @export doRun_prc
#' @export doRun_rej
#' @export doSimulation
#' @export doSimulationForPlotting
#' @export fitContinuous.hacked
#' @export getSimulationSplits
#' @export HPD
#' @export interparticleDistance
#' @export mutateState
#' @export pairings
#' @export PairwiseESS
#' @export PairwiseKS
#' @export parallelSimulation
#' @export parentOffspringPlots
#' @export plotPosteriors
#' @export plotPrior
#' @export pullFromPrior
#' @export simulateData
#' @export summarizeTaxonStates
#' @export summaryStatsLong
#' @export ThreeD.ABCplots
#' @export sumStatNames
#' @export getUnivariatePriorCurve
#' @export getUnivariatePosteriorCurve
#' @export plotUnivariatePosteriorVsPrior
#' @export PLSRejection
#' @export returnPLSModel
#' @export PLSTransform
#' @export GetBMRatePrior
#' @export doSimulationWithPossibleExtinction
#' @export getTaxonDFWithPossibleExtinction
#' @export createAbcTaxonFromHeightsRow)

#'
NULL
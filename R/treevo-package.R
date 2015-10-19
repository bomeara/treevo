

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
#' char<-doSimulation(
#' 	splits=getSimulationSplits(phy), 
#' 	intrinsicFn=brownianIntrinsic, 
#' 	extrinsicFn=nullExtrinsic, 
#' 	startingValues=c(30), 
#' 	intrinsicValues=c(.01), 
#' 	extrinsicValues=c(0), 
#' 	timeStep=0.001)
#' 
NULL




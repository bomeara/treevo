#' Simulate data for initial TreEvo analysis
#'
#' This function pulls parameters from prior distributions and simulates
#' continuous characters
#'
#' Used by TreEvo doRun function to calculate simulations for doRun_prc and
#' doRun_rej.
#'

# @param taxon.df starting object from getTaxonDFWithPossibleExtinction

# @param phy Tree (Phylogenetic tree in phylo format)

#' @param phy A phylogenetic tree, in package \code{ape}'s \code{phylo} format.

#' @param startingPriorsValues Matrix with ncol=number of states (characters)
#' at root and nrow=2 (two parameters to pass to prior distribution)

#' @param intrinsicPriorsValues Matrix with ncol=number of states (characters)
#' at root and nrow=2 (two parameters to pass to prior distribution)

#' @param extrinsicPriorsValues Matrix with ncol=number of states (characters)
#' at root and nrow=2 (two parameters to pass to prior distribution)

#' @param startingPriorsFns Vector containing names of prior distributions to
#' use for root states: can be one of fixed, uniform, normal, lognormal, gamma,
#' exponential

#' @param intrinsicPriorsFns Vector containing names of prior distributions to
#' use for root states: can be one of fixed, uniform, normal, lognormal, gamma,
#' exponential

#' @param extrinsicPriorsFns Vector containing names of prior distributions to
#' use for root states: can be one of fixed, uniform, normal, lognormal, gamma,
#' exponential

#' @param freevector A vector (length=number of parameters) of free (T) and
#' fixed (F) parameters

#' @param timeStep This value corresponds to the number of discrete time steps

#' @param intrinsicFn Name of intrinsic function characters should be simulated
#' under (as used by doSimulation)

#' @param extrinsicFn Name of extrinsic function characters should be simulated
#' under (as used by doSimulation)

#' @param giveUpAttempts Value for when to stop the analysis if NAs are present

#' @param niter.brown Number of random starts for BM model (min of 2)

#' @param niter.lambda Number of random starts for lambda model (min of 2)

#' @param niter.delta Number of random starts for delta model (min of 2)

#' @param niter.OU Number of random starts for OU model (min of 2)

#' @param niter.white Number of random starts for white model (min of 2)

#' @param verbose If TRUE, chat about how the sim is going

#' @return Returns matrix of trueFreeValues and summary statistics for
#' simulations

#' @author Brian O'Meara and Barb Banbury

# @references O'Meara and Banbury, unpublished

#' @examples
#'
#' #Continuous character simulation under Brownian motion
# library(ape)
#' phy<-rcoal(20)
#' char<-doSimulationWithPossibleExtinction(
#' 	taxon.df=getTaxonDFWithPossibleExtinction(phy),
#' 	intrinsicFn=brownianIntrinsic,
#' 	extrinsicFn=nullExtrinsic,
#' 	startingValues=c(30),
#' 	intrinsicValues=c(.01),
#' 	extrinsicValues=c(0),
#' 	timeStep=0.001
#'	)
#' 

#' @name simulateData
#' @rdname simulateData
#' @export
simulateData<-function(taxon.df, phy, startingPriorsValues, intrinsicPriorsValues, extrinsicPriorsValues, startingPriorsFns, intrinsicPriorsFns, extrinsicPriorsFns, freevector, timeStep, intrinsicFn, extrinsicFn,giveUpAttempts=10, niter.brown=25, niter.lambda=25, niter.delta=25, niter.OU=25, niter.white=25, verbose=FALSE) {

	taxon.df <- getTaxonDFWithPossibleExtinction(phy)

	simTrueAndStats<-rep(NA,10) #no particular reason for it to be 10 wide
	n.attempts<-0
	while (length(which(is.na(simTrueAndStats)))>0) {
	    n.attempts<-n.attempts+1
	    if (n.attempts>giveUpAttempts) {
	    	stop("Error: keep getting NA in the output of simulateData")
	    }
		trueStarting<-rep(NaN, dim(startingPriorsValues)[2])
		trueIntrinsic<-rep(NaN, dim(intrinsicPriorsValues)[2])
		trueExtrinsic<-rep(NaN, dim(extrinsicPriorsValues)[2])
		for (j in 1:dim(startingPriorsValues)[2]) {
			trueStarting[j]=pullFromPrior(startingPriorsValues[,j],startingPriorsFns[j])
		}
		for (j in 1:dim(intrinsicPriorsValues)[2]) {
			trueIntrinsic[j]=pullFromPrior(intrinsicPriorsValues[,j],intrinsicPriorsFns[j])
		}
		for (j in 1:dim(extrinsicPriorsValues)[2]) {
			trueExtrinsic[j]=pullFromPrior(extrinsicPriorsValues[,j],extrinsicPriorsFns[j])
		}
		trueInitial<-c(trueStarting, trueIntrinsic, trueExtrinsic)
		trueFreeValues<-trueInitial[freevector]

		cat(".")
		simTraits<-doSimulationWithPossibleExtinction(taxon.df=taxon.df, intrinsicFn=intrinsicFn, extrinsicFn=extrinsicFn, startingValues=trueStarting, intrinsicValues=trueIntrinsic, extrinsicValues=trueExtrinsic, timeStep=timeStep, verbose=verbose)
		simSumStats<-summaryStatsLong(phy, simTraits, niter.brown, niter.lambda, niter.delta, niter.OU, niter.white)
		simTrueAndStats <-c(trueFreeValues, simSumStats)
	}
	if(n.attempts>1) {
		warning(paste("Had to run simulateData()",n.attempts,"times to get results with no NA. This could bias results if runs with certain parameters failed more often and this happens in many attempted simulations"))
	}
	return(simTrueAndStats)
}

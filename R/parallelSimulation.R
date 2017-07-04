#' Simulate data for initial TreEvo analysis
#'
#' This is a wrapper function for simulateData that allows for multithreading
#' and checkpointing
#'
#'
#' @param nrepSim Number of simulations
#' @param coreLimit Number of cores to be used
#' @param taxon.df object from getTaxonDFWithPossibleExtinction
#' @param phy Tree (Phylogenetic tree in phylo format)
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

#' @param multicore Whether to use multicore, default is FALSE. If TRUE, one of
#' two suggested packages must be installed, either 'doMC' (for UNIX systems) or
#' 'doParallel' (for Windows), which are used to activate multithreading.
#' If neither package is installed, this function will fail if multicore=TRUE.

#' @param checkpointFile Optional file name for checkpointing simulations
#' @param checkpointFreq Saving frequency for checkpointing
#' @param niter.brown Number of random starts for BM model (min of 2)
#' @param niter.lambda Number of random starts for lambda model (min of 2)
#' @param niter.delta Number of random starts for delta model (min of 2)
#' @param niter.OU Number of random starts for OU model (min of 2)
#' @param niter.white Number of random starts for white model (min of 2)
#' @return Returns matrix of trueFreeValues and summary statistics for
#' simulations
#' @author Brian O'Meara and Barb Banbury
#' @references O'Meara and Banbury, unpublished

parallelSimulation<-function(nrepSim, coreLimit, taxon.df, phy, startingPriorsValues, intrinsicPriorsValues, extrinsicPriorsValues, startingPriorsFns, intrinsicPriorsFns, extrinsicPriorsFns, freevector, timeStep, intrinsicFn, extrinsicFn, multicore,checkpointFile=NULL,checkpointFreq=24, niter.brown=25, niter.lambda=25, niter.delta=25, niter.OU=25, niter.white=25) {
	#library(doMC, quietly=T)
	#library(foreach, quietly=T)
	cores=1
	if (multicore) {
		if (is.na(coreLimit)){
			registerMulticoreEnv()
			getDoParWorkers()->cores
		}else{
			registerMulticoreEnv(coreLimit)
			coreLimit->cores
			}
		}
	cat(paste("Using", cores, "core(s) for simulations \n\n"))
	if (nrepSim %%cores != 0) {
		warning("The simulation is most efficient if the number of nrepSim is a multiple of the number of cores")
	}

	cat("Doing simulations: ")
	if (is.null(checkpointFile)) {
		trueFreeValuesANDSummaryValues<-foreach(1:nrepSim, .combine=rbind) %dopar% simulateData(taxon.df, phy, startingPriorsValues, intrinsicPriorsValues, extrinsicPriorsValues, startingPriorsFns, intrinsicPriorsFns, extrinsicPriorsFns, freevector, timeStep, intrinsicFn, extrinsicFn)
	}
	else {
		checkpointFileName<-paste(checkpointFile,".trueFreeValuesANDSummaryValues.Rsave",sep="")
		trueFreeValuesANDSummaryValues<-c()
		checkpointFreqAdjusted<-max(cores*round(checkpointFreq/cores),1)
		numberSimsInCheckpointRuns<-checkpointFreqAdjusted * floor(nrepSim/checkpointFreqAdjusted)
		numberLoops<-floor(numberSimsInCheckpointRuns/checkpointFreqAdjusted)
		numberSimsPerLoop<-numberSimsInCheckpointRuns/numberLoops
		numberSimsAfterLastCheckpoint<-nrepSim - numberSimsInCheckpointRuns
		if (checkpointFreqAdjusted != checkpointFreq ) {
			warning(paste("Checkpoint frequency adjusted from",checkpointFreq,"to",checkpointFreqAdjusted,"to reduce the wasted time on unused cores"))
		}
		for (rep in sequence(numberLoops)) {
			trueFreeValuesANDSummaryValues<-rbind(trueFreeValuesANDSummaryValues, foreach(1:numberSimsPerLoop, .combine=rbind) %dopar% simulateData(taxon.df, phy, startingPriorsValues, intrinsicPriorsValues, extrinsicPriorsValues, startingPriorsFns, intrinsicPriorsFns, extrinsicPriorsFns, freevector, timeStep, intrinsicFn, extrinsicFn))
			save(trueFreeValuesANDSummaryValues,file=checkpointFileName)
			print(paste("Just finished",dim(trueFreeValuesANDSummaryValues)[1],"of",nrepSim,"simulations; progress so far saved in",checkpointFileName))
		}
		trueFreeValuesANDSummaryValues<-rbind(trueFreeValuesANDSummaryValues, foreach(1:numberSimsAfterLastCheckpoint, .combine=rbind) %dopar% simulateData(taxon.df, phy, startingPriorsValues, intrinsicPriorsValues, extrinsicPriorsValues, startingPriorsFns, intrinsicPriorsFns, extrinsicPriorsFns, freevector, timeStep, intrinsicFn, extrinsicFn))
	}
	return(trueFreeValuesANDSummaryValues)
}

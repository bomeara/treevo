#' Simulate data for initial TreEvo analysis
#'
#' The \code{simulateWithPriors} function pulls parameters from prior distributions and conducts a single simulation of
#' continuous trait evolution (using \code{\link{doSimulation}} functions), returning useful summary statistics for ABC.
#' \code{parallelSimulateWithPriors} is a wrapper function for \code{simulateWithPriors} that allows for multithreading
#' and checkpointing. This family of functions is mostly used as internal components, generating simulations
#' within ABC analyses using the \code{\link{doRun}} functions. See \emph{Note} below.
#'

#' @note 
#' The \code{\link{simulateWithPriors}} functions are effectively the engine that powers the \code{\link{doRun}}
#' functions, while the \code{\link{doSimulation}} functions are the pistons within the \code{\link{simulateWithPriors}} engine.
#' In general, most users will just drive the car - they will just use \code{\link{doRun}}, but some users may
#' want to use \code{\link{simulateWithPriors}} or \code{\link{doSimulation}} functions to do various simulations.

#' @inheritParams doSimulation


#' @param startingPriorsValues Matrix with number of columns equal to the number of states (characters)
#' at root and number of rows equal to two (representing two parameters to pass to prior distribution).

#' @param intrinsicPriorsValues Matrix with number of columns equal to the number of parameters to pass
#' to the intrinsic function and nrow=2 (two parameters to pass to prior distribution)

#' @param extrinsicPriorsValues Matrix with number of columns equal to the number of parameters to pass
#' to the extrinsic function and nrow=2 (two parameters to pass to prior distribution)


#' @param startingPriorsFns Vector containing names of prior distributions to
#' use for root states: can be one of \code{"fixed"}, \code{"uniform"}, \code{"normal"},
#" \code{"lognormal"}, \code{"gamma"}, \code{"exponential"}.

#' @param intrinsicPriorsFns Vector containing names of prior distributions to
#' use for intrinsic function parameters: can be one of \code{"fixed"}, \code{"uniform"}, \code{"normal"},
#" \code{"lognormal"}, \code{"gamma"}, \code{"exponential"}.

#' @param extrinsicPriorsFns Vector containing names of prior distributions to
#' use for extrinsic function parameters: can be one of \code{"fixed"}, \code{"uniform"}, \code{"normal"},
#" \code{"lognormal"}, \code{"gamma"}, \code{"exponential"}.

#' @param freevector A logical vector (with length equal to the number of parameters), indicating free (\code{TRUE}) and
#' fixed (\code{FALSE}) parameters.

#' @param giveUpAttempts Value for when to stop the analysis if \code{NA} values are present.

#' @param nrepSim Number of replicated simulations to run.

#' @param coreLimit Maximum number of cores to be used.

#' @param multicore Whether to use multicore, default is \code{FALSE}. If \code{TRUE}, one of
#' two suggested packages must be installed, either \code{doMC} (for UNIX systems) or
#' \code{doParallel} (for Windows), which are used to activate multithreading.
#' If neither package is installed, this function will fail if \code{multicore=TRUE}.

#' @param niter.brown Number of random starts for the Brownian Motion (BM) model (minimum of 2).

#' @param niter.lambda Number of random starts for the lambda model (minimum of 2).

#' @param niter.delta Number of random starts for the delta model (minimum of 2).

#' @param niter.OU Number of random starts for the Ornstein-Uhlenbeck (OU) model (minimum of 2).

#' @param niter.white Number of random starts for the white noise model (minimum of 2).


#' @param checks If \code{TRUE}, checks inputs for consistency. This activity is skipped (\code{checks = FALSE})
#' when run in parallel by \code{parallelSimulateWithPriors}, and instead is only checked once. This
#' argument also controls whether \code{simulateWithPriors} assigns \code{freevector} as an attribute to the
#' output produced.

#' @param checkpointFile Optional file name for checkpointing simulations

#' @param checkpointFreq Saving frequency for checkpointing


#' @return Function \code{simulateWithPriors} returns a vector of \code{trueFreeValues},
#' the true generating parameters used in the simulation
#' (a set of values as long as the number of freely varying parameters), concatenated with a set of summary statistics for
#' the simulation. Function \code{parallelSimulateWithPriors} returns a matrix of such vectors bound
#' together, with each row representing a different simulation. By default,
#' both functions also assign a logical vector named \code{freevector}, indicating the total number of 
#' parameters and which parameters are freely-varying (have \code{TRUE} values), as an attribute of
#' the output.

#' @author Brian O'Meara and Barb Banbury

# @references O'Meara and Banbury, unpublished

#' @examples
#'
#' simPhy<-rcoal(30)
#' 
#' # example simulation
#' 	
#' simData<-simulateWithPriors(phy=simPhy, 
#'   intrinsicFn=brownianIntrinsic,
#'   extrinsicFn=nullExtrinsic,
#'   startingPriorsFns="normal",
#'   startingPriorsValues=matrix(c(mean(simChar[,1]), sd(simChar[,1]))),
#'   intrinsicPriorsFns=c("exponential"),
#'   intrinsicPriorsValues=matrix(c(10, 10), nrow=2, byrow=FALSE),
#'   extrinsicPriorsFns=c("fixed"),
#'   extrinsicPriorsValues=matrix(c(0, 0), nrow=2, byrow=FALSE),
#' 	 timeStep=0.0001,
#' 	 freevector=NULL, 	
#' 	 giveUpAttempts=10, 
#' 	 verbose=TRUE,
#' 	 niter.brown=25, niter.lambda=25, niter.delta=25, niter.OU=25, niter.white=25) 
#' 
#' simData
#' 
#' simDataParallel<-parallelSimulateWithPriors( 
#'   nrepSim=10, multicore=FALSE, coreLimit=1, 
#'   phy=simPhy,
#'   intrinsicFn=brownianIntrinsic,
#'   extrinsicFn=nullExtrinsic,
#'   startingPriorsFns="normal",
#'   startingPriorsValues=matrix(c(mean(simChar[,1]), sd(simChar[,1]))),
#'   intrinsicPriorsFns=c("exponential"),
#'   intrinsicPriorsValues=matrix(c(10, 10), nrow=2, byrow=FALSE),
#'   extrinsicPriorsFns=c("fixed"),
#'   extrinsicPriorsValues=matrix(c(0, 0), nrow=2, byrow=FALSE), 
#'   timeStep=0.0001,
#'   checkpointFile=NULL, checkpointFreq=24,
#'   verbose=FALSE,
#'   freevector=NULL, taxon.df=NULL,
#'   niter.brown=25, niter.lambda=25, niter.delta=25, niter.OU=25, niter.white=25) 
#' 
#' simDataParallel
#' 


#' @name simulateWithPriors
#' @rdname simulateWithPriors
#' @export
simulateWithPriors<-function(
	phy=NULL, intrinsicFn, extrinsicFn, startingPriorsFns, startingPriorsValues, 
	intrinsicPriorsFns, intrinsicPriorsValues, extrinsicPriorsFns, extrinsicPriorsValues, timeStep,
	giveUpAttempts=10, verbose=FALSE, checks=TRUE, taxon.df=NULL, freevector=NULL, 
	niter.brown=25, niter.lambda=25, niter.delta=25, niter.OU=25, niter.white=25) {

	if(is.null(taxon.df)){
		taxon.df <- getTaxonDFWithPossibleExtinction(phy)
		}
	
	
	# checks
	if(checks){
		checkNiter(niter.brown=niter.brown, niter.lambda=niter.lambda,
			niter.delta=niter.delta, niter.OU=niter.OU, niter.white=niter.white)
		}
		
	if(is.null(freevector)){
		freevector<-getFreeVector(startingPriorsFns=startingPriorsFns, startingPriorsValues=startingPriorsValues, 
					intrinsicPriorsFns=intrinsicPriorsFns, intrinsicPriorsValues=intrinsicPriorsValues,
					extrinsicPriorsFns=extrinsicPriorsFns, extrinsicPriorsValues=extrinsicPriorsValues)
		}
			
	simTrueAndStats<-rep(NA,10) #no particular reason for it to be 10 wide
	n.attempts<-0
	while (length(which(is.na(simTrueAndStats)))>0) {
	    n.attempts<-n.attempts+1
	    if (n.attempts>giveUpAttempts) {
	    	stop("Error: keep getting NA in the output of simulateWithPriors")
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
		simTraits<-doSimulationWithPossibleExtinction(phy=phy, taxon.df=taxon.df, intrinsicFn=intrinsicFn, extrinsicFn=extrinsicFn, 
			startingValues=trueStarting, intrinsicValues=trueIntrinsic, extrinsicValues=trueExtrinsic, timeStep=timeStep, verbose=verbose)
		simSumStats<-summaryStatsLong(phy=phy, traits=simTraits, 
			niter.brown=niter.brown, niter.lambda=niter.lambda, niter.delta=niter.delta,
			niter.OU=niter.OU, niter.white=niter.white)
		simTrueAndStats <-c(trueFreeValues, simSumStats)
	}
	if(n.attempts>1) {
		warning(paste("Had to run simulateWithPriors()",n.attempts,"times to get results with no NA. This could bias results if runs with certain parameters failed more often and this happens in many attempted simulations"))
	}
	if(!checks){
		attr(trueFreeValuesANDSummaryValues,"freevector")<-freevector
		}
	return(simTrueAndStats)
}




#' @rdname simulateWithPriors
#' @export
parallelSimulateWithPriors<-function(
	nrepSim, multicore, coreLimit,  
	phy, 
	intrinsicFn, extrinsicFn, startingPriorsFns, startingPriorsValues,
	intrinsicPriorsFns, intrinsicPriorsValues, extrinsicPriorsFns, extrinsicPriorsValues, timeStep,
	checkpointFile=NULL, checkpointFreq=24, verbose=FALSE, freevector=NULL, taxon.df=NULL, giveUpAttempts=10, 
	niter.brown=25, niter.lambda=25, niter.delta=25, niter.OU=25, niter.white=25) {
	
	#library(doMC, quietly=T)
	#library(foreach, quietly=T)
	
	# checks
	checkNiter(niter.brown=niter.brown, niter.lambda=niter.lambda,
		niter.delta=niter.delta, niter.OU=niter.OU, niter.white=niter.white)

	if(is.null(freevector)){
		freevector<-getFreeVector(startingPriorsFns=startingPriorsFns, startingPriorsValues=startingPriorsValues, 
					intrinsicPriorsFns=intrinsicPriorsFns, intrinsicPriorsValues=intrinsicPriorsValues,
					extrinsicPriorsFns=extrinsicPriorsFns, extrinsicPriorsValues=extrinsicPriorsValues)
		}
	
	if(is.null(taxon.df)){
		taxon.df <- getTaxonDFWithPossibleExtinction(phy)
		}

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
		trueFreeValuesANDSummaryValues<-foreach(1:nrepSim, .combine=rbind) %dopar% simulateWithPriors(phy=phy, taxon.df=taxon.df,
			startingPriorsValues=startingPriorsValues, intrinsicPriorsValues=intrinsicPriorsValues, extrinsicPriorsValues=extrinsicPriorsValues,
			startingPriorsFns=startingPriorsFns, intrinsicPriorsFns=intrinsicPriorsFns, extrinsicPriorsFns=extrinsicPriorsFns, giveUpAttempts=giveUpAttempts, 
			freevector=freevector, timeStep=timeStep, intrinsicFn=intrinsicFn, extrinsicFn=extrinsicFn, verbose=verbose, checks=FALSE)
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
			trueFreeValuesANDSummaryValues<-rbind(trueFreeValuesANDSummaryValues, 
				foreach(1:numberSimsPerLoop, .combine=rbind) %dopar% simulateWithPriors(phy=phy,  taxon.df=taxon.df,
					startingPriorsValues=startingPriorsValues, intrinsicPriorsValues=intrinsicPriorsValues, extrinsicPriorsValues=extrinsicPriorsValues,
					startingPriorsFns=startingPriorsFns, intrinsicPriorsFns=intrinsicPriorsFns, extrinsicPriorsFns=extrinsicPriorsFns,  giveUpAttempts=giveUpAttempts, 
					freevector=freevector, timeStep=timeStep, intrinsicFn=intrinsicFn, extrinsicFn=extrinsicFn, verbose=verbose, checks=FALSE))
			save(trueFreeValuesANDSummaryValues,file=checkpointFileName)
			print(paste("Just finished",dim(trueFreeValuesANDSummaryValues)[1],"of",nrepSim,"simulations; progress so far saved in",checkpointFileName))
		}
		trueFreeValuesANDSummaryValues<-rbind(trueFreeValuesANDSummaryValues, 
			foreach(1:numberSimsAfterLastCheckpoint, .combine=rbind) %dopar% simulateWithPriors(phy=phy,  taxon.df=taxon.df,
				startingPriorsValues=startingPriorsValues, intrinsicPriorsValues=intrinsicPriorsValues, extrinsicPriorsValues=extrinsicPriorsValues,
				startingPriorsFns=startingPriorsFns, intrinsicPriorsFns=intrinsicPriorsFns, extrinsicPriorsFns=extrinsicPriorsFns,  giveUpAttempts=giveUpAttempts, 
				freevector=freevector, timeStep=timeStep, intrinsicFn=intrinsicFn, extrinsicFn=extrinsicFn, verbose=verbose, checks=FALSE))
	}
	attr(trueFreeValuesANDSummaryValues,"freevector")<-freevector
	return(trueFreeValuesANDSummaryValues)
}

checkNiter<-function(niter.brown=25, niter.lambda=25, niter.delta=25, niter.OU=25, niter.white=25){
	if(niter.brown<2){stop("niter.brown must be at least 2")}
	if(niter.lambda<2){stop("niter.lambda must be at least 2")}
	if(niter.delta<2){stop("niter.delta must be at least 2")}
	if(niter.OU<2){stop("niter.OU must be at least 2")}
	if(niter.white<2){stop("niter.white must be at least 2")}
	}
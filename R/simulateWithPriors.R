#' Simulate data for initial TreEvo analysis
#' 
#' The \code{simulateWithPriors} function pulls parameters from prior
#' distributions and conducts a single simulation of continuous trait evolution
#' (using the \code{\link{doSimulation}} function),
#' returning useful summary statistics for ABC.
#' \code{parallelSimulateWithPriors} is a wrapper function for
#' \code{simulateWithPriors} that allows for multithreading and
#' checkpointing. This family of functions is mostly used as
#' internal components, generating simulations within ABC analyses
#' using the \code{\link{doRun}} functions. See \emph{Note} below.
#' 

#' @note
#' The \code{\link{simulateWithPriors}} functions are effectively the
#' engine that powers the  \code{\link{doRun}} functions, while the
#' \code{\link{doSimulation}} function is the pistons within the
#' \code{\link{simulateWithPriors}} engine. In general, most users
#' will just drive the car - they will just use \code{\link{doRun}},
#' but some users may want to use \code{\link{simulateWithPriors}}
#' or \code{\link{doSimulation}} to do various simulations.

#' @inheritParams doSimulation

#' @param startingPriorsValues A list of the same length as
#' the number of prior distributions specified in
#' \code{startingPriorsFns} (for starting values, this should
#' be one prior function specified for each trait
#' - thus one for most univariate trait analyses), with each
#' element of the list a vector the same length
#' as the appropriate number of parameters for that prior distribution
#' (1 for \code{"fixed"}, 2 for \code{"uniform"}, 2 for \code{"normal"},
#' 2 for \code{"lognormal"}, 2 for \code{"gamma"}, 1 for \code{"exponential"}).

#' @param intrinsicPriorsValues A list of the same length
#' as the number of prior distributions specified
#' in \code{intrinsicPriorsFns} (one prior function specified
#' for each parameter in the intrinsic model),
#' with each element of the list a vector the same length#
#' as the appropriate number of parameters for that prior distribution
#' (1 for \code{"fixed"}, 2 for \code{"uniform"}, 2 for \code{"normal"},
#' 2 for \code{"lognormal"}, 2 for \code{"gamma"}, 1 for \code{"exponential"}).

#' @param extrinsicPriorsValues A list of the same length as
#' the number of prior distributions specified in \code{extrinsicPriorsFns} 
#'(one prior function specified for each parameter in the extrinsic model),
#' with each element of the list a vector the same length as the appropriate
#' number of parameters for that prior distribution
#' (1 for \code{"fixed"}, 2 for \code{"uniform"}, 2 for \code{"normal"},
#' 2 for \code{"lognormal"}, 2 for \code{"gamma"}, 1 for \code{"exponential"}).

#' @param startingPriorsFns Vector containing names of prior distributions to
#' use for root states: can be one of
#' \code{"fixed"}, \code{"uniform"}, \code{"normal"}, 
#' \code{"lognormal"}, \code{"gamma"}, \code{"exponential"}.

#' @param intrinsicPriorsFns Vector containing names of prior distributions to
#' use for intrinsic function parameters: can be one of
#' \code{"fixed"}, \code{"uniform"}, \code{"normal"}, 
#' \code{"lognormal"}, \code{"gamma"}, \code{"exponential"}.

#' @param extrinsicPriorsFns Vector containing names of prior distributions to
#' use for extrinsic function parameters: can be one of
#' \code{"fixed"}, \code{"uniform"}, \code{"normal"}, 
#' \code{"lognormal"}, \code{"gamma"}, \code{"exponential"}.

#' @param freevector A logical vector
#' (with length equal to the number of parameters),
#' indicating free (\code{TRUE}) and
#' fixed (\code{FALSE}) parameters.

#' @param giveUpAttempts Value for when to stop the
#' analysis if \code{NA} values are present.

#' @param nrepSim Number of replicated simulations to run.

#' @param coreLimit Maximum number of cores to be used.

#' @param multicore Whether to use multicore, default is \code{FALSE}.
#' If \code{TRUE}, one of two suggested packages must be installed,
#' either \code{doMC} (for UNIX systems) or
#' \code{doParallel} (for Windows), which are used to activate multithreading.
#' If neither package is installed, this function will fail if \code{multicore = TRUE}.

# @param niter.brown Number of random starts for the Brownian Motion (BM) model (minimum of 2).
# @param niter.lambda Number of random starts for the lambda model (minimum of 2).
# @param niter.delta Number of random starts for the delta model (minimum of 2).
# @param niter.OU Number of random starts for the Ornstein-Uhlenbeck (OU) model (minimum of 2).
# @param niter.white Number of random starts for the white noise model (minimum of 2).

#' @param checks If \code{TRUE}, checks inputs for consistency.
#' This activity is skipped (\code{checks = FALSE})
#' when run in parallel by \code{parallelSimulateWithPriors},
#' and instead is only checked once. This
#' argument also controls whether \code{simulateWithPriors}
#' assigns \code{freevector} as an attribute to the
#' output produced.

#' @param checkpointFile Optional file name for checkpointing simulations

#' @param checkpointFreq Saving frequency for checkpointing

#' @param verbose If \code{TRUE}, gives messages about how
#' the simulation is progessing via \code{message}.

#' @param verboseNested Should looped runs of \code{simulateWithPriors} be verbose?

#' @return Function \code{simulateWithPriors} returns a
#' vector of \code{trueFreeValues}, the true generating parameters
#' used in the simulation (a set of values as long as
#' the number of freely varying parameters), concatenated with
#' a set of summary statistics for the simulation.
#' 
#' Function \code{parallelSimulateWithPriors} returns a matrix
#' of such vectors bound together, with each row representing
#' a different simulation. 
#' 
#' By default, both functions also assign a logical vector named
#' \code{freevector}, indicating the total number of parameters and
#' which parameters are freely-varying (have \code{TRUE} values),
#' as an attribute of the output.

#' @author Brian O'Meara and Barb Banbury

# @references O'Meara and Banbury, unpublished

#' @examples
#
#' \donttest{
#
#' set.seed(1)
#' tree <- rcoal(20)
#' # get realistic edge lengths
#' tree$edge.length <- tree$edge.length*20
#' 
#' # example simulation
#' 
#' # NOTE: the example analyses involve too few simulations,
#'     # as well as overly coarse time-units...
#'     # ...all for the sake of examples that reasonably test the functions
#'     
#' simData <- simulateWithPriors(
#'   phy = tree, 
#'   intrinsicFn = brownianIntrinsic, 
#'   extrinsicFn = nullExtrinsic, 
#'   startingPriorsFns = "normal", 
#'   startingPriorsValues = list(
#'       c(mean(simChar[, 1]), sd(simChar[, 1]))), 
#'   intrinsicPriorsFns = c("exponential"), 
#'   intrinsicPriorsValues = list(10), 
#'   extrinsicPriorsFns = c("fixed"), 
#'   extrinsicPriorsValues = list(0), 
#'   generation.time = 100000, 
#'   freevector = NULL,     
#'   giveUpAttempts = 10, 
#'   verbose = TRUE)
#' 
#' simData
#' 
#' simDataParallel <- parallelSimulateWithPriors(
#'   nrepSim = 2, 
#'   multicore = FALSE, 
#'   coreLimit = 1, 
#'   phy = tree, 
#'   intrinsicFn = brownianIntrinsic, 
#'   extrinsicFn = nullExtrinsic, 
#'   startingPriorsFns = "normal", 
#'   startingPriorsValues = list(
#'      c(mean(simChar[, 1]), sd(simChar[, 1]))), 
#'   intrinsicPriorsFns = c("exponential"), 
#'   intrinsicPriorsValues = list(10), 
#'   extrinsicPriorsFns = c("fixed"), 
#'   extrinsicPriorsValues = list(0), 
#'   generation.time = 100000, 
#'   checkpointFile = NULL, checkpointFreq = 24, 
#'   verbose = TRUE, 
#'   freevector = NULL, 
#'   taxonDF = NULL)
#' 
#' simDataParallel
#' 
#' }
#' 



#' @name simulateWithPriors
#' @rdname simulateWithPriors
#' @export
simulateWithPriors <- function(
		phy = NULL, 
		intrinsicFn, 
		extrinsicFn, 
		startingPriorsFns, 
		startingPriorsValues, 
		intrinsicPriorsFns, 
		intrinsicPriorsValues, 
		extrinsicPriorsFns, 
		extrinsicPriorsValues, 
		generation.time = 1000, 
		TreeYears = max(branching.times(phy)) * 1e6, 
		timeStep = NULL, 
		giveUpAttempts = 10, 
		verbose = FALSE, 
		checks = TRUE, 
		taxonDF = NULL, 
		freevector = NULL
		#, niter.brown = 25, niter.lambda = 25, 
		# niter.delta = 25, niter.OU = 25, niter.white = 25
		) {
	########################################
	#
    if(is.null(taxonDF)){
        taxonDF <- getTaxonDFWithPossibleExtinction(phy)
        }
        
    if(is.null(timeStep)){
        timeStep <- generation.time/TreeYears
        }
    
    # checks
    if(checks){
        #checkNiter(
		#	niter.brown = niter.brown, niter.lambda = niter.lambda, 
        #	niter.delta = niter.delta, niter.OU = niter.OU, niter.white = niter.white
		#	)
		#
        # check TimeStep
        numberofsteps <- max(taxonDF$endTime)/timeStep
        mininterval <- min(taxonDF$endTime - taxonDF$startTime)
        #
        if (floor(mininterval/timeStep)<50 & floor(mininterval/timeStep) >= 3) {
            warning(paste0(
				"You have only ", floor(mininterval/timeStep), 
                " timeSteps on the shortest branch in this dataset\n",
				"  but should probably have a lot more if you expect change on this branch.\n",
				"Please consider decreasing timeStep to no more than ", 
                signif(mininterval/50, 2)
				))
            }
        if (floor(mininterval/timeStep)<3) {
            warning(paste0(
				"You have only ", floor(mininterval/timeStep), 
                " timeSteps on the shortest branch in this dataset\n",
				"  but should probably have a lot more if you expect change on this branch.\n",
				"Please consider decreasing timeStep to no more than ", 
                signif(mininterval/50, 2), " or at the very least ", 
				signif(mininterval/3, 2)))
            #    timeStep <- mininterval/3
            }
        }
        
    if(is.null(freevector)){
        freevector <- getFreeVector(
			startingPriorsFns = startingPriorsFns, 
			startingPriorsValues = startingPriorsValues, 
			intrinsicPriorsFns = intrinsicPriorsFns, 
			intrinsicPriorsValues = intrinsicPriorsValues, 
			extrinsicPriorsFns = extrinsicPriorsFns, 
			extrinsicPriorsValues = extrinsicPriorsValues
			)
        }
            
	#no particular reason for it to be 10 wide
    simTrueAndStats <- rep(NA, 10) 
    n.attempts <- 0
    while (length(which(is.na(simTrueAndStats)))>0) {
        n.attempts <- n.attempts+1
        if (n.attempts>giveUpAttempts) {
            stop(
				"Error: keep getting NA in the output of simulateWithPriors"
				)
        }
        trueStarting <- rep(NaN, length(startingPriorsValues))
        trueIntrinsic <- rep(NaN, length(intrinsicPriorsValues))
        trueExtrinsic <- rep(NaN, length(extrinsicPriorsValues))
        for (j in 1:length(startingPriorsValues)) {
            trueStarting[j] <- pullFromPrior(
				startingPriorsValues[[j]], startingPriorsFns[j]
				)
			}
        for (j in 1:length(intrinsicPriorsValues)) {
            trueIntrinsic[j] <- pullFromPrior(
				intrinsicPriorsValues[[j]], intrinsicPriorsFns[j]
				)
			}
        for (j in 1:length(extrinsicPriorsValues)) {
            trueExtrinsic[j] <- pullFromPrior(
				extrinsicPriorsValues[[j]], extrinsicPriorsFns[j]
				)
			}
        trueInitial <- c(trueStarting, trueIntrinsic, trueExtrinsic)
        trueFreeValues <- trueInitial[freevector]
        #
        #message(".")
        #
        simTraits <- doSimulationInternal(
            taxonDF = taxonDF, 
			timeStep = timeStep, 
            intrinsicFn = intrinsicFn,
			extrinsicFn = extrinsicFn, 
            startingValues = trueStarting, 
			intrinsicValues = trueIntrinsic, 
			extrinsicValues = trueExtrinsic
			)
		#
        simSumStats <- summaryStatsLong(
			phy = phy, traits = simTraits
            #, niter.brown = niter.brown, niter.lambda = niter.lambda,
			# niter.delta = niter.delta, 
            # niter.OU = niter.OU, niter.white = niter.white
            )
        simTrueAndStats  <- c(trueFreeValues, simSumStats)
        }
    if(n.attempts>1) {
        warning(
            paste0(
				"Had to run simulateWithPriors ", n.attempts, 
				" times to get results with no NA.\n",
				"This could bias results if runs with certain parameters failed more often,\n",
				"  and is repeated in many attempted simulations."
            )
        )
    }
    if(checks){
        attr(simTrueAndStats, "freevector") <- freevector
        }
    return(simTrueAndStats)
	}




#' @rdname simulateWithPriors
#' @export
parallelSimulateWithPriors <- function(
    nrepSim, 
	multicore, 
	coreLimit, 
	phy, 
    intrinsicFn, 
	extrinsicFn, 
	startingPriorsFns, 
	startingPriorsValues, 
    intrinsicPriorsFns, 
	intrinsicPriorsValues, 
	extrinsicPriorsFns, 
	extrinsicPriorsValues, 
    generation.time = 1000, 
	TreeYears = max(branching.times(phy)) * 1e6, 
	timeStep = NULL, 
	#timeStep = 1e-04, 
    checkpointFile = NULL, 
	checkpointFreq = 24, 
	verbose = TRUE, 
	checkTimeStep = TRUE, 
    verboseNested = FALSE, 
	freevector = NULL, 
	taxonDF = NULL, 
	giveUpAttempts = 10
    #, niter.brown = 25, niter.lambda = 25, 
	# niter.delta = 25, niter.OU = 25, niter.white = 25
    ) {
    #
    #library(doMC, quietly = TRUE)
    #library(foreach, quietly = TRUE)
    #
    # #get values for arguments that can be NULL if not provided
    #
    if(is.null(timeStep)){
        timeStep <- generation.time/TreeYears
        }
    #
    if(is.null(freevector)){
        freevector <- getFreeVector(
			startingPriorsFns = startingPriorsFns, 
			startingPriorsValues = startingPriorsValues, 
			intrinsicPriorsFns = intrinsicPriorsFns, 
			intrinsicPriorsValues = intrinsicPriorsValues, 
			extrinsicPriorsFns = extrinsicPriorsFns, 
			extrinsicPriorsValues = extrinsicPriorsValues
			)
        }
    #
    if(is.null(taxonDF)){
        taxonDF <- getTaxonDFWithPossibleExtinction(phy)
        }
    #
    #################################################
    # check TimeStep
    if(checkTimeStep){
        numberofsteps <- max(taxonDF$endTime)/timeStep
        mininterval <- min(taxonDF$endTime - taxonDF$startTime)
        #
        if (floor(mininterval/timeStep)<50 & floor(mininterval/timeStep) >= 3) {
            message(paste0(
				"You have only ", floor(mininterval/timeStep), 
				" timeSteps on the shortest branch in this dataset\n",
				"  but should probably have a lot more if you expect change on this branch.\n",
				"Please consider decreasing timeStep to no more than ", 
				signif(mininterval/50, 2)))
        }
        if (floor(mininterval/timeStep)<3) {
            message(paste0(
				"You have only ", floor(mininterval/timeStep), 
				" timeSteps on the shortest branch in this dataset\n",
				"  but should probably have a lot more if you expect change on this branch.\n",
				"Please consider decreasing timeStep to no more than ", 
				signif(mininterval/50, 2), " or at the very least ", 
				signif(mininterval/3, 2)
				# timeStep <- mininterval/3
				))
            }
        }
    #    
    # multicore
    # set up multicore
    cluster <- setupMulticore(
		multicore, 
		nSim = nrepSim, 
		coreLimit = coreLimit
		)
    #
    # verbosity
    nCores <- attr(cluster, "nCores")
    if(verbose){
        message(paste("Using", nCores, "core(s) for simulations \n\n"))
        if (nrepSim %%nCores  !=  0) {
            warning(paste0(
				"The simulation is most efficient if the number of nrepSim\n",
				" is a multiple of the number of nCores"
				))
            }
        message("Doing simulations: ")
        }
    #
    # now run simulations
    if (is.null(checkpointFile)) {
        # no checkpointFile to generate!
        repSimFE <- foreach(1:nrepSim, .combine = rbind)
        trueFreeValuesANDSummaryValues <- (    #makeQuiet(
            repSimFE %dopar% simulateWithPriors(
                phy = phy, 
				taxonDF = taxonDF, 
                startingPriorsValues = startingPriorsValues, 
				intrinsicPriorsValues = intrinsicPriorsValues, 
				extrinsicPriorsValues = extrinsicPriorsValues, 
                startingPriorsFns = startingPriorsFns, 
				intrinsicPriorsFns = intrinsicPriorsFns, 
				extrinsicPriorsFns = extrinsicPriorsFns, 
                giveUpAttempts = giveUpAttempts, 
				freevector = freevector, 
				timeStep = timeStep, 
                intrinsicFn = intrinsicFn, 
				extrinsicFn = extrinsicFn, 
				verbose = verboseNested, 
				checks = FALSE
				)
            #)
            )
    }else{
        checkpointFileName <- paste(checkpointFile, 
			".trueFreeValuesANDSummaryValues.Rsave", sep = "")
        trueFreeValuesANDSummaryValues <- c()
        checkpointFreqAdjusted <- max(nCores*round(checkpointFreq/nCores), 1)
        numberSimsInCheckpointRuns <- checkpointFreqAdjusted * floor(
			nrepSim/checkpointFreqAdjusted
			)
        numberLoops <- floor(numberSimsInCheckpointRuns/checkpointFreqAdjusted)
        numberSimsPerLoop <- numberSimsInCheckpointRuns/numberLoops
        numberSimsAfterLastCheckpoint <- nrepSim - numberSimsInCheckpointRuns
        if (checkpointFreqAdjusted  !=  checkpointFreq ) {
            warning(paste(
				"Checkpoint frequency adjusted from", 
				checkpointFreq, "to", 
				checkpointFreqAdjusted, 
				"to reduce the wasted time on unused nCores"))
            }
        for (rep in sequence(numberLoops)) {
            trueFreeValuesANDSummaryValues <- rbind(trueFreeValuesANDSummaryValues, 
                foreach(1:numberSimsPerLoop, .combine = rbind) %dopar% simulateWithPriors(
					phy = phy, 
					taxonDF = taxonDF, 
                    startingPriorsValues = startingPriorsValues, 
					intrinsicPriorsValues = intrinsicPriorsValues, 
					extrinsicPriorsValues = extrinsicPriorsValues, 
                    startingPriorsFns = startingPriorsFns, 
					intrinsicPriorsFns = intrinsicPriorsFns, 
					extrinsicPriorsFns = extrinsicPriorsFns, 
					giveUpAttempts = giveUpAttempts, 
                    freevector = freevector, 
					timeStep = timeStep, 
					intrinsicFn = intrinsicFn, 
					extrinsicFn = extrinsicFn, 
					verbose = verboseNested, 
					checks = FALSE
					)
				)
            save(
				trueFreeValuesANDSummaryValues, 
				file = checkpointFileName
				)
            message(paste(
				"Just finished", 
				dim(trueFreeValuesANDSummaryValues)[1], "of", 
				nrepSim, "simulations; progress so far saved in", 
				checkpointFileName
				))
            }
        trueFreeValuesANDSummaryValues <- rbind(trueFreeValuesANDSummaryValues, 
            foreach(1:numberSimsAfterLastCheckpoint, .combine = rbind) %dopar% simulateWithPriors(
				phy = phy, 
				taxonDF = taxonDF, 
                startingPriorsValues = startingPriorsValues, 
				intrinsicPriorsValues = intrinsicPriorsValues, 
				extrinsicPriorsValues = extrinsicPriorsValues, 
                startingPriorsFns = startingPriorsFns, 
				intrinsicPriorsFns = intrinsicPriorsFns, 
				extrinsicPriorsFns = extrinsicPriorsFns, 
				giveUpAttempts = giveUpAttempts, 
                freevector = freevector, 
				timeStep = timeStep, 
				intrinsicFn = intrinsicFn,
				extrinsicFn = extrinsicFn, 
				verbose = verboseNested, 
				checks = FALSE
				)
			)
		}
    # stop multicore processes
    stopMulticore(cluster)
    #
    attr(trueFreeValuesANDSummaryValues, "freevector") <- freevector
	#
    return(trueFreeValuesANDSummaryValues)
    }

#checkNiter <- function(
#  niter.brown = 25, niter.lambda = 25,
#  niter.delta = 25, niter.OU = 25, niter.white = 25
#  ){
#
#    if(niter.brown<2){stop("niter.brown must be at least 2")}
#    if(niter.lambda<2){stop("niter.lambda must be at least 2")}
#    if(niter.delta<2){stop("niter.delta must be at least 2")}
#    if(niter.OU<2){stop("niter.OU must be at least 2")}
#    if(niter.white<2){stop("niter.white must be at least 2")}
#
#    }



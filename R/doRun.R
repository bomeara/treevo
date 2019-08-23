#' Run Approximate Bayesian Computation for Phylogenetic Comparative Methods
#' 
#' The \code{doRun} functions are the main interface for
#' \code{TreEvo} users to do Approximate Bayesian Computation (ABC) analysis, 
#' effectively wrapping the \code{\link{simulateWithPriors}} functions
#' to perform simulations, which themselves are wrappers for
#' the \code{\link{doSimulation}} function. The two current \code{doRun}
#' functions are \code{doRun_prc} which applies a partial
#' rejection control (PRC) ABC analysis over multiple generations of simulations, 
#' and \code{doRun_rej} which performs a full rejection ('rej') ABC analysis.
#' 
#' Both \code{doRun} functions take an input phylogeny (\code{phy}), 
#' observed trait data set (\code{traits}), models
#' (\code{intrinsicFn}, \code{extrinsicFn}), and priors
#' (\code{startingPriorsValues}, \code{startingPriorsFns},
#' \code{intrinsicPriorsValues}, 
#' \code{intrinsicPriorsFns}, \code{extrinsicPriorsValues},
#' \code{extrinsicPriorsFns}). Pulling from the priors, it
#' simulates an initial set of simulations (\code{nInitialSims}).
#' This set of simulations is boxcox transformed , and a PLS regression (see
#' \code{\link{methodsPLS}}) is performed for each free parameter
#' to determine the most informative summary statistics (using
#' \code{variance.cutoff}). The euclidean distance is calculated
#' between each each initial simulation's most informative summary
#' statistics and the input observed data.
#' 
#' Following that step, the approach of the two functions diverge drastically.
#' 
#' For function \code{doRun_rej}, those simulations whose
#' most informative summary statistics fall below
#' \code{abcTolerance} are kept as accepted 'particles'
#' (simulations runs), describing the
#' posterior distribution of parameters. No additional
#' generations of simulations are performed.
#' 
#' Coversely, function \code{doRun_prc} performs an ABC-PRC analysis, 
#' which follows a much more complicated algorithm with additional
#' generations. In PRC, tolerance is set based on
#' \code{epsilonProportion} and \code{epsilonMultiplier}, and the analysis
#' begins with generation 1 and continues through a
#' number of generations equal to \code{nStepsPRC}.
#' Single simulation runs ('particles') are accepted if the distance of the
#' most informative summary stats to the original data
#' is less than this tolerance. A generation
#' is complete when enough particles (\code{numParticles}) have been accepted. These
#' particles make up the distribution for the next generation. The accepted
#' particles from the final generation describe the posterior distributions of
#' parameters.
#' 


#' @inheritParams doSimulation
#' @inheritParams simulateWithPriors
#' @inheritParams methodsPLS

#' @param traits Data matrix with rownames identical to \code{phy@tip.label}. 
#' If a vector, \code{traits} will be coerced to a matrix, with element names as rownames.

# @param startingValuesGuess Optional guess of starting values.

# @param intrinsicValuesGuess Optional guess of intrinsic values.

# @param extrinsicValuesGuess Optional guess of extrinsic values.

# TreeYears = 1000000 if tree length is 1 million of years, 1000 if 1 thousand years, etc.

#' @param numParticles Number of accepted particles per PRC generation.

#' @param standardDevFactor Standard deviation for mutating states each time a new particle is generated in a PRC generation.

#' @param nInitialSims Number of initial simulations used to calibrate particle rejection control algorithm.
#' If not given, this will be a function of the number of freely varying parameters.

# @param plot If \code{TRUE}, plots distance of each simulation.

#' @param epsilonProportion Sets tolerance for initial simulations.

#' @param epsilonMultiplier Sets tolerance on subsequent PRC generations.

#' @param nStepsPRC Number of PRC generations to run.

#' @param nRuns Number of independent PRC runs to be performed,
#' each consisting of independent sets of
#' initial simulations and PRC generations. Note that runs are
#' run \emph{sequentially}, and not in parallel, 
#' as the generation of particles within each run
#' makes use of the multicore functionality.
#' If \code{nRuns} is greater than 1, the output from
#' \code{doRun_prc} will be a list object composed of
#' multiple output lists, as described.

#' @param nInitialSimsPerParam If \code{nInitialSims}
#' is not given by the user, the number of initial simulations
#' performed to calibrate the particle rejection algorithm
#' will instead be the number of free parameters
#' in the model multiplied by \code{nInitialSimsPerParam}.

#' @param jobName Optional job name.

#' @param stopRule If \code{TRUE}, an analysis will be terminated prior to set number
#' of generations if the ratio of all parameters from one generation to the
#' next falls below the \code{stopValue}.

#' @param stopValue Threshold value for terminating an analysis prior to
#' \code{nStepesPRC}.

#' @param maxAttempts Maximum number attempts made
#' in \code{while} loop within the PRC algorithm.
#' If reached, the algorithm will terminate with an error message.
#' By default, this is infinite, and thus there is no
#' effective limit without user intervention.

#' @param abcTolerance Proportion of accepted simulations.

#' @param checkpointFile Optional file name for checkpointing simulations.

#' @param checkpointFreq Saving frequency for checkpointing.

#' @param savesims Option to save individual simulations, output to a .Rdata file.

#' @param saveData Option to save various
#' run information during the analysis, 
#' including summary statistics from analyses,
#' output to external .Rdata and .txt files.


#' @param verboseParticles If \code{TRUE} (the default),
#' a large amount of information about parameter estimates
#' and acceptance of particles is output to
#' console via \code{message} as \code{doRun_prc} runs.

#' @param diagnosticPRCmode If \code{TRUE} 
#'(\emph{not} the default), the function will be very
#' noisy about characteristics of the PRC algorithm as it runs.

#' @param multicoreSuppress This argument suppresses the
#' multicore code and will apply a plain vanilla serial \code{for()} loop instead
#' of \code{doPar}. Mainly to be used for developer diagnostic purposes.

# param niter.goal Adjust number of starting points for package \code{geiger}
# to return the best parameter estimates this number of times on average.

#' @return
#' The output of these two functions are lists,
#' composed of multiple objects, 
#' which differ slightly in their content among the
#' two functions. For \code{doRun_prc}, the output is:

#' \describe{
#' \item{input.data}{Input variables: 
#' \code{jobName}, \code{nTaxa} (number of taxa),
#' \code{nInitialSims}, \code{generation.time}, \code{TreeYears},
#' \code{timeStep}, \code{totalGenerations},
#' \code{epsilonProportion}, \code{epsilonMultiplier}, 
#' \code{nRuns} (number of runs), \code{nStepsPRC},
#' \code{numParticles}, \code{standardDevFactor}.}

#' \item{priorList}{List of prior distributions. This is used
#' for doing post-analysis comparisons between prior and
#' posterior distributions, such as with function \code{plotPosteriors}.}

#' \item{particleDataFrame}{DataFrame with information from each 
#' simulation, including generation, attempt, \code{id}, \code{parentid},
#' distance, weight, and parameter states.}

#' \item{toleranceVector}{Final tolerance vector.}

#' \item{phy}{Input phylogeny.}

#' \item{traits}{Input traits.} 

#' \item{simTime}{Processor time for initial simulations.}

#' \item{time.per.gen}{Processor time for subsequent generations.}



# \item{credibleInt}{Credible Interval calculation for
# each free parameter of the final generation.}

#' \item{postSummary}{Summarizes the posterior distribution from
#' the final generation for all free parameters, giving the mean,
#' standard deviation and Highest Posterior
#' Density (at a 0.8 alpha) for each parameter.}
#' }
#' 
#' If \code{nRuns} is greater than 1, the output from
#' \code{doRun_prc} will be a list object composed of
#' multiple output lists, as described.
#' 

#' For \code{doRun_rej}, the output is:

#' \describe{
#' \item{input.data}{Input variables: 
#' \code{jobName}, \code{nTaxa} (number of taxa),
#' \code{nInitialSims}, \code{generation.time}, \code{TreeYears},
#' \code{timeStep}, \code{totalGenerations},
#' \code{epsilonProportion}, \code{epsilonMultiplier}, 
#' \code{nRuns} (number of runs), \code{nStepsPRC},
#' \code{numParticles}, \code{standardDevFactor}.
#' }

#' \item{priorList}{List of prior distributions. This is
#' used for doing post-analysis comparisons between prior
#' and posterior distributions, such as
#' with function \code{plotPosteriors}.}

#' \item{phy}{Input phylogeny.}

#' \item{traits}{Input traits.}

#' \item{trueFreeValuesANDSummaryValues}{Parameter
#' estimates and summary stats from all sims.}

#' \item{simTime}{Processor time for simulations.}

# \item{whichVip}{Matrix of vip summary statistics for each free parameter.}

#' \item{abcDistancesRaw}{Euclidean distances
#' for each simulation and free parameter.}

#' \item{particleDataFrame}{DataFrame with information from each
#' simulation, including generation, attempt, id,
#' parentid, distance, weight, and parameter states.}

# \item{credibleInt}{Credible Interval calculation for each
# free parameter of the final generation.}

#' \item{postSummary}{Summarizes the posterior distribution
#' from the final generation for all free parameters,
#' giving the mean, standard deviation and Highest Posterior
#' Density (at a 0.8 alpha) for each parameter.}

#' \item{parMeansList}{A list of parameter means for both
#' fixed and unfixed parameters, sorted into a list of
#' three vectors (starting, intrinsic, extrinsic).}
#' }
#' 
#' 

#' @author Brian O'Meara and Barb Banbury

#' @references Sisson et al. 2007, Wegmann et
#' al. 2009

# @references O'Meara and Banbury, unpublished;
# Sisson et al. 2007, Wegmann et al. 2009

#' @examples
#' 
#' \donttest{
#' set.seed(1)
#' data(simRunExample)
#' 
#' # NOTE: the example analyses below sample too few particles, 
#'     # over too few steps, with too few starting simulations
#'     # - all for the sake of examples that demonstrate these
#'     # within a reasonable time-frame
#' 
#' ## Please set these values to more realistic levels for your analyses!
#' 
#' resultsPRC <- doRun_prc(
#'   phy = simPhyExample, 
#'   traits = simCharExample, 
#'   intrinsicFn = brownianIntrinsic, 
#'   extrinsicFn = nullExtrinsic, 
#'   startingPriorsFns = "normal", 
#'   startingPriorsValues = list(c(mean(simCharExample[, 1]), 
#'          sd(simCharExample[, 1]))), 
#'   intrinsicPriorsFns = c("exponential"), 
#'   intrinsicPriorsValues = list(10), 
#'   extrinsicPriorsFns = c("fixed"), 
#'   extrinsicPriorsValues = list(0), 
#'   generation.time = 10000, 
#'   nRuns = 2, 
#'   nStepsPRC = 3, 
#'   numParticles = 20, 
#'   nInitialSimsPerParam = 10, 
#'   jobName = "examplerun_prc", 
#'   stopRule = FALSE, 
#'   multicore = FALSE, 
#'   coreLimit = 1
#'   )
#' 
#' resultsPRC
#' 
#' #one should make sure priors are uniform with doRun_rej!
#' 
#' resultsRej <- doRun_rej(
#'     phy = simPhyExample, 
#'     traits = simCharExample, 
#'     intrinsicFn = brownianIntrinsic, 
#'     extrinsicFn = nullExtrinsic, 
#'     startingPriorsFns = "normal", 
#'     startingPriorsValues = list(c(mean(simCharExample[, 1]),
#'          sd(simCharExample[, 1]))), 
#'     intrinsicPriorsFns = c("exponential"), 
#grep for normal in pkg
#'     intrinsicPriorsValues = list(10), 
#'     extrinsicPriorsFns = c("fixed"), 
#'     extrinsicPriorsValues = list(0), 
#'     generation.time = 10000, 
#'     jobName = "examplerun_rej", 
#'     abcTolerance = 0.05, 
#'     multicore = FALSE, 
#'     coreLimit = 1
#'     )
#' 
#' resultsRej
#' }
#' 


# ARE each particle-gather interval of the PRC algorithm a GENERATION OR a STEP? PLEASE CLARIFY
    # what a weird comment - generations are the time increments (in abs time) used for evolving traits along the tree

##old, from before DWB: This seems to be working if partialResults does not exist.  If checkpoint = TRUE, then run fails.
    # DWB - I guess I should test that... sigh...



    
#' @name doRun
#' @rdname doRun
#' @export
doRun_prc <- function(
    phy, traits, intrinsicFn, extrinsicFn, 
	startingPriorsValues, startingPriorsFns, 
    intrinsicPriorsValues, intrinsicPriorsFns, 
	extrinsicPriorsValues, extrinsicPriorsFns, 
    #
    # OLD COMMENTING (before DWB):    
    #the doRun_prc function takes input from the user and then 
		# automatically guesses optimal parameters, though user overriding is also possible.
    #the guesses are used to do simulations near the expected region.
		# If omitted, they are set to the midpoint of the input parameter matrices
    #
    # DWB :so it seems like the above guess parameters are for the 'override' of this feature
		# So... I just got rid of them, its too confusing..
    #
    #startingValuesGuess = c(), intrinsicValuesGuess = c(), extrinsicValuesGuess = c(), 
    #
    numParticles = 300, 
    nStepsPRC = 5, 
    nRuns = 2, 
    nInitialSims = NULL, 
    nInitialSimsPerParam = 100, 
    #
    generation.time = 1000, 
    TreeYears = max(branching.times(phy)) * 1e6, 
    #
    multicore = FALSE, 
    coreLimit = NA, 
    multicoreSuppress = FALSE, 
    #
    standardDevFactor = 0.20, 
    epsilonProportion = 0.70, 
    epsilonMultiplier = 0.70, 
    #
    validation = "CV", 
    scale = TRUE, 
    variance.cutoff = 95, 
    #niter.goal = 5, 
    #
    stopRule = FALSE, 
    stopValue = 0.05, 
    maxAttempts = Inf, 
    diagnosticPRCmode = FALSE, 
    #
    jobName = NA, 
    saveData = FALSE, 
    verboseParticles = TRUE
    ){                        
    #
    #
    functionStartTime <- proc.time()[[3]]
    #	
    if (!ape::is.binary.phylo(phy)) {
        warning("Tree is not fully dichotomous, this may cause issues!")
        }
	#
	if(is.null(phy$edge.length)){
		stop("Input phylogeny lacks necessary edge lengths")
		}
    #    
    if(!is.numeric(maxAttempts)){
        stop("maxAttempts must be numeric")
        }
	#
	if(!is.numeric(nRuns)){
		stop("nRuns must be numeric")
		}
	if(nRuns < 1){
		stop("nRuns must 1 or greater")
		}
    #
    timeStep <- generation.time/TreeYears
    message(paste0(
		"The effective timeStep for this tree will be ",
		signif(timeStep, 2), ",\n",
		"     as a proportion of tree height (root to furthest tip)..."))
    #
    edgesRescaled <- phy$edge.length/max(node.depth.edgelength(phy))
    message("Rescaling edge lengths relative to maximum tip-to-root distance...")
	#
    #
    minEdgeRescaled <- min(edgesRescaled)
    minEdgeRescaledNZ <- min(edgesRescaled[edgesRescaled>0])
    #
    if(minEdgeRescaled == 0){
        message(
			"The smallest edge length on the input tree is ZERO LENGTH..."
			)
        message(paste0(
			"  Edge with smallest rescaled NON-ZERO length on tree\n",
			"      is ",
			signif(minEdgeRescaledNZ, 2), 
			" as a proportion of tip-to-root distance"
			))
        message(paste0(
			"   (This is ", signif(minEdgeRescaledNZ*TreeYears, 2),
			" in the same TreeYears units\n",
			"      as used for the input generation.time ( = ",
			generation.time, "))"
			))
    }else{
        message(paste0(
			"The smallest edge length on the input tree is ", 
			signif(minEdgeRescaled, 2)
			))
        message(paste0(
			"   (This is ", signif(minEdgeRescaled*TreeYears, 2), 
			" in the same TreeYears units\n",
			"      as used for the input generation.time ( = ",
			generation.time, "))"
			))
        }
    #
    if(max(edgesRescaled) < timeStep) {
        stop(paste0(
			"Tree has *NO* rescaled branches longer than generation.time/TreeYears \n",
			"     thus *NO* simulated evolutionary change can occur!"
			))
        }
    #
    if(minEdgeRescaled < timeStep) {
        warning(paste0(
			"Tree has rescaled branches shorter than generation.time/TreeYears\n",
			"    No trait evolution will be simulated on these branches and\n",
			"    thus ML summary statistic functions may fail or produce unexpected results."
			))
        #message("Tree has zero or nearly zero length branches")
        }
    #
    totalGenerations <- sum(sapply(edgesRescaled, function(x) floor(x/timeStep)))
    message(paste0( 
		"Given generation time, a total of ", round(totalGenerations), 
		" generations are\n", "     expected to occur over this tree"
		))
	#####################################
	#
	# save input data for use later
    input.data <- rbind(
		jobName = jobName, 
		nTaxa = Ntip(phy), 
        nInitialSims = nInitialSims, 
		nInitialSimsPerParam = nInitialSimsPerParam, 
        generation.time = generation.time, 
		TreeYears = TreeYears, 
		timeStep = timeStep, 
        totalGenerations = totalGenerations, 
		epsilonProportion = epsilonProportion, 
        epsilonMultiplier = epsilonMultiplier, 
		nRuns = nRuns, 
		nStepsPRC = nStepsPRC, 
        numParticles = numParticles, 
		standardDevFactor = standardDevFactor
		)
	#
	####################################
    #
    #splits <- getSimulationSplits(phy) #initialize this info
    taxonDF <- getTaxonDFWithPossibleExtinction(phy)
    #
    # get freeVector
    freeVector <- getFreeVector(
		startingPriorsFns = startingPriorsFns, 
		startingPriorsValues = startingPriorsValues, 
        intrinsicPriorsFns = intrinsicPriorsFns, 
		intrinsicPriorsValues = intrinsicPriorsValues, 
        extrinsicPriorsFns = extrinsicPriorsFns, 
		extrinsicPriorsValues = extrinsicPriorsValues
		)
    numberParametersTotal <- length(freeVector)
    numberParametersFree <- sum(freeVector)
    if(numberParametersFree<1){
      stop(paste0(
		"No freely varying parameters found; analysis cannot continue.\n",
		"    Check prior functions and values."
		))
      }
    namesParFree <- names(freeVector)[freeVector]
    #
    # get prior list
    priorList <- getPriorList(
        startingPriorsValues = startingPriorsValues, 
        intrinsicPriorsValues = intrinsicPriorsValues, 
        extrinsicPriorsValues = extrinsicPriorsValues, 
        startingPriorsFns = startingPriorsFns, 
        intrinsicPriorsFns = intrinsicPriorsFns, 
        extrinsicPriorsFns = extrinsicPriorsFns, 
        numberParametersTotal = numberParametersTotal
        )
    ##    
    #initialize weighted mean sd matrices
    weightedMeanParam <- matrix(
		nrow = nStepsPRC, 
        ncol = numberParametersFree
		)
    param.stdev <- matrix(nrow = nStepsPRC, ncol = numberParametersFree)
    #
	colnames(weightedMeanParam) <- colnames(param.stdev) <- namesParFree 
	genRowNames <- paste0("Gen ", c(1: nStepsPRC), sep = "") 
    rownames(weightedMeanParam) <- rownames(param.stdev) <- genRowNames
    #
    if (is.null(nInitialSims)) {
        #modified from 1000, which is rather computationally abusive
		nInitialSims <- nInitialSimsPerParam*numberParametersFree    
		#
        message(paste0(
			"Number of initial simulations to be performed (nInitialSims) not given\n", 
            "    Number of initial simulations will thus be the\n",
			"    product of the number of free parameters multiplied by ",
            nInitialSimsPerParam, "\n",
            "    ( = ", nInitialSims, " initial simulations)\n"))
        }
    #
    # get summary values for observed data
		# INCREDIBLY!
		#
		# this is really the **ONLY** time when doRun ever really looks at input traits
		#
    originalSummaryValues <- summaryStatsLong(
        phy = phy, 
        traits = traits
        )
    #
    results <- list()
    class(results) <- c("multiRun_doRun_prc", class(results))
    #
    for(runID in 1:nRuns){
        if(nRuns>1){
            message(paste0("Beginning PRC Run ", runID, 
				" out of a total of ", nRuns, "..."))
            }
        # INITIAL SIMULATIONS
        initialSimsRes <- initialSimsPRC(
            nrepSim = nInitialSims, 
            phy = phy, 
            taxonDF = taxonDF, 
            startingPriorsValues = startingPriorsValues, 
            intrinsicPriorsValues = intrinsicPriorsValues, 
            extrinsicPriorsValues = extrinsicPriorsValues, 
            startingPriorsFns = startingPriorsFns, 
            intrinsicPriorsFns = intrinsicPriorsFns, 
            extrinsicPriorsFns = extrinsicPriorsFns, 
            freevector = freeVector, 
            timeStep = timeStep, 
            intrinsicFn = intrinsicFn, 
            extrinsicFn = extrinsicFn, 
            nStepsPRC = nStepsPRC, 
            coreLimit = coreLimit, 
            multicore = multicore, 
            numberParametersFree = numberParametersFree, 
            originalSummaryValues = originalSummaryValues, 
            epsilonProportion = epsilonProportion, 
            epsilonMultiplier = epsilonMultiplier, 
            validation = validation, 
            variance.cutoff = variance.cutoff, 
            saveData = saveData, 
            jobName = jobName
            )        
        #
        # pull the PLS model list out
        pls.model.list <- initialSimsRes$pls.model.list
        #
        #------------------ ABC-PRC (Start) ------------------
        message("Beginning partial rejection control algorithm...")
        #
        nameVector <- c(
			"generation", "attempt",
			"id", "parentid", "distance",
			"weight", namesParFree)
        #
        # save before initiating PRC procedure
        if(saveData){
            save.image(
				file = paste0("WS", jobName, ".Rdata", sep = "")
				)
            }
        #
        # here I removed the initial loop (DWB, December 2017)
        #
        time.per.gen <- numeric(length = nStepsPRC)            #time.per.Step <- 
        for(dataGenerationStep in 1:nStepsPRC) {
            #    
            #dataGenerationStep <- dataGenerationStep+1
            if(verboseParticles){
                message("\n", "STARTING DATA GENERATION STEP ", dataGenerationStep, "\n")
                }
            #
            if(dataGenerationStep == 1){
                particleDataFrame <- data.frame()        
                prevGenParticleList <- NULL            
            }else{
				##
                #prevGenParticleWeights <- sapply(prevGenParticleList, function(x) x$weight)
				#
                #particleWeights <- numeric(length = numParticles)				
                #stores weights for each particle.
                    #Initially, assume infinite number of possible particles 
					# (so might not apply in discrete case)
				#
                # why pull particle weights and particle list as this stupid seperate vector
					# why not just get from particleDataFrame ?
					# okay let's do that
                prevGenParticleList <- currParticleList
                }
            #stores parameters in model for each particle
            particleParameters <- matrix(
				nrow = numParticles, 
                ncol = (length(startingPriorsValues) +  
					length(intrinsicPriorsValues) + 
					length(extrinsicPriorsValues))
				)
            #
            weightScaling = 0
            particleDistance = rep(NA, numParticles)
            particle <- 1
            attempts <- 0        
            # acceptance rate of particles - reset to 0.5
            expParticleAcceptanceRate <- 0.5			
            #
            # create a list where we'll save information for each particle for this generation
            currParticleList <- list()
            #
            if(verboseParticles){
                message(
					"Successes ", 
                    "  Attempts ", 
                    "  Distance   ", 
                    paste0(namesParFree, "  ", collapse = "")
                    #"  Exp. # of Attempts Req. ", 
                    #"Params." #\n
                    )
                }        
            #message("Beginning partial rejection control algorithm...")
            #
            start.time <- particleStartTime <- proc.time()[[3]]
            #
            # set nAcceptedParticles to 0
            nAcceptedParticles <- nFailedParticles <- 0
            #
            # get tolerance value
            toleranceValue <- initialSimsRes$toleranceVector[dataGenerationStep]
            #
            while (nAcceptedParticles <= numParticles) {
                #
                if(attempts>maxAttempts){
                    stop("maxAttempts exceeded in while() loop of PRC algorithm")
                    }
                #
                # check how many particles are needed
                nParticlesNeeded <- (numParticles-nAcceptedParticles)
                # siulate that number times the apparent rate of success
                nSim <- nParticlesNeeded/expParticleAcceptanceRate
                # round up so its a multiple of the number of cores
                nSim <- boostNsim(nSims = nSim, nCores = coreLimit)
                # repeat until you have enough particles
                newParticleList <- list()
                #print(nSim)
                if(!multicoreSuppress){
                    newParticleList <- simParticlePRCParallel(
                        nSim = nSim, 
                        multicore = multicore, 
                        coreLimit = coreLimit, 
                        phy = phy, 
                        taxonDF = taxonDF, 
                        timeStep = timeStep, 
                        intrinsicFn = intrinsicFn, 
                        extrinsicFn = extrinsicFn, 
                        startingPriorsValues = startingPriorsValues, 
                        intrinsicPriorsValues = intrinsicPriorsValues, 
                        extrinsicPriorsValues = extrinsicPriorsValues, 
                        startingPriorsFns = startingPriorsFns, 
                        intrinsicPriorsFns = intrinsicPriorsFns, 
                        extrinsicPriorsFns = extrinsicPriorsFns, 
                        originalSummaryValues = originalSummaryValues, 
                        pls.model.list = pls.model.list, 
                        toleranceValue = toleranceValue, 
                        prevGenParticleList = prevGenParticleList, 
                        standardDevFactor = standardDevFactor, 
                        numParticles = numParticles
                        )
                }else{
                    newParticleList <- simParticlePRCNotParallel(
                        nSim = nSim, 
                        multicore = multicore, 
                        coreLimit = coreLimit, 
                        phy = phy, 
                        taxonDF = taxonDF, 
                        timeStep = timeStep, 
                        intrinsicFn = intrinsicFn, 
                        extrinsicFn = extrinsicFn, 
                        startingPriorsValues = startingPriorsValues, 
                        intrinsicPriorsValues = intrinsicPriorsValues, 
                        extrinsicPriorsValues = extrinsicPriorsValues, 
                        startingPriorsFns = startingPriorsFns, 
                        intrinsicPriorsFns = intrinsicPriorsFns, 
                        extrinsicPriorsFns = extrinsicPriorsFns, 
                        originalSummaryValues = originalSummaryValues, 
                        pls.model.list = pls.model.list, 
                        toleranceValue = toleranceValue, 
                        prevGenParticleList = prevGenParticleList, 
                        standardDevFactor = standardDevFactor, 
                        numParticles = numParticles
                        )            
                    }
                #print("ended")
                #print(newParticleList)
                # will need to count successful particles
                acceptedParticles <- !is.na(newParticleList)
                nAcceptedNew <- sum(acceptedParticles)
                nFailedNew <- nSim-nAcceptedNew
                nAcceptedParticles <- nAcceptedParticles+nAcceptedNew
                nFailedParticles <- nFailedParticles+nFailedNew
                # updated expParticleAcceptanceRate for properly calculate number of necc nSim
                    # but only if its going to be finite, otherwise go with a 0.1 rate
                if(nAcceptedParticles>0){
                    expParticleAcceptanceRate <- nAcceptedParticles/(nAcceptedParticles+nFailedParticles)
                }else{
                    expParticleAcceptanceRate <- 0.1
                    }
                #
                # now, if there are new particles, let's probably clean newParticleList
                    # and add to other saved particles
                if(nAcceptedNew>1){
                    # remove NAs (failed particles)
                    newParticleList <- newParticleList[acceptedParticles]
                    #
                    # remove particles that go beyond the number needed
                    #if(length(res)>nParticlesNeeded){
                    #    res <- res[1:particlesNeeded]}
                    #       is this really a bad thing??
                    #
                    #change particle count value
                    #nAcceptedParticles <- nAcceptedParticles+length(newParticleList)
                    # save successful particles to currParticleList
                    currParticleList <- append(currParticleList, newParticleList)
                    #if( plotrix::listDepth(currParticleList)>2){
                    #    print(newParticleList)
                    #    print(currParticleList)
                    #    stop()
                    #    }
                    # updated number of attemped particles so far
                    attempts <- attempts+nSim
                    }
                }
            #print(currParticleList)
            #
            # append particle IDs to each accepted particle    
            for(i in 1:length(currParticleList)){
                currParticleList[[i]]$id <- i
                parVector <- c(currParticleList[[i]]$startingValues, 
                    currParticleList[[i]]$intrinsicValues, currParticleList[[i]]$extrinsicValues)
                # keep only free parameters
                parVector <- parVector[freeVector]
                vectorForDataFrame <- c(dataGenerationStep, attempts, currParticleList[[i]]$id, currParticleList[[i]]$parentid, 
                    currParticleList[[i]]$distance, currParticleList[[i]]$weight, parVector)
                #    
                #print(i)
                #print(currParticleList[[i]])
                #print(currParticleList[[i]]$distance)
                #print(signif(currParticleList[[i]]$distance, 2))
                #
                if(diagnosticPRCmode){
                    message(
						"\n\nlength of vectorForDataFrame = ", 
							length(vectorForDataFrame), 
                        "\n", "length of startingValues = ", 
							length(currParticleList[[i]]$startingValues), 
                        "\nlength of intrinsicValues = ", 
							length(currParticleList[[i]]$intrinsicValues), 
                        "\nlength of extrinsicValues = ", 
							length(currParticleList[[i]]$extrinsicValues), 
                        "\ndistance = ", 
							currParticleList[[i]]$distance, 
                        "\nweight = ", 
							currParticleList[[i]]$weight, 
						"\n", vectorForDataFrame, "\n")
                    }
                #
                #NOTE THAT WEIGHTS AREN'T NORMALIZED IN THIS DATAFRAME
                particleDataFrame <- rbind(particleDataFrame, vectorForDataFrame)
                #
                if(verboseParticles){
                    message(
                        paste(
                            i, "         ", attempts, "       "
                            #floor(nAcceptedParticles*attempts/i), 
							#"                    ", 
                            , signif(currParticleList[[i]]$distance, 2), "      "
                            , if(!identical(startingPriorsFns, "fixed")){
                                paste0(signif(
									currParticleList[[i]]$startingValues, 2), 
									collapse = "       ")
                                }, "       "
                            , if(!identical(intrinsicPriorsFns, "fixed")){
                                paste0(signif(
									currParticleList[[i]]$intrinsicValues, 2), 
									collapse = "        ")
                                }, "        "
                            , if(!identical(extrinsicPriorsFns, "fixed")){
                                paste0(signif(
									currParticleList[[i]]$extrinsicValues, 2), 
									collapse = "        ")
                                }, "        "
                                                    
                            #, signif(currParticleList[[i]]$distance, 2)
                            )
                        )
                    }    
                }            
            #
            if(dataGenerationStep == 1){
                names(particleDataFrame) <- nameVector
                }
            #
            #rejects.gen.one <- (dim(subset(particleDataFrame, particleDataFrame$id<0))[1])/(dim(subset(particleDataFrame, ))[1])
            #rejects <- c()
            #
            
            # normalizing the weights for this generation to the sum of the weights
            #
            particlesThisGen <- which(particleDataFrame$generation == dataGenerationStep)
            particleDataFrame[particlesThisGen, ]$weight <- (
                (particleDataFrame[particlesThisGen, ]$weight) / 
                (sum(particleDataFrame[particlesThisGen, ]$weight))
                )
            #    
            timePRCStep <- proc.time()[[3]]-start.time
            time.per.gen[dataGenerationStep] <- timePRCStep
            #rejects.per.gen <- (dim(subset(particleDataFrame, 
				# particleDataFrame$id<0))[1])/(
                # dim(subset(particleDataFrame[which(
					# particleDataFrame$generation == dataGenerationStep), ], ))[1])
            #
            #rejects <- c(rejects, rejects.per.gen)
            #
            # create subsets, then calculate param stdev and weighted means
            #
            sub1 <- subset(particleDataFrame, particleDataFrame$generation == dataGenerationStep)
            sub2 <- subset(sub1, sub1$id>0)
            for (i in 1:numberParametersFree){
                param.stdev[dataGenerationStep, i] <- c(sd(sub2[, 6+i]))
                weightedMeanParam[dataGenerationStep, i] <- weighted.mean(sub2[, 6+i], sub2[, 6])
                #c(mean(subset(particleDataFrame, X3>0)[,
					# 7:dim(particleDataFrame)[2]])/subset(particleDataFrame, X3>0)[, 6])
                }        
            #
            if (stopRule){
				####################################################
				#this will stop the PRC from running out to max number of generations
					# if all params are below stopValue
				#################################################
                FF <- rep(1, dim(weightedMeanParam)[2])
                for (check.weightedMeanParam in 1:length(FF)){
                    #
                    if (
						is.na(
							(abs(weightedMeanParam[dataGenerationStep, check.weightedMeanParam] 
							- weightedMeanParam[dataGenerationStep-1, check.weightedMeanParam])
							/ mean(weightedMeanParam[dataGenerationStep, check.weightedMeanParam], 
								weightedMeanParam[dataGenerationStep-1, check.weightedMeanParam]
							))  <=  stopValue
						#this && is here to make sure any NAs are
							# from fixed params and not miscalculations.
						) && mean(
                            weightedMeanParam[dataGenerationStep, check.weightedMeanParam],
							weightedMeanParam[dataGenerationStep-1, check.weightedMeanParam]
							)  ==  0
					){
                        FF[check.weightedMeanParam] <- 0
                    }else{
                        stopValueWeightTest <- abs(
							weightedMeanParam[dataGenerationStep, check.weightedMeanParam]
                            -weightedMeanParam[dataGenerationStep-1, check.weightedMeanParam])
                        stopValueWeightTest <- stopValueWeightTest/
                                mean(weightedMeanParam[dataGenerationStep, check.weightedMeanParam], 
                                    weightedMeanParam[dataGenerationStep-1, check.weightedMeanParam])
                        if (stopValueWeightTest  <=  stopValue){
                            FF[check.weightedMeanParam] <- 0
                            }
                        }
                    #message(FF)
                    }
                if (sum(FF) == 0){
                    message("\n\n\nweightedMeanParam is < ", stopValue, 
						"Analysis is being terminated at", dataGenerationStep
                        , "instead of continuing to ", nStepsPRC, "\n\n\n")
                    dataGenerationStep <- nStepsPRC
                    }
                }    # end if-stopRule
            #
            if(saveData){
                save.image(file = paste0("WS", jobName, ".Rdata", sep = ""))
                #
                prcResults <- list()
                prcResults$input.data <- input.data
                prcResults$priorList <- priorList
                prcResults$particleDataFrame <- particleDataFrame
                names(prcResults$particleDataFrame) <- nameVector
                prcResults$toleranceVector <- initialSimsRes$toleranceVector
                prcResults$phy <- phy
                prcResults$traits <- traits
                prcResults$simTime <- initialSimsRes$simTime
                prcResults$time.per.gen <- time.per.gen
                #
                save(prcResults, 
					file = paste0("partialResults", jobName, ".txt", sep = ""))
                }
            #
            } #for (dataGenerationStep in 1:nStepsPRC) bracket
        #
        # ######### end of PRC particle collection algorithm loop ##############
        #
        names(particleDataFrame) <- nameVector
        #
        particleTime <- proc.time()[[3]]-particleStartTime
        message(paste0(
			"Collection of simulation particles under PRC completed in ",
			signif(particleTime,2), " seconds..."))
        #---------------------- ABC-PRC (End) --------------------------------
        #
        time3 <- proc.time()[[3]]
        genTimes <- c(time.per.gen, time3)
        #
        # save them to prcResults (this hasn't been done yet if save.data = FALSE)
        prcResults <- list()
		#
        prcResults$input.data <- input.data
		#
		# composition of input.data as of 05/07/19
			# jobName = jobName, 
			# nTaxa = Ntip(phy), 
			# nInitialSims = nInitialSims, 
			# nInitialSimsPerParam = nInitialSimsPerParam, 
			# generation.time = generation.time, 
			# TreeYears = TreeYears, 
			# timeStep = timeStep, 
			# totalGenerations = totalGenerations, 
			# epsilonProportion = epsilonProportion, 
			# epsilonMultiplier = epsilonMultiplier, 
			# nRuns = nRuns, 
			# nStepsPRC = nStepsPRC, 
			# numParticles = numParticles, 
			# standardDevFactor = standardDevFactor
		#
        prcResults$phy <- phy
        prcResults$traits <- traits
		prcResults$generation.time <- generation.time
		prcResults$intrinsicFn <- intrinsicFn
		prcResults$extrinsicFn <- extrinsicFn
		#
        prcResults$priorList <- priorList		
		prcResults$freeVector <- freeVector
		#
        prcResults$particleDataFrame <- particleDataFrame
        #names(prcResults$particleDataFrame) <- nameVector
        prcResults$toleranceVector <- initialSimsRes$toleranceVector
		#
        prcResults$simTime <- initialSimsRes$simTime
        prcResults$time.per.gen <- genTimes
        #######################################
        if(nrow(particleDataFrame)>2){
            prcResults$postSummary  <- summarizePosterior(
				particleDataFrame, 
				verboseMultimodal = FALSE,
				# avoid fatal errors due to flat posteriors...
				stopIfFlat = FALSE
				)
			prcResults$parMeansList  <- getSummaryMeans(
				doRunOutput = prcResults)
        }else{
          warning(paste0(
			"Posterior Summaries were not calculated as number",
			" of accepted particles was less than 2"
			))
          }
        ######################################
        if(nRuns<2){
            results <- prcResults
        }else{
            results[[runID]] <- prcResults
            message("Run finished...\n...\n\n")
            }
        }
	#
	if(nRuns>1){
		
		}
    #
    if(multicore){
        registerMulticoreEnv(nCore = 1)
        }
    #
    functionTime <- proc.time()[[3]]-functionStartTime
    message(paste0("Function completed in ", 
		signif(functionTime, 2), " seconds."))
    #
    #message(prcResults)
    return(results)
    }


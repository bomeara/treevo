#' @rdname doRun
#' @export
doRun_rej <- function(
        phy, traits, 
        intrinsicFn, extrinsicFn, 
        startingPriorsValues, startingPriorsFns, 
        intrinsicPriorsValues, intrinsicPriorsFns, 
        extrinsicPriorsValues, extrinsicPriorsFns, 
        #startingValuesGuess = c(), intrinsicValuesGuess = c(), extrinsicValuesGuess = c(), 
        generation.time = 1000, 
        TreeYears = max(branching.times(phy)) * 1e6, 
        multicore = FALSE, 
        coreLimit = NA, 
        validation = "CV", 
        scale = TRUE, 
        #
        nInitialSims = NULL, 
        nInitialSimsPerParam = 100, 
        #
        variance.cutoff = 95, 
        #niter.goal = 5, 
        standardDevFactor = 0.20, 
        #
        jobName = NA, abcTolerance = 0.1, 
        checkpointFile = NULL, checkpointFreq = 24,
        savesims = FALSE) {    
    ####################################################
    
    
    if (!ape::is.binary.phylo(phy)) {
        warning("Tree is not fully dichotomous, this may lead to issues")
    }
    startTime <- proc.time()[[3]]

    timeStep <- generation.time/TreeYears
    message(paste0("The effective timeStep for this tree will be ", signif(timeStep, 2), 
        ", as a proportion of tree height (root to furthest tip)..."))

    edgesRescaled <- phy$edge.length/max(node.depth.edgelength(phy))
    message("Rescaling edge lengths relative to maximum tip-to-root distance...")
    
    if(max(edgesRescaled) < timeStep) {
        stop("Tree has *NO* rescaled branches longer than generation.time/TreeYears, no simulated evol change can occur!")
    }
    if(min(edgesRescaled) < timeStep) {
        warning("Tree has rescaled branches shorter than generation.time/TreeYears; no evol change can be assigned to these, and ML summary stat functions may fail!")
        #message("Tree has zero or nearly zero length branches")
    }
    
    totalGenerations <- sum(sapply(edgesRescaled, function(x) floor(x/timeStep)))
    message("Given generation time, a total of ", round(totalGenerations), " generations will occur over this tree")
        
    #splits <- getSimulationSplits(phy) #initialize this info
    taxonDF <- getTaxonDFWithPossibleExtinction(phy = phy)

    # get freevector
    freevector <- getFreeVector(startingPriorsFns = startingPriorsFns, 
                        startingPriorsValues = startingPriorsValues, 
                        intrinsicPriorsFns = intrinsicPriorsFns, 
                        intrinsicPriorsValues = intrinsicPriorsValues, 
                        extrinsicPriorsFns = extrinsicPriorsFns, 
                        extrinsicPriorsValues = extrinsicPriorsValues
                        )
	#
    numberParametersTotal <- length(freevector)
    numberParametersFree <- sum(freevector)
	#
    if(numberParametersFree<1){
        stop("No freely varying parameters found; analysis cannot continue. Check prior functions and values")
        }
    namesParFree <- names(freevector)[freevector]
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
    #
    #initialize guesses, if needed
    #if (length(startingValuesGuess) == 0) { #if no user guesses, try pulling a value from the prior
    #    startingValuesGuess <- rep(NA, length(startingPriorsFns))
    #    for (i in 1:length(startingPriorsFns)) {
    #        startingValuesGuess[i] <- pullFromPrior(startingPriorsValues[[i]], startingPriorsFns[i])
    #    }
    #}
    #if (length(intrinsicValuesGuess) == 0) { #if no user guesses, try pulling a value from the prior
    #    intrinsicValuesGuess <- rep(NA, length(intrinsicPriorsFns))
    #    for (i in 1:length(intrinsicPriorsFns)) {
    #        intrinsicValuesGuess[i] <- pullFromPrior(intrinsicPriorsValues[[i]], intrinsicPriorsFns[i])
    #    }
    #}
    #if (length(extrinsicValuesGuess) == 0) { #if no user guesses, try pulling a value from the prior
    #    extrinsicValuesGuess <- rep(NA, length(extrinsicPriorsFns))
    #    for (i in 1:length(extrinsicPriorsFns)) {
    #        extrinsicValuesGuess[i] <- pullFromPrior(extrinsicPriorsValues[[i]], extrinsicPriorsFns[i])
    #    }
    #}
    #
    if (is.null(nInitialSims)) {
        nInitialSims <- nInitialSimsPerParam*numberParametersFree    
        message(paste0("Number of initial simulations to be performed (nInitialSims) not given\n", 
            "Number of initial simulations will instead be product of the number of free parameters multiplied by ",
            nInitialSimsPerParam, 
            "\n( = ", nInitialSims, " initial simulations)"))
        }
    #
    #Used to be multiple tries where nrepSim = nInitialSims*((2^try)/2).
        #If initial simulations are not enough, and we need to try again then new analysis will double number of initial simulations
    nrepSim <- nInitialSims
    #message(paste0("Number of initial simulations set to ", nrepSim))     #, "\n"
    if(!is.null(checkpointFile)) {
        save(list = ls(), file = paste0(checkpointFile, ".intialsettings.Rsave", sep = ""))
        }
    #
    #Figure out how many iterations to use for optimization in Geiger.
    #it actually runs faster without checking for cores. And we parallelize elsewhere
    #niter.brown.g <- getBM(phy = phy, dat = traits, niterN = 100, niter.goal = niter.goal)$niter.g
    #niter.lambda.g <- getLambda(phy = phy, dat = traits, niterN = 100, niter.goal = niter.goal)$niter.g
    #niter.delta.g <- getDelta(phy = phy, dat = traits, niterN = 100, niter.goal = niter.goal)$niter.g
    #niter.OU.g <- getOU(phy = phy, dat = traits, niterN = 100, niter.goal = niter.goal)$niter.g
    #niter.white.g <- getWhite(phy = phy, dat = traits, niterN = 100, niter.goal = niter.goal)$niter.g
    #
    #message(paste0("Setting number of starting points for Geiger optimization to ", 
    #    paste0("\n", niter.brown.g, " for Brownian motion"), 
    #    paste0("\n", niter.lambda.g, " for lambda"), 
    #    paste0("\n", niter.delta.g, " for delta"), 
    #    paste0("\n", niter.OU.g, " for OU"), 
    #    paste0("\n", niter.white.g, " for white noise")))
    #
    trueFreeValuesANDSummaryValues <- parallelSimulateWithPriors(
        nrepSim = nrepSim, coreLimit = coreLimit, phy = phy, taxonDF = taxonDF, 
        startingPriorsValues = startingPriorsValues,
        intrinsicPriorsValues = intrinsicPriorsValues,
        extrinsicPriorsValues = extrinsicPriorsValues, 
        startingPriorsFns = startingPriorsFns,
        intrinsicPriorsFns = intrinsicPriorsFns,
        extrinsicPriorsFns = extrinsicPriorsFns, 
        freevector = freevector, timeStep = timeStep,
        intrinsicFn = intrinsicFn, extrinsicFn = extrinsicFn, 
        multicore = multicore,
        checkpointFile = checkpointFile,
        checkpointFreq = checkpointFreq, 
        )
    #
    #message("\n\n")
    simTime <- proc.time()[[3]]-startTime
    message(paste0("Initial simulations took ", round(simTime, digits = 3), " seconds")) #, "\n"

    #separate the simulation results: true values and the summary values
    trueFreeValuesMatrix <- trueFreeValuesANDSummaryValues[, 1:numberParametersFree]
    summaryValuesMatrix <- trueFreeValuesANDSummaryValues[, -1:-numberParametersFree]

    if (savesims){
        save(trueFreeValuesMatrix, summaryValuesMatrix, simTime, file = paste0("sims", jobName, ".Rdata", sep = ""))
        }

    res <- makeQuiet(PLSRejection(
        summaryValuesMatrix = summaryValuesMatrix, 
        trueFreeValuesMatrix = trueFreeValuesMatrix, 
        phy = phy, 
        traits = traits, 
        abcTolerance = abcTolerance
        ))
    # number of accepted simulations
    nAcceptedSims<-nrow(res$particleDataFrame)
    # print acceptance rate
    message(paste0("Of ",nInitialSims,"initial simulations, ",nAcceptedSims,
        "were accepted (",nAcceptedSims/nInitialSims," Particle Acceptance Rate)"))
    #
    rejectionResults <- list()
    # save input data, conditions
    Ntips <- Ntip(phy)
    input.data <- rbind(jobName, Ntips,
        generation.time, TreeYears, timeStep, totalGenerations,
        nInitialSims, nAcceptedSims, standardDevFactor, abcTolerance)
    #
    rejectionResults$input.data <- input.data
    rejectionResults$priorList <- priorList
    rejectionResults$phy <- phy
    rejectionResults$traits <- traits
    rejectionResults$trueFreeValuesANDSummaryValues <- trueFreeValuesANDSummaryValues
    rejectionResults$SimTime <- simTime
    rejectionResults$abcDistances <- res$abcDistances
    rejectionResults$particleDataFrame <- res$particleDataFrame
    #rejectionResults$credibleInt <- credibleInt(res$particleDataFrame)
    if(nAcceptedSims>2){
        rejectionResults$postSummary  <- summarizePosterior(
			res$particleDataFrame, 
			verboseMultimodal = FALSE,
			stopIfFlat = FALSE
			)
    }else{
        warning(
			"Posterior Summaries were not calculated as the number of accepted particles was less than 2")
        }
    #
    return(rejectionResults)
    }

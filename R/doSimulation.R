#' Discrete-Time Character Simulation
#' 
#' The function \code{doSimulation} evolves continuous characters under a discrete time process.
#' These functions are mainly used as internal components, generating simulations
#' within ABC analyses using the \code{\link{doRun}} functions. See \emph{Note} below.
#' 
#' The phylogenetic tree used is rescaled such that the
#' distance from the root to the furthest tip is rescaled to equal 1 time-unit, 
#' and it is this rescaled edge lengths to with arguments
#' \code{timeStep} refers to. Typically, this will be determined though
#' as a ratio of \code{TreeYears} (which is the number
#' of calender years constituing the root-to-furthest-tip distance, and
#' is determined by default as if the user had supplied
#' a tree with edge lengths in time-units of 1 time-unit = 1 million years), and
#' \code{generation.time}, which gives the length of
#' \code{timeSteps} in calender years (e.g. \code{generation.time = 1000} means
#' an evolutionary change in trait values every 1000 years).
#' Note that the real number of trait change events simulated may be less
#' less, because simulations start ver at each branching
#' node, but intervals between changes should be so fine that this
#' should be negligible (related, the results should be
#' independent of your choice for \code{generation.time} or \code{timeStep}).
#' We recommend that the effective \code{timeStep} should be 
#' as short as is computationally possible.
# 
# When \code{saveHistory} is \code{TRUE}, processor time will increase quite a bit.
#
#' \code{SaveRealParams} is useful for tracking the \emph{real} true values if simulating
#' data to test the performance of ABC analyses.  
#' It is not useful for ABC analyses of empirical data.
#' 

#' @note
#' The \code{\link{simulateWithPriors}} functions are effectively
#' the engine that powers the \code{\link{doRun}}
#' functions, while the \code{\link{doSimulation}} functions are
#' the pistons within the \code{\link{simulateWithPriors}} engine.
#' In general, most users will just drive the car - they will just
#' use \code{\link{doRun}}, but some users may
#' want to use \code{\link{simulateWithPriors}} or
#' \code{\link{doSimulation}} functions to do various simulations.


#' @param phy A phylogenetic tree, in package \code{ape}'s \code{phylo} format.

#' @param taxonDF A data.frame containing data on nodes (both
#' tips and internal nodes) output by various internal functions.
#' Can be supplied as input to spead up repeated calculations, but by default is
#' \code{NULL}, which instead forces a calculation from input \code{phy}.

#' @param intrinsicFn Name of (previously-defined) function that governs how
#' traits evolve within a lineage, regardless of the trait values of other taxa.

#' @param extrinsicFn Name of (previously-defined) function that governs how
#' traits evolve within a lineage, based on their own
#' ('internal') trait vlaue and the trait values
#' of other taxa.

#' @param startingValues State at the root.

#' @param intrinsicValues Vector of values corresponding to the parameters of the
#' intrinsic model.

#' @param extrinsicValues Vector of values corresponding to the parameters of the
#' extrinsic model.

#' @param generation.time The number of years per generation.
#' This sets the coarseness of the simulation; if it's set to 1000, 
#' for example, the population's trait values change every 1000
#' calender years. Note that this is in calender years (see description
#' for argument \code{TreeYears}), and not in millions of years
#' (as is typical for dated trees in macroevolutionary studies).
#' Thus, if a branch is 1 million-year time-unit long, and a
#' user applies the default \code{generation.time = 1000}, 
#' then 1000 evolutionary changes will be simulated along that branch.
#' See documentation for \code{\link{doSimulation}} for further details.

#' @param TreeYears The amount of calender time from the root
#' to the furthest tip. Most trees in macroevolutionary studies are dated with
#' branch lengths in units of millions of years, and thus the
#' default for this argument is \code{max(branching.times(phy)) * 1e6}.
#' If your tree has the most recent tip at time zero (i.e.,
#' the modern day), this would be the same as the root age of the tree. If your
#' branch lengths are not in millions of years, you should
#' alter this argument. Otherwise, leave this argument alone.
#' See documentation for \code{\link{doSimulation}} for further details.

#' @param timeStep This value corresponds to the length
#' of intervals between discrete evolutionary events ('generations')
#' simulated along branches, relative to a rescaled tree
#' where the root to furthest tip distance is 1. For example, 
#' \code{timeStep = 0.01} of would mean 100 (i.e., 1 / 0.01)
#' evolutionary changes would be expected to occur from
#' the root to the furthest tip. (Note that the real number
#' simulated will be much less, because simulations start
#' over at each branching node.) Ideally, \code{timeStep}
#' (or its effective value, via other arguments) should be
#' as short as is computationally possible.
#' Typically \code{NULL} by default and
#' determined internally as follows: \code{timeStep = generation.time / TreeYears}.
#' Can be provided a value as an alternative to using arguments \code{generation.time}
#' and \code{TreeYears}, which would then be overridden.
#' See documentation for \code{\link{doSimulation}} for further details.

#' @param checkTimeStep If \code{TRUE}, warnings will be issued if \code{TimeStep} is too short.

#' @param returnAll If \code{TRUE}, the output returned is a
#' \code{data.frame} containing the values at each node from the simulation.


# @param maxAttempts How many attempts should be tried if a
# run produces an \code{NA} result? If \code{maxAttempts}
# is reached without producing a non-\code{NA} result, the simulation is terminated.

# @param saveHistory If \code{TRUE}, saves the character history throughout the simulation.
# When \code{saveHistory} is \code{TRUE}, processor time will increase quite a bit.

# @param saveRealParams Saves \code{intrinsicValues} and \code{extrinsicValues} as both a
# matrix and a vector to an external .Rdata file.

# @param jobName Optional name for the job.

# @param verbose If \code{TRUE}, gives messages about
# how the simulation is progessing via \code{message}.

# @param reject.NaN If \code{TRUE}, stop run if any simulated value is \code{NaN}.

# @param plot Will create a new interactive window that plots character values
# throughout the history of the tree.

# @param savePlot Saves the character tree using \code{jobName}.


#' @return If \code{returnAll = FALSE} (the default),
#' this function returns a data frame of species character (tip)
#' values in the tree, as a single column named 'states'
#' with rownames reflecting the taxon names given in 
#' \code{phy$tip.label}.
#' If \code{returnAll = TRUE}, the raw \code{data.frame}
#' from the simulation will instead be returned.

# OLD RETURN INFO
# column headings \code{taxonid}
# (representing the index for the corresponding tip label
# for that taxon, as given in \code{phy$tip.label)),
# \code{taxonname}, \code{taxontimesincespeciation}
# (the time to the most recent divergence event for that lineage),
# and \code{statesmatrix} (the simulated trait data).


#' @author Brian O'Meara and Barb Banbury

#' @examples
#' 
#' \donttest{
#' 
#' set.seed(1)
#' tree <- rcoal(20)
#' # get realistic edge lengths
#' tree$edge.length <- tree$edge.length*20
#' 
#' #Simple Brownian motion
#'
#' char <- doSimulation(
#'     phy = tree, 
#'     generation.time = 100000, 
#'     intrinsicFn = brownianIntrinsic, 
#'     extrinsicFn = nullExtrinsic, 
#'     startingValues = c(10), #root state
#'     intrinsicValues = c(0.01), 
#'     extrinsicValues = c(0))
#' 
#' #Character displacement model with minimum bound
#' 
#' char <- doSimulation(
#'     phy = tree, 
#'     generation.time = 100000, 
#'     intrinsicFn = boundaryMinIntrinsic, 
#'     extrinsicFn = ExponentiallyDecayingPushExtrinsic, 
#'     startingValues = c(10), #root state
#'     intrinsicValues = c(0.05, 10, 0.01), 
#'     extrinsicValues = c(0, .1, .25))
#' 
#' }

#' @rdname doSimulation
#' @export
doSimulation <- function(
	phy = NULL, 
	intrinsicFn, extrinsicFn, 
	startingValues, intrinsicValues, extrinsicValues, 
    generation.time = 1000, 
	TreeYears = max(branching.times(phy)) * 1e6, 
	timeStep = NULL, 
	returnAll = FALSE, 
	taxonDF = NULL, 
	checkTimeStep = TRUE
    #maxAttempts = 100, saveHistory = FALSE, plot = FALSE, savePlot = FALSE, 
    #reject.NaN = TRUE, saveRealParams = FALSE, jobName = "", verbose = FALSE, 
	) {
    #
    #
    if(is.null(timeStep)){
        timeStep <- generation.time/TreeYears
        }
    #
    if(is.null(taxonDF)){
		if(is.null(phy$edge.length)){
			stop("Input phylogeny lacks necessary edge lengths")
			}
		#
        taxonDF <- getTaxonDFWithPossibleExtinction(phy)
        }
    #
    if(is.null(taxonDF) & is.null(phy)){
        stop("phy or taxonDF must be provided as input")
        }
    #
    #if (saveRealParams){
    #    RealParams <- vector("list", 2)
    #    names(RealParams) <- c("matrix", "vector")
    #    RealParams$vector <- c(startingValues, intrinsicValues, extrinsicValues)
    #    maxLength <- (max(length(startingValues), length(intrinsicValues), length(extrinsicValues)))
    #    RealParams$matrix <- matrix(ncol = maxLength, nrow = 3)
    #    rownames(RealParams$matrix) <- c("startingValues", "intrinsicFn", "extrinsicFn")
    #    RealParams$matrix[1, ] <- c(startingValues, rep(NA, maxLength-length(startingValues)))
    #    RealParams$matrix[2, ] <- c(intrinsicValues, rep(NA, maxLength-length(intrinsicValues)))
    #    RealParams$matrix[3, ] <- c(extrinsicValues, rep(NA, maxLength-length(extrinsicValues)))
    #    save(RealParams, file = paste0("RealParams", jobName, ".Rdata", sep = ""))
    #    }
    #
    #if (plot || savePlot || saveHistory) {
    #    startVector <- c()
    #    endVector <- c()
    #    startTime <- c()
    #    endTime <- c()
    #    }
    #
    if(checkTimeStep){
        #
        numberofsteps <- max(taxonDF$endTime)/timeStep
        mininterval <- min(taxonDF$endTime - taxonDF$startTime)
        #if (numberofsteps<1000) {
            #warning(paste0("You have only ", numberofsteps, " but should probably have a lot more. Please consider 
            #decreasing timeStep to no more than ", taxonDF[1, 1]/1000))
        #    }
        if (floor(mininterval/timeStep)<50 & floor(mininterval/timeStep) >= 3) {
            warning(paste0("You have only ", floor(mininterval/timeStep), 
                " timeSteps on the shortest branch in this dataset\n",
				" but should probably have a lot more if you expect change on this branch.\n",
				" Please consider decreasing timeStep to no more than ", 
                signif(mininterval/50, 2)))
            }
        if (floor(mininterval/timeStep)<3) {
            warning(paste0("You have only ", floor(mininterval/timeStep), 
                " timeSteps on the shortest branch in this dataset\n",
				" but should probably have a lot more if you expect change on this branch.\n",
				" Please consider decreasing timeStep to no more than ", 
                    signif(mininterval/50, 2), 
				" or at the very least ", signif(mininterval/3, 2)))
        #    timeStep <- mininterval/3
            }
        }
    #
    final.result.df <- doSimulationInternal(
        taxonDF = taxonDF, timeStep = timeStep, 
        intrinsicFn = intrinsicFn, extrinsicFn = extrinsicFn, 
        startingValues = startingValues, 
		intrinsicValues = intrinsicValues, extrinsicValues = extrinsicValues
		)
    
    #if(plot){
    #    #dev.new()
    #    plot(x = c(min(c(startVector, endVector)), 
    #            max(c(startVector, endVector))), 
    #         y = c(0, max(c(startTime, endTime))), 
    #         type = "n", 
    #         ylab = "Time", xlab = "Trait value", main = "", bty = "n")
    #    for (i in 1:length(startVector)){
    #        lines(x = c(startVector[i], endVector[i]), 
    #              y = max(c(startTime, endTime)) - c(startTime[i], endTime[i])
    #              )
    #        }
    #    }
    #if (savePlot) {
    #    pdf(paste0("SimTree", jobName, ".pdf", sep = ""))    
    #    plot(x = c(min(c(startVector, endVector)), max(c(startVector, endVector))), 
    #         y = c(0, max(c(startTime, endTime))), 
    #         type = "n", ylab = "Time", xlab = "Trait value", 
    #         main = "", bty = "n")
    #    for (i in 1:length(startVector)) {
    #        lines(x = c(startVector[i], endVector[i]), 
    #        y = max(c(startTime, endTime)) - c(startTime[i], endTime[i]))
    #        }
    #    dev.off()
    #    }

    return(final.result.df)
    }
    
    
# 
doSimulationInternal <- function(
    taxonDF, timeStep, intrinsicFn, extrinsicFn, 
    startingValues, intrinsicValues, extrinsicValues){
    #
    # just the meat of doSimulation, no checks, no nothing
    #
    ####################################################################
    #initial setup
    #
    numberofsteps <- max(taxonDF$endTime)/timeStep
    mininterval <- min(taxonDF$endTime - taxonDF$startTime)
    #
    depthfrompresent = max(taxonDF$endTime)
    heightfromroot = 0
    # which element of DF is states
    whichStatesCol <- which(names(taxonDF) == "states")
    #
    taxonDF[[whichStatesCol]][which(taxonDF$startTime == 0)] <- startingValues
    #    
    taxonDF.previous <- taxonDF
    #
    taxonID <- taxonDF$id
    taxonStartTime <- taxonDF$startTime
    taxonEndTime <- taxonDF$endTime
    taxonTerminal <- taxonDF$terminal
    taxonAnc <- taxonDF$ancestorId
    while(depthfrompresent>0) {
        #if(reject.NaN) {
        #    taxonDF.previous <- taxonDF
        #    }
        #
        depth.start <- depthfrompresent
        depth.end <- depthfrompresent - timeStep
        height.start <- heightfromroot
        height.end <- heightfromroot + timeStep
        ids.alive.at.start <- taxonID[which(taxonStartTime  <=  height.start & taxonEndTime > height.start)]
        ids.alive.at.end <-  taxonID[which(taxonEndTime > height.end & taxonStartTime  <=  height.end)]
        ids.only.alive.in.interval <- taxonID[which(taxonStartTime > height.start & taxonEndTime < height.end)]
        ids.changing.status <-  c(ids.alive.at.start[!(ids.alive.at.start  %fin% ids.alive.at.end)], ids.only.alive.in.interval)
        ids.speciating <- c(taxonID[which((taxonID %fin% ids.changing.status) & (!taxonTerminal))], ids.only.alive.in.interval)
        #
        aliveRows <- which(taxonID %fin% ids.alive.at.start)
        #
        #if(any(is.na(aliveRows))){
        #    stop("some aliveRows are NA")
        #    }
        #
        taxonStates <- taxonDF[[whichStatesCol]]
        currentStates <- taxonStates[aliveRows]
        #
        #if(any(is.na(currentStates))) {
        #    message(paste0("currentStates ", currentStates))
        #    message(paste0("taxonID %fin% ids.alive.at.start ", paste0(taxonID %fin% ids.alive.at.start, collapse = " ")))
        #    message(c(height.start, height.end))
        #    message(taxonDF[taxonID %fin% ids.alive.at.start, ])
        #    stop("there are NAs in currentStates! How?? Something is very wrong")
        #    }
        #
        # get vector of ancestors
        ancestors <- match(taxonAnc[aliveRows], taxonID)    
        #
        #first evolve in this interval, then speciate
        for (whichTaxon in 1:length(aliveRows)) {
            # find match within aliveRows for match to currentStates
            taxonIndex <- aliveRows[whichTaxon]
            taxonState <- taxonStates[taxonIndex]
            # check if the ancestor is NA
            if(is.na(taxonState)) {
                #whichAncestor <-  which(taxonID == taxonAnc[taxonIndex])
                taxonState <- taxonStates[ancestors[whichTaxon]]
                if(is.na(taxonState)){
                    stop("A taxon's ancestor has an NA character")
                    }
                currentStates[whichTaxon] <- taxonState
                }                
            #
            newState <- taxonState + 
                intrinsicFn(params = intrinsicValues, states = currentStates[whichTaxon], 
                    timefrompresent = depthfrompresent) + 
                extrinsicFn(params = extrinsicValues, selfstates = currentStates[whichTaxon], 
                    otherstates = currentStates[-whichTaxon], timefrompresent = depthfrompresent)
            #
            ##
            ##if(is.na(newState)) {    # what happens if I change this to a stop? oh that's no good
            ##    stop("A simulation run produced a state of NA - something is probably very wrong")
            ##    if(any(is.na(currentStates))) {
            ##        stop(paste0("there are NAs in currentStates! How?? Something is very wrong\n", 
            ##            "currentStates ", currentStates))
            ##        }
            ##    attempt.count <- 0
            ##    while(is.na(newState) & attempt.count  <=  maxAttempts) {
            ##        old = taxonState
            ##        #check
            ##        if(is.na(old)){
            ##            stop("Ancestral state for a simulated character state is NA - something is very wrong")
            ##            }
            ##        intrinsic.displacement = intrinsicFn(params = intrinsicValues, states = currentStates[taxonIndex], 
            ##            timefrompresent  = depthfrompresent)
            ##        #check
            ##        if(is.na(intrinsic.displacement)){
            ##            stop("The intrinsicFn is returning NAs; something terrible has happened")
            ##            }
            ##        extrinsic.displacement = extrinsicFn(params = extrinsicValues, selfstates = currentStates[taxonIndex], 
            ##            otherstates = currentStates[-taxonIndex], timefrompresent  = depthfrompresent)
            ##        #check
            ##        if(is.na(extrinsic.displacement)){
            ##            stop(
            ##                paste("The extrinsicFn is returning NAs; something terrible has happened:", 
            ##                    "\ntaxonIndex", taxonIndex, 
            ##                    "\naliveRows", aliveRows, 
            ##                    "\ncurrentStates", currentStates, 
            ##                    "\nparams:", extrinsicValues, 
            ##                    "\nselfstates:", currentStates[taxonIndex], 
            ##                    "\notherstates:", currentStates[-taxonIndex], 
            ##                    "\ntimefrompresent:", depthfrompresent
            ##                    #, "\nextFun:", extrinsicFn
            ##                    )
            ##                )
            ##            }
            ##        newState <- old + intrinsic.displacement + extrinsic.displacement
            ##        warning(paste0("Attempt ", attempt.count, " led to using old value of ", old, " intrinsicFn return of ", intrinsic.displacement, " and extrinsicFn return of ", extrinsic.displacement))
            ##        #message(paste0("For diagnostic purposes: IntrinsicValues ", intrinsicValues))
            ##        attempt.count <- attempt.count+1
            ##        }
            ##    if(is.na(newState) & attempt.count>maxAttempts) {
            ##        if(is.na(extrinsic.displacement)){
            ##            message(paste0(ls(), collapse = ", "))
            ##            #message(str(aliveRows))
            ##            message(paste0("taxonIndex ", taxonIndex, "\n", 
            ##                            "aliveRows ", paste0(aliveRows, collapse = ", "), "\n", 
            ##                            "length(aliveRows) ", length(aliveRows), "\n", 
            ##                            "1:length(aliveRows)", paste(1:length(aliveRows)), collapse = ", "), "\n", 
            ##                            "currentStates ", paste(currentStates, collapse = ", "), "\n", 
            ##                            "params ", extrinsicValues, "\n", 
            ##                            "selfstates ", currentStates[taxonIndex], "\n", 
            ##                            "otherstates ", paste(
            ##                                currentStates[-taxonIndex], collapse = " "), "\n", 
            ##                            "timefrompresent ", depthfrompresent, "\n"))
            ##            }
            ##        stop(paste0(
            ##            "Simulating with these parameters resulted in problematic results despite ", maxAttempts, " attempts", 
            ##            "\nFor one example, taxonDF$states[taxonIndex] was ", 
            ##            taxonStates[taxonIndex], ", for which intrinsicFn returned ", 
            ##            intrinsicFn(params = intrinsicValues, states = currentStates[taxonIndex], 
            ##                timefrompresent  = depthfrompresent)
            ##            , "\nand extrinsicFn returned ", 
            ##            extrinsicFn(params = extrinsicValues, 
            ##                selfstates = currentStates[taxonIndex], otherstates = currentStates[-taxonIndex], 
            ##                timefrompresent  = depthfrompresent
            ##                )
            ##            , " with currentStates[taxonIndex] = ", currentStates[taxonIndex])
            ##            )
            ##        }
            ##    if(is.na(newState)) {
            ##        stop("where are these NA newStates coming from?? Something is very wrong")
            ##        }                    
            ##    }
            #
            #    if (plot || savePlot || saveHistory) {
            #        startVector <- append(startVector, taxa[[i]]$states)
            #        endVector  <- append(endVector, newvalues)
            #        startTime  <- append(startTime, timefrompresent+timeStep)
            #        endTime  <- append(endTime, timefrompresent)
            #        if(saveHistory){
            #            save(startVector, endVector, startTime, endTime, file = paste0("savedHistory", jobName, ".Rdata", sep = ""))
            #        }
            #    }
            #
            #
            taxonStates[taxonIndex] <- newState
            }
        #
        # now speciate and pass one state to descendant
        if(length(ids.speciating)>0) {
            for (speciating.taxonIndex in 1:length(ids.speciating)) {
                ancestor.row <- which(taxonID == ids.speciating[speciating.taxonIndex])
                descendant.rows <- which(taxonAnc == taxonID[ancestor.row])
                taxonStates[descendant.rows] <- taxonStates[ancestor.row]
                }
            }    
        #
        depthfrompresent <- depth.end
        heightfromroot <- height.end
        #
        if(any(is.nan(taxonStates[aliveRows])) | any(is.na(taxonStates[aliveRows]))) {
                stop(paste0("A simulation run produced a state of NA or NaN - something is probably very wrong\n", 
                    "Here are the states for currently living taxa: ",     paste(taxonStates[aliveRows], collapse = " ")))
                }
        #                
        #if(verbose) {
        #    message(paste0("now at height", height.end, "finishing at", max(taxonEndTime)))
        #    #message(taxonDF)
        #    }
        #if(reject.NaN) {
        #    if(any(is.nan(taxonStates))) {
        #        save(list = ls(), file = "ErrorRun.rda")
        #        stop(paste0("There was an NaN generated. See saved objects in ", getwd(), "/ErrorRun.rda", sep = ""))
        #        }
        #    }
        #
        # finally update the dataframe
        taxonDF[[whichStatesCol]] <- taxonStates
        }
    #if(returnAll) {
    #    return(taxonDF)
    #    }
    final.result <- subset(taxonDF, taxonTerminal == TRUE)
    final.result.df <- data.frame(states = final.result[[whichStatesCol]])
    rownames(final.result.df) <- final.result$name
    return(final.result.df)
    }    





#   Create a data frame of taxon states
#
#   This function creates a data frame of taxon states while simulating
#   characters with doSimulation and doSimulationsForPlotting TreEvo functions
#
#   Used by TreEvo doSimulation and doSimulation functions to
#   summarize a list of objects into a data frame of taxon values
#
#   @param taxa a list of objects
#   @return Returns a data frame of taxon values
#   @author Brian O'Meara and Barb Banbury
# @references O'Meara and Banbury, unpublished

# @name summarizeTaxonStates
# @rdname summarizeTaxonStates
# @export
summarizeTaxonStates <- function(taxa) {
#message("in summarizeTaxonStates")
    statesvector <- c()
    taxonid <- c()
    taxonname <- c()
    taxontimesincespeciation <- c()
    for (i in 1:length(taxa)) {
        statesvector <- c(statesvector, taxa[[i]]$states)
        taxonid <- c(taxonid, taxa[[i]]$id )
        taxonname <- c(taxonname, taxa[[i]]$name )
        taxontimesincespeciation <- c(taxontimesincespeciation, taxa[[i]]$timeSinceSpeciation)
    }
    statesmatrix <- matrix(statesvector, ncol = length(taxa[[1]]$states), byrow = TRUE) #each row represents one taxon
    taxonframe <- data.frame(taxonid, taxonname, taxontimesincespeciation, statesmatrix)
    return(taxonframe)
}


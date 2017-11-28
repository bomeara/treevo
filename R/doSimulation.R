#' Discrete-Time Character Simulation
#'
#' The \code{doSimulation} family of functions evolve continuous characters under a discrete time process.
#' These functions are mainly used as internal components, generating simulations
#' within ABC analyses using the \code{\link{doRun}} functions. See \emph{Note} below.
#'
#' The phylogenetic tree used is rescaled such that the distance from the root to the furthest tip is rescaled to equal 1 time-unit,
#' and it is this rescaled edge lengths to with arguments \code{timeStep} refers to. Typically, this will be determined though
#' as a ratio of \code{TreeYears} (which is the number of calender years constituing the root-to-furthest-tip distance, and
#' is determined by default as if the user had supplied a tree with edge lengths in time-units of 1 time-unit = 1 million years), and
#' \code{generation.time}, which gives the length of \code{timeSteps} in calender years (e.g. \code{generation.time = 1000} means
#' an evolutionary change in trait values every 1000 years). Note that the real number of trait change events simulated may be less
#' less, because simulations start ver at each branching node, but intervals between changes should be so fine that this
#' should be negligible (related, the results should be independent of your choice for \code{generation.time} or \code{timeStep}).
#' We recommend that the effective \code{timeStep} should be  as short as is computationally possible.

#' 
#' When \code{saveHistory} is \code{TRUE}, processor time will increase quite a bit.
#' \code{SaveRealParams} is useful for tracking the \emph{real} true values if simulating
#' data to test the performance of ABC analyses.  It is not useful for ABC analyses of empirical data.
#'

#' @note 
#' The \code{\link{simulateWithPriors}} functions are effectively the engine that powers the \code{\link{doRun}}
#' functions, while the \code{\link{doSimulation}} functions are the pistons within the \code{\link{simulateWithPriors}} engine.
#' In general, most users will just drive the car - they will just use \code{\link{doRun}}, but some users may
#' want to use \code{\link{simulateWithPriors}} or \code{\link{doSimulation}} functions to do various simulations.


#' @param phy A phylogenetic tree, in package \code{ape}'s \code{phylo} format.

#' @param taxon.df A data.frame containing data on nodes (both tips and internal nodes) output by various internal functions.
#' Can be supplied as input to spead up repeated calculations, but by default is
#' \code{NULL}, which instead forces a calculation from input \code{phy}.

#' @param intrinsicFn Name of (previously-defined) function that governs how
#' traits evolve within a lineage, regardless of the trait values of other taxa.

#' @param extrinsicFn Name of (previously-defined) function that governs how
#' traits evolve within a lineage, based on their own ('internal') trait vlaue and the trait values
#' of other taxa.

#' @param startingValues State at the root.

#' @param intrinsicValues Vector of values corresponding to the parameters of the
#' intrinsic model.

#' @param extrinsicValues Vector of values corresponding to the parameters of the
#' extrinsic model.

#' @param maxAttempts How many attempts should be tried if a run produces an \code{NA} result? If \code{maxAttempts}
#' is reached without producing a non-\code{NA} result, the simulation is terminated.

#' @param generation.time The number of years per generation. This sets the coarseness of the simulation; if it's set to 1000, 
#' for example, the population's trait values change every 1000 calender years. Note that this is in calender years (see description
#' for argument \code{TreeYears}), and not in millions of years (as is typical for dated trees in macroevolutionary studies).
#' Thus, if a branch is 1 million-year time-unit long, and a user applies the default \code{generation.time = 1000},
#' then 1000 evolutionary changes will be simulated along that branch.
#' See documentation for \code{\link{doSimulation}} for further details.

#' @param TreeYears The amount of calender time from the root to the furthest tip. Most trees in macroevolutionary studies are dated with
#' branch lengths in units of millions of years, and thus the default for this argument is \code{max(branching.times(phy)) * 1e6}.
#' If your tree has the most recent tip at time zero (i.e., the modern day), this would be the same as the root age of the tree. If your
#' branch lengths are not in millions of years, you should alter this argument. Otherwise, leave this argument alone.
#' See documentation for \code{\link{doSimulation}} for further details.

#' @param timeStep This value corresponds to the lenght of intervals between discrete evolutionary events ('generations')
#' simulated along branches, relative to a rescaled tree where the root to furthest tip distance is 1. For example,
#' \code{timeStep = 0.01} of would mean 100 (i.e., 1 / 0.01) evolutionary changes would be expected to occur from
#' the root to the furthest tip. (Note that the real number simulated will be much less, because simulations start
#' over at each branching node.) Ideally, \code{timeStep} (or its effective value, via other arguments) should be 
#' as short as is computationally possible.
#' Typically \code{NULL} by default and
#' determined internally as follows: \code{timeStep = generation.time / TreeYears}.
#' Can be provided a value as an alternative to using arguments \code{generation.time}
#' and \code{TreeYears}, which would then be overridden. 
#' See documentation for \code{\link{doSimulation}} for further details.

#' @param checkTimeStep If \code{TRUE}, warnings will be issued if \code{TimeStep} is too short.

#' @param saveHistory If \code{TRUE}, saves the character history throughout the simulation.
#' When \code{saveHistory} is \code{TRUE}, processor time will increase quite a bit.

#' @param saveRealParams Saves \code{intrinsicValues} and \code{extrinsicValues} as both a
#' matrix and a vector to an external .Rdata file.

#' @param jobName Optional name for the job.

#' @param returnAll If \code{TRUE}, the output returned is a \code{data.frame} containing the values at each node from the simulation.

#' @param verbose If \code{TRUE}, gives messages about how the simulation is progessing via \code{print}.

#' @param reject.NaN If \code{TRUE}, stop run if any simulated value is \code{NaN}.

#' @param plot Will create a new interactive window that plots character values
#' throughout the history of the tree.

#' @param savePlot Saves the character tree using \code{jobName}.



#' @return If \code{returnAll = FALSE} (the default), this function returns a data frame of species character (tip)
#' values in the tree, with column headings \code{taxonid} (representing the index for the corresponding tip label
#" for that taxon, as given in \code{phy$tip.label)), \code{taxonname}, \code{taxontimesincespeciation} (the time to
#' the most recent divergence event for that lineage), and \code{statesmatrix} (the simulated trait data).
#' If \code{returnAll = TRUE}, the raw \code{data.frame} from the simulation will instead be returned.

#' @author Brian O'Meara and Barb Banbury

#' @author Brian O'Meara and Barb Banbury

# @references O'Meara and Banbury, unpublished
# @keywords doSimulation

#' @examples
# 
#' \donttest{
#
#' tree<-rcoal(30)
#'
#' #Simple Brownian motion
#' char<-doSimulation(
#' 	phy=tree,
#' 	intrinsicFn=brownianIntrinsic,
#' 	extrinsicFn=nullExtrinsic,
#' 	startingValues=c(10), #root state
#' 	intrinsicValues=c(0.01),
#' 	extrinsicValues=c(0),
# 	timeStep=0.0001,
#' 	saveHistory=FALSE)
#'
#' #Character displacement model with minimum bound
#' char<-doSimulation(
#' 	phy=tree,
#' 	intrinsicFn=boundaryMinIntrinsic,
#' 	extrinsicFn=ExponentiallyDecayingPushExtrinsic,
#' 	startingValues=c(10), #root state
#' 	intrinsicValues=c(0.05, 10, 0.01),
#' 	extrinsicValues=c(0, .1, .25),
# 	timeStep=0.001,
#' 	saveHistory=FALSE)
#'
#' #Simple Brownian motion
#' char<-doSimulationForPlotting(
#' 	phy=tree,
#' 	intrinsicFn=brownianIntrinsic,
#' 	extrinsicFn=nullExtrinsic,
#' 	startingValues=c(10), #root state
#' 	intrinsicValues=c(0.01),
#' 	extrinsicValues=c(0),
# 	timeStep=0.0001,
#' 	plot=FALSE,
#' 	saveHistory=FALSE)
#' 
#' 
#' #Character displacement model with minimum bound
#' char<-doSimulationForPlotting(
#' 	phy=tree,
#' 	intrinsicFn=boundaryMinIntrinsic,
#' 	extrinsicFn=ExponentiallyDecayingPushExtrinsic,
#' 	startingValues=c(10), #root state
#' 	intrinsicValues=c(0.05, 10, 0.01),
#' 	extrinsicValues=c(0, .1, .25),
#' 	plot=TRUE,
#' 	saveHistory=FALSE)
#' 
#' 
#' # with extinction
#'
#' #Simple Brownian motion
#' char<-doSimulationWithPossibleExtinction(
#' 	phy=tree,
#' 	intrinsicFn=brownianIntrinsic,
#' 	extrinsicFn=nullExtrinsic,
#' 	startingValues=c(10), #root state
#' 	intrinsicValues=c(0.01),
#' 	extrinsicValues=c(0),
#' 	saveHistory=FALSE)
#
#' }


#' @name doSimulation
#' @rdname doSimulation
#' @export
doSimulation<-function(phy=NULL, intrinsicFn, extrinsicFn, startingValues, intrinsicValues, extrinsicValues, 
	generation.time=1000, TreeYears=max(branching.times(phy)) * 1e6, 
	timeStep=NULL, saveHistory=FALSE, saveRealParams=FALSE, jobName="", taxon.df = NULL) {

	if(!is.ultrametric(phy)){
		stop("phy must be ultrametric for function doSimulation")
		}
	
	if(is.null(timeStep)){
		timeStep<-generation.time/TreeYears
		}
	
	if(is.null(taxon.df)){
		taxon.df<-getSimulationSplits(phy)
		}
		
	#
	if (saveRealParams){
		RealParams<-vector("list", 2)
		names(RealParams)<-c("matrix", "vector")
		RealParams$vector<-c(startingValues, intrinsicValues, extrinsicValues)
		maxLength<-(max(length(startingValues), length(intrinsicValues), length(extrinsicValues)))
		RealParams$matrix<-matrix(ncol=maxLength, nrow=3)
		rownames(RealParams$matrix)<-c("startingValues", "intrinsicFn", "extrinsicFn")
		RealParams$matrix[1,]<-c(startingValues, rep(NA, maxLength-length(startingValues)))
		RealParams$matrix[2,]<-c(intrinsicValues, rep(NA, maxLength-length(intrinsicValues)))
		RealParams$matrix[3,]<-c(extrinsicValues, rep(NA, maxLength-length(extrinsicValues)))
		save(RealParams, file=paste0("RealParams", jobName, ".Rdata", sep=""))
	}
	if (saveHistory) {
		startVector<-c()
		endVector<-c()
		startTime<-c()
		endTime<-c()
	}
	numberofsteps<-floor(taxon.df[1, 1]/timeStep)
	mininterval<-min(taxon.df[1:(dim(taxon.df)[1]-1), 1]-taxon.df[2:(dim(taxon.df)[1]), 1])

	#if (numberofsteps<1000) {
		#warning(paste0("You have only ", numberofsteps, " but should probably have a lot more. Please consider decreasing timeStep to no more than ", taxon.df[1, 1]/1000))
	#}
	#if (floor(mininterval/timeStep)<50 {
		#warning(paste0("You have only ", floor(mininterval/timeStep),
		#" timeSteps on the shortest branch in this dataset but should probably have a lot more if you expect change on this branch. Please consider decreasing timeStep to no more than ",
		#	signif(mininterval/50,2)))
	#}

	#initial setup
	timefrompresent=taxon.df[1, 1]
	taxa<-list(abctaxon(id=taxon.df[1, 3], states=startingValues), abctaxon(id=taxon.df[1, 4], states=startingValues))
	taxon.df<-taxon.df[2:dim(taxon.df)[1], ] #pop off top value

	#start running
	while(timefrompresent>0) {
		#print(timefrompresent)
		#speciation if needed
		while ((timefrompresent-timeStep)<=taxon.df[1, 1]) { #do speciation. changed from if to while to deal with effectively polytomies
			originallength<-length(taxa)
			taxontodelete<-Inf
			originallength<-length(taxa)
			taxontodelete<-Inf
			for (i in 1:originallength) { #need to merge this from my branch still -DG
				if (taxa[[i]]$id==taxon.df[1, 2]) {
					taxontodelete<-i
					taxa[[originallength+1]] <- taxa[[i]]
					taxa[[originallength+2]] <- taxa[[i]]
					taxa[[originallength+1]]$id<-taxon.df[1, 3]
					taxa[[originallength+1]]$timeSinceSpeciation<-0
					taxa[[originallength+2]]$id<-taxon.df[1, 4]
					taxa[[originallength+2]]$timeSinceSpeciation<-0
				}
			}
			#message("taxontodelete = ", taxontodelete)
			taxa<-taxa[-1*taxontodelete]
			if(dim(taxon.df)[1]>1) {
				taxon.df<-taxon.df[2:(dim(taxon.df)[1]), ] #pop off top value
			}
			else {
				taxon.df[1, ]<-c(-1, 0, 0, 0)
			}
			#print("------------------- speciation -------------------")
			#print(taxa)
			#summarizeTaxonStates(taxa)
		}
		#trait evolution step
		otherstatefn<-function(x){
			taxa[[x]]$states
		}

		otherMatrix<-function(i){
			taxvec<-c(1:length(taxa))
			taxvec<-taxvec[-which(taxvec==i)]
			#
			# NOTE this step is slow. Figure out way to make it faster. taxa is a list of abctaxon objects, so taxa$state won't work
			otherstatesvector<-sapply(taxvec,otherstatefn) 
			#
			otherstatesmatrix<-matrix(otherstatesvector, ncol=length(taxa[[i]]$states), byrow=TRUE) #each row represents one taxon
			newvalues<-taxa[[i]]$states+intrinsicFn(params=intrinsicValues, states=taxa[[i]]$states, timefrompresent =timefrompresent)+extrinsicFn(params=extrinsicValues, selfstates=taxa[[i]]$states, otherstates=otherstatesmatrix, timefrompresent =timefrompresent)
			taxa[[i]]$nextstates<-newvalues
			if (saveHistory) {
				startVector<-append(startVector, taxa[[i]]$states)
				endVector <-append(endVector, newvalues)
				startTime <-append(startTime, timefrompresent+timeStep)
				endTime <-append(endTime, timefrompresent)
				save(startVector, endVector, startTime, endTime, file=paste0("savedHistory", jobName, ".Rdata", sep=""))
			}


			return(taxa[[i]])
		}
		taxvec<-c(1:length(taxa))
		taxa<-lapply(taxvec,otherMatrix)

		stateNextState<-function(i){
			i$states<-i$nextstates
			return(i)

		}
		taxa<-lapply(taxa,stateNextState)
		#print("------------------- step -------------------")
		#print(taxa)
		#summarizeTaxonStates(taxa)

		timefrompresent<-timefrompresent-timeStep
		timeSinceSp<-function(i) {
			i$timeSinceSpeciation<-i$timeSinceSpeciation+timeStep
			return(i)
		}
		taxa<-lapply(taxa,timeSinceSp)
	}
	return(summarizeTaxonStates(taxa))
}

#' @rdname doSimulation
#' @export
doSimulationForPlotting<-function(phy=NULL, intrinsicFn, extrinsicFn, startingValues, intrinsicValues, 
	extrinsicValues, 	generation.time=1000, TreeYears=max(branching.times(phy)) * 1e6, 
	timeStep=NULL, plot=FALSE, savePlot=FALSE, saveHistory=FALSE, saveRealParams=FALSE, jobName="", taxon.df=NULL) {

	if(!is.ultrametric(phy)){
		stop("phy must be ultrametric for function doSimulationForPlotting")
		}	
	
	if(is.null(timeStep)){
		timeStep<-generation.time/TreeYears
		}	
	
	if(is.null(taxon.df)){
		taxon.df<-getSimulationSplits(phy)
		}

	if (saveRealParams){
		RealParams<-vector("list", 2)
		names(RealParams)<-c("matrix", "vector")	
		RealParams$vector<-c(startingValues, intrinsicValues, extrinsicValues)
		maxLength<-(max(length(startingValues), length(intrinsicValues), length(extrinsicValues)))
		RealParams$matrix<-matrix(ncol=maxLength, nrow=3)
		rownames(RealParams$matrix)<-c("startingValues", "intrinsicFn", "extrinsicFn")
		RealParams$matrix[1,]<-c(startingValues, rep(NA, maxLength-length(startingValues)))
		RealParams$matrix[2,]<-c(intrinsicValues, rep(NA, maxLength-length(intrinsicValues)))
		RealParams$matrix[3,]<-c(extrinsicValues, rep(NA, maxLength-length(extrinsicValues)))
		save(RealParams, file=paste0("RealParams", jobName, ".Rdata", sep=""))
	}
	if (plot || savePlot || saveHistory) {
		startVector<-c()
		endVector<-c()
		startTime<-c()
		endTime<-c()
	}
		numberofsteps<-floor(taxon.df[1, 1]/timeStep)
		mininterval<-min(taxon.df[1:(dim(taxon.df)[1]-1), 1]-taxon.df[2:(dim(taxon.df)[1]), 1])
		if (numberofsteps<1000) {
			#warning(paste0("You have only ", numberofsteps, " but should probably have a lot more. Please consider decreasing timeStep to no more than ", taxon.df[1, 1]/1000))
		}
		if (floor(mininterval/timeStep)<50) {
			#warning(paste0("You have only ", floor(mininterval/timeStep), " timeSteps on the shortest branch in this dataset but should probably have a lot more. Please consider decreasing timeStep to no more than ", signif(mininterval/50)))
		}
	#initial setup
		timefrompresent=taxon.df[1, 1]
		taxa<-list(abctaxon(id=taxon.df[1, 3], states=startingValues), abctaxon(id=taxon.df[1, 4], states=startingValues))
		taxon.df<-taxon.df[2:dim(taxon.df)[1], ] #pop off top value
		
	#start running
		while(timefrompresent>0) {
	#print(timefrompresent)
	#speciation if needed
			while ((timefrompresent-timeStep)<=taxon.df[1, 1]) { #do speciation. changed from if to while to deal with effectively polytomies
				originallength<-length(taxa)
				taxontodelete<-Inf
				for (i in 1:originallength) {
					if (taxa[[i]]$id==taxon.df[1, 2]) {
						taxontodelete<-i
						taxa[[originallength+1]] <- taxa[[i]]
						taxa[[originallength+2]] <- taxa[[i]]
						taxa[[originallength+1]]$id<-taxon.df[1, 3]
						taxa[[originallength+1]]$timeSinceSpeciation<-0
						taxa[[originallength+2]]$id<-taxon.df[1, 4]
						taxa[[originallength+2]]$timeSinceSpeciation<-0
					}
				}
	#message("taxontodelete = ", taxontodelete)
				taxa<-taxa[-1*taxontodelete]
				if(dim(taxon.df)[1]>1) {
					taxon.df<-taxon.df[2:(dim(taxon.df)[1]), ] #pop off top value
				}
				else {
					taxon.df[1, ]<-c(-1, 0, 0, 0)
				}
	#print("------------------- speciation -------------------")
	#print(taxa)
	#summarizeTaxonStates(taxa)
			}
	#trait evolution step
			for (i in 1:length(taxa)) {
				otherstatesvector<-c()
				for (j in 1:length(taxa)) {
					if (j!=i) {
						otherstatesvector<-c(otherstatesvector, taxa[[j]]$states)
					}
				}
	#print(taxa)
	#print(length(otherstatesvector))
				otherstatesmatrix<-matrix(otherstatesvector, ncol=length(taxa[[i]]$states), byrow=TRUE) #each row represents one taxon
				newvalues<-taxa[[i]]$states+intrinsicFn(params=intrinsicValues, states=taxa[[i]]$states, timefrompresent =timefrompresent)+extrinsicFn(params=extrinsicValues, selfstates=taxa[[i]]$states, otherstates=otherstatesmatrix, timefrompresent =timefrompresent)
				taxa[[i]]$nextstates<-newvalues
				if (plot || savePlot || saveHistory) {
					startVector<-append(startVector, taxa[[i]]$states)
					endVector <-append(endVector, newvalues)
					startTime <-append(startTime, timefrompresent+timeStep)
					endTime <-append(endTime, timefrompresent)
					if(saveHistory){
						save(startVector, endVector, startTime, endTime, file=paste0("savedHistory", jobName, ".Rdata", sep=""))
					}
				}
				
			}
			for (i in 1:length(taxa)) {

				taxa[[i]]$states<-taxa[[i]]$nextstates

				
			}
	#print("------------------- step -------------------")
	#print(taxa)
	#summarizeTaxonStates(taxa)
			
			timefrompresent<-timefrompresent-timeStep
			for (i in 1:length(taxa)) {
				taxa[[i]]$timeSinceSpeciation<-taxa[[i]]$timeSinceSpeciation+timeStep
			}
		}
		if (plot) {
			#dev.new()
			plot(x=c(min(c(startVector, endVector)), max(c(startVector, endVector))), y=c(0, max(c(startTime, endTime))),
				type="n", ylab="Time", xlab="Trait value", main="", bty="n")
			for (i in 1:length(startVector)) {
				lines(x=c(startVector[i], endVector[i]), y=max(c(startTime, endTime)) - c(startTime[i], endTime[i]))
			}
		}
		if (savePlot) {
			pdf(paste0("SimTree", jobName, ".pdf", sep=""))	
			plot(x=c(min(c(startVector, endVector)), max(c(startVector, endVector))), y=c(0, max(c(startTime, endTime))),
				type="n", ylab="Time", xlab="Trait value", main="", bty="n")
			for (i in 1:length(startVector)) {
				lines(x=c(startVector[i], endVector[i]), y=max(c(startTime, endTime)) - c(startTime[i], endTime[i]))
			}
			dev.off()
		}
		return(summarizeTaxonStates(taxa))
	}




#' @rdname doSimulation
#' @export
doSimulationWithPossibleExtinction<-function(phy=NULL, intrinsicFn, extrinsicFn, startingValues, intrinsicValues, extrinsicValues,
	generation.time=1000, TreeYears=max(branching.times(phy)) * 1e6, 
	timeStep=NULL, saveHistory=FALSE, saveRealParams=FALSE, jobName="", maxAttempts = 100,
	returnAll = FALSE, verbose=FALSE, reject.NaN=TRUE, taxon.df=NULL, checkTimeStep=TRUE) {
	#
	
	if(is.null(timeStep)){
		timeStep<-generation.time/TreeYears
		}
	
	if(is.null(taxon.df)){
		taxon.df <- getTaxonDFWithPossibleExtinction(phy)
		}
		
	if(is.null(taxon.df) & is.null(phy)){
		stop("phy or taxon.df must be provided as input")
		}
	#
	if (saveRealParams){
		RealParams<-vector("list", 2)
		names(RealParams)<-c("matrix", "vector")
		RealParams$vector<-c(startingValues, intrinsicValues, extrinsicValues)
		maxLength<-(max(length(startingValues), length(intrinsicValues), length(extrinsicValues)))
		RealParams$matrix<-matrix(ncol=maxLength, nrow=3)
		rownames(RealParams$matrix)<-c("startingValues", "intrinsicFn", "extrinsicFn")
		RealParams$matrix[1,]<-c(startingValues, rep(NA, maxLength-length(startingValues)))
		RealParams$matrix[2,]<-c(intrinsicValues, rep(NA, maxLength-length(intrinsicValues)))
		RealParams$matrix[3,]<-c(extrinsicValues, rep(NA, maxLength-length(extrinsicValues)))
		save(RealParams, file=paste0("RealParams", jobName, ".Rdata", sep=""))
		}
	#
	if (saveHistory) {
		startVector<-c()
		endVector<-c()
		startTime<-c()
		endTime<-c()
		}
	#
	numberofsteps<-max(taxon.df$endTime)/timeStep
	mininterval<-min(taxon.df$endTime - taxon.df$startTime)
	#
	if(checkTimeStep){
		#if (numberofsteps<1000) {
			#warning(paste0("You have only ", numberofsteps, " but should probably have a lot more. Please consider decreasing timeStep to no more than ", taxon.df[1, 1]/1000))
		#	}
		if (floor(mininterval/timeStep)<50 & floor(mininterval/timeStep)>=3) {
			warning(paste0("You have only ", floor(mininterval/timeStep), 
				" timeSteps on the shortest branch in this dataset but should probably have a lot more if you expect change on this branch. Please consider decreasing timeStep to no more than ",
				signif(mininterval/50,2)))
			}
		if (floor(mininterval/timeStep)<3) {
			warning(paste0("You have only ", floor(mininterval/timeStep), 
				" timeSteps on the shortest branch in this dataset but should probably have a lot more if you expect change on this branch. Please consider decreasing timeStep to no more than ", 
					signif(mininterval/50,2)," or at the very least ", signif(mininterval/3,2)))
		#	timeStep <- mininterval/3
			}
		}
	#
	#initial setup
	depthfrompresent = max(taxon.df$endTime)
	heightfromroot = 0
	taxon.df$states[which(taxon.df$startTime==0)] <- startingValues
	taxon.df.previous <- taxon.df
	while(depthfrompresent>0) {
		if(reject.NaN) {
			taxon.df.previous <- taxon.df
			}
		#
		depth.start <- depthfrompresent
		depth.end <- depthfrompresent - timeStep
		height.start <- heightfromroot
		height.end <- heightfromroot + timeStep
		ids.alive.at.start <- taxon.df$id[which(taxon.df$startTime <= height.start & taxon.df$endTime > height.start)]
		ids.alive.at.end <-  taxon.df$id[which(taxon.df$endTime > height.end & taxon.df$startTime <= height.end)]
		ids.only.alive.in.interval <- taxon.df$id[which(taxon.df$startTime > height.start & taxon.df$endTime < height.end)]
		ids.changing.status <-  c(ids.alive.at.start[!(ids.alive.at.start  %in% ids.alive.at.end)], ids.only.alive.in.interval)
		ids.speciating <- c(taxon.df$id[which((taxon.df$id %in% ids.changing.status) & (!taxon.df$terminal))], ids.only.alive.in.interval)
		alive.rows <- which(taxon.df$id %in% ids.alive.at.start)
		current.states <- taxon.df$states[alive.rows]
		#if(any(is.na(current.states))) {
		#	#message(paste0("current.states ",current.states))
		#	#message(paste0("taxon.df$id %in% ids.alive.at.start ",paste0(taxon.df$id %in% ids.alive.at.start,collapse=" ")))
		#	print(c(height.start,height.end))
		#	print(taxon.df[taxon.df$id %in% ids.alive.at.start,])
		#	stop("there are NAs in current.states! How?? Something is very wrong")
		#	}
		#first evolve in this interval, then speciate
		for (taxon.index in sequence(length(alive.rows))) {
			if(is.na(taxon.df$states[alive.rows[taxon.index]])) {
				taxon.df$states[alive.rows[taxon.index]] <- taxon.df$states[which(taxon.df$id==taxon.df$ancestorId[alive.rows[taxon.index]])]
				if(is.na(taxon.df$states[alive.rows[taxon.index]])){
					stop("A taxon's ancestor has an NA character")
					}
				current.states[taxon.index]<-taxon.df$states[alive.rows[taxon.index]]
				}
			#
			new.state <- taxon.df$states[alive.rows[taxon.index]] + intrinsicFn(params=intrinsicValues,
				states=current.states[taxon.index], timefrompresent =depthfrompresent)+extrinsicFn(params=extrinsicValues,
				selfstates=current.states[taxon.index], otherstates=current.states[-taxon.index], timefrompresent =depthfrompresent)
			#
			if(is.na(new.state)) {
				warning("A simulation run produced a state of NA - something is probably very wrong")
				attempt.count=0
				while(is.na(new.state) & attempt.count <= maxAttempts) {
					old = taxon.df$states[alive.rows[taxon.index]]
					intrinsic.displacement = intrinsicFn(params=intrinsicValues, states=current.states[taxon.index],
						timefrompresent =depthfrompresent)
					extrinsic.displacement = extrinsicFn(params=extrinsicValues, selfstates=current.states[taxon.index],
						otherstates=current.states[-taxon.index], timefrompresent =depthfrompresent)
					#if(is.na(intrinsic.displacement)){
					#	stop("The intrinsicFn is returning NAs; something terrible has happened")
					#	}
					#if(is.na(extrinsic.displacement)){
					#	stop("The extrinsicFn is returning NAs; something terrible has happened")
					#	}
					new.state <- old + intrinsic.displacement + extrinsic.displacement
					warning(paste0("Attempt ", attempt.count, " led to using old value of ", old, " intrinsicFn return of ",intrinsic.displacement, " and extrinsicFn return of ", extrinsic.displacement))
					#message(paste0("For diagnostic purposes: IntrinsicValues ",intrinsicValues)) 
					attempt.count<-attempt.count+1
					}
				if(is.na(new.state) & attempt.count>maxAttempts) {
					if(is.na(extrinsic.displacement)){
						message(paste0(ls(),collapse=", "))
						#message(str(alive.rows))
						message(paste0("taxon.index ",taxon.index,"\n",
										"alive.rows ",paste0(alive.rows, collapse=", "),"\n",
										"length(alive.rows) ",length(alive.rows),"\n",
										"sequence(length(alive.rows))", paste(sequence(length(alive.rows)),collapse=", "), "\n",
										"current.states ",paste(current.states,collapse=", "),"\n",
										"params ",extrinsicValues,"\n",
										"selfstates ",current.states[taxon.index],"\n",
										"otherstates ",paste(
											current.states[-taxon.index],collapse=" "),"\n",
										"timefrompresent ",depthfrompresent,"\n"))
						}
					stop(paste0(
						"Simulating with these parameters resulted in problematic results despite ", maxAttempts, " attempts",
						"\nFor one example, taxon.df$states[alive.rows[taxon.index]] was ",
						taxon.df$states[alive.rows[taxon.index]], ", for which intrinsicFn returned ",
						intrinsicFn(params=intrinsicValues, states=current.states[taxon.index], 
							timefrompresent =depthfrompresent)
						, "\nand extrinsicFn returned ", 
						extrinsicFn(params=extrinsicValues, 
							selfstates=current.states[taxon.index], otherstates=current.states[-taxon.index], 
							timefrompresent =depthfrompresent
							)
						," with current.states[taxon.index] = ", current.states[taxon.index])
						)
					}
				}
			if(is.na(new.state)) {
				stop("where are these NA new.states coming from?? Something is very wrong")
				}
			taxon.df$states[alive.rows[taxon.index]] <- new.state
		}
		if(length(ids.speciating)>0) {
			for (speciating.taxon.index in sequence(length(ids.speciating))) {
				ancestor.row <- which(taxon.df$id==ids.speciating[speciating.taxon.index])
				descendant.rows <- which(taxon.df$ancestorId==taxon.df$id[ancestor.row])
				taxon.df$states[descendant.rows] <- taxon.df$states[ancestor.row]
			}
		}
		depthfrompresent <- depth.end
		heightfromroot <- height.end
		if(verbose) {
			print(paste0("now at height", height.end, "finishing at", max(taxon.df$endTime)))
			print(taxon.df)
		}
		if(reject.NaN) {
			if(any(is.nan(taxon.df$states))) {
				save(list=ls(), file="ErrorRun.rda")
				stop(paste0("There was an NaN generated. See saved objects in ", getwd(), "/ErrorRun.rda", sep=""))
			}
		}
	}
	if(returnAll) {
		return(taxon.df)
	}
	final.results <- subset(taxon.df, taxon.df$terminal==TRUE)
	final.result.df <- data.frame(states=final.results$states)
	rownames(final.result.df) <- final.results$name

	return(final.result.df)
}

#   Create a data frame of taxon states
#   
#   This function creates a data frame of taxon states while simulating
#   characters with doSimulation and doSimulationsForPlotting TreEvo functions
#   
#   Used by TreEvo doSimulation and doSimulationForPlotting functions to
#   summarize a list of objects into a data frame of taxon values
#   
#   @param taxa a list of objects
#   @return Returns a data frame of taxon values
#   @author Brian O'Meara and Barb Banbury
# @references O'Meara and Banbury, unpublished

# @name summarizeTaxonStates
# @rdname summarizeTaxonStates
# @export
summarizeTaxonStates<-function(taxa) {
#print("in summarizeTaxonStates")
	statesvector<-c()
	taxonid<-c()
	taxonname<-c()
	taxontimesincespeciation<-c()
	for (i in 1:length(taxa)) {
		statesvector<-c(statesvector, taxa[[i]]$states)
		taxonid<-c(taxonid, taxa[[i]]$id )
		taxonname<-c(taxonname, taxa[[i]]$name )
		taxontimesincespeciation<-c(taxontimesincespeciation, taxa[[i]]$timeSinceSpeciation)
	}
	statesmatrix<-matrix(statesvector, ncol=length(taxa[[1]]$states), byrow=TRUE) #each row represents one taxon
	taxonframe<-data.frame(taxonid, taxonname, taxontimesincespeciation, statesmatrix)
	return(taxonframe)
}


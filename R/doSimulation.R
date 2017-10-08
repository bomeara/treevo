#' Discrete-Time Character Simulation
#'
#' The \code{doSimulation} family of functions evolve continuous characters under a discrete time process.
#' These functions are mainly used as internal components, generating simulations
#' within ABC analyses using the \code{\link{doRun}} functions. See \emph{Note} below.
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


#' @param timeStep This value corresponds to the number of discrete time steps
#' on the shortest branch.

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
#' 
#' \donttest{
#'
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
#' 	timeStep=0.0001,
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
#' 	timeStep=0.001,
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
#' 	timeStep=0.0001,
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
#' 	timeStep=0.001,
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
#' 	timeStep=0.0001,
#' 	saveHistory=FALSE)
#'
#' }


#' @name doSimulation
#' @rdname doSimulation
#' @export
doSimulation<-function(phy=NULL, intrinsicFn, extrinsicFn, startingValues, intrinsicValues, extrinsicValues, 
	timeStep, saveHistory=FALSE, saveRealParams=FALSE, jobName="", taxon.df = NULL) {
	
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
		save(RealParams, file=paste("RealParams", jobName, ".Rdata", sep=""))
	}
	if (saveHistory) {
		startVector<-c()
		endVector<-c()
		startTime<-c()
		endTime<-c()
	}
	numberofsteps<-floor(taxon.df[1, 1]/timeStep)
	mininterval<-min(taxon.df[1:(dim(taxon.df)[1]-1), 1]-taxon.df[2:(dim(taxon.df)[1]), 1])
	if (numberofsteps<1000) {
		#warning(paste("You have only ", numberofsteps, " but should probably have a lot more. I would suggest decreasing timeStep to no more than ", taxon.df[1, 1]/1000))
	}
	if (floor(mininterval/timeStep)<50) {
		#warning(paste("You have only ", floor(mininterval/timeStep), " on the shortest interval but should probably have a lot more if you expect change on this branch. I would suggest decreasing timeStep to no more than ", mininterval/50))
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
			#cat("taxontodelete = ", taxontodelete)
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
				save(startVector, endVector, startTime, endTime, file=paste("savedHistory", jobName, ".Rdata", sep=""))
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
	extrinsicValues, timeStep, plot=FALSE, savePlot=FALSE, saveHistory=FALSE, saveRealParams=FALSE, jobName="", taxon.df=NULL) {

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
		save(RealParams, file=paste("RealParams", jobName, ".Rdata", sep=""))
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
			#warning(paste("You have only ", numberofsteps, " but should probably have a lot more. I would suggest decreasing timeStep to no more than ", taxon.df[1, 1]/1000))
		}
		if (floor(mininterval/timeStep)<50) {
			#warning(paste("You have only ", floor(mininterval/timeStep), " on the shortest interval but should probably have a lot more. I would suggest decreasing timeStep to no more than ", mininterval/50))
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
	#cat("taxontodelete = ", taxontodelete)
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
						save(startVector, endVector, startTime, endTime, file=paste("savedHistory", jobName, ".Rdata", sep=""))
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
			pdf(paste("SimTree", jobName, ".pdf", sep=""))	
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
	timeStep, saveHistory=FALSE, saveRealParams=FALSE, jobName="", returnAll = FALSE, verbose=FALSE, reject.NaN=TRUE, taxon.df=NULL, checkTimeStep=TRUE) {
	
	if(is.null(taxon.df)){
		taxon.df <- getTaxonDFWithPossibleExtinction(phy)
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
		save(RealParams, file=paste("RealParams", jobName, ".Rdata", sep=""))
	}
	if (saveHistory) {
		startVector<-c()
		endVector<-c()
		startTime<-c()
		endTime<-c()
	}
	numberofsteps<-max(taxon.df$endTime)/timeStep
	mininterval<-min(taxon.df$endTime - taxon.df$startTime)
	#
	if(checkTimeStep){
		if (numberofsteps<1000) {
			#warning(paste("You have only ", numberofsteps, " but should probably have a lot more. I would suggest decreasing timeStep to no more than ", taxon.df[1, 1]/1000))
			}
		if (floor(mininterval/timeStep)<50) {
			warning(paste("You have only ", floor(mininterval/timeStep), " on the shortest interval but should probably have a lot more if you expect change on this branch. I would suggest decreasing timeStep to no more than ", mininterval/50))
			}
		if (floor(mininterval/timeStep)<3) {
			warning(paste("You have only ", floor(mininterval/timeStep), " on the shortest interval but should probably have a lot more if you expect change on this branch. I would suggest decreasing timeStep to no more than ", mininterval/50,"we would suggest at least", mininterval/3))
		#	timeStep <- mininterval/3
			}
		}
	#initial setup
	depthfrompresent = max(taxon.df$endTime)
	heightfromroot = 0
	taxon.df$states[which(taxon.df$startTime==0)] <- startingValues
	taxon.df.previous <- taxon.df
	while(depthfrompresent>0) {
		if(reject.NaN) {
			taxon.df.previous <- taxon.df
		}
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
		#first evolve in this interval, then speciate
		for (taxon.index in sequence(length(alive.rows))) {
			if(is.na(taxon.df$states[alive.rows[taxon.index]])) {
				taxon.df$states[alive.rows[taxon.index]] <- taxon.df$states[which(taxon.df$id==taxon.df$ancestorId[alive.rows[taxon.index]])]
			}
			new.state <- taxon.df$states[alive.rows[taxon.index]] + intrinsicFn(params=intrinsicValues, states=current.states[taxon.index],
				timefrompresent =depthfrompresent)+extrinsicFn(params=extrinsicValues, selfstates=current.states[taxon.index], otherstates=current.states[-taxon.index], timefrompresent =depthfrompresent)
			if(is.na(new.state)) {
				warning("got bad sim")
				attempt.count=0
				while(is.na(new.state) & attempt.count < 100) {
					old = taxon.df$states[alive.rows[taxon.index]]
					intrinsic.displacement = intrinsicFn(params=intrinsicValues, states=current.states[taxon.index], timefrompresent =depthfrompresent)
					extrinsic.displacement = extrinsicFn(params=extrinsicValues, selfstates=current.states[taxon.index], otherstates=current.states[-taxon.index], timefrompresent =depthfrompresent)
					new.state <- old + intrinsic.displacement + extrinsic.displacement
					warning(paste("Attempt ", attempt.count, "led to using old value of", old, "intrinsicFn return of ",intrinsic.displacement, "and extrinsicFn return of ", extrinsic.displacement))
					print("IntrinsicValues")
					print(intrinsicValues)
				}
				if(is.na(new.state) & attempt.count==100) {
					stop(paste(
						"Simulating with these parameters resulted in problematic results; for one example, taxon.df$states[alive.rows[taxon.index]] was ",
						taxon.df$states[alive.rows[taxon.index]], "intrinsicFn returned ",
						intrinsicFn(params=intrinsicValues, states=current.states[taxon.index],
						timefrompresent =depthfrompresent), "and extrinsicFn returned", extrinsicFn(params=extrinsicValues,
						selfstates=current.states[taxon.index], otherstates=current.states[-taxon.index], timefrompresent =depthfrompresent)))
				}
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
			print(paste("now at height", height.end, "finishing at", max(taxon.df$endTime)))
			print(taxon.df)
		}
		if(reject.NaN) {
			if(any(is.nan(taxon.df$states))) {
				save(list=ls(), file="ErrorRun.rda")
				stop(paste("There was an NaN generated. See saved objects in ", getwd(), "/ErrorRun.rda", sep=""))
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


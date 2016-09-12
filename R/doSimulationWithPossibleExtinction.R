#' Discrete Time Character Simulation
#'
#' This function evolves continuous characters in a discrete time process.
#'
#' When saveHistory is TRUE, processor time will increase quite a bit.
#' SaveRealParams is useful for tracking the "real" true values if simulating
#' data for abc runs.  It is not useful for empirical abc runs.
#'
#' @param taxon.df Output from the function getTaxonDFWithPossibleExtinction; is a data frame
#' with info on all the taxa (including internal ones)
#' @param intrinsicFn Name of intrinsic function characters should be simulated
#' under
#' @param extrinsicFn Name of extrinsic function characters should be simulated
#' under
#' @param startingValues State at the root
#' @param intrinsicValues Vector of values corresponding to the params of the
#' intrinsic model
#' @param extrinsicValues Vector of values corresponding to the params of the
#' extrinsic model
#' @param timeStep This value corresponds to the number of discrete time steps
#' on the shortest branch
#' @param saveHistory Saves the character history throughout the simulation
#' @param saveRealParams Saves intrinsicValues and extrinsicValues as both a
#' matrix and a vector
#' @param jobName Optional name for the job
#' @return A data frame of species character (tip) values in the tree.
#' @author Brian O'Meara and Barb Banbury
#' @references O'Meara and Banbury, unpublished
#' @keywords doSimulation
#' @examples
#'
#'
#' phy<-rcoal(30)
#'
#' #Simple Brownian motion
#' char<-doSimulation(
#' 	splits=getSimulationSplits(phy),
#' 	intrinsicFn=brownianIntrinsic,
#' 	extrinsicFn=nullExtrinsic,
#' 	startingValues=c(10), #root state
#' 	intrinsicValues=c(0.01),
#' 	extrinsicValues=c(0),
#' 	timeStep=0.0001,
#' 	saveHistory=FALSE)
#'
#'
#' #Character displacement model with minimum bound
#' char<-doSimulation(
#' 	splits=getSimulationSplits(phy),
#' 	intrinsicFn=boundaryMinIntrinsic,
#' 	extrinsicFn=ExponentiallyDecayingPush,
#' 	startingValues=c(10), #root state
#' 	intrinsicValues=c(0.05, 10, 0.01),
#' 	extrinsicValues=c(0, .1, .25),
#' 	timeStep=0.001,
#' 	saveHistory=FALSE)
#'
doSimulationWithPossibleExtinction<-function(taxon.df, intrinsicFn, extrinsicFn, startingValues, intrinsicValues, extrinsicValues, timeStep, saveHistory=FALSE, saveRealParams=FALSE, jobName="") {
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
	if (numberofsteps<1000) {
		#warning(paste("You have only ", numberofsteps, " but should probably have a lot more. I would suggest decreasing timeStep to no more than ", splits[1, 1]/1000))
	}
	if (floor(mininterval/timeStep)<50) {
		warning(paste("You have only ", floor(mininterval/timeStep), " on the shortest interval but should probably have a lot more if you expect change on this branch. I would suggest decreasing timeStep to no more than ", mininterval/50))
	}
	if (floor(mininterval/timeStep)<3) {
		warning(paste("You have only ", floor(mininterval/timeStep), " on the shortest interval but should probably have a lot more if you expect change on this branch. I would suggest decreasing timeStep to no more than ", mininterval/50, "but we are automatically adjusting it to", mininterval/3))
		timeStep <- mininterval/3
	}
	#initial setup
	depthfrompresent = max(taxon.df$endTime)
	heightfromroot = 0
	taxon.df$states[which(taxon.df$startTime==0)] <- startingValues

	while(depthfrompresent>0) {
		depth.start <- depthfrompresent
		depth.end <- depthfrompresent - timeStep
		height.start <- heightfromroot
		height.end <- heightfromroot + timeStep
		ids.alive.at.start <- taxon.df$id[which(taxon.df$startTime <= height.start & taxon.df$endTime > height.start)]
		ids.alive.at.end <-  taxon.df$id[which(taxon.df$endTime > height.end & taxon.df$startTime <= height.end)]
		ids.changing.status <-  ids.alive.at.start[!(ids.alive.at.start  %in% ids.alive.at.end)]
		ids.speciating <- taxon.df$id[which((taxon.df$id %in% ids.changing.status) & (!taxon.df$terminal))]
		alive.rows <- which(taxon.df$id %in% ids.alive.at.start)
		current.states <- taxon.df$states[alive.rows]
		#first evolve in this interval, then speciate
		for (taxon.index in sequence(length(alive.rows))) {
			taxon.df$states[alive.rows[taxon.index]] <- taxon.df$states[alive.rows[taxon.index]] + intrinsicFn(params=intrinsicValues, states=current.states[taxon.index], timefrompresent =depthfrompresent)+extrinsicFn(params=extrinsicValues, selfstates=current.states[taxon.index], otherstates=curent.states[-taxon.index], timefrompresent =depthfrompresent)
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
	}
	final.results <- subset(taxon.df, terminal==TRUE)
	final.result.df <- data.frame(states=final.results$states)
	rownames(final.result.df) <- final.results$name

	return(final.result.df)
}

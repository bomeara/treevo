#' Discrete Time Character Simulation
#' 
#' This function evolves continuous characters in a discrete time process
#' 
#' When saveHistory is TRUE, processor time will increase quite a bit.
#' SaveRealParams is useful for tracking the "real" true values if simulating
#' data for abc runs.  It is not useful for empirical abc runs.
#' 

#' @param splits Output from the function getSimulationSplits; is a data frame
#' of branching times, ancestor and descendant vectors
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

#' @param plot Will create a new interactive window that plots character values
#' throughout the history of the tree

#' @param savePlot Saves the character tree using jobName

#' @param saveHistory Saves the character history throughout the simulation

#' @param saveRealParams Saves intrinsicValues and extrinsicValues as both a
#' matrix and a vector

#' @param jobName Optional name for the job

#' @return A data frame of species character (tip) values in the tree.

#' @author Brian O'Meara and Barb Banbury

# @references O'Meara and Banbury, unpublished
# @keywords doSimulation doSimulationForPlotting


#' @examples
#' 
#' 
#' phy<-rcoal(30)
#' 
#' #Simple Brownian motion
#' char<-doSimulationForPlotting(
#' 	splits=getSimulationSplits(phy), 
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
#' 	splits=getSimulationSplits(phy), 
#' 	intrinsicFn=boundaryMinIntrinsic,
#' 	extrinsicFn=ExponentiallyDecayingPushExtrinsic,
#' 	startingValues=c(10), #root state
#' 	intrinsicValues=c(0.05, 10, 0.01),
#' 	extrinsicValues=c(0, .1, .25),
#' 	timeStep=0.001,
#' 	plot=TRUE,
#' 	saveHistory=FALSE)
#' 

#' @name doSimulationForPlotting
#' @rdname doSimulationForPlotting
#' @export
doSimulationForPlotting<-function(splits, intrinsicFn, extrinsicFn, startingValues, intrinsicValues, extrinsicValues, timeStep, plot=FALSE, savePlot=FALSE, saveHistory=FALSE, saveRealParams=FALSE, jobName="") {
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
	numberofsteps<-floor(splits[1, 1]/timeStep)
	mininterval<-min(splits[1:(dim(splits)[1]-1), 1]-splits[2:(dim(splits)[1]), 1])
	if (numberofsteps<1000) {
		#warning(paste("You have only ", numberofsteps, " but should probably have a lot more. I would suggest decreasing timeStep to no more than ", splits[1, 1]/1000))
	}
	if (floor(mininterval/timeStep)<50) {
		#warning(paste("You have only ", floor(mininterval/timeStep), " on the shortest interval but should probably have a lot more. I would suggest decreasing timeStep to no more than ", mininterval/50))
	}
#initial setup
	timefrompresent=splits[1, 1]
	taxa<-list(abctaxon(id=splits[1, 3], states=startingValues), abctaxon(id=splits[1, 4], states=startingValues))
	splits<-splits[2:dim(splits)[1], ] #pop off top value
	
#start running
	while(timefrompresent>0) {
#print(timefrompresent)
#speciation if needed
		while ((timefrompresent-timeStep)<=splits[1, 1]) { #do speciation. changed from if to while to deal with effectively polytomies
			originallength<-length(taxa)
			taxontodelete<-Inf
			for (i in 1:originallength) {
				if (taxa[[i]]$id==splits[1, 2]) {
					taxontodelete<-i
					taxa[[originallength+1]] <- taxa[[i]]
					taxa[[originallength+2]] <- taxa[[i]]
					taxa[[originallength+1]]$id<-splits[1, 3]
					taxa[[originallength+1]]$timeSinceSpeciation<-0
					taxa[[originallength+2]]$id<-splits[1, 4]
					taxa[[originallength+2]]$timeSinceSpeciation<-0
				}
			}
#cat("taxontodelete = ", taxontodelete)
			taxa<-taxa[-1*taxontodelete]
			if(dim(splits)[1]>1) {
				splits<-splits[2:(dim(splits)[1]), ] #pop off top value
			}
			else {
				splits[1, ]<-c(-1, 0, 0, 0)
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
		dev.new()
		plot(x=c(min(c(startVector, endVector)), max(c(startVector, endVector))), y=c(0, max(c(startTime, endTime))), type="n", ylab="Time", xlab="Trait value", main="", bty="n")
		for (i in 1:length(startVector)) {
			lines(x=c(startVector[i], endVector[i]), y=max(c(startTime, endTime)) - c(startTime[i], endTime[i]))
		}
	}
	if (savePlot) {
		pdf(paste("SimTree", jobName, ".pdf", sep=""))	
		plot(x=c(min(c(startVector, endVector)), max(c(startVector, endVector))), y=c(0, max(c(startTime, endTime))), type="n", ylab="Time", xlab="Trait value", main="", bty="n")
		for (i in 1:length(startVector)) {
			lines(x=c(startVector[i], endVector[i]), y=max(c(startTime, endTime)) - c(startTime[i], endTime[i]))
		}
		dev.off()
	}
	return(summarizeTaxonStates(taxa))
}

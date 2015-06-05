
parallelSimulation<-function(nrepSim, coreLimit, splits, phy, startingPriorsValues, intrinsicPriorsValues, extrinsicPriorsValues, startingPriorsFns, intrinsicPriorsFns, extrinsicPriorsFns, freevector, timeStep, intrinsicFn, extrinsicFn, multicore,checkpointFile=NULL,checkpointFreq=24) {
	#library(doMC, quietly=T)
	#library(foreach, quietly=T)
	cores=1
	if (multicore) {
		if (is.na(coreLimit)){
			registerDoMC()
			getDoParWorkers()->cores
		}
		else {
			registerDoMC(coreLimit)
			coreLimit->cores
		}
	}
	cat(paste("Using", cores, "core(s) for simulations \n\n"))
	if (nrepSim %%cores != 0) {
		warning("The simulation is most efficient if the number of nrepSim is a multiple of the number of cores")
	}

	cat("Doing simulations: ")
	if (is.null(checkpointFile)) {
		trueFreeValuesANDSummaryValues<-foreach(1:nrepSim, .combine=rbind) %dopar% simulateData(splits, phy, startingPriorsValues, intrinsicPriorsValues, extrinsicPriorsValues, startingPriorsFns, intrinsicPriorsFns, extrinsicPriorsFns, freevector, timeStep, intrinsicFn, extrinsicFn)
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
			trueFreeValuesANDSummaryValues<-rbind(trueFreeValuesANDSummaryValues, foreach(1:numberSimsPerLoop, .combine=rbind) %dopar% simulateData(splits, phy, startingPriorsValues, intrinsicPriorsValues, extrinsicPriorsValues, startingPriorsFns, intrinsicPriorsFns, extrinsicPriorsFns, freevector, timeStep, intrinsicFn, extrinsicFn))
			save(trueFreeValuesANDSummaryValues,file=checkpointFileName)
			print(paste("Just finished",dim(trueFreeValuesANDSummaryValues)[1],"of",nrepSim,"simulations; progress so far saved in",checkpointFileName))
		}
		trueFreeValuesANDSummaryValues<-rbind(trueFreeValuesANDSummaryValues, foreach(1:numberSimsAfterLastCheckpoint, .combine=rbind) %dopar% simulateData(splits, phy, startingPriorsValues, intrinsicPriorsValues, extrinsicPriorsValues, startingPriorsFns, intrinsicPriorsFns, extrinsicPriorsFns, freevector, timeStep, intrinsicFn, extrinsicFn))
	}
	return(trueFreeValuesANDSummaryValues) 
}

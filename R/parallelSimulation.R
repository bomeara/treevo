
parallelSimulation<-function(nrepSim, coreLimit, startingPriorsValues, intrinsicPriorsValues, extrinsicPriorsValues, startingPriorsFns, intrinsicPriorsFns, extrinsicPriorsFns, trueFreeValues, freevector, timeStep, intrinsicFn, extrinsicFn, multicore, jobName, filename) {
	library(doMC, quietly=T)
	library(foreach, quietly=T)
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
	cat(paste("Using", cores, "core(s) for initial simulations \n\n"))

	trueFreeValuesANDSummaryValues<-foreach(1:nrepSim, .combine=rbind) %dopar% simulateData(startingPriorsValues, intrinsicPriorsValues, extrinsicPriorsValues, startingPriorsFns, intrinsicPriorsFns, extrinsicPriorsFns, trueFreeValues, freevector, timeStep, intrinsicFn, extrinsicFn)
	save(trueFreeValuesANDSummaryValues,file=filename,compress=TRUE)
	save(trueFreeValuesANDSummaryValues, file=paste("trueFreeValuesANDSummaryValues", jobName, ".Rdata", sep=""))
}

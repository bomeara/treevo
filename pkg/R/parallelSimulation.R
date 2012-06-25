
parallelSimulation<-function(nrepSim, coreLimit, splits, phy, startingPriorsValues, intrinsicPriorsValues, extrinsicPriorsValues, startingPriorsFns, intrinsicPriorsFns, extrinsicPriorsFns, freevector, timeStep, intrinsicFn, extrinsicFn, multicore) {
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
	cat("Doing simulations: ")
	
	trueFreeValuesANDSummaryValues<-foreach(1:nrepSim, .combine=rbind) %dopar% simulateData(splits, phy, startingPriorsValues, intrinsicPriorsValues, extrinsicPriorsValues, startingPriorsFns, intrinsicPriorsFns, extrinsicPriorsFns, freevector, timeStep, intrinsicFn, extrinsicFn)
}

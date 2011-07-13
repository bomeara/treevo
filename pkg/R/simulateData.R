simulateData<-function(simIndex, nrepSim, startingPriorsValues, intrinsicPriorsValues, extrinsicPriorsValues, startingPriorsFns, intrinsicPriorsFns, extrinsicPriorsFns, trueFreeValues, freevector, timeStep, intrinsicFn, extrinsicFn) {
	simIndex<-simIndex+1
	print(simIndex)
	trueStarting<-rep(NaN, dim(startingPriorsValues)[2])
	trueIntrinsic<-rep(NaN, dim(intrinsicPriorsValues)[2])
	trueExtrinsic<-rep(NaN, dim(extrinsicPriorsValues)[2])
	for (j in 1:dim(startingPriorsValues)[2]) {
		trueStarting[j]=pullFromPrior(startingPriorsValues[,j],startingPriorsFns[j])
	}
	for (j in 1:dim(intrinsicPriorsValues)[2]) {
		trueIntrinsic[j]=pullFromPrior(intrinsicPriorsValues[,j],intrinsicPriorsFns[j])
	}
	for (j in 1:dim(extrinsicPriorsValues)[2]) {
		trueExtrinsic[j]=pullFromPrior(extrinsicPriorsValues[,j],extrinsicPriorsFns[j])
	}
	trueInitial<-c(trueStarting, trueIntrinsic, trueExtrinsic)
	trueFreeValues<-rbind(trueFreeValues, trueInitial[freevector])

	cat("Now doing simulation rep ",simIndex," of ",nrepSim,"\n",sep="")
	simdata<-convertTaxonFrameToGeigerData(doSimulation(splits=splits, intrinsicFn=intrinsicFn, extrinsicFn=extrinsicFn, startingStates=trueStarting, intrinsicValues=trueIntrinsic, extrinsicValues=trueExtrinsic, timeStep=timeStep), phy)
	
	simSummaryStats <-summaryStatsLong(phy, simdata, jobName=jobName)
	print(simSummaryStats)
	return(simSummaryStats)
}

		summaryValues<-rbind(summaryValues, summaryStatsLong(phy, simdata, jobName=jobName))

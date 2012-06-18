simulateDataS3<-function(startingPriorsValues, intrinsicPriorsValues, extrinsicPriorsValues, startingPriorsFns, intrinsicPriorsFns, extrinsicPriorsFns, trueFreeValues, freevector, timeStep, intrinsicFn, extrinsicFn, jobName) {
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

	cat("~")
	simdata<-summaryStatsLong(phy, convertTaxonFrameToGeigerData(doSimulationS3(splits=splits, intrinsicFn=intrinsicFn, extrinsicFn=extrinsicFn, startingStates=trueStarting, intrinsicValues=trueIntrinsic, extrinsicValues=trueExtrinsic, timeStep=timeStep), phy), jobName=jobName)
	
	simSummaryStats <-c(trueFreeValues, simdata)
	#print(simSummaryStats)
	#print(rbind(simSummaryStats, simSummaryStats))
	#save(simSummaryStats,file=paste("RunningSimulations",jobName,".Rdata",sep=""))
	#print(simSummaryStats)
	return(simSummaryStats)
}


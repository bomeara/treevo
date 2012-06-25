simulateData<-function(splits, phy, startingPriorsValues, intrinsicPriorsValues, extrinsicPriorsValues, startingPriorsFns, intrinsicPriorsFns, extrinsicPriorsFns, freevector, timeStep, intrinsicFn, extrinsicFn) {
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
	trueFreeValues<-trueInitial[freevector]

	cat(".")
	simTraits<-convertTaxonFrameToGeigerData(doSimulation(splits=splits, intrinsicFn=intrinsicFn, extrinsicFn=extrinsicFn, startingValues=trueStarting, intrinsicValues=trueIntrinsic, extrinsicValues=trueExtrinsic, timeStep=timeStep), phy)
	simSumStats<-summaryStatsLong(phy, simTraits)
	simTrueAndStats <-c(trueFreeValues, simSumStats)
	return(simTrueAndStats)
}


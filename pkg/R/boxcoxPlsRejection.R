boxcoxPlsRejection2<-function(summaryValuesMatrix, trueFreeValuesMatrix, phy, traits, vipthresh, abcTolerance, verbose=TRUE){
	#boxcox transform summary values
	if (verbose) {
		print("Starting boxcox transformation")
	}
	boxcoxEstimates<-boxcoxEstimation(summaryValuesMatrix)
	boxcoxAddition<-boxcoxEstimates$boxcoxAddition
	boxcoxLambda<-boxcoxEstimates$boxcoxLambda
	boxcoxSummaryValuesMatrix<-boxcoxEstimates$boxcoxSummaryValuesMatrix
	boxcoxOriginalSummaryStats<-boxcoxTransformation(summaryStatsLong(phy=phy, data=traits), boxcoxAddition, boxcoxLambda)

	if (verbose) {
		print("Starting to identify important summary stats")
	}
	library("mixOmics")
	getVipSingleColumn<-function(trueFreeValuesColumn, boxcoxSummaryValuesMatrix) {
		return(vip(pls(X=boxcoxSummaryValuesMatrix, Y=trueFreeValuesColumn, ncomp=1)))
	}
	
	#note this has vip for summary val 5, for true param 1, in row 5, column 1, and so forth
	allVip<-apply(trueFreeValuesMatrix, 2, getVipSingleColumn, boxcoxSummaryValuesMatrix)
	
	#this will have which summary vals have importance above the threshold
	whichVip<-(allVip>vipthresh)
	rownames(whichVip)<-sumStatNames(phy)
	
	save(trueFreeValuesMatrix, whichVip, boxcoxOriginalSummaryStats, boxcoxSummaryValuesMatrix, file="BarbsDistObjects.Rdata")
	if (verbose) {
		print("Calculating distances")
	}

	calculatedDist<-abcDistance(trueFreeValuesMatrix, whichVip, boxcoxOriginalSummaryStats, boxcoxSummaryValuesMatrix) 
	abcDistances<-calculatedDist$abcDistances
	
	#abcDistancesRaw<-matrix(nrow=dim(trueFreeValuesMatrix)[1], ncol=dim(trueFreeValuesMatrix)[2])  #rep(0,dim(trueFreeValuesMatrix)[1]) #will hold the distances for each particle
	#abcDistancesRawTotal<-vector(length=dim(trueFreeValuesMatrix)[1])
	#now we go true parameter by true parameter, using only the summary values with enough importance for each
	#we get a distance for each particle from the observed particle for each true param
	#then get a euclidean distance for all of these
	#it's like getting delta X, then delta Y, and figuring out distance from the origin using
	#sqrt((X-0)^2 + (Y-0)^2)
	#for (freeParamIndex in sequence(dim(trueFreeValuesMatrix)[2])) {
	#	abcDistancesRaw[,freeParamIndex]<-abc(target=boxcoxOriginalSummaryStats[whichVip[,freeParamIndex]], param=trueFreeValuesMatrix[,freeParamIndex], sumstat= boxcoxSummaryValuesMatrix[,whichVip[,freeParamIndex]], tol=1, method=abcMethod)$dist^2 #because we're doing euclidean distance, from observed, which has value 0, 0, 0, etc.
	#}
	#abcDistancesRawTotal<-apply(abcDistancesRaw, 1, sum)
	#abcDistances<-sqrt(abcDistancesRawTotal) #Euclid rules.
	if (verbose) {
		print("Getting accepted particles")
	}


	acceptedParticles<-trueFreeValuesMatrix[which(abcDistances<=quantile(abcDistances, prob=abcTolerance)), ] #here's where we diy abc
	acceptedDistances<-abcDistances[which(abcDistances<=quantile(abcDistances, prob=abcTolerance))]
	
	particleDataFrame<-data.frame(cbind(rep(1, dim(acceptedParticles)[1]), as.vector(which(abcDistances<=quantile(abcDistances, prob=abcTolerance))), seq(1:dim(acceptedParticles)[1]), rep(0, dim(acceptedParticles)[1]), acceptedDistances, rep(1, dim(acceptedParticles)[1]), acceptedParticles))
	colnames(particleDataFrame)<-c("generation", "attempt", "id", "parentid", "distance", "weight",  paste("param", seq(dim(trueFreeValuesMatrix)[2])))

	return(list(particleDataFrame=particleDataFrame, abcDistancesRaw=abcDistancesRaw, whichVip=whichVip, boxcoxSummaryValuesMatrix=boxcoxSummaryValuesMatrix, boxcoxOriginalSummaryStats=boxcoxOriginalSummaryStats))
}
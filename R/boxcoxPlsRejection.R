boxcoxPlsRejection<-function(summaryValuesMatrix, trueFreeValuesMatrix, phy, traits, vipthresh, abcTolerance){
	#boxcox transform summary values
	boxcoxEstimates<-boxcoxEstimation(summaryValuesMatrix)
	boxcoxAddition<-boxcoxEstimates$boxcoxAddition
	boxcoxLambda<-boxcoxEstimates$boxcoxLambda
	boxcoxSummaryValuesMatrix<-boxcoxEstimates$boxcoxSummaryValuesMatrix
	boxcoxOriginalSummaryStats<-boxcoxTransformation(summaryStatsLong(phy=phy, data=traits), boxcoxAddition, boxcoxLambda)

	library("mixOmics")
	getVipSingleColumn<-function(trueFreeValuesColumn, boxcoxSummaryValuesMatrix) {
		return(vip(pls(X=boxcoxSummaryValuesMatrix, Y=trueFreeValuesColumn, ncomp=1)))
	}
	
	#note this has vip for summary val 5, for true param 1, in row 5, column 1, and so forth
	allVip<-apply(trueFreeValuesMatrix, 2, getVipSingleColumn, boxcoxSummaryValuesMatrix)
	
	#this will have which summary vals have importance above the threshold
	whichVip<-(allVip>vipthresh)
	rownames(whichVip)<-sumStatNames(phy)
	
	abcDistancesRaw<-matrix(nrow=dim(trueFreeValuesMatrix)[1], ncol=dim(trueFreeValuesMatrix)[2])  #rep(0,dim(trueFreeValuesMatrix)[1]) #will hold the distances for each particle
	#abcDistancesRawTotal<-vector(length=dim(trueFreeValuesMatrix)[1])
	#now we go true parameter by true parameter, using only the summary values with enough importance for each
	#we get a distance for each particle from the observed particle for each true param
	#then get a euclidean distance for all of these
	#it's like getting delta X, then delta Y, and figuring out distance from the origin using
	#sqrt((X-0)^2 + (Y-0)^2)
	for (freeParamIndex in sequence(dim(trueFreeValuesMatrix)[2])) {
		abcDistancesRaw[,freeParamIndex]<-abc(target=boxcoxOriginalSummaryStats[whichVip[,freeParamIndex]], param=trueFreeValuesMatrix[,freeParamIndex], sumstat= boxcoxSummaryValuesMatrix[,whichVip[,freeParamIndex]], tol=1, method="rejection")$dist^2 #because we're doing euclidean distance, from observed, which has value 0, 0, 0, etc.
	}
	abcDistancesRawTotal<-apply(abcDistancesRaw, 1, sum)
	abcDistances<-sqrt(abcDistancesRawTotal) #Euclid rules.

	abcResults<-vector("list")
	abcResults$unadj.values<-trueFreeValuesMatrix[which(abcDistances<=quantile(abcDistances, prob=abcTolerance)), ] #here's where we diy abc
	abcResults$dist<-abcDistances[which(abcDistances<=quantile(abcDistances, prob=abcTolerance))]
	
	particleDataFrame<-data.frame(cbind(rep(1, dim(abcResults$unadj.values)[1]), as.vector(which(abcDistances<=quantile(abcDistances, prob=abcTolerance))), seq(1:dim(abcResults$unadj.values)[1]), rep(0, dim(abcResults$unadj.values)[1]), abcResults$dist, rep(1, dim(abcResults$unadj.values)[1]), abcResults$unadj.values  ))
	colnames(particleDataFrame)<-c("generation", "attempt", "id", "parentid", "distance", "weight",  paste("param", seq(dim(trueFreeValuesMatrix)[2])))

	return(list(particleDataFrame=particleDataFrame, abcDistancesRaw=abcDistancesRaw, whichVip=whichVip))
}
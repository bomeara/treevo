
#remember, when doing with PRC, set
#boxcoxSummaryValuesMatrix<-matrix(boxcoxParticleSummaryStats,nrow=1)
abcDistance<-function(trueFreeValuesMatrix,whichVip,boxcoxOriginalSummaryStats,boxcoxSummaryValuesMatrix) {
	abcDistancesRaw<-matrix(nrow=dim(trueFreeValuesMatrix)[1], ncol=dim(trueFreeValuesMatrix)[2])  
	#rep(0,dim(trueFreeValuesMatrix)[1]) #will hold the distances for each particle
	#abcDistancesRawTotal<-vector(length=dim(trueFreeValuesMatrix)[1])
	#now we go true parameter by true parameter, using only the summary values with enough importance for each
	#we get a distance for each particle from the observed particle for each true param
	#then get a euclidean distance for all of these
	#it's like getting delta X, then delta Y, and figuring out distance from the origin using
	#sqrt((X-0)^2 + (Y-0)^2)
	distanceByRow<-function(x,boxcoxOriginalSummaryStats) {
		return(dist(matrix(c(x,boxcoxOriginalSummaryStats),byrow=TRUE,nrow=2))[1])
	}
	for (freeParamIndex in sequence(dim(trueFreeValuesMatrix)[2])) {
		abcDistancesRaw[,freeParamIndex]<-(apply(boxcoxSummaryValuesMatrix,1,distanceByRow,boxcoxOriginalSummaryStats=boxcoxOriginalSummaryStats))^2 #need to square for euclidean distance, finished being calculated below
	}
	abcDistancesRawTotal<-apply(abcDistancesRaw, 1, sum)
	abcDistances<-sqrt(abcDistancesRawTotal) #Euclid rules.
	return(list(abcDistances=abcDistances, abcDistancesRaw=abcDistancesRaw)) #this will be a single value if we've passed it a single particle to compare with observed
}
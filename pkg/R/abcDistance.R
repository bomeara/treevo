
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
    numParticles<-dim(trueFreeValuesMatrix)[1]
    if (numParticles==1) { #because abc() fails with just one particle,make up data, then only return the real data (due to the sequence(numParticles) at the end
    	boxcoxSummaryValuesMatrix<-rbind(boxcoxSummaryValuesMatrix,boxcoxOriginalSummaryStats)
    	trueFreeValuesMatrix<-rbind(trueFreeValuesMatrix,rnorm(dim(trueFreeValuesMatrix)[2]))
    }
	for (freeParamIndex in sequence(dim(trueFreeValuesMatrix)[2])) {
		abcDistancesRaw[,freeParamIndex]<-abc(target=boxcoxOriginalSummaryStats[whichVip[,freeParamIndex]], param=trueFreeValuesMatrix[,freeParamIndex], sumstat= boxcoxSummaryValuesMatrix[,whichVip[,freeParamIndex]], tol=1, method=abcMethod)$dist^2 #because we're doing euclidean distance, from observed, which has value 0, 0, 0, etc.
	}
	abcDistancesRawTotal<-apply(abcDistancesRaw, 1, sum)
	abcDistances<-sqrt(abcDistancesRawTotal) #Euclid rules.
	return(list(abcDistances=abcDistances[sequence(numParticles),], abcDistancesRaw=abcDistancesRaw[sequence(numParticles),])) #this will be a single value if we've passed it a single particle to compare with observed
}
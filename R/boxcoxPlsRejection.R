boxcoxPlsRejection<-function(summaryValuesMatrix, trueFreeValuesMatrix, phy, traits, vipthresh, abcTolerance, verbose=TRUE){
	#uses library(pls)
	#first, find number of components using cross validation
	pls.model <- plsr(trueFreeValuesMatrix~summaryValuesMatrix,validation="CV",scale=TRUE) #scaling is important
	explained.variance <-cumsum(sort(attr(scores(pls.model),"explvar"),decreasing=TRUE))
	ncomp.final<-min(c(as.numeric(which(explained.variance>=95)[1]),length(explained.variance)),na.rm=TRUE) #min is to deal with case of never explaining >95%
	pls.model.final <- plsr(trueFreeValuesMatrix~summaryValuesMatrix,ncomp=ncomp.final, validation="none",scale=TRUE) #now rerun with the ideal number of components
	original.summary.stats<-summaryStatsLong(phy=phy, data=traits)
	original.summary.stats.transformed<-predict(pls.model.final,matrix(original.summary.stats,nrow=1),type="scores")
	simulated.summary.stats.transformed<-predict(pls.model.final,summaryValuesMatrix,type="scores")
 	
	
	calculatedDist<-abcDistance(trueFreeValuesMatrix,matrix(rep(TRUE,ncomp.final),nrow=1), original.summary.stats.transformed, simulated.summary.stats.transformed)$abcDistancesRaw[,1] #later, fix abcDistance 

#open question: does the above procedure work best or is it better to do PLS for each parameter separately?


	if (verbose) {
		print("Getting accepted particles")
	}


	acceptedParticles<-trueFreeValuesMatrix[which(abcDistances<=quantile(abcDistances, prob=abcTolerance)), ] #here's where we diy abc
	acceptedDistances<-abcDistances[which(abcDistances<=quantile(abcDistances, prob=abcTolerance))]
	
	particleDataFrame<-data.frame(cbind(rep(1, dim(acceptedParticles)[1]), as.vector(which(abcDistances<=quantile(abcDistances, prob=abcTolerance))), seq(1:dim(acceptedParticles)[1]), rep(0, dim(acceptedParticles)[1]), acceptedDistances, rep(1, dim(acceptedParticles)[1]), acceptedParticles))
	colnames(particleDataFrame)<-c("generation", "attempt", "id", "parentid", "distance", "weight",  paste("param", seq(dim(trueFreeValuesMatrix)[2])))

	return(list(particleDataFrame=particleDataFrame, calculatedDist=calculatedDist,  boxcoxSummaryValuesMatrix=boxcoxSummaryValuesMatrix, boxcoxOriginalSummaryStats=boxcoxOriginalSummaryStats))
}
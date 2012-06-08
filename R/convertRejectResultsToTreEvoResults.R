convertabcResultsToTreEvoResults<-function(abcResults){
	#wrapper function to convert results from rejection analyses in abc() to our particleDataFrame results for diagnostic or plotting functions
	#could think about subsetting only those particles/simulations athat are accepted so that we don't have to do that in all the different functions downstream (could add object to do all gens or just last gen too)
	if (class(abcResults) == "data.frame"){ #if TreEvo results are input
		particleDataFrame<-abcResults
		print(paste("Already particleDataFrame"))
	}
	if (class(abcResults) == "matrix"){ #if the abc$unadj.values are input
		particleDataFrame<-data.frame(cbind(rep(1, dim(abcResults)[1]), seq(1:dim(abcResults)[1]), seq(1:dim(abcResults)[1]), rep(0, dim(abcResults)[1]), rep(0, dim(abcResults)[1]), rep(1, dim(abcResults)[1]), abcResults))
		colnames(particleDataFrame)<-c("generation", "attempt", "id", "parentid", "distance", "weight")
	}
	if (class(abcResults) == "abc"){ #if the entire abc result is input
		particleDataFrame<-data.frame(cbind(rep(1, dim(abcResults$unadj.values)[1]), as.vector(which(abcResults$region)), seq(1:dim(abcResults$unadj.values)[1]), rep(0, dim(abcResults$unadj.values)[1]), abcResults$dist, rep(1, dim(abcResults$unadj.values)[1]), abcResults$unadj.values  ))
	
		colnames(particleDataFrame)<-c("generation", "attempt", "id", "parentid", "distance", "weight",  abcResults$names$parameter.names)
	}
	
	return(particleDataFrame)
}


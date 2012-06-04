convertRejectResultsToTreEvoResults<-function(RejectResults){
	#wrapper function to convert results from rejection analyses in abc() to our particleDataFrame results for diagnostic or plotting functions
	#could think about subsetting only those particles/simulations athat are accepted so that we don't have to do that in all the different functions downstream (could add object to do all gens or just last gen too)
	if (class(RejectResults) == "data.frame"){ #if TreEvo results are input
		particleDataFrame<-RejectResults
		print(paste("Already particleDataFrame"))
	}
	if (class(RejectResults) == "matrix"){ #if the abc$unadj.values are input
		particleDataFrame<-data.frame(cbind(rep(1, dim(RejectResults)[1]), seq(1:dim(RejectResults)[1]), seq(1:dim(RejectResults)[1]), rep(0, dim(RejectResults)[1]), rep(0, dim(RejectResults)[1]), rep(1, dim(RejectResults)[1]), RejectResults))
		colnames(particleDataFrame)<-c("generation", "attempt", "id", "parentid", "distance", "weight")
	}
	if (class(RejectResults) == "abc"){ #if the entire abc result is input
		particleDataFrame<-data.frame(cbind(rep(1, dim(RejectResults$unadj.values)[1]), as.vector(which(a$region)), seq(1:dim(RejectResults$unadj.values)[1]), rep(0, dim(RejectResults$unadj.values)[1]), RejectResults$dist, rep(1, dim(RejectResults$unadj.values)[1]), RejectResults$unadj.values  ))
	
		colnames(particleDataFrame)<-c("generation", "attempt", "id", "parentid", "distance", "weight",  a$names$parameter.names)
	}
	
	return(particleDataFrame)
}

as.vector(which(a$region))
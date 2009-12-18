#models
#note that these work for univariate, but need to be generalized for multivariate
#otherstates has one row per taxon, one column per state
#states is a vector for each taxon, with length=nchar

nullExtrinsic<-function(params,selfstates,otherstates) {
	newstates<-0*selfstates
	return(newstates)
}

nullIntrinsic<-function(params,states) {
	newstates<-0*states
	return(newstates)
}



brownianIntrinsic<-function(params,states) {
	newstates<-rnorm(n=length(states),mean=states,sd=params)
	return(newstates)
}

thresholdIntrinsic<-function(params, states) {
	#params[2] is min, params[3] is max. params[2] could be 0 or -Inf, for example
	newstates<-min(max(rnorm(n=length(states),mean=states,sd=params[1]),params[2]),params[3])
	for (i in length(newstates)) {
		newstates[i]<-max(newstates[i],params[2]) #whichever is larger, the current value or the min value
		newstates[i]<-min(newstates[i],params[3])
	}
	return(newstates)

}

autoregressiveIntrinsic<-function(params,states) { #a discrete time OU, same sd, mean, and attraction for all chars
	sd<-params[1]
	attractor<-params[2]
	attraction<-params[3]	#in this model, this should be between zero and one
	newstates<-rnorm(n=length(states),mean=attraction*states + attractor,sd=params)
	return(newstates)
}

nearestNeighborDisplacementExtrinsic<-function(params,selfstates,otherstates) { #this is set up for one character only right now
	repulsor<-otherstates[which.min(abs(otherstates-selfstates)),1] #the ,1 makes it for one char only
	#add things here
	
	return(newstates)
}



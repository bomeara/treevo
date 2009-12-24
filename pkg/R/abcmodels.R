#models
#note that these work for univariate, but need to be generalized for multivariate
#otherstates has one row per taxon, one column per state
#states is a vector for each taxon, with length=nchar

nullExtrinsic<-function(params,selfstates,otherstates, timefrompresent) {
	newstates<-0*selfstates
	return(newstates)
}

nullIntrinsic<-function(params,states, timefrompresent) {
	newstates<-0*states
	return(newstates)
}



brownianIntrinsic<-function(params,states, timefrompresent) {
	newstates<-rnorm(n=length(states),mean=states,sd=params)
	return(newstates)
}

boundaryIntrinsic<-function(params, states, timefrompresent) {
	#params[2] is min, params[3] is max. params[2] could be 0 or -Inf, for example
	newstates<-min(max(rnorm(n=length(states),mean=states,sd=params[1]),params[2]),params[3])
	for (i in length(newstates)) {
		newstates[i]<-max(newstates[i],params[2]) #whichever is larger, the current value or the min value
		newstates[i]<-min(newstates[i],params[3])
	}
	return(newstates)

}

boundaryMinIntrinsic<-function(params, states, timefrompresent) {
	#params[2] is min
	newstates<-max(rnorm(n=length(states),mean=states,sd=params[1]),params[2])
	for (i in length(newstates)) {
		newstates[i]<-max(newstates[i],params[2]) #whichever is larger, the current value or the min value
	}
	return(newstates)

}

autoregressiveIntrinsic<-function(params,states, timefrompresent) { #a discrete time OU, same sd, mean, and attraction for all chars
	sd<-params[1]
	attractor<-params[2]
	attraction<-params[3]	#in this model, this should be between zero and one
	newstates<-rnorm(n=length(states),mean=attraction*states + attractor,sd=sd)
	return(newstates)
}

nearestNeighborDisplacementExtrinsic<-function(params,selfstates,otherstates, timefrompresent) { #this is set up for one character only right now
	repulsorTaxon<-otherstates[which.min(abs(otherstates-selfstates)),1] #the ,1 makes it for one char only
	repulsorValue<-otherstates[repulsorTaxon]
	sd<-params[1]
	springK<-params[2]
	maxforce<-params[3]
	newstates<-rnorm(n=1,mean=selfstates[1]+sign(repulsorValue-selfstates[1])*min(abs(springK/((selfstates[1]-repulsorValue)*(selfstates[1]-repulsorValue))),maxforce),sd=sd)
	return(newstates)
}

everyoneDisplacementExtrinsic<-function(params,selfstates,otherstates, timefrompresent) { #this is set up for one character only right now
	sd<-params[1]
	springK <-params[2]
	maxforce<-params[3]
	netforce<-0
	for (i in 1:length(otherstates)) {
			netforce<-netforce+sign(otherstates[i]-selfstates[1])*min(abs(springK/((selfstates[1]-otherstates[i])*(selfstates[1]-otherstates[i]))),maxforce)
	}
	newstates<-rnorm(n=1,mean=selfstates[1]+ netforce,sd=sd)
	return(newstates)
}


autoregressiveIntrinsicTimeSlices<-function(params,states, timefrompresent) { #a discrete time OU, differing mean, sigma, and attaction with time
	#params=[sd1, attractor1, attraction1, timethreshold1, sd2, attractor2, attraction2, timethreshold2, ...]
	#time is time before present (i.e., 65 could be 65 MYA). The last time threshold should be 0, one before that is the end of the previous epoch, etc.
	numTimeSlices<-length(params)/4
	sd<-params[1]
	attractor<-params[2]
	attraction<-params[3]	#in this model, this should be between zero and one
	previousThresholdTime<-Inf
	for (slice in 0:(numTimeSlices-1)) {
		thresholdTime<-params[4+4*slice]
		if (thresholdTime >= timefrompresent) {
			if (thresholdTime<previousThresholdTime) {
				sd<-params[1+4*slice]
				attractor<-params[2+4*slice]
				attraction<-params[3+4*slice]
			}		
		}	
		previousThresholdTime<-thresholdTime
	}
	newstates<-rnorm(n=length(states),mean=attraction*states + attractor,sd=sd)
	return(newstates)
}

autoregressiveIntrinsicTimeSlicesConstantMean<-function(params,states, timefrompresent) { #a discrete time OU, constant mean, differing sigma, and differing attaction with time
	#params=[sd1, attraction1, timethreshold1, sd2, attraction2, timethreshold2, ..., attractor]
	#time is time before present (i.e., 65 could be 65 MYA). The last time threshold should be 0, one before that is the end of the previous epoch, etc.
	numTimeSlices<-(length(params)-1)/3
	sd<-params[1]
	attractor<-params[length(params)]
	attraction<-params[2]	#in this model, this should be between zero and one
	previousThresholdTime<-Inf
	for (slice in 0:(numTimeSlices-1)) {
		thresholdTime<-params[3+3*slice]
		if (thresholdTime >= timefrompresent) {
			if (thresholdTime<previousThresholdTime) {
				sd<-params[1+3*slice]
				attraction<-params[2+3*slice]
			}		
		}	
		previousThresholdTime<-thresholdTime
	}
	newstates<-rnorm(n=length(states),mean=attraction*states + attractor,sd=sd)
	return(newstates)
}



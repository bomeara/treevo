#models

brownianIntrinsic<-function(params,states) {
	newstates<-rnorm(n=length(states),mean=states,sd=params)
	return(newstates)
}

brownianExtrinsic<-function(params,selfstates,otherstates) {
	newstates<-0*selfstates
	return(newstates)
}

autoregressiveIntrinsic<-function(params,states) { #a discrete time OU
	sd<-params[1]
	attractor<-params[2]
	attraction<-params[3]	#in this model, this should be between zero and one
	newstates<-rnorm(n=length(states),mean=attraction*states + attractor,sd=params)
	return(newstates)
}

nullExtrinsic-function(params,selfstates,otherstates) {
	newstates<-0*selfstates
	return(newstates)
}

nullIntrinsic-function(params,states) {
	newstates<-0*states
	return(newstates)
}


#extrinsic models
#note that these work for univariate, but need to be generalized for multivariate
#otherstates has one row per taxon, one column per state
#states is a vector for each taxon, with length=nchar

nullExtrinsic<-function(params,selfstates,otherstates, timefrompresent) {
	newdisplacement<-0*selfstates
	return(newdisplacement)
}


nearestNeighborDisplacementExtrinsic<-function(params,selfstates,otherstates, timefrompresent) { 
	#params[1] is sd, params[2] is springK, params[3] is maxforce
	repulsorTaxon<-which.min(abs(otherstates-selfstates))
	repulsorValue<-otherstates[repulsorTaxon]
	sd<-params[1]
	springK<-params[2]
	maxforce<-params[3]
	localsign<-sign(selfstates[1]- repulsorValue)
	#print(abs((selfstates[1]-repulsorValue)))
	if(localsign==0) { 
		localsign=sign(rnorm(n=1))	
	}
	newdisplacement<-rnorm(n=1,mean=localsign*min(c(abs(springK/((selfstates[1]-repulsorValue)*(selfstates[1]-repulsorValue))),maxforce),na.rm=TRUE),sd=sd)
	return(newdisplacement)
}

ExponentiallyDecayingPush<-function(params,selfstates,otherstates, timefrompresent) { 
	#params[1] is sd, params[2] is maxForce when character difference = 0, params[3] is half distance (the phenotypic distance at which repulsion is half maxForce)
	repulsorTaxon<-which.min(abs(otherstates-selfstates))
	repulsorValue<-otherstates[repulsorTaxon]
	sd<-params[1]
	maxForce<-params[2]
	halfDistance<-params[3] #is like half life
	rate<-log(2, base = exp(1))/ halfDistance
	localsign<-sign(selfstates[1]- repulsorValue)
	if(localsign==0) {  #to deal with case of identical values
		localsign=sign(rnorm(n=1))	
	}
	newdisplacement<-rnorm(n=1,mean=maxForce*localsign*exp(-1*rate*abs((selfstates[1]-repulsorValue))),sd=sd)
	return(newdisplacement)
}

everyoneDisplacementExtrinsic<-function(params,selfstates,otherstates, timefrompresent) { #this is set up for one character only right now
	#params[1] is sd, params[2] is springK, params[3] is maxforce
	sd<-params[1]
	springK <-params[2]
	maxforce<-params[3]
	netforce<-0
	for (i in 1:length(otherstates)) {
			localsign<-sign(selfstates[1]-otherstates[i])
			if(localsign==0) {
				localsign=sign(rnorm(n=1))	
			}
			netforce<-netforce+localsign*min(c(abs(springK/((selfstates[1]-otherstates[i])*(selfstates[1]-otherstates[i]))),maxforce),na.rm=TRUE)
	}
	newdisplacement<-rnorm(n=1,mean=netforce,sd=sd)
	return(newdisplacement)
}


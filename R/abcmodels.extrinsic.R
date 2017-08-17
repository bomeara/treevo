#extrinsic models
#note that these work for univariate, but need to be generalized for multivariate
#otherstates has one row per taxon, one column per state
#states is a vector for each taxon, with length=nchar



#' Extrinsic Character Evolution Models
#' 
#' Functions describing various models of 'extrinsive' evolution (i.e. evolutionary processes
#' dependent on factors extrinsic to the evolving lineage, such as environmental change, or
#' other evolving lineages that interact with the lineage in question (competitors, predators, etc).
#'
#' The following extrinsic models are:
#'



#' 
#' This function describes a model of no extrinsic character evolution
#' 
#' 
#' @param params describes input paramaters for the model
#' @param selfstates vector of states for each taxon
#' @param otherstates matrix of character states, one row per taxon and once
#' column per state
#' @param timefrompresent which time slice in the tree
#' @return A matrix of values representing character displacement from a single
#' time step in the tree.
#' @author Brian O'Meara and Barb Banbury
# @references O'Meara and Banbury, unpublished
# @keywords nullExtrinsic extrinsic


#' @name extrinsicModels
#' @rdname extrinsicModels
#' @export
nullExtrinsic<-function(params,selfstates,otherstates, timefrompresent) {
	newdisplacement<-0*selfstates
	return(newdisplacement)
}




#' \code{nearestNeighborDisplacementExtrinsic} describes a model of extrinsic character evolution.  Character
#' values of a focal taxon depend on values of closest relatives on the tree


#' 
#' 
#' @param params describes input paramaters for the model.
#' \code{nearestNeighborDisplacementExtrinsic} params = sd, springK, maximum
#' force
#' @param selfstates vector of states for each taxon
#' @param otherstates matrix of character states, one row per taxon and once
#' column per state
#' @param timefrompresent which time slice in the tree
#' @return A matrix of values representing character displacement from a single
#' time step in the tree.
#' @author Brian O'Meara and Barb Banbury
# @references O'Meara and Banbury, unpublished
# @keywords nearestNeighborDisplacementExtrinsic extrinsic

#' @rdname extrinsicModels
#' @export
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



#' Extrinsic Character Evolution Models
#' 
#' \code{ExponentiallyDecayingPush} describes a model of extrinsic character evolution.  Character
#' values of a focal taxon pushes away harder from other taxa with like values;
#' "push" exponentially decays as the values become less similar
#' 
#' 
#' @param params describes input paramaters for the model.
#' \code{ExponentiallyDecayingPush} params = sd, maximum force, half distance
#' @param selfstates vector of states for each taxon
#' @param otherstates matrix of character states, one row per taxon and once
#' column per state
#' @param timefrompresent which time slice in the tree
#' @return A matrix of values representing character displacement from a single
#' time step in the tree.
#' @author Brian O'Meara and Barb Banbury
# @references O'Meara and Banbury, unpublished
# @keywords ExponentiallyDecayingPush extrinsic

#' @rdname extrinsicModels
#' @export
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



#' \code{nearestNeighborDisplacementExtrinsic} describes a model of extrinsic character evolution.  Character
#' values of a focal taxon depend on values of all relatives on the tree

#' 
#' 
#' @param params describes input paramaters for the model.
#' \code{nearestNeighborDisplacementExtrinsic} params = sd, springK, maximum
#' force
#' @param selfstates vector of states for each taxon
#' @param otherstates matrix of character states, one row per taxon and once
#' column per state
#' @param timefrompresent which time slice in the tree
#' @return A matrix of values representing character displacement from a single
#' time step in the tree.
#' @author Brian O'Meara and Barb Banbury
# @references O'Meara and Banbury, unpublished
# @keywords everyoneDisplacementExtrinsic extrinsic

#' @rdname extrinsicModels
#' @export
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



#Extra functions for calculating Exponential Decay Push priors
GetExpPushPriors<-function(numSteps, phy, data) {
  #returns a matrix with exponential rates for the three Exponential push priors
  timeStep<-1/numSteps  #out of doRun_rej code
  sd<-GetBMRatePrior(phy, data, timeStep)  #new TreEvo function
  data.sort<-sort(data[,1])
  data.min.diff<-min(abs(data.sort[1:(length(data.sort)-1)]-data.sort[2:(length(data.sort))]))
  data.max.diff<-abs(max(data.sort)-min(data.sort))
  half.distance.prior<-EstimateRate(data.min.diff)
  max.push.prior<-EstimateMaxForce(data.max.diff, halfDistance=qexp(0.5, half.distance.prior), numSteps)
  return(matrix(c(sd, sd, half.distance.prior, half.distance.prior, max.push.prior, max.push.prior), nrow=2))  #each one is doubled so they will exp rate param will fit into intrinsicPriorValues matrix in doRun_rej
}
qexpOptimality<-function(x, q, diff) {
  return(abs(qexp(q,x)-diff))
}
EstimateRate<-function(diff, quantile=0.95) {
  result<-optimize(f=qexpOptimality, interval = c(0, 100000), q=quantile, diff=diff)
  return(result$minimum)
}
OneStepPush<-function(maxForce, halfDistance, state) {
  rate<-log(2, base=exp(1))/ halfDistance
  push<-maxForce*2*exp(-1*rate*abs(state))
  return(push)
}
MultiStepPush<-function(maxForce, halfDistance, numSteps) {
  state<-0
  for (i in sequence(numSteps)) {
    state<-OneStepPush(maxForce, halfDistance, state)+state
  }
  return(state)
}
ForceOptimality<-function(x, halfDistance, numSteps, diff) {
  return(abs(diff - MultiStepPush(x, halfDistance, numSteps)))
}
EstimateMaxForce<-function(diff, halfDistance, numSteps) {
  result<-optimize(f=ForceOptimality, interval = c(0, 100000), halfDistance=halfDistance, numSteps = numSteps, diff=diff)
  return(result$minimum)
}

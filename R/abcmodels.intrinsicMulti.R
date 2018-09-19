#' Single Evolutionary Regime Model with Multiple Optima
#' 
#' This model describes an evolutionary process where multiple optima exist. Lineages are attracted
#' to an optima in their vicinity (this is handled as a stochastic process where which attract affects
#' which population at a given point in time is weighted with respect to the proximity of a given
#' population's trait value to the given optima), but as a lineage approaches the attractor that
#' controls it, it may experience less directional evolution in the direction of the optima,
#' and overall show more variation. This framework allows for the model to describe the evolution
#' with respect to multiple optima simultaneously, without needing to treat individual transitions
#' between optima (or macroevolutionary 'regimes') as a separate parameter. Thus, the optima exist
#' within a single evolutionary regime, with only a scaling parameter present that allows more
#' or less frequent transitions between optima by altering the impact of proximity to optima.


#' @details
#' \code{autoregressiveIntrinsic} describes a model of intrinsic character evolution. New
#' character values are generated after one time step via a discrete-time OU
#' process.
#' The input parameters for this model are:
#' \code{autoregressiveIntrinsic} with \code{params = sd (sigma), attractor (character mean), attraction (alpha)}


#' @details

#' @inheritParams abcmodels.intrinsic


#' @return
#' A vector of values representing character displacement of that lineage over a single time step.
 


#' @author David W. Bapst


#' @seealso
#' An alternative approach in \code{TreEvo} to estimating a macroevolutionary landscape with multiple optima
#' is the Fokker–Planck–Kolmogorov (FPK) model, which can be fit with \code{\link{landscapeFPK_Intrinsic}}.
#' This model does not assume a set number of optima nor that they have similar attractor strength, but
#' parameters may be difficult to interpret in isolation, and fitting this model may be slower with \code{TreEvo}
#' due to necessary linear algebra transformations. 
#' Other intrinsic models are described at \code{\link{abcmodels.intrinsic}}.

#' @examples
#' 
#' # three optima model, with strong attraction	
#' set.seed(1)
#' params<-c(
#' 	sigma=0.1,
#' 	alpha=0.7,
#' 	rho=1,
#' 	theta=c(-20,20,50)
#' 	)	
#' 	
#' multiOptimaIntrinsic(params=params, states=0, timefrompresent=NA)
#' 
#' # simulate n time-steps, repeat many times, plot results
#' repeatSimSteps<-function(params,trait=0,nSteps){
#' 	for(i in 1:nSteps){
#' 	# add to original trait value to get new trait value
#' 		trait<-trait+multiOptimaIntrinsic(
#' 			params=params, states=trait, timefrompresent=NA)
#' 		}
#' 	trait
#' 	}
#' repSim<-replicate(300,repeatSimSteps(params,trait=0,100))
#' hist(repSim,main="Simulated Trait Values",breaks=20)
#' 
#' 
#' # same model above, with more switching between optima
#' set.seed(1)
#' params<-c(
#' 	sigma=0.1,
#' 	alpha=0.7,
#' 	rho=0.5,
#' 	theta=c(-20,20,50)
#' 	)	
#' 	
#' multiOptimaIntrinsic(params=params, states=0, timefrompresent=NA)
#' 
#' # simulate n time-steps, repeat many times, plot results
#' repeatSimSteps<-function(params,trait=0,nSteps){
#' 	for(i in 1:nSteps){
#' 	# add to original trait value to get new trait value
#' 		trait<-trait+multiOptimaIntrinsic(
#' 			params=params, states=trait, timefrompresent=NA)
#' 		}
#' 	trait
#' 	}
#' repSim<-replicate(300,repeatSimSteps(params,trait=0,100))
#' hist(repSim,main="Simulated Trait Values",breaks=20)
#' 


# 08-06-18 moser_multi-optima-single evolutionary-regime-model.R 
# multi optima single evolutionary regime model
# MOSER?

#' @name multiOptimaIntrinsic
#' @rdname multiOptimaIntrinsic
#' @export
multiOptimaIntrinsic <- function(params, states, timefrompresent) {
    #a discrete time OU with multiple optima in the same regime 
		# with equal attraction (alpha) to all optima (theta 1:N)
	# breakdown of params:
		# params[1] is dispersion (sigma)
		# params[2] is alpha (strength of attraction to an optima)
		# params[3] is rho, an exponent scaling the weighting of distance to optima
			# this parameter will control switching optima
		# params[4:n] describes theta values
			# n-2 = N # of optima describe by this model
	# In this model, optima represent fixed trait values conveying adaptive benefit
		# the proximity of a population to an optima makes it more likely to be under that regime
		# a point equidistant between multiple regimes may be drawn to any
	# the draw to any specific optima is inverse to distance from optima
	# thus a lineage at an optima may show large variance as it circles the plateau
		# then suddenly feel drawn to another optima, and show sudden, giant shifts toward that optima
	# this all seems realistic...
	#
	sigma<-params[1]
	alpha<-params[2]
	rho <- params[3]
	theta<-params[-(1:3)]
	#
	# measure distances to theta
	# convert to probabilistic weights
		# raised to the power of rho - scaling parameter
	thetaWeights<-(1/abs(theta-states))^rho
	# rescale so sum to 1, as probabilities
	thetaWeights<-thetaWeights/sum(thetaWeights)
	# sample a theta
	theta<-sample(theta,1,prob=thetaWeights)
	# now 
	#subtract current states because we want displacement
    newdisplacement <- rpgm::rpgm.rnorm(n = length(states), mean = (theta-states)*alpha, sd = sd) 
    return(newdisplacement)
    }

	
	
	
	
	

multiOptimaIntrinsicMaxBoudary <- function(params, states, timefrompresent) {
    #a discrete time OU with multiple optima in the same regime 
		# with equal attraction (alpha) to all optima (theta 1:N)
		# and each regime having its own max trait value
	# breakdown of params:
		# params[1] is dispersion (sigma)
		# params[2] is alpha (strength of attraction to an optima)
		# params[3] is rho, an exponent scaling the weighting of distance to optima
			# this parameter will control switching optima
		# params[4:n] describes theta values
			# n-2 = N # of optima describe by this model
		# params 
		
		
		
	# In this model, optima represent fixed trait values conveying adaptive benefit
		# the proximity of a population to an optima makes it more likely to be under that regime
		# a point equidistant between multiple regimes may be drawn to any
	# the draw to any specific optima is inverse to distance from optima
	# thus a lineage at an optima may show large variance as it circles the plateau
		# then suddenly feel drawn to another optima, and show sudden, giant shifts toward that optima
	# this all seems realistic...
	#
	sigma<-params[1]
	alpha<-params[2]
	rho <- params[3]
	theta<-params[-(1:3)]
	#
	# measure distances to theta
	# convert to probabilistic weights
		# raised to the power of rho - scaling parameter
	thetaWeights<-(1/abs(theta-states))^rho
	# rescale so sum to 1, as probabilities
	thetaWeights<-thetaWeights/sum(thetaWeights)
	# sample a theta
	theta<-sample(theta,1,prob=thetaWeights)
	# now 
	#subtract current states because we want displacement
    newdisplacement <- rpgm::rpgm.rnorm(n = length(states), mean = (theta-states)*alpha, sd = sd) 
    return(newdisplacement)
    }
	
maxBoundaryAutoregressiveIntrinsic <- function(params, states, timefrompresent) {
    #a discrete time OU, same sd, mean, and attraction for all chars
    #params[1] is sd (sigma), 
		# params[2] is attractor (ie. character mean), params[3] is attraction (ie. alpha), 
		# params[4] is max bound
    sd <- params[1]
    attractor <- params[2]
    attraction <- params[3]    #in this model, this should be between zero and one
    minBound <- params[4]
    newdisplacement <- rpgm::rpgm.rnorm(
		n = length(states), 
		mean = (attractor-states)*attraction, 
		sd = sd) #subtract current states because we want displacement
    #message(newdisplacement)
	for (i in length(newdisplacement)) {
        newstate <- newdisplacement[i]+states[i]
        if (newstate>params[2]) { #newstate less than min
            newdisplacement[i] <- params[2]-states[i] 
			#so, rather than go below the minimum, this moves the new state to the maximum
        }
    }
    return(newdisplacement)
}



	
    #a discrete time OU, same sd, mean, and attraction for all chars
    #params[1] is sd (sigma), 
		# params[2] is attractor (ie. character mean), params[3] is attraction (ie. alpha), 
		# params[4] is max bound	
	#
	# September 2018
	# Make two more models for Aquilegia
	# 1) trait values have three optima on gradient, with some rate of switching to next-largest optima, cannot reverse
	# 2) trait values evolve in three regimes with  successive upper bounds on gradient (so only two upper-bounds, highest regime has no bounds) with some rate of switching to next-largest regime, cannot reverse
	# addendum - maybe make rate of switching to next optima dependent on trait value?

	
	
	
	# September 2018
	# Make two more models for Aquilegia
	# 1) trait values have three optima on gradient, with some rate of switching to next-largest optima, cannot reverse
	# 2) trait values evolve in three regimes with  successive upper bounds on gradient (so only two upper-bounds, highest regime has no bounds) with some rate of switching to next-largest regime, cannot reverse
	# addendum - maybe make rate of switching to next optima dependent on trait value?

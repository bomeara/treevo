#' Intrinsic Trait Evolution Model for a Macroevolutionary Landscape on a Fokker-Planck-Kolmogorov Potential Surface (Boucher et al., 2018)
#' 
#' This function describes a discrete-time Fokker-Planck-Kolmogorov Potential Surface model
#' (the FPK model for short), described for trait macroevolution by Boucher et al. (2018),
#' and included here in \code{TreEvo} for use as an intrinsic trait model. This model can
#' be used to describe a complex landscape of hills and valleys for a bounded univariate trait space.


#' @details
#' The FPK model is a four parameter model for describing the evolution of a trait without
#' reference to the traits of other taxa (i.e. an intrinsic model in the terminology of
#' package \code{TreEvo}). Three of these parameters are used to describe the shape of the
#' landscape, which dictates a lineage's overall deterministic evolutionary trajectory on
#' the landscape, while the fourth parameter (sigma) is a dispersion parameter that
#' describes the rate of unpredictable, stochastic change. 
#' 
#' The bounds on the trait space are treated as two additional parameters, but these are
#' intended in the model as described by Boucher et al. to be treated as nuisance parameters.
#' Boucher et al. intentionally fix these bounds at a far distance beyond the range of
#' observed trait values, so in order to ensure they have as little effect on the evolutionary trajectories as possible.
#' 
#' The discrete time Fokker-Planck-Kolmogorov model used here describes a landscape specific
#' to a population at a particular time, with a particular trait value. This landscape
#' represents a potential surface, where height corresponds to a tendency to change in that
#' direction. This allows for multiple optima, and assigns different heights to those optima,
#' which represent the current attraction for a population to move toward that trait value.
#' The shape of this potential surface is the macroevolutionary landscape. which Boucher
#' et al describe the shape of using a fourth-order polynomial with the cubic component removed:
#' 
#' \deqn{V(x) = ax^4 + bx^2 + cx}
#' 
#' Where x is the trait values within the bounded interval - for simplicity, these are
#' internally rescaled to sit within the arbitrary interval \code{(-1.5 : 1.5)}. The parameters
#' thus that describe the landscape shape are the coefficients \emph{a}, \emph{b}, and
#' \emph{c}. To calculate the landscape, the trait space is discretized into fine intervals,
#' and the potential calculated for each interval. The scale of this discretization can be controlled by the user.
#' 
#' Note that if the landscape has a single peak, or if the potential surface is flat, the
#' model effectively collapses to Brownian Motion, or Ornstein-Uhlenbeck with a single
#' optima. To (roughly) quote Boucher et al: \sQuote{Finally, note that both BM and the
#' OU model are special cases of the FPK model: BM corresponds to V(x)=0 and OU to V(x)=((alpha/sigma^2)*x^2)-((2*alpha*theta/(sigma^2))*x).}
#' 
#' 
			
#' @inheritParams intrinsicModels

#' @param	grainScaleFPK To calculate the potential-surface landscape, the trait
#' space is discretized into fine intervals, and the potential calculated for each
#' interval. The scale of this discretization can be controlled with this argument,
#' which specifies the number of intervals used (default is 1000 intervals).

#' @param	traitName The name given to the trait, mainly for use in the
#' macroevolutionary landscape plot.

#' @param	plotLandscape If \code{TRUE}, the estimated macroevolutionary
#' landscape is plotted by the function.

#' @param	traitData A set of trait data to calculate distant bounds from for use with this model.

#' @return
#' A vector of values representing character displacement of that lineage over a single time step.

#' @seealso 
#' An alternative approach in \code{TreEvo} to estimating a macroevolutionary landscape
#' with multiple optima is the intrinsic model \code{\link{multiOptimaIntrinsic}}, 
#' however this model requires an \emph{a priori} choice on the number of optima
#' and the assumption that the optima have similar attractor strength. 
#' Other intrinsic models are described at \code{\link{intrinsicModels}}.
		
#' @author David W. Bapst, loosely based on studying the code for function \code{Sim_FPK} from package \code{BBMV}.

#' @references		
#' Boucher, F. C., V. Demery, E. Conti, L. J. Harmon, and J. Uyeda. 2018.
#' A General Model for Estimating Macroevolutionary Landscapes. \emph{Systematic Biology}
#' 67(2):304-319.
						
			
#' @examples
#' 
#' set.seed(444)
#' traitData<-rnorm(100,0,1)
#' # need traits to calculate bounds
#' bounds<-getTraitBoundsFPK(traitData)
#' # pick a value at random
#' trait<-0
#' 
#' # two peak symmetric landscape example
#' params<-c(
#' 	a=2,
#' 	b=-4,
#' 	c=0,
#' 	sigma=1,
#' 	bounds)
#' 
#' plot_landscapeFPK_model(params)
#' 
#' # simulate under this model - simulated trait DIVERGENCE
#' landscapeFPK_Intrinsic(params=params, states=trait, timefrompresent=NULL)
#' 
#' # simulate n time-steps, repeat many times, plot results
#' repeatSimSteps<-function(params,trait,nSteps){
#' 	for(i in 1:nSteps){
#' 	# add to original trait value to get new trait value
#' 		trait<-trait+landscapeFPK_Intrinsic(
#' 			params=params, states=trait, timefrompresent=NULL)
#' 			}
#' 	trait
#' 	}
#' repSim<-replicate(30,repeatSimSteps(params,trait,20))
#' hist(repSim,main="Simulated Trait Values")
#' 
#' 
#' # uneven two peak symmetric landscape example
#' params<-c(
#' 	a=2,
#' 	b=-4,
#' 	c=0.3,
#' 	sigma=1,
#' 	bounds)
#' plot_landscapeFPK_model(params)
#' 
#' # simulate under this model - simulated trait DIVERGENCE
#' landscapeFPK_Intrinsic(params=params, states=trait, timefrompresent=NULL)
#' 
#' repSim<-replicate(30,repeatSimSteps(params,trait,20))
#' hist(repSim,main="Simulated Trait Values")
#' 
			

			
# all of the following only need to be run
		# when the parameters of FPK are changed
	# this *could* be pre-calculated for a single run with lexical scoping	
	#	

	
	
#' @name landscapeFPK_Intrinsic
#' @rdname landscapeFPK_Intrinsic
#' @export
landscapeFPK_Intrinsic <- function(params, states, timefrompresent,	
		grainScaleFPK = 100	# sim controls
		) {
	#
	#a discrete time Fokker-Planck-Kolmogorov model (FPK)
		# V(x)=ax4+bx2+cx 
	# parameters: a,b,c, sigma
		# dependent on former trait value
	#
	# describes a potential surface where height corresponds to
		# tendency to change in that direction
	# allows for multiple optima, different heights to those optima
		#also collapses to BM and OU1
	#
	# From Boucher et al:
	# Finally, note that both BM and the OU model are special cases of the FPK
	# model: BM corresponds to V(x)=0 and OU to
	# V(x)=((alpha/sigma^2)*x^2)-((2*alpha*theta/(sigma^2))*x)
	#
	# following code is loosely based on function Sim_FPK from package BBMV
	#
	# all of the following only need to be run
		# when the parameters of FPK are changed
	# this *could* be pre-calculated for a single run with lexical scoping	
	#	
	# landscape descriptor function
	# over the arbitrary interval (-1.5 : 1.5)
	arbSequence<-seq(from=-1.5,to=1.5,
		length.out=grainScaleFPK)
	#
	# get bounds from params
	bounds <- params[5:6]
	#
	# translate  to original trait scale
	origSequence<-seq(from=bounds[1],to=bounds[2],
		length.out=grainScaleFPK)
	origIntLength<-abs((bounds[2]-bounds[1])/(grainScaleFPK-1))
	# # potentialVector is numeric vector representing the potential
		# length = grainScaleFPK
	# V(x)=ax4+bx2+cx 
	potentialVector<-potentialFunFPK(
		x=arbSequence,
		a=params[1],b=params[2],c=params[3])	
	#
	# Coefficient of Diffusion of the Model
	dCoeff <- log((params[4])^2/2)   # log((sigma^2)/2)
	#
	# Transition matrix describing prob of evolving between two sites in the trait grid in an infinitesimal time step.	
	# Create and diagonalize the transition matrix that has been discretized
	# returns: the transition matrix going forward in time, for simulating traits only
	#
	# make empty matrix
	expD <- tranMatrix <- matrix(0,grainScaleFPK,grainScaleFPK)
	#assign values not on outer rows/columns
	for (i in 1:(grainScaleFPK)){
		if(i>1){
			tranMatrix[i-1,i] <- exp((potentialVector[i]-potentialVector[i-1])/2)
			}
		if(i<grainScaleFPK){
			tranMatrix[i+1,i] <- exp((potentialVector[i]-potentialVector[i+1])/2)
			}
		# rate of staying in place is negative sum of neighboring cells
		neighbors<-c(ifelse(i>1,tranMatrix[i-1,i],0),
				ifelse(i<grainScaleFPK,tranMatrix[i+1,i],0)
				)
		tranMatrix[i,i] <- (-sum(neighbors))
		}
	# eigenvalues and eigenvectors of transition matrix
		# take only real components
	eigTranMatrix <- lapply(eigen(tranMatrix),Re)
	# 			
	# solve the eigenvectors
	solvedEigenvectors <- solve(eigTranMatrix$vectors,tol = 1e-30) 
	#
	# get expected dispersion
	# scale expected dispersion to original trait scale
		# squared distance between points in resolution of trait scale
		# (tau from Boucher et al.'s original code)
	origScaler <- origIntLength^2
	# assign dispersion to diagonal of expD
	diag(expD) <- exp(exp(dCoeff)/origScaler*eigTranMatrix$values)
	# previous time-dep version from Boucher et al's code
		# diag(expD) <- exp(t*diag_expD)
	#
	# take dot product of expD, eigenvectors and solved eigenvectors
	# get matrix of potential for future trait values
		# given potentialVector alone, not yet considering previous values
	potentialMatrix <- eigTranMatrix$vectors%*%expD%*%solvedEigenvectors
	#
	###############################################
	#######################################################
	#
	# need a vector, length = grainScaleFPK
		# with zeroes in intervals far from current trait value
	# and with density=1 distributed in interval=origIntLength
		# centered around the original trait value
	intDensity<-getTraitIntervalDensityFPK(
		trait=states,
		origIntLength=origIntLength,
		origSequence=origSequence,
		grainScaleFPK=grainScaleFPK)
	#
	# take product to get potential with respect to position
	probDivergence <- potentialMatrix %*% intDensity
	# round all up to zero at least
	probDivergence[probDivergence<0] <-	0
	# convert potential to a probability
	probDivergence<-probDivergence / sum(probDivergence)
	#					
	# sample from this probability distribution
		# to get divergence over a time-step
	newTraitPosition <- sample(
		x=origSequence,
		size=1,
		prob=probDivergence)
	# subtract the current trait position so to get divergence
	newDisplacement<-newTraitPosition-states
	return(newDisplacement)
	}
	

#' @rdname landscapeFPK_Intrinsic
#' @export				
getTraitBoundsFPK<-function(traitData){
	# get actualistic bounds on function
	bounds <- c(
		min(traitData)-((max(traitData) - min(traitData))/2),
		max(traitData)+((max(traitData) - min(traitData))/2)
		)
	return(bounds)
	}

	
#' @rdname landscapeFPK_Intrinsic
#' @export	
plot_landscapeFPK_model<-function(params,
		grainScaleFPK=1000,
		traitName="Trait",
		plotLandscape=TRUE
		){
	#
	# get bounds from params
	bounds <- params[5:6]
	# get trait sequence
	origSequence<-seq(from=bounds[1],to=bounds[2],
		length.out=grainScaleFPK)
	# V(x)=ax4+bx2+cx 
	potentialVector<-potentialFunFPK(
		a=params[1],b=params[2],c=params[3],
		x=seq(from=-1.5,to=1.5,
			length.out=grainScaleFPK))
	#
	if(plotLandscape){
		#potentialVector<-1*exp(-potentialVector)
		# equation from Boucher's BBMV tutorial (??)
		potentialVector<-exp(-potentialVector)/sum(exp(-potentialVector)
				*((bounds[2]-bounds[1])/grainScaleFPK))
		yLabel<-"Macroevolutionary Landscape (N*exp(-V))"	
	}else{
		potentialVector<-potentialVector/max(potentialVector)
		yLabel<-"Evolutionary Potential (Rescaled to Max Potential)"
		}
	#
	plot(origSequence,potentialVector,type="l",
		xlab=paste0(traitName," (Original Scale)"),
		ylab=yLabel)
	}
	

getTraitIntervalDensityFPK<-function(trait,origIntLength,
		origSequence,grainScaleFPK){
	#### example dataset for testing
	# grainScaleFPK<-100
	# origIntLength<-0.23
	# origSequence<-(-20:(-20+grainScaleFPK))*origIntLength
	# trait<-(-1.83)
	# trait<-max(origSequence)
	# trait<-min(origSequence)
	####################################################
	traitRange<-c(trait-origIntLength/2,
		trait+origIntLength/2)
	#
	intDensity<-rep(NA,grainScaleFPK)
	intDensity[2:(grainScaleFPK-1)]<-sapply(2:(grainScaleFPK-1), function(i) max(0,
		min(origSequence[i+1],traitRange[2])-max(origSequence[i],traitRange[1])))
	#
	# special calculations for first and last
	if(traitRange[2]<origSequence[2]){
		intDensity[1]<-origIntLength
	}else{
		intDensity[1]<-0
		}
	#
	if(traitRange[1]>origSequence[grainScaleFPK-1]){
		intDensity[grainScaleFPK]<-origIntLength
	}else{
		intDensity[grainScaleFPK]<-0
		}
	#
	#if(length(intDensity)!=grainScaleFPK){
	#	stop("intDensity is not calculated with correct length")
	#	}
	#
	#sum(intDensity)==origIntLength
	return(intDensity)
	}



# equation for getting potential under FPK	
potentialFunFPK<-function(x,a,b,c){
	# V(x)=ax4+bx2+cx 
	Vres <- (a*(x^4))+(b*(x^2))+(c*x)
	return(Vres)
	}

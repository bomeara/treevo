


# need function that simply outputs a list with those parameters
	# (median?) expectations from the last MCMC generation
# this function would be run with doRun such that doRun would include with output
# these parameter estimates would be given as a list
	# split into 3 vectors: starting/intrinsic/extrinsic parameters
	# formatted for immediate use as parameter estimates for doSimulation
	# with matching intrinsic/extrinsic functions
	
# internal function

# this could then be used for easily generating new simulations from the median param values



getSummaryMeans<-function(doRunOutput, freeVector){
	# this tells us which parameters are fixed
	
	
	
	# need to sort into intrinsic, extrinsic, starting (and which are fixed)
		# this should be indicated easily by other aspects of the output
	
	
	# this contains means for unfixed parameters
	sumPostOutput<-doRunOutput$postSummary	
	# need to get out both fixed and unfixed parameters
	
	
	sumPostOutput 
	


	startingPar <- 
	intrinsicPar <- 
	extrinsicPar <- 
	# output list of vectors for starting, intrinsic, extrinsic
	res <- list(
		starting = startingPar,
		intrinsic = intrinsicPar,
		extrinsic = extrinsicPar,
		)
	return(res)
	}
	
	
# a function that would take a tree, models and means of parameters
  # and then simulate data using those means as fixed priors
  
#

function(){
	
	
	char <- doSimulation(
		phy = tree, 
		generation.time = 100000, 
		intrinsicFn = brownianIntrinsic, 
		extrinsicFn = nullExtrinsic, 
		startingValues = c(10), #root state
		intrinsicValues = c(0.01), 
		extrinsicValues = c(0)
		)
	res(char)
	}

  
 
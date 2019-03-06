


# need function that simply outputs a list with those parameters
	# (median?) expectations from the last MCMC generation
# this function would be run with doRun such that doRun would include with output
# these parameter estimates would be given as a list
	# split into 3 vectors: starting/intrinsic/extrinsic parameters
	# formatted for immediate use as parameter estimates for doSimulation
	# with matching intrinsic/extrinsic functions
	
# internal function

# this could then be used for easily generating new simulations from the median param values



getSummaryMeans<-function(doRunOutput){
	# this contains means for unfixed parameters
	sumPostOutput<-doRunOutput$postSummary
	
	
	# need to sort into intrinsic, extrinsic, starting
		# this should be indicated easily by other aspects of the outut

	
	
	# need to get out both fixed and unfixed parameters
	
	sumPostOutput 
	



	# output list of vectors for starting, intrinsic, extrinsic


	}
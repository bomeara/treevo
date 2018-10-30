


# need function that simply outputs a list with those parameters
	# (median?) expectations from the last MCMC generation
# this function would be run with doRun such that doRun would include with output
# these parameter estimates would be given as a list
	# split into 3 vectors: starting/intrinsic/extrinsic parameters
	# formatted for immediate use as parameter estimates for doSimulation
	# with matching intrinsic/extrinsic functions
	
# internal function



getSummaryMeans<-function(doRunOutput){

	sumPostOutput<-doRunOutput$postSummary
	sumPostOutput # this contains means for unfixed parameters
	# need to sort into intrinsic, extrinsic, starting

# need to get out both fixed and unfixed parameters


	# output list of vectors for starting, intrinsic, extrinsic


	}
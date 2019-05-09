# internal.postRunSimulateFromMeans.R

# postRunSimulateFromMeans
# internal function

# a function that would take a tree, models and means of parameters
	# and then simulate data using those means as fixed priors

postRunSimulateFromMeans <- function(prcOut, verbose = TRUE){
	# a function that would take a tree, models and means of parameters
		# and then simulate data using those means as fixed priors
	# 
	# needs to be output from a single Run
	if ( inherits(prcOut, "multiRun_doRun_prc") ) {
		if(verbose){
			message("Multiple runs of doRun_prc found, using means from the first run")
			}
        prcOut <- prcOut[[1]]
        }
	#
	startingValues <- prcOut$parMeansList$starting
	intrinsicValues <- prcOut$parMeansList$intrinsic
	extrinsicValues <- prcOut$parMeansList$extrinsic
	#
	char <- doSimulation(
		phy = prcOut$phy, 
		generation.time = prcOut$generation.time, 
		intrinsicFn = prcOut$intrinsicFn, 
		extrinsicFn = prcOut$extrinsicFn, 
		startingValues = startingValues, 
		intrinsicValues = intrinsicValues, 
		extrinsicValues = extrinsicValues
		)
	res(char)
	}




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
	# get the data from doRun output
	#
	# this tells us which parameters are fixed
	freeVector  <- doRunOutput$freeVector
	# this contains means for unfixed parameters
		# from the post-analysis summary function
	postSummary <- doRunOutput$postSummary
	# this has the values for fixed parameters
	priorList   <- doRunOutput$priorList
	#
	summaryValues <- numeric()
	#
	# get the summary values (means, or fixed values)
	for(i in 1:length(freeVector)){
		# need to sort out which are free and which are fixed
		if(freeVector[i]){
			# IF IS FREE
			matchPar <- which(names(postSummary) == names(freeVector[i]))
			summaryValues[i] <- postSummary[[matchPar]]$mean
		}else{
			# IF IS FIXED
			summaryValues[i] <- priorList[[i]]$params
			}
		}
	#
	# name the vector the parameter names
	names(summaryValues) <- names(freeVector)
	# need to sort into intrinsic, extrinsic, starting 
		# this should be indicated easily by name
		# eg starting_ intrinsic_ extrinsic_
	parTypeMat <- startsWithMultiple(test = names(freeVector),
		prefixes = c("starting","intrinsic","extrinsic")
	# now test
	# all row sums should be identical to one
	if(any(rowSums(parTypeMat) != 1){
		stop("Some parameter labels from freeVector do not uniquely match parameter types")
		}
	# now sort into parameter types

	startingPar  <- summaryValues[parTypeMat[,"starting"]]
	intrinsicPar <- summaryValues[parTypeMat[,"intrinsic"]]
	extrinsicPar <- summaryValues[parTypeMat[,"extrinsic"]]	
	#
	# output list of vectors for starting, intrinsic, extrinsic
	res <- list(
		starting = startingPar,
		intrinsic = intrinsicPar,
		extrinsic = extrinsicPar,
		)
	return(res)
	}

		
startsWithMultiple <- function(test, prefixes){
	# constructs a matrix of logicals indicating which labels match which prefixes
		# when multiples of each exist
	#
	# example
	# startsWithMultiple(c("hellooo","bye"),prefixes=c("hello","bye","hell","b","bell"))
	#
	res <- sapply(prefixes, function(x) 
		 startsWith(test, prefix = x)
		 )
	rownames(res) <- test
	colnames(res) <- prefixes
	return(res)
	}	
	
	
	
	
	
	
	
	
	
	

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

  
 
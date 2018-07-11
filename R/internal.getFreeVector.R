# internal function for calculating freevector

getFreeVector <- function(startingPriorsFns, startingPriorsValues,
						intrinsicPriorsFns, intrinsicPriorsValues, 
						extrinsicPriorsFns, extrinsicPriorsValues){
	#checks
	if(!is.matrix(startingPriorsValues)){
		stop("startingPriorsValues must be a matrix, with columns representing separate prior distributions parameters, and rows representing parameters for those priors")
		}
	if(!is.matrix(intrinsicPriorsValues)){
		stop("intrinsicPriorsValues must be a matrix, with columns representing separate prior distributions parameters, and rows representing parameters for those priors")
		}	
	if(!is.matrix(extrinsicPriorsValues)){
		stop("extrinsicPriorsValues must be a matrix, with columns representing separate prior distributions parameters, and rows representing parameters for those priors")
		}	
	#
	if(dim(startingPriorsValues)[2] != length(startingPriorsFns)){
		stop("startingPriorsValues must have the same number of columns as functions names in startingPriorsFns")
		}
	if(dim(intrinsicPriorsValues)[2] != length(intrinsicPriorsFns)){
		stop("intrinsicPriorsValues must have the same number of columns as functions names in intrinsicPriorsFns")
		}
	if(dim(extrinsicPriorsValues)[2] != length(extrinsicPriorsFns)){
		stop("extrinsicPriorsValues must have the same number of columns as functions names in extrinsicPriorsFns")
		}
	#					
	#figure out number of free params
	#numberParametersTotal <- dim(startingPriorsValues)[2] +  dim(intrinsicPriorsValues)[2] + dim(extrinsicPriorsValues)[2]
	#numberParametersFree <- numberParametersTotal
	#freevariables <- matrix(data = NA, nrow = 2, ncol = 0)
	#numberParametersStarting <- 0
	#numberParametersIntrinsic <- 0
	#numberParametersExtrinsic <- 0				
	#titlevector <- c()
	freevector <- logical()
	#
	#Calculate freevector
	for (i in 1:dim(startingPriorsValues)[2]) {
		if(dim(startingPriorsValues)[2])
		priorFn <- match.arg(arg = startingPriorsFns[i],
			choices = c("fixed", "uniform", "normal", "lognormal", "gamma", "exponential"),
			several.ok = FALSE)
		#
		if (priorFn == "fixed") {
			#numberParametersFree <- numberParametersFree-1
			parFree <- FALSE
		}else{
			#numberParametersStarting <- numberParametersStarting+1
			#freevariables <- cbind(freevariables, startingPriorsValues[, i])
			#titlevector  <- c(titlevector, paste("Starting", numberParametersStarting))
			parFree <- TRUE
			}
		names(parFree) <- paste0("starting_",i)
		freevector <- c(freevector, parFree)
		}
	for (j in 1:dim(intrinsicPriorsValues)[2]) {
		priorFn <- match.arg(arg = intrinsicPriorsFns[j],
			choices = c("fixed", "uniform", "normal", "lognormal", "gamma", "exponential"),
			several.ok = FALSE)
		#
		if (priorFn == "fixed") {
			#numberParametersFree <- numberParametersFree-1
			parFree <- FALSE
		}else{
			#numberParametersIntrinsic <- numberParametersIntrinsic+1
			#freevariables <- cbind(freevariables, intrinsicPriorsValues[, j])
			#titlevector  <- c(titlevector, paste("Intrinsic", numberParametersIntrinsic))
			parFree <- TRUE
			}
		names(parFree) <- paste0("intrinsic_",j)
		freevector <- c(freevector, parFree)
		}
	for (k in 1:dim(extrinsicPriorsValues)[2]) {
		priorFn <- match.arg(arg = extrinsicPriorsFns[k],
			choices = c("fixed", "uniform", "normal", "lognormal", "gamma", "exponential"),
			several.ok = FALSE)
		#
		if (priorFn == "fixed") {
			#numberParametersFree <- numberParametersFree-1
			parFree <- FALSE
		}else{
			#numberParametersExtrinsic <- numberParametersExtrinsic+1
			#freevariables <- cbind(freevariables, extrinsicPriorsValues[, k])
			#titlevector  <- c(titlevector, paste("Extrinsic", numberParametersExtrinsic))
			parFree <- TRUE
			}
		names(parFree) <- paste0("extrinsic_",k)
		freevector <- c(freevector, parFree)
		}
		
	#
	res <- freevector
	return(res)
	}


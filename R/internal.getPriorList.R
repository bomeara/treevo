# function for getting prior matrix for doRun_prc
getPriorList <- function(
    startingPriorsValues, 
    intrinsicPriorsValues, 
    extrinsicPriorsValues, 
    startingPriorsFns, 
    intrinsicPriorsFns, 
    extrinsicPriorsFns, 
    numberParametersTotal){
    #
    #
    #combing prior functions
	priorFuns <- c(startingPriorsFns, 
		intrinsicPriorsFns, extrinsicPriorsFns)
    #
	# combine prior values lists
	priorList <- c(startingPriorsValues, 
            intrinsicPriorsValues, 
            extrinsicPriorsValues
            )
	if(length(priorFuns)!=length(priorList)){
		stop("the number of priors doesn't match the number of priors given parameters as input")
		}
	# append function name to each element of priorList
	priorList <- lapply(1:length(priorList), function(x) 
			list(fun = priorFuns[[x]], params=priorList[[x]]))
	# create vector of names for priorList
    for (a in 1:length(startingPriorsValues)){
        namesForPriorList <- c(paste0("starting_", a, sep = ""))
        }
    #
    for (b in 1:length(intrinsicPriorsValues)){
        namesForPriorList <- append(namesForPriorList, paste0("intrinsic_", b, sep = ""))
        }
    #
    #message(extrinsicPriorsValues)
    for (c in 1:length(extrinsicPriorsValues)){
        namesForPriorList  <- append(namesForPriorList, paste0("extrinsic_", c, sep = ""))
        }	
	#
    names(priorList) <- namesForPriorList
    #rownames(priorList) <- c("shape", "value1", "value2")
    return(priorList)
    }

# internal function for calculating freevector

getFreeVector<-function(startingPriorsFns, startingPriorsValues, 
						intrinsicPriorsFns, intrinsicPriorsValues, extrinsicPriorsFns, extrinsicPriorsValues){
	
	#figure out number of free params
	#numberParametersTotal<-dim(startingPriorsValues)[2] +  dim(intrinsicPriorsValues)[2] + dim(extrinsicPriorsValues)[2]
	#numberParametersFree<-numberParametersTotal
	#freevariables<-matrix(data=NA, nrow=2, ncol=0)
	#numberParametersStarting<-0
	#numberParametersIntrinsic<-0
	#numberParametersExtrinsic<-0				
	#titlevector<-c()
	freevector<-c()
	
	#Calculate freevector
	for (i in 1:dim(startingPriorsValues)[2]) {
		priorFn<-match.arg(arg=startingPriorsFns[i],choices=c("fixed", "uniform", "normal", "lognormal", "gamma", "exponential"),several.ok=FALSE)
		if (priorFn=="fixed") {
			#numberParametersFree<-numberParametersFree-1
			freevector<-c(freevector, FALSE)
		}
		else {
			#numberParametersStarting<-numberParametersStarting+1
			#freevariables<-cbind(freevariables, startingPriorsValues[, i])
			#titlevector <-c(titlevector, paste("Starting", numberParametersStarting))
			freevector<-c(freevector, TRUE)
		}
	}
	for (i in 1:dim(intrinsicPriorsValues)[2]) {
		priorFn<-match.arg(arg=intrinsicPriorsFns[i],choices=c("fixed", "uniform", "normal", "lognormal", "gamma", "exponential"),several.ok=FALSE)
		if (priorFn=="fixed") {
			#numberParametersFree<-numberParametersFree-1
			freevector<-c(freevector, FALSE)
		}
		else {
			#numberParametersIntrinsic<-numberParametersIntrinsic+1
			freevariables<-cbind(freevariables, intrinsicPriorsValues[, i])
			#titlevector <-c(titlevector, paste("Intrinsic", numberParametersIntrinsic))
			freevector<-c(freevector, TRUE)
		}
	}
	for (i in 1:dim(extrinsicPriorsValues)[2]) {
		priorFn<-match.arg(arg=extrinsicPriorsFns[i],choices=c("fixed", "uniform", "normal", "lognormal", "gamma", "exponential"),several.ok=FALSE)
		if (priorFn=="fixed") {
			#numberParametersFree<-numberParametersFree-1
			freevector<-c(freevector, FALSE)
		}
		else {
			#numberParametersExtrinsic<-numberParametersExtrinsic+1
			#freevariables<-cbind(freevariables, extrinsicPriorsValues[, i])
			#titlevector <-c(titlevector, paste("Extrinsic", numberParametersExtrinsic))
			freevector<-c(freevector, TRUE)
			}
		}
		
	#
	res<-freevector
	return(res)
	}
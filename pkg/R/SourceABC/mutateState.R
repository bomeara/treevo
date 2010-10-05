mutateState<-function(startingState, standardDevFactor, priorValues, priorFn) {
	newState<-NA
	minBound=-Inf
	maxBound=Inf
	validNewState<-FALSE  #was lowercase, but not recognised 
	priorFn<-match.arg(arg=priorFn,choices=c("fixed", "uniform", "normal", "lognormal", "gamma", "exponential"),several.ok=FALSE);
	if (priorFn=="fixed" || priorFn=="uniform") {
		minBound<-min(priorValues)
		maxBound<-max(priorValues)
	}
	else if (priorFn=="lognormal" || priorFn=="gamma" || priorFn=="exponential") {
		minBound<-0
	}
	
	sdToUse<-standardDevFactor
	if (priorFn=="fixed") {
		sdToUse<-0
	}
	else if (priorFn=="uniform") {
		sdToUse<-standardDevFactor*(abs(max(priorValues)-min(priorValues)))
			if (sdToUse<0){
				print(paste("priorFN=", priorFn, "standardDevFactor=", standardDevFactor, "range(priorValues)=", range(priorValues)))
			}
	}
	else if (priorFn=="normal") {
		sdToUse<-standardDevFactor*priorValues[2]
			if (sdToUse<0){
				print(paste("priorFN=", priorFn, "standardDevFactor=", standardDevFactor, "range(priorValues)=", range(priorValues)))
			}
	}
	else if (priorFn=="lognormal") {
		sdToUse<-standardDevFactor*priorValues[2]
			if (sdToUse<0){
				print(paste("priorFN=", priorFn, "standardDevFactor=", standardDevFactor, "range(priorValues)=", range(priorValues)))
			}
	}
	else if (priorFn=="gamma") {
		sdToUse<-standardDevFactor*sqrt(priorValues[1]*priorValues[2]*priorValues[2])
			if (sdToUse<0){
				print(paste("priorFN=", priorFn, "standardDevFactor=", standardDevFactor, "range(priorValues)=", range(priorValues)))
			}
	}
	else if (priorFn=="exponential") {
		sdToUse<-standardDevFactor/priorValues[1]
			if (sdToUse<0){
				print(paste("priorFN=", priorFn, "standardDevFactor=", standardDevFactor, "range(priorValues)=", range(priorValues)))
			}
	}
	else {
		stop(priorFn," was not a recognized prior function")
	}

	while(!validNewState) {
		newState<-rnorm(n=1, mean=startingState, sd=sdToUse)
		validNewState<-TRUE
		if(is.na(newState)) {
			print(paste("MUTATESTATE_ERROR: newState = ",newState," sdToUse=",sdToUse," startingState=",startingState," priorFn=",priorFn," startingState=",startingState," priorValues=\n",sep=""))
			print(priorValues)
		}
		if (newState<minBound){
			validNewState<-FALSE
			}
		if (newState>maxBound){
			validNewState<-FALSE
			}	
		if (!validNewState)	{
			cat("newState ",newState," does not fit into one of the bounds (", minBound, "--", maxBound, ")\n")
		}	
	}
	
newState
}

meanFromPrior<-function(priorValues, priorFn) {
	#fixed, uniform, normal, lognormal, gamma, exponential
	x<-NA
	priorFn<-match.arg(arg=priorFn,choices=c("fixed", "uniform", "normal", "lognormal", "gamma", "exponential"),several.ok=FALSE);
	if (priorFn=="fixed") {
		x<-priorValues[1]
	}
	else if (priorFn=="uniform") {
		x<-mean(priorValues)
	}
	else if (priorFn=="normal") {
		x<-priorValues[1]
	}
	else if (priorFn=="lognormal") {
		x<- priorValues[1]
	}
	else if (priorFn=="gamma") {
		x<-priorValues[1]*priorValues[2]
	}
	else if (priorFn=="exponential") {
		x<-1.0/priorValues[1]
	}
	else {
		stop(priorFn," was not a recognized prior function")
	}
	x
}


pullFromPrior<-function(priorValues, priorFn) {
	#fixed, uniform, normal, lognormal, gamma
	x<-NA
	priorFn<-match.arg(arg=priorFn,choices=c("fixed", "uniform", "normal", "lognormal", "gamma", "exponential"),several.ok=FALSE);
	if (priorFn=="fixed") {
		x<-priorValues[1]
	}
	else if (priorFn=="uniform") {
		x<-runif(n=1,min=min(priorValues), max=max(priorValues))
	}
	else if (priorFn=="normal") {
		x<-rnorm(n=1, mean=priorValues[1], sd=priorValues[2])
	}
	else if (priorFn=="lognormal") {
		x<-rlnorm(n=1, meanlog = priorValues[1], sdlog = priorValues[2])
	}
	else if (priorFn=="gamma") {
		x<-rgamma(n=1, shape=priorValues[1], scale = priorValues[2])
	}
	else if (priorFn=="exponential") {
		x<-rexp(n=1, rate=priorValues[1])
	}
	else {
		stop(priorFn," was not a recognized prior function")
	}
	x
}

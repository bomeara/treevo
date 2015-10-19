#' getUnivariatePriorCurve
#' 
#' This function pulls random values from the prior
#' 
#' 
#' @param priorValues Variables needed to describe the shape of the
#' distribution.  uniform distributions=c(min, max); normal
#' distributions=c(mean, sd); lognormal distributions=c(mean, sd); gamma
#' distributions=c(shape, scale); exponential distributions=c(rate)
#' @param priorFn Prior shape; either "fixed", "uniform", "normal",
#' "lognormal", "gamma", "exponential"
#' @param nPoints Number of points
#' @param from Lower bound, if any
#' @param to Upper bound, if any
#' @param prob Probability content of HPD
#' @return Returns a list of x and y density coordinates, mean, and lower and
#' upper HPD
#' @author Barb Banbury and Brian O'Meara
#' @references O'Meara and Banbury, unpublished
#' @keywords getUnivariatePriorCurve
getUnivariatePriorCurve<-function(priorValues, priorFn, nPoints=100000, from=NULL, to=NULL, prob=0.95) {
	samples<-replicate(nPoints,pullFromPrior(priorValues, priorFn))
	if (is.null(from)) {
		from<-min(samples)
	}
	if (is.null(to)) {
		to<-max(samples)
	}
	result<-density(samples,from=from, to=to)
	hpd.result<-HPDinterval(as.mcmc(samples), prob)
  if (priorFn=="uniform") {
    result<-list(x=result$x,y=rep(max(result$y),length(result$x))) #kludge so uniform looks uniform 
  }
	return(list(x=result$x,y=result$y, mean=mean(samples), lower=hpd.result[1,1], upper=hpd.result[1,2]))
}



#' getUnivariatePosteriorCurve
#' 
#' This function pulls values from the posterior distribution
#' 
#' 
#' @param acceptedValues Vector of accepted particle values
#' @param from Lower bound, if any
#' @param to Upper bound, if any
#' @param prob Probability content of HPD
#' @return Returns a list of x and y density coordinates, mean, and lower and
#' upper HPD
#' @author Barb Banbury and Brian O'Meara
#' @references O'Meara and Banbury, unpublished
#' @keywords getUnivariatePosteriorCurve
getUnivariatePosteriorCurve<-function(acceptedValues, from=NULL, to=NULL, prob=0.95) {
	if (is.null(from)) {
		from<-min(acceptedValues)
	}
	if (is.null(to)) {
		to<-max(acceptedValues)
	}
	result<-density(acceptedValues,from=from, to=to)
	hpd.result<-HPDinterval(as.mcmc(acceptedValues), prob)
	return(list(x=result$x, y=result$y, mean=mean(acceptedValues), lower=hpd.result[1,1], upper=hpd.result[1,2]))
}



#' plotUnivariatePosteriorVsPrior
#' 
#' This function plots density distributions from the prior and posterior
#' 
#' 
#' @param posteriorCurve Kernal density estimates for the posterior
#' @param priorCurve Kernal density estimates for the prior
#' @param label X label for plot
#' @param trueValue True parameter value, if any
#' @param prob Probability content of HPD
#' @return Returns a plot
#' @author Barb Banbury and Brian O'Meara
#' @references O'Meara and Banbury, unpublished
#' @keywords plotUnivariatePosteriorVsPrior
plotUnivariatePosteriorVsPrior<-function(posteriorCurve, priorCurve, label="parameter", trueValue=NULL, prob=0.95) {
	plot(x=range(c(posteriorCurve$x, priorCurve$x)), y=range(c(posteriorCurve$y, priorCurve$y)), type="n", xlab=label, ylab="", bty="n", yaxt="n")
	polygon(x=c(priorCurve$x, max(priorCurve$x), priorCurve$x[1]), y=c(priorCurve$y, 0, 0), col=rgb(0,0,0,0.3),border=rgb(0,0,0,0.3))
	lines(x=rep(priorCurve$mean,2),y=c(0,max(c(priorCurve$y, posteriorCurve$y))),col="gray")
	lines(x=rep(priorCurve$lower,2),y=c(0,max(c(priorCurve$y, posteriorCurve$y))),col="gray",lty="dotted")
	lines(x=rep(priorCurve$upper,2),y=c(0,max(c(priorCurve$y, posteriorCurve$y))),col="gray",lty="dotted")
	
	polygon(x=c(posteriorCurve$x, max(posteriorCurve$x), posteriorCurve$x[1]), y=c(posteriorCurve$y,0, 0), col=rgb(1,0,0,0.3),border=rgb(1,0,0,0.3))
	lines(x=rep(posteriorCurve$mean,2),y=c(0,max(c(priorCurve$y, posteriorCurve$y))),col="red")
	lines(x=rep(posteriorCurve$lower,2),y=c(0,max(c(priorCurve$y, posteriorCurve$y))),col="red",lty="dotted")
	lines(x=rep(posteriorCurve$upper,2),y=c(0,max(c(priorCurve$y, posteriorCurve$y))),col="red",lty="dotted")
	if (!is.null(trueValue)) {
		lines(x=rep(trueValue,2),y=c(0,max(c(priorCurve$y, posteriorCurve$y))),col="blue")
	}	
}

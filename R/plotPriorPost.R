#' Plotting, Summarizing and Comparing the Prior and Posterior Distributions
#' 
#' Assorted functions for visualizing and summarizing the prior and posterior probability distributions associated with ABC analyses.
#' 
#' Function \code{plotPrior} visualizes the shape of various prior probability distributions available in TreEvo ABC analyses, while
#' function \code{plotUnivariatePosteriorVsPrior} plots the univariate density distributions from the prior and posterior against each other for comparison.
#' Function \code{getUnivariatePriorCurve} returns density coordinates and summary statistics from a selected prior probability distribution, while
#' function \code{getUnivariatePosteriorCurve} returns density coordinates and summary statistics from the posterior distribution of an ABC analysis.
#' The summaries from \code{getUnivariatePriorCurve} and \code{getUnivariatePosteriorCurve} are used as the input for \code{plotUnivariatePosteriorVsPrior}.
#' 

#' @param priorFn Prior Shape of the distribution; one of either "fixed", "uniform", "normal",
#' "lognormal", "gamma", or "exponential".

#' @param priorVariables Variables needed to describe the shape of the
#' distribution, dependent on \code{priorFn}:
#' \describe{

#' \item{priorFn = "uniform"}{priorVariables = c(min, max)}

#' \item{priorFn = "normal"}{priorVariables = c(mean, standard deviation)}

#' \item{priorFn = "lognormal"}{priorVariables = c(mean, standard deviation)}

#' \item{priorFn = "gamma"}{priorVariables = c(shape, scale)}

#' \item{priorFn = "exponential"}{priorVariables = c(rate)}

# \item{priorFn = "fixed"}{priorVariables = c(x)}

#' }

#' @param plotQuants If \code{TRUE}, plots line segments at the quantiles

#' @param plotLegend If \code{TRUE}, plots legend box with quantile values

#' @param nPoints Number of points to draw.

#' @param from Lower bound, if any.

#' @param to Upper bound, if any.

#' @param prob Probability content of the highest posterior density (HPD).

#' @param acceptedValues Vector of accepted particle values.

#' @param posteriorCurve Kernal density estimates for the posterior distribution from \code{getUnivariatePosteriorCurve}.

#' @param priorCurve Kernal density estimates for the prior distribution from \code{getUnivariatePriorCurve}.

#' @param label X label for plot.

#' @param trueValue True parameter value, if any such exists and is known (usually only true of simulations).

#' @return
#' \code{plotPrior} and \code{plotUnivariatePosteriorVsPrior} produce plots of the respective distributions (see above).
#' 
#' \code{getUnivariatePriorCurve} and \code{getUnivariatePosteriorCurve} returns a list of x and y density coordinates, mean, and lower and
#' upper highest posterior density (HPD), for their respective distribution.
#' 

#' @seealso
#' \code{\link{plotPosteriors}} Plots multiple posteriors against their priors and potential known values.

#' @author Brian O'Meara and Barb Banbury

# @references O'Meara and Banbury, unpublished
# @keywords getUnivariatePosteriorCurve
# @references O'Meara and Banbury, unpublished
# @keywords plotPrior prior
# @references O'Meara and Banbury, unpublished
# @keywords getUnivariatePriorCurve
# @references O'Meara and Banbury, unpublished
# @keywords plotUnivariatePosteriorVsPrior



#' @examples
#' 
#' data(simRunExample)
#' 
#' # examples with plotPrior
#' 
#' 	plotPrior(priorFn="exponential", priorVariables=c(10))
#' 	
#' 	plotPrior(priorFn="normal", priorVariables=c(1,2))
#' 	
#' 	plotPrior(priorFn="gamma", priorVariables=c(2, .2), plotQuants=FALSE, plotLegend=FALSE)
#' 
#' # examples of getting density coordinates and summary statistics from distributions
#' 
#' priorKernal<-getUnivariatePriorCurve(priorFn="normal", priorVariables=c(28,2),
#' 	nPoints=100000, from=NULL, to=NULL, prob=0.95)
#' 
#' postKernal<-getUnivariatePosteriorCurve(acceptedValues=results[[1]]$particleDataFrame$starting_1,
#' 	from=NULL, to=NULL, prob=0.95)
#' 
#' priorKernal
#' postKernal
#' 
#' # let's compare this (supposed) prior against the posterior in a plot
#' 
#' plotUnivariatePosteriorVsPrior(posteriorCurve=postKernal, priorCurve=priorKernal,
#' 	label="parameter", trueValue=NULL, prob=0.95)
#' 
#' # cool!
#' 

# priorVariables depend on priorFn.  uniform=c(min, max); normal=c(mean, stdev);
# lognormal=c(mean, stdev), gamma=c(shape, scale), exponential=c(rate)
#
# plotPrior("normal", c(1,2), plotQuants=TRUE)
# plotPrior("exponential", c(1), plotQuants=FALSE)

#' @name plotPriorPost
#' @rdname plotPriorPost
#' @export
plotPrior<-function(
	priorFn=match.arg(arg=priorFn,choices=c("fixed", "uniform", "normal", 
		"lognormal", "gamma", "exponential"), several.ok=FALSE),
	priorVariables, plotQuants=TRUE, plotLegend=TRUE){
	#plot.new()
	x<-NA
	quant<-c(0.01, 0.05, 0.25, .50, 0.75, 0.95, 0.99)
	quant.value<-vector()
	mm<-vector()
	if (priorFn=="fixed") {
		cat ("can't plot fixed prior")
	}else{
		if (priorFn=="uniform") {
			min=priorVariables[1]
			max=priorVariables[2]
			x<-runif(1000, min, max)
			curve(dunif(x), xlim=c(min-(.3*(max-min)), max+(.3*(max-min))), ylim=c(0, 1.5),
				type="n", xlab="", ylab=paste("Density (parameters: min=", min, "; max =", max, ")", sep=""))
			#title(priorFn)
			rect(min, 0, max, 1, col=rgb(0, 0, 0, .2))
			for (i in 1:length(quant)) {
				quant.value[i]<-qunif(quant[i], min, max)
				
				if (plotQuants) {
					#mm[i]<-dunif(quant.value[i], min, max)
					segments(quant.value[i], 0, quant.value[i], 1)
					segments(qunif(.5, min, max), 0, qunif(.5, min, max), 1, lwd=2)
				}
			}
		}else{
			if (priorFn=="normal") {
				mean=priorVariables[1]
				stdev=priorVariables[2]
				x<-rnorm(1000, mean, stdev)
				poly<-curve(dnorm(x, mean, stdev), from=min(x), to=max(x), xlab="",
					ylab=paste("Density (parameters: mean=", mean, "; stdev =",stdev, ")", sep=""))
				poly$x<-c(min(poly$x),poly$x, max(poly$x))
				poly$y<-c(min(poly$y), poly$y, min(poly$y))
				polygon(poly, col=rgb(0, 0, 0, 0.3))  #to include 95%?
				for (i in 1:length(quant)) {
					quant.value[i]<-qnorm(quant[i], mean, stdev)
					if (plotQuants) {
						mm[i]<-dnorm(quant.value[i], mean, stdev)
						segments(quant.value[i], min(poly$y), quant.value[i], mm[i])
						segments(qnorm(.5, mean, stdev), min(poly$y), qnorm(.5, mean, stdev),
							dnorm(qnorm(.5, mean, stdev), mean, stdev), lwd=2)
					}
				}
			}else{
				if (priorFn=="lognormal") {  #messed up quant lines and polygon
					mean=priorVariables[1]
					stdev=priorVariables[2]
					x<-rlnorm(1000, mean, stdev)
					poly<-curve(dlnorm(x, mean, stdev), from=0, to=qlnorm(0.99, mean, stdev),
						xlab="", ylab=paste("Density (parameters: mean=", mean, "; stdev =",stdev, ")", sep=""))
					poly$x<-c(poly$x, max(poly$x))
					poly$y<-c(poly$y, min(poly$y))
					polygon(poly, col=rgb(0, 0, 0, 0.3))  #to include 95%?
					for (i in 1:length(quant)) {
						quant.value[i]<-qlnorm(quant[i], mean, stdev)
						if(plotQuants) {
							mm[i]<-dlnorm(quant.value[i], mean, stdev)
							segments(quant.value[i], min(poly$y), quant.value[i], mm[i])
							segments(qlnorm(.5, mean, stdev), min(poly$y), qlnorm(.5, mean, stdev),
								dlnorm(qlnorm(.5, mean, stdev), mean, stdev), lwd=2)
						}
					}
				}else{
					if (priorFn=="gamma") {
						shape=priorVariables[1]
						scale=priorVariables[2]
						x<-rgamma(1000, shape, scale)
						poly<-curve(dgamma(x, shape, scale), from=0, to=qgamma(0.99, shape, scale),
							xlab="", ylab=paste("Density (parameters: shape=", shape, "; scale =",scale, ")", sep=""))
						poly$x<-c(0, poly$x, max(poly$x), 0)
						poly$y<-c(0, poly$y, min(poly$y), 0)
						#title(priorFn)
						polygon(poly, col=rgb(0, 0, 0, 0.3))
						for (i in 1:length(quant)) {
							quant.value[i]<-qgamma(quant[i], shape, scale)
							if (plotQuants) {
								mm[i]<-dgamma(quant.value[i], shape, scale)
								segments(quant.value[i], min(poly$y), quant.value[i], mm[i])
								segments(qgamma(.5, shape, scale), min(poly$y), qgamma(.5, shape, scale),
									dgamma(qgamma(.5, shape, scale), shape, scale), lwd=2)
							}
						}
					}else{
						if (priorFn=="exponential") {
							rate=priorVariables[1]
							x<-rexp(1000, rate)
							poly<-curve(dexp(x, rate), from=0, to=qexp(0.99, rate), xlab="",
								ylab=paste("Density (parameters: rate=", rate, ")", sep=""))
							poly$x<-c(poly$x, 0)
							poly$y<-c(poly$y, min(poly$y))
							polygon(poly, col=rgb(0, 0, 0, 0.3))
							for (i in 1:length(quant)) {
								quant.value[i]<-qexp(quant[i], rate)
								if (plotQuants) {
									mm[i]<-dexp(quant.value[i], rate)
									segments(quant.value[i], min(poly$y), quant.value[i], mm[i])
									segments(qexp(.5, rate), min(poly$y), qexp(.5, rate), dexp(qexp(.5, rate), rate), lwd=2)
									}
								}
							}
						}
					}
				}
			}
		}
	results<-data.frame(cbind(quant, quant.value))
	if (plotLegend){
		legend("topright", legend=paste(c(quant, signif(quant.value, digits=3))), title="Quantiles", ncol=2, bty="n")
		}
	return(results)			
	}

#' @rdname plotPriorPost
#' @export
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




#' @rdname plotPriorPost
#' @export
getUnivariatePriorCurve<-function(priorFn, priorVariables,
		nPoints=100000, from=NULL, to=NULL, prob=0.95) {
	#
	samples<-replicate(nPoints,pullFromPrior(priorVariables, priorFn))
	if (is.null(from)) {
		from<-min(samples)
	}
	if (is.null(to)) {
		to<-max(samples)
	}
	result<-density(samples,from=from, to=to)
	hpd.result<-coda::HPDinterval(coda::as.mcmc(samples), prob)
  if (priorFn=="uniform") {
    result<-list(x=result$x,y=rep(max(result$y),length(result$x))) #kludge so uniform looks uniform
  }
	return(list(x=result$x,y=result$y, mean=mean(samples), lower=hpd.result[1,1], upper=hpd.result[1,2]))
}




#' @rdname plotPriorPost
#' @export
getUnivariatePosteriorCurve<-function(acceptedValues, from=NULL, to=NULL, prob=0.95) {
	if (is.null(from)) {
		from<-min(acceptedValues)
	}
	if (is.null(to)) {
		to<-max(acceptedValues)
	}
	result<-density(acceptedValues,from=from, to=to)
	hpd.result<-coda::HPDinterval(coda::as.mcmc(acceptedValues), prob)
	return(list(x=result$x, y=result$y, mean=mean(acceptedValues), lower=hpd.result[1,1], upper=hpd.result[1,2]))
}







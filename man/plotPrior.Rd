\name{plotPrior}

\title{Plot prior distribution}

\description{This function plots the prior of several distributions}

\usage{
plotPrior(priorFn=match.arg(arg=priorFn,choices=c("fixed", "uniform", "normal", "lognormal", "gamma", "exponential"),several.ok=FALSE), priorVariables, plotQuants=TRUE, plotLegend=TRUE)
}

\arguments{
  \item{priorFn}{Either "fixed", "uniform", "normal", "lognormal", "gamma", "exponential"}
  \item{priorVariables}{One or two variables needed to describe the shape of the distribution.  uniform distributions=c(min, max); normal distributions=c(mean, sd); lognormal distributions=c(mean, sd); gamma distributions=c(shape, scale); exponential distributions=c(rate)}
  \item{plotQuants}{Plots line segments at the quantiles}
  \item{plotLegend}{Plots legend box with quantile values}
}
\details{
	This functions is to help visualize the shape of the priors for TreEvo abc analyses.  
}
\value{
	Returns a distribution plot
}
\references{O'Meara and Banbury, unpublished}
\author{Brian O'Meara and Barb Banbury}
\seealso{
}
\examples{
	plotPrior(priorFn="exponential", priorVariables=c(10))
	
	plotPrior("normal", c(1,2))
	
	plotPrior("gamma", c(2, .2), plotQuants=F, plotLegend=F)
}
\keyword{}


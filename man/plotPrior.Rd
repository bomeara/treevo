\name{plotPrior}
\alias{plotPrior}


\title{Plot prior distribution}

\description{This function plots the distribution of priors}

\usage{
plotPrior(priorFn=match.arg(arg=priorFn,choices=c("fixed", "uniform", "normal", "lognormal", "gamma", "exponential"),several.ok=FALSE), priorVariables, plotQuants=TRUE, plotLegend=TRUE)
}

\arguments{
  \item{priorFn}{Prior shape; either "fixed", "uniform", "normal", "lognormal", "gamma", "exponential"}
  \item{priorVariables}{Variables needed to describe the shape of the distribution.  uniform distributions=c(min, max); normal distributions=c(mean, sd); lognormal distributions=c(mean, sd); gamma distributions=c(shape, scale); exponential distributions=c(rate)}
  \item{plotQuants}{If TRUE, plots line segments at the quantiles}
  \item{plotLegend}{If TRUE, plots legend box with quantile values}
}
\details{
	This functions is to help visualize the shape of the priors for TreEvo abc analyses.  
}
\value{
	Returns a distribution plot
}
\references{O'Meara and Banbury, unpublished}
\author{Brian O'Meara and Barb Banbury}

\examples{
	plotPrior(priorFn="exponential", priorVariables=c(10))
	
	plotPrior("normal", c(1,2))
	
	plotPrior("gamma", c(2, .2), plotQuants=FALSE, plotLegend=FALSE)
}
\keyword{plotPrior}
\keyword{prior}



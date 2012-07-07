\name{getUnivariatePriorCurve}
\alias{getUnivariatePriorCurve}


\title{getUnivariatePriorCurve}

\description{This function pulls random values from the prior}

\usage{getUnivariatePriorCurve(priorValues, priorFn, nPoints=100000, from=NULL, to=NULL, prob=0.95)}

\arguments{
  \item{priorValues}{Variables needed to describe the shape of the distribution.  uniform distributions=c(min, max); normal distributions=c(mean, sd); lognormal distributions=c(mean, sd); gamma distributions=c(shape, scale); exponential distributions=c(rate)}
  \item{priorFn}{Prior shape; either "fixed", "uniform", "normal", "lognormal", "gamma", "exponential"}
  \item{from}{Lower bound, if any}
  \item{to}{Upper bound, if any}
  \item{prob}{Probability content of HPD}
}
\value{
	Returns a list of x and y density coordinates, mean, and lower and upper HPD  
}

\references{O'Meara and Banbury, unpublished}
\author{Barb Banbury and Brian O'Meara}


\keyword{getUnivariatePriorCurve}


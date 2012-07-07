\name{plotUnivariatePosteriorVsPrior}
\alias{plotUnivariatePosteriorVsPrior}


\title{plotUnivariatePosteriorVsPrior}

\description{This function plots density distributions from the prior and posterior}

\usage{plotUnivariatePosteriorVsPrior(posteriorCurve, priorCurve, label="parameter", trueValue=NULL, prob=0.95)}

\arguments{
  \item{posteriorCurve}{Kernal density estimates for the posterior}
  \item{priorCurve}{Kernal density estimates for the prior}
  \item{label}{X label for plot}
  \item{trueValue}{True parameter value, if any}
  \item{prob}{Probability content of HPD}
}
\value{
	Returns a plot
}

\references{O'Meara and Banbury, unpublished}
\author{Barb Banbury and Brian O'Meara}


\keyword{plotUnivariatePosteriorVsPrior}


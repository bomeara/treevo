\name{getUnivariatePosteriorCurve}
\alias{getUnivariatePosteriorCurve}


\title{getUnivariatePosteriorCurve}

\description{This function pulls values from the posterior distribution}

\usage{getUnivariatePosteriorCurve(acceptedValues, from=NULL, to=NULL, prob=0.95)}

\arguments{
  \item{acceptedValues}{Vector of accepted particle values}
  \item{from}{Lower bound, if any}
  \item{to}{Upper bound, if any}
  \item{prob}{Probability content of HPD}
}
\value{
	Returns a list of x and y density coordinates, mean, and lower and upper HPD  
}

\references{O'Meara and Banbury, unpublished}
\author{Barb Banbury and Brian O'Meara}


\keyword{getUnivariatePosteriorCurve}


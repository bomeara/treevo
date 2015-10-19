\name{boxcoxTransformation}
\alias{boxcoxTransformation}

\title{Box-Cox Transformation}

\description{This function Box-Cox transforms summary values from a single simulation}

\usage{
boxcoxTransformation(summaryValues, boxcoxAddition, boxcoxLambda)
}

\arguments{
  \item{summaryValues}{Vector of summary statistics}
  \item{boxcoxAddition}{Vector of boxcox additions from boxcoxEstimation}
  \item{boxcoxLambda}{Vector of boxcox lambda values from boxcoxEstimation}
}
\value{
	Returns a vector of Box-Cox transformed summary statistics from a single observation.
}
\references{O'Meara and Banbury, unpublished; Bates et al. 2009}
\author{Brian O'Meara and Barb Banbury}

\keyword{boxcoxTransformation}
\keyword{Box-Cox}


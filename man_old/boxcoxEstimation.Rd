\name{boxcoxEstimation}
\alias{boxcoxEstimation}

\title{Box-Cox Estimation}

\description{This function Box-Cox transforms summary values}

\usage{
boxcoxEstimation(summaryValuesMatrix)
}

\arguments{
  \item{summaryValuesMatrix}{Matrix of summary statistics from simulations}
}
\value{
	Returns a matrix of Box-Cox transformed summary statistics with the same dimensions as summaryValuesMatrix.
}
\references{O'Meara and Banbury, unpublished; Bates et al. 2009}
\author{Brian O'Meara and Barb Banbury}

\examples{
#data(res)

#boxcoxEstimation(summaryValuesMatrix)
}

\keyword{boxcoxEstimation}
\keyword{Box-Cox}


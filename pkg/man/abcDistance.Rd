\name{abcDistance}
\alias{abcDistance}

\title{Calculate ABC Distances}

\description{
	This function uses results from pls model to calculate the distance from each sim to the original data.
}

\usage{
abcDistance(summaryValuesMatrix, originalSummaryValues, pls.model.list)
}
\arguments{
  \item{summaryValuesMatrix}{Matrix of summary statistics from simulations}
  \item{originalSummaryValues}{Original summary stats}
  \item{pls.model.list}{indexing of a list of pls models for different params}
}
\details{
This function runs a PLS regression on each free param in the model, unlike a true multivariate PLS regression.  For ABC, this seems to result in much better results, without one parameter dominating the combined variance. 
}	
\value{
	Returns euclidean distance of each simulation's summary values to the original summary stats.
}	
\references{O'Meara and Banbury, unpublished}
\author{Brian O'Meara and Barb Banbury}

\examples{
#abcDistance(summaryValuesMatrix, originalSummaryValues, pls.model.list)
}
\keyword{abcDistance}

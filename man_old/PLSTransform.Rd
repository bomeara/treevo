\name{PLSTransform}
\alias{PLSTransform}

\title{PLS Transform}

\description{Uses results from the PLS model to transform summary values}

\usage{
PLSTransform(summaryValuesMatrix, pls.model)
}
\arguments{
  \item{summaryValuesMatrix}{Matrix of summary statistics from simulations}
  \item{pls.model}{results from returnPLSModel}
}
\details{
	This function uses the scores from the pls model to transform the summary statistics. 
}	
\value{
	Returns transformed summary statistics
}	
\references{O'Meara and Banbury, unpublished}
\author{Brian O'Meara and Barb Banbury}

\examples{
#PLSTransform(summaryValuesMatrix, pls.model)
}
\keyword{PLSTransform}
\keyword{PLS}

\name{returnPLSModel}
\alias{returnPLSModel}

\title{return PLS Model}

\description{Returns a PLS Model}

\usage{
returnPLSModel(trueFreeValuesMatrix, summaryValuesMatrix, validation="CV", scale=TRUE, variance.cutoff=95)
}
\arguments{
  \item{trueFreeValuesMatrix}{Matrix of true free values from simulations}
  \item{summaryValuesMatrix}{Matrix of summary statistics from simulations}
  \item{validation}{Cross Validation procedure for abc}
  \item{scale}{scale for pls.model.list}
  \item{variance.cutoff}{variance cutoff for pls.model.list}
}
\details{
This function runs a PLS regression on each free param in the model, unlike a true multivariate PLS regression.  For ABC, this seems to result in much better results, without one parameter dominating the combined variance. 
}	
\value{
	Returns a pls model	
}	
\references{O'Meara and Banbury, unpublished}
\author{Brian O'Meara and Barb Banbury}

\examples{
#returnPLSModel(trueFreeValuesMatrix,summaryValuesMatrix, validation="CV", scale=scale, variance.cutoff=variance.cutoff)
}
\keyword{returnPLSModel}
\keyword{PLSTransform}
\keyword{PLS}

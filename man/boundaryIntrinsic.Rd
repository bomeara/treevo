\name{abcmodels.intrinsic.R}
\alias{boundaryIntrinsic}

\title{Intrinsic Character Evolution Models}

\description{This function describes a model of intrinsic character evolution. Character change is restricted above minimum and below maximum thresholds}

\usage{
boundaryIntrinsic(params, states, timefrompresent)
}

\arguments{
  \item{params}{describes input paramaters for the model. \code{boundaryIntrinsic} params = sd, minimum, maximum}   
  \item{states}{vector of states for each taxon}
  \item{timefrompresent}{which time slice in the tree}
}
\details{
}
\value{
	A matrix of values representing character displacement from a single time step in the tree.
	}
\references{O'Meara and Banbury, unpublished}
\author{Brian O'Meara and Barb Banbury}
\seealso{
}
\examples{

}
\keyword{}


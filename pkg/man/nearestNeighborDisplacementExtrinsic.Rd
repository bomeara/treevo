\name{nearestNeighborDisplacementExtrinsic}
\alias{nearestNeighborDisplacementExtrinsic}

\title{Extrinsic Character Evolution Models}

\description{This function describes a model of extrinsic character evolution.  Character values of a focal taxon depend on values of closest relatives on the tree}

\usage{
nearestNeighborDisplacementExtrinsic(params,selfstates,otherstates, timefrompresent)
}

\arguments{
  \item{params}{describes input paramaters for the model. \code{nearestNeighborDisplacementExtrinsic} params = sd, springK, maximum force}  \item{selfstates}{vector of states for each taxon}
  \item{otherstates}{matrix of character states, one row per taxon and once column per state}
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


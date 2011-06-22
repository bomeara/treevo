\name{abcmodels.extrinsic.R}
\alias{nullExtrinsic}
\alias{nearestNeighborDisplacementExtrinsic}
\alias{ExponentiallyDecayingPush}
\alias{everyoneDisplacementExtrinsic}
\title{Extrinsic Character Evolution Models}

\description{These functions describe models of character evolution by extrinsic forces.  Each function generates new character values after one time step: \code{nullExtrinsic}; no extrinsic character change; \code{nearestNeighborDisplacementExtrinsic}: character value pushes away from closest relative on the tree; \code{ExponentiallyDecayingPush}: new character value pushes away harder from taxa with like values; and \code{everyoneDisplacementExtrinsic}: character values push away from all other taxa}

\usage{
nullExtrinsic(params,selfstates,otherstates, timefrompresent)
nearestNeighborDisplacementExtrinsic(params,selfstates,otherstates, timefrompresent)
ExponentiallyDecayingPush(params,selfstates,otherstates, timefrompresent)
everyoneDisplacementExtrinsic(params,selfstates,otherstates, timefrompresent)
}

\arguments{
  \item{params}{describes input paramaters for the model. \code{nearestNeighborDisplacementExtrinsic} params = sd, springK, maximum force; \code{ExponentiallyDecayingPush} params = sd, maximum force, half distance; and \code{everyoneDisplacementExtrinsic} sd, springK, maximum force}
  \item{selfstates}{vector of states for each taxon}
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


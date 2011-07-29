\name{everyoneDisplacementExtrinsic}
\alias{everyoneDisplacementExtrinsic}
\title{Extrinsic Character Evolution Models}

\description{This function describes a model of extrinsic character evolution.  Character values push away from all other taxa}

\usage{
everyoneDisplacementExtrinsic(params,selfstates,otherstates, timefrompresent)
}

\arguments{
  \item{params}{describes input paramaters for the model. \code{everyoneDisplacementExtrinsic} sd, springK, maximum force}
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


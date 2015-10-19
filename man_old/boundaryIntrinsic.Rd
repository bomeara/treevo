\name{boundaryIntrinsic}
\alias{boundaryIntrinsic}

\title{Intrinsic Character Evolution Models}

\description{This function describes a model of intrinsic character evolution. Character change is restricted above a minimum and below a maximum threshold}

\usage{
boundaryIntrinsic(params, states, timefrompresent)
}

\arguments{
  \item{params}{describes input paramaters for the model. \code{boundaryMinIntrinsic} params = sd, minimum, maximum}   
  \item{states}{vector of states for each taxon}
  \item{timefrompresent}{which time slice in the tree}
}

\value{
	A matrix of values representing character displacement from a single time step in the tree.
}

\references{O'Meara and Banbury, unpublished}
\author{Brian O'Meara and Barb Banbury}

\keyword{boundaryIntrinsic}
\keyword{intrinsic}


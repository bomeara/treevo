\name{nullIntrinsic}
\alias{nullIntrinsic}

\title{Intrinsic Character Evolution Models}

\description{This function describes a model of no intrinsic character change}

\usage{
nullIntrinsic(params,states, timefrompresent)
}

\arguments{
  \item{params}{describes input paramaters for the model}
  \item{states}{vector of states for each taxon}
  \item{timefrompresent}{which time slice in the tree}
}

\value{
	A matrix of values representing character displacement from a single time step in the tree.
}

\references{O'Meara and Banbury, unpublished}
\author{Brian O'Meara and Barb Banbury}

\keyword{nullIntrinsic}
\keyword{intrinsic}


\name{brownianIntrinsic}
\alias{brownianIntrinsic}

\title{Intrinsic Character Evolution Models}

\description{This function describes a model of intrinsic character evolution via Brownian motion.}

\usage{
brownianIntrinsic(params, states, timefrompresent)
}

\arguments{
  \item{params}{describes input paramaters for the model. \code{boundaryIntrinsic} params = sd}   
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


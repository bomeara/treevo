\name{autoregressiveIntrinsicTimeSlices}
\alias{autoregressiveIntrinsicTimeSlices}

\title{Intrinsic Character Evolution Models}

\description{This function describes a model of intrinsic character evolution. New character values are generated after one time step via a discrete-time OU process with differing means, sigma, and attraction over time}

\usage{
autoregressiveIntrinsicTimeSlices(params,states, timefrompresent)
}

\arguments{
  \item{params}{describes input paramaters for the model. \code{autoregressiveIntrinsicTimeSlices} params = sd-1 (sigma-1), attractor-1 (character mean-1), attraction-1 (alpha-1), time threshold-1, sd-2 (sigma-2), attractor-2 (character mean-2), attraction-2 (alpha-2), time threshold-2}   
  \item{states}{vector of states for each taxon}
  \item{timefrompresent}{which time slice in the tree}
}
\details{
	In the TimeSlices models, time threshold units are in time before present (i.e., 65 could be 65 MYA). The last time threshold should be 0.
}
\value{
	A matrix of values representing character displacement from a single time step in the tree.
}

\references{O'Meara and Banbury, unpublished}
\author{Brian O'Meara and Barb Banbury}

\keyword{autoregressiveIntrinsicTimeSlices}
\keyword{intrinsic}


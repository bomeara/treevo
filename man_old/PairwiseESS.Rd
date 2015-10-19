\name{PairwiseESS}
\alias{PairwiseESS}
\title{Pairwise Effective Sample Size}
\description{
	This function calculates ESS on results.  Performs the best when results are from multiple runs.   
}
\usage{PairwiseESS(particleDataFrame)}
\arguments{
  \item{particleDataFrame}{particleDataFrame can be a single data frame or a list of data frames}
}
\value{Returns a matrix with ESS values of all pairwise runs}

\references{O'Meara and Banbury, unpublished}
\author{Brian O'Meara and Barb Banbury}

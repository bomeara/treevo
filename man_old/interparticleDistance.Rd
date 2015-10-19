\name{interparticleDistance}
\alias{interparticleDistance}
\title{Calculate interparticle distance}
\description{
	This function calculates interparticle distance between two vectors     
}
\usage{interparticleDistance(x,y,abs=TRUE)}
\arguments{
  \item{x}{Vector of parameter values}
  \item{y}{Vector of parameter values}
  \item{abs}{If TRUE, uses absolute value for distance}
}
\value{Returns a matrix of distances between pairwise points in the vectors}
\details{IPD is used by several internal TreEvo function to compare two distributions.}

\references{O'Meara and Banbury, unpublished}
\author{Brian O'Meara and Barb Banbury}

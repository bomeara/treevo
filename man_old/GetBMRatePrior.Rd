\name{GetBMRatePrior}
\alias{GetBMRatePrior}

\title{Get BM Rate Prior}

\description{This function automatically calculates prior distributions for BM model of evolution}

\usage{
GetBMRatePrior(phy, traits, timeStep)
}
\arguments{
  \item{phy}{Tree (Phylogenetic tree in phylo format)}
  \item{traits}{data matrix with rownames equal to phy}
  \item{timeStep}{time in a single iteration of the discrete-time simulation}
}
\details{
	Returns a matrix of prior values that can be used in the doRun functions.  Builds on functions in geiger to estimate distribution.  
}
\value{
	Returns a matrix of prior values
}
\references{O'Meara and Banbury, unpublished}
\author{Brian O'Meara and Barb Banbury}

\examples{
#GetBMRatePrior(phy, traits, timeStep)

}
\keyword{GetBMRatePrior}


\name{PLSRejection}
\alias{PLSRejection}

\title{PLS Rejection}

\description{This function automatically calculates prior distributions for BM model of evolution}

\usage{
PLSRejection(summaryValuesMatrix, trueFreeValuesMatrix, phy, traits, abcTolerance, verbose=TRUE, validation="CV", scale=TRUE, variance.cutoff=95)
}
\arguments{
  \item{summaryValuesMatrix}{Matrix of summary statistics from simulations}
  \item{trueFreeValuesMatrix}{Matrix of true free values from simulations}
  \item{phy}{Tree (Phylogenetic tree in phylo format)}
  \item{traits}{data matrix with rownames equal to phy}
  \item{abcTolerance}{Proportion of accepted simulations}
  \item{verbose}{option to print progress to screen}
  \item{validation}{Cross Validation procedure for abc}
  \item{scale}{scale for pls.model.list}
  \item{variance.cutoff}{variance cutoff for pls.model.list}
}
\details{
This function performs the abc-rejection analysis using an input simulation data. Particles are accepted is they fall sufficiently close to the target data (within the tolerance). Distances are calculated using abcDistance. 
}	
\value{
	Returns a list of the particle data frame and abc distances.
}
\references{O'Meara and Banbury, unpublished}
\author{Brian O'Meara and Barb Banbury}

\examples{
#PLSRejection(summaryValuesMatrix, trueFreeValuesMatrix, phy, traits, abcTolerance)
}
\keyword{PLSRejection}
\keyword{doRun}
\keyword{doRun_rej}
\keyword{abc}

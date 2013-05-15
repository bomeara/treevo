\name{summaryStatsLong}
\alias{summaryStatsLong}

\title{Calculate summary statistics}

\description{This function creates a vector of summary statistics for TreEvo analysis}

\usage{
summaryStatsLong(phy, traits)
}

\arguments{
  \item{phy}{Tree (Phylogenetic tree in phylo format)}
  \item{traits}{data matrix with rownames equal to phy}
}
\details{
	Caluclates 17 summary statistics from the fitContinuous.hacked function (brown.lnl, brown.beta, brown.aic, lambda.lnl, lambda.beta, lambda.lambda, lambda.aic, delta.lnl, delta.beta, delta.delta, delta.aic, ou.lnl, ou.beta, ou.alpha, ou.aic, white.lnl, white.aic), plus  raw.mean, raw.max, raw.min, raw.var, raw.median, and all tip character values, phylogenetic independent contrasts, ancestral state reconstruction values, and the range of ancestral state reconstruction confidence interval. 
}
\value{
	Returns a vector of summary statistics
}
\references{O'Meara and Banbury, unpublished}
\author{Brian O'Meara and Barb Banbury}

\examples{
#summaryStatsLong(phy, char)
}
\keyword{summaryStatsLong}


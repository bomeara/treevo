\name{summaryStatsLong}

\title{Calculate summary statistics}

\description{This function creates a vector of summary statistics, including 22 summary statistics and raw tip data}

\usage{
summaryStatsLong(phy, data, todo=c(), jobName="")
}

\arguments{
  \item{phy}{Tree (ape phylo format)}
  \item{data}{Data frame of character data}
  \item{todo}{A vector (length = 22 + number of taxa) representing the value of the summary statistic (1 = valuable, 0 = not explanatory)}
  \item{jobName}{Optional name for the job}
}
\details{
	Caluclates 22 summary statistics from fitContinuous.hacked function (brown.lnl, brown.beta, brown.aic, lambda.lnl, lambda.beta, lambda.lambda, lambda.aic, delta.lnl, delta.beta, delta.delta, delta.aic, ou.lnl, ou.beta, ou.alpha, ou.aic, white.lnl, white.aic), plus  raw.mean, raw.max, raw.min, raw.var, raw.median, and all tip character values.  
}
\value{
	Returns a vector of summary statistics
}
\references{O'Meara and Banbury, unpublished}
\author{Brian O'Meara and Barb Banbury}
\seealso{
}
\examples{

}
\keyword{}


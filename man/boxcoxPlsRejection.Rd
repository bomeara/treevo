\name{boxcoxPlsRejection}
\alias{boxcoxPlsRejection}
\title{Approximate Bayesian computation for comparative methods-Rejection}
\description{
	Analyzes simulations for rejection 
}
\usage{
boxcoxPlsRejection(summaryValuesMatrix, trueFreeValuesMatrix, phy, traits, vipthresh, abcTolerance, verbose=TRUE) 
}
\arguments{
  \item{summaryValuesMatrix}{Matrix of summary statistics from simulations}
  \item{trueFreeValuesMatrix}{Matrix of parameter values from simulations}
  \item{phy}{Tree (Phylogenetic tree in phylo format)}
  \item{traits}{data matrix with rownames equal to phy}
  \item{vipthresh}{Threshold value for pls regression of vip summary stats}
  \item{abcTolerance}{Proportion of accepted simulations}
  \item{verbose}{Optional screen updates}

}

\details{
	This function is an internal function used by doRun_rej that performs an abc-rejection analysis using simulations, but can also be used independently using simulations as input. Simulations are boxcox transformed and a pls regression is performed for each free parameter to determine the most informative summary statistics (vipthresh). The euclidean distance of each initial simulation's vip summary stats to the input data is calculated, and those that fall under the abcTolerance are kept as accepted particles. These describe the posterior distributions of parameters.  
	
}

\value{

particleDataFrame=particleDataFrame, abcDistancesRaw=abcDistancesRaw, whichVip=whichVip, boxcoxSummaryValuesMatrix=boxcoxSummaryValuesMatrix, boxcoxOriginalSummaryStats=boxcoxOriginalSummaryStats

  \item{particleDataFrame}{DataFrame with information from each simulation, including generation, attempt, id, parentid, distance, weight, and parameter states}
  \item{abcDistancesRaw}{Euclidean distances for each simulation and free parameter}
  \item{whichVip}{Matrix of vip summary statistics for each free parameter}
  \item{boxcoxSummaryValuesMatrix}{Boxcox transformed summary statistics for simulations}
  \item{boxcoxOriginalSummaryStats}{Boxcox transformed summary statistics for original data}
  
}

\references{O'Meara and Banbury, unpublished; Sisson et al. 2007, Wegmann et al. 2009}
\author{Brian O'Meara and Barb Banbury}
\examples{
#data(res)

#boxcoxPlsRejection(summaryValuesMatrix=res$summaryValuesMatrix, trueFreeValuesMatrix=res$trueFreeValuesMatrix, phy=res$phy, traits=res$traits, vipthresh=0.8, abcTolerance=0.05)


}
\keyword{boxcoxPlsRejection}
\keyword{abc}


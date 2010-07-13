\name{doRun}
\alias{doRun}
\title{Start approximate Bayesian computation comparative methods run}
\description{
	Starts the run
}
\usage{
doRun(phy, traits, intrinsicFn, extrinsicFn, summaryFns=c(rawValuesSummaryStats, geigerUnivariateSummaryStats2), startingPriorsValues, startingPriorsFns, intrinsicPriorsValues, intrinsicPriorsFns, extrinsicPriorsValues, extrinsicPriorsFns, startingStatesGuess=c(), intrinsicValuesGuess=c(), extrinsicValuesGuess=c(), timeStep, toleranceVector=c(), numParticles=1000, standardDevFactor=0.05, nrepSim=100, plot=T, vipthresh=0.8, epsilonProportion=0.2, epsilonMultiplier=0.5, nStepsPRC=4) 
}
\arguments{
  \item{phy}{Tree (ape phylo format)}
  \item{traits}{Traits}
  \item{intrinsicFn}{name of (previously-defined) function that governs how traits evolve within a lineage, regardless of the states of other taxa}
  \item{extrinsicFn}{name of (previously-defined) function that governs how traits evolve within a lineage, based on the internal state and the states of other taxa}
  \item{summaryFns}{choice of which summary stats functions to use}
  \item{startingPriorsValues}{Matrix with ncol=number of states (characters) at root and nrow=2 (two parameters to pass to prior distribution)}
  \item{startionPriorsFns}{Vector containing names of prior distributions to use for root states: can be one of fixed, uniform, normal, lognormal, gamma (note that gamma with alpha=1 is the exponential distribution)}
  \item{intrinsicPriorsFns}{Matrix with ncol=number of parameters to pass to the intrinsic function and nrow=2 (two parameters to pass to prior distribution)}
  \item{intrinsicPriorsFns}{Vector containing names of prior distributions to use for intrinsic function parameters: can be one of fixed, uniform, normal, lognormal, gamma (note that gamma with alpha=1 is the exponential distribution)}
  \item{extrinsicPriorsFns}{Matrix with ncol=number of parameters to pass to the extrinsic function and nrow=2 (two parameters to pass to prior distribution)}
  \item{extrinsicPriorsFns}{Vector containing names of prior distributions to use for extrinsic function parameters: can be one of fixed, uniform, normal, lognormal, gamma (note that gamma with alpha=1 is the exponential distribution)}
  \item{startingStatesGuess}{Optional guess of starting value}
}
\details{
	This is awesome but not ready to use
}
\value{
	This has some value
	}
\references{O'Meara and Banbury, unpublished}
\author{Brian O'Meara}
\seealso{
}
\examples{
}
\keyword{arith}


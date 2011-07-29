\name{simulateData}
\alias{simulateData}


\title{Simulate data for initial TreEvo analysis}

\description{This function retreives a random value from prior distribution}

\usage{
simulateData(nrepSim, startingPriorsValues, intrinsicPriorsValues, extrinsicPriorsValues, startingPriorsFns, intrinsicPriorsFns, extrinsicPriorsFns, trueFreeValues, freevector, timeStep, intrinsicFn, extrinsicFn)
}

\arguments{
  \item{nrepSim}{number of initial simulations}
  \item{startingPriorsValues}{Matrix with ncol=number of states (characters) at root and nrow=2 (two parameters to pass to prior distribution)}
  \item{intrinsicPriorsValues}{Matrix with ncol=number of states (characters) at root and nrow=2 (two parameters to pass to prior distribution)}
  \item{extrinsicPriorsValues}{Matrix with ncol=number of states (characters) at root and nrow=2 (two parameters to pass to prior distribution)}
  \item{startingPriorsFns}{Vector containing names of prior distributions to use for root states: can be one of fixed, uniform, normal, lognormal, gamma, exponential}
  \item{intrinsicPriorsFns}{Vector containing names of prior distributions to use for root states: can be one of fixed, uniform, normal, lognormal, gamma, exponential}
  \item{extrinsicPriorsFns}{Vector containing names of prior distributions to use for root states: can be one of fixed, uniform, normal, lognormal, gamma, exponential}
  \item{trueFreeValues}{Initially will be an empty matrix, values are added as simulations run.  These are the parameter values pulled from the prior from each simulation.}
  \item{freevector}{A vector (length=number of parameters) of free (T) and fixed (F) parameters}
  \item{timeStep}{This value corresponds to the number of discrete time steps on the shortest branch}
  \item{intrinsicFn}{Name of intrinsic function characters should be simulated under (as used by doSimulation)}
  \item{extrinsicFn}{Name of extrinsic function characters should be simulated under (as used by doSimulation)}
}
\details{
	Used by TreEvo doRun function to calculate initial simulations
}
\value{
	Returns trueFreeValues and summary stats for simulations, both contained in an object called simSummaryStats
}
\references{O'Meara and Banbury, unpublished}
\author{Brian O'Meara and Barb Banbury}
\seealso{
}
\examples{

}
\keyword{}


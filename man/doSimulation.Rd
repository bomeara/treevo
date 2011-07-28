\name{doSimulation}

\title{Discrete Time Character Simulation}

\description{This function evolves continuous characters in a discrete time process.}

\usage{
doSimulation(splits, intrinsicFn, extrinsicFn, startingStates, intrinsicValues, extrinsicValues, timeStep, saveRealParams=FALSE, jobName="")
}

\arguments{
  \item{splits}{Output from the function getSimulationSplits; is a data frame of branching times, ancestor and descendant vectors}
  \item{intrinsicFn}{Name of intrinsic function characters should be simulated under}
  \item{extrinsicFn}{Name of extrinsic function characters should be simulated under}
  \item{startingStates}{State at the root}
  \item{intrinsicValues}{Vector of values corresponding to the params of the intrinsic model}
  \item{extrinsicValues}{Vector of values corresponding to the params of the extrinsic model}
  \item{timeStep}{This value corresponds to the number of discrete time steps on the shortest branch}
  \item{saveRealParams}{Saves intrinsicValues and extrinsicValues as both a matrix and a vector}
  \item{jobName}{Optional name for the job}
}
\details{
	SaveRealParams is useful for tracking the "real" true values if simulating data for abc runs.  It is not useful for emperical abc runs.  
}
\value{
	A data frame of species character (tip) values in the tree.
	}
\references{O'Meara and Banbury, unpublished}
\author{Brian O'Meara and Barb Banbury}
\seealso{
}
\examples{

}
\keyword{}


\name{parallelSimulation}
\alias{parallelSimulation}


\title{Simulate data for initial TreEvo analysis}

\description{This is a wrapper function for simulateData that allows for multithreading and checkpointing}

\usage{
parallelSimulation(nrepSim, coreLimit, splits, phy, startingPriorsValues, intrinsicPriorsValues, extrinsicPriorsValues, startingPriorsFns, intrinsicPriorsFns, extrinsicPriorsFns, freevector, timeStep, intrinsicFn, extrinsicFn, multicore,checkpointFile=NULL,checkpointFreq=24)
}

\arguments{
  \item{nrepSim}{Number of simulations}
  \item{coreLimit}{Number of cores to be used}
  \item{splits}{simulation splits (from getSimulationSplits}
  \item{phy}{Tree (Phylogenetic tree in phylo format)}
  \item{startingPriorsValues}{Matrix with ncol=number of states (characters) at root and nrow=2 (two parameters to pass to prior distribution)}
  \item{intrinsicPriorsValues}{Matrix with ncol=number of states (characters) at root and nrow=2 (two parameters to pass to prior distribution)}
  \item{extrinsicPriorsValues}{Matrix with ncol=number of states (characters) at root and nrow=2 (two parameters to pass to prior distribution)}
  \item{startingPriorsFns}{Vector containing names of prior distributions to use for root states: can be one of fixed, uniform, normal, lognormal, gamma, exponential}
  \item{intrinsicPriorsFns}{Vector containing names of prior distributions to use for root states: can be one of fixed, uniform, normal, lognormal, gamma, exponential}
  \item{extrinsicPriorsFns}{Vector containing names of prior distributions to use for root states: can be one of fixed, uniform, normal, lognormal, gamma, exponential}
  \item{freevector}{A vector (length=number of parameters) of free (T) and fixed (F) parameters}
  \item{timeStep}{This value corresponds to the number of discrete time steps}
  \item{intrinsicFn}{Name of intrinsic function characters should be simulated under (as used by doSimulation)}
  \item{extrinsicFn}{Name of extrinsic function characters should be simulated under (as used by doSimulation)}
  \item{multicore}{Whether to use multicore, default is FALSE}
  \item{checkpointFile}{Value for when to stop the analysis if NAs are present}
  \item{checkpointFile}{Optional file name for checkpointing simulations}
  \item{checkpointFreq}{Saving frequency for checkpointing}
}

\value{
	Returns matrix of trueFreeValues and summary statistics for simulations
}
\references{O'Meara and Banbury, unpublished}
\author{Brian O'Meara and Barb Banbury}


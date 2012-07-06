\name{getSimulationSplits}
\alias{getSimulationSplits}

\title{Simulation Splits}

\description{This function retrieves branch lengths and ancestor-decendant relationships from a tree}

\usage{
getSimulationSplits(phy)
}

\arguments{
  \item{phy}{Tree (ape phylo format)}
}
\details{
	This function is used by other TreEvo functions for internal calculations. 
}
\value{
A data frame of branching times, ancestor and descendant vectors
}
\references{O'Meara and Banbury, unpublished}
\author{Brian O'Meara and Barb Banbury}


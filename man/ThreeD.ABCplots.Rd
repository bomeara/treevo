\name{ThreeD.ABCplots}
\alias{ThreeD.ABCplots}

\title{Plot posterior density distribution for each generation in 3d}
\description{Plot posterior density distribution for each generation in 3d plot window}
\usage{
ThreeD.ABCplots(a, parameter, show.particles="none", plot.parent=FALSE, realParam=FALSE, realParamValues=NA) 
}
\arguments{
  \item{a}{$particleDataFrame output from doRun}
  \item{parameter}{column number of parameter of interest from $particleDataFrame}
  \item{show.particles}{option to show particles on 3d plot as "none" or as a function of "weights" or "distance"}
  \item{plot.parent}{option to plot lines on the floor of the 3d plot to show particle parantage}
  \item{realParam}{option to display real parameter value as a solid line, also must give actual value for this (realParamValues).  Note: this should only be done with simulated data where real param values are recorded}
  \item{realParamValues}{Value for (realParam)}
}
\details{
	This opens a new interactive 3d plotting window and plots the posterior density distribution of accepted particles from each generation.  Several options are available to add to the plot: plotting particles by weight or distance, plotting particle parantage, and plotting the real parameter values (if known).
}

\references{O'Meara and Banbury, unpublished}
\author{Barb Banbury}
\seealso{ThreeD.ABCplotsWPrior}
\examples{}
\keyword{}


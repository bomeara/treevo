\name{CredInt}
\alias{CredInt}
\title{Bayesian Credible Interval}
\description{
	This function calculates credible interval for each free parameter     
}
\usage{CredInt(particleDataFrame, percent=0.95)}
\arguments{
  \item{particleDataFrame}{particleDataFrame output from doRun}
  \item{percent}{Probability content of HPD}
  
}
\value{Returns a matrix with weighted mean, sd, upper and lower credible intervals for each free parameter}

\references{O'Meara and Banbury, unpublished}
\author{Brian O'Meara and Barb Banbury}

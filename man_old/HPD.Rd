\name{HPD}
\alias{HPD}
\title{Highest Posterior Density}
\description{
	This function calculates HPD for each free parameter     
}
\usage{HPD(particleDataFrame, percent=0.95, returnData=F)}
\arguments{
  \item{particleDataFrame}{particleDataFrame output from doRun}
  \item{percent}{Probability content of HPD}
  \item{returnData}{Option to return data that falls within HPD interval}
  
}
\value{Returns a matrix with weighted mean, sd, upper and lower HPD for each free parameter}

\references{O'Meara and Banbury, unpublished}
\author{Brian O'Meara and Barb Banbury}

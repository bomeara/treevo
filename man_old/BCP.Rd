\name{BCP}
\alias{BCP}
\title{Bayesian Coverage Probability}
\description{
	This function calculates coverage probability for a list of HPDs     
}
\usage{BCP(RealParam, HPD, verbose=F)}
\arguments{
  \item{RealParam}{Real parameter values}
  \item{HPD}{List of HPD from doRun_rej or doRun_prc results}
  \item{verbose}{Commented screen output}

}
\value{Returns a value for each free parameter that describes the percentage the real value falls within the HPD}
\details{Only for use with simulated data to test models}

\references{O'Meara and Banbury, unpublished}
\author{Brian O'Meara and Barb Banbury}

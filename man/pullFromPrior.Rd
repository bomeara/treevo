\name{pullFromPrior}
\alias{pullFromPrior}


\title{Pull value from prior}

\description{This function retreives random value from variable prior functions}

\usage{
pullFromPrior(priorValues, priorFn)
}

\arguments{
  \item{priorValues}{Matrix with ncol=number of states (characters) at root and nrow=2 (two parameters to pass to prior distribution)}
  \item{priorFn}{Prior shape}
}
\value{
	Returns mean
}
\references{O'Meara and Banbury, unpublished}
\author{Brian O'Meara and Barb Banbury}

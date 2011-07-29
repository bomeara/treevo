\name{pullFromPrior}
\alias{pullFromPrior}


\title{Get Value from Prior}

\description{This function retreives a random value from prior distribution}

\usage{
pullFromPrior(priorValues, priorFn)
}

\arguments{
  \item{priorValues}{Matrix with ncol=number of states (characters) at root and nrow=2 (two parameters to pass to prior distribution)}
  \item{priorFn}{Prior shape}
}
\details{
	Used by TreEvo functions
}
\value{
	Returns value from prior
}
\references{O'Meara and Banbury, unpublished}
\author{Brian O'Meara and Barb Banbury}
\seealso{
}
\examples{

}
\keyword{}


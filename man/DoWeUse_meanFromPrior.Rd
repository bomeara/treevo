\name{meanFromPrior}
\alias{meanFromPrior}


\title{Get Mean Value from Prior}

\description{This function retreives a mean value from variable prior functions}

\usage{
meanFromPrior(priorValues, priorFn)
}

\arguments{
  \item{priorValues}{Matrix with ncol=number of states (characters) at root and nrow=2 (two parameters to pass to prior distribution)}
  \item{priorFn}{Prior shape}
}
\details{
	Used by other TreEvo functions
}
\value{
	Returns mean
}
\references{O'Meara and Banbury, unpublished}
\author{Brian O'Meara and Barb Banbury}
\seealso{
}
\examples{

}
\keyword{}


\name{parentOffspringPlots}
\alias{parentOffspringPlots}


\title{Plotting parent-offspring particles}

\description{This function uses the particleDataFrame from TreEvo results and plots parent-offspring particles from generation to generation.}

\usage{
parentOffspringPlots(x)
}

\arguments{
  \item{x}{particleDataFrame from TreEvo results}
}
\details{
	Each parameter is plotted twice for parent-offspring relationships through the generations.  In the top row, particles are plotted as a measure of distance to the observed data; the farther away the particle the bigger the circle.  Gray circles indicate rejected particles.  In the bottom row, particles are plotted as a measure of their weights; larger circles are closer to observed data and therefore carry more weight in the analysis.  Grayed out stars represent rejected particles.
}
\value{
	Returns a layout of parent-offspring plots
}
\references{O'Meara and Banbury, unpublished}
\author{Brian O'Meara and Barb Banbury}
\seealso{
}
\examples{
	parentOffspringPlots(a$particleDataFrame)
}
\keyword{}


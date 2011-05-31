\name{boxcoxplsSummary.R}
\alias{boxcoxplsSummary.R}
\title{Box-Cox Transformations + PLS}
\description{Performs Box-Cox transformation and partial least squares (PLS) on summaryValues }
\usage{
boxcoxplsSummary.R(todo, summaryValues, prunedPlsResult, boxcoxLambda, boxcoxAddition)
}
\arguments{
  \item{todo}{vector the same length as summaryValues, describing importance (1=important, 0=not important)}
  \item{summaryValues}{summary statistics}
  \item{prunedPlsResult}{summary statistics found to be important via PLS}
  \item{boxcoxLambda}{Vector of best lambda for each summaryValue ($lambda from powerTransform used in Box-Cox transformation)}
  \item{boxcoxAddition}{low value addition}
}
\details{
	Given the Box-Cox transformation variable lambda, this function transforms the summary statistics and uses them to predict...
}
\value{
	boxcoxplsSummary returns predictResult...
	}
\references{Bates et al., 2009, Wegmann et al. (2009), Dejean et al., 2009, O'Meara and Banbury, unpublished}
\author{Brian O'Meara}
\seealso{
}
\examples{

}
\keyword{}


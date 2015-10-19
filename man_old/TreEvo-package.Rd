\name{TreEvo-package}
\alias{TreEvo-package}
\alias{treevo}
\docType{package}
\title{TreEvo--abc for comparative methods}
\description{
	A package for using approximate Bayesian computation for comparative methods in trait evolution 
}
\details{
\tabular{ll}{
Package: \tab TreEvo\cr
Type: \tab Package\cr
Version: \tab 0.3.3\cr
Date: \tab 2012-07-02\cr
License: \tab GPL\cr
}

}
\author{
Brian O'Meara, Barb L. Banbury 

Maintainer: Barb Banbury <darwinthesun@gmail.com>
}

\examples{
#Continuous character simulation under Brownian motion
library(ape)
phy<-rcoal(20)
char<-doSimulation(
	splits=getSimulationSplits(phy), 
	intrinsicFn=brownianIntrinsic, 
	extrinsicFn=nullExtrinsic, 
	startingValues=c(30), 
	intrinsicValues=c(.01), 
	extrinsicValues=c(0), 
	timeStep=0.001)
}

\keyword{treevo}
\keyword{abc}


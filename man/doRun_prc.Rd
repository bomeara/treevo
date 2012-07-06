\name{doRun_prc}
\alias{doRun_prc}
\title{Approximate Bayesian computation for comparative methods-Partial Rejection Control}
\description{
	Starts the abc-prc run
}
\usage{
doRun_prc(phy, traits, intrinsicFn, extrinsicFn, startingPriorsValues, startingPriorsFns, intrinsicPriorsValues, intrinsicPriorsFns, extrinsicPriorsValues, extrinsicPriorsFns, startingValuesGuess=c(), intrinsicValuesGuess=c(), extrinsicValuesGuess=c(), TreeYears=1e+04, numParticles=300, standardDevFactor=0.20, StartSims=300, plot=FALSE, epsilonProportion=0.7, epsilonMultiplier=0.7, nStepsPRC=5, jobName=NA, stopRule=FALSE, stopValue=0.05, vipthresh=0.8, multicore=FALSE, coreLimit=NA) 
}
\arguments{
  \item{phy}{Tree (Phylogenetic tree in phylo format)}
  \item{traits}{data matrix with rownames equal to phy}
  \item{intrinsicFn}{Name of (previously-defined) function that governs how traits evolve within a lineage, regardless of the states of other taxa}
  \item{extrinsicFn}{Name of (previously-defined) function that governs how traits evolve within a lineage, based on the internal state and the states of other taxa}
  \item{startingPriorsValues}{Matrix with ncol=number of states (characters) at root and nrow=2 (two parameters to pass to prior distribution)}
  \item{startingPriorsFns}{Vector containing names of prior distributions to use for root states: can be one of fixed, uniform, normal, lognormal, gamma, exponential}
  \item{intrinsicPriorsValues}{Matrix with ncol=number of parameters to pass to the intrinsic function and nrow=2 (two parameters to pass to prior distribution)}
  \item{intrinsicPriorsFns}{Vector containing names of prior distributions to use for intrinsic function parameters: can be one of fixed, uniform, normal, lognormal, gamma, exponential}
  \item{extrinsicPriorsValues}{Matrix with ncol=number of parameters to pass to the extrinsic function and nrow=2 (two parameters to pass to prior distribution)}
  \item{extrinsicPriorsFns}{Vector containing names of prior distributions to use for extrinsic function parameters: can be one of fixed, uniform, normal, lognormal, gamma, exponential}
  \item{startingValuesGuess}{Optional guess of starting values}
  \item{intrinsicValuesGuess}{Optional guess of intrinsic values}
  \item{extrinsicValuesGuess}{Optional guess of extrinsic values}
  \item{TreeYears}{Unit length of phy}
  \item{numParticles}{Number of accepted particles per generation}
  \item{standardDevFactor}{Standard deviation for mutating states}
  \item{StartSims}{Number of initial simulations}
  \item{plot}{If TRUE, plots distance of each simulation}
  \item{epsilonProportion}{Sets tolerance for initial simulations}
  \item{epsilonMultiplier}{Sets tolerance on subsequent generations}
  \item{nStepsPRC}{Number of generations}
  \item{jobName}{Optional job name}
  \item{stopRule}{If TRUE, an analysis will be terminated prior to set number of generations if the ratio of all parameters from one generation to the next falls below the stopValue}
  \item{stopValue}{Threshold value for terminating an analysis prior to nStpesPRC}
  \item{vipthresh}{Threshold value for pls regression of vip summary stats}
  \item{multicore}{If TRUE, initial simulations will be split among coreLimit nodes}
  \item{coreLimit}{Number of cores for initial simulations}

}

\details{
	This function performs an abc-prc analysis using an input phylogeny (phy), character data set (traits), models (intrinsicFn, extrinsicFn), and priors (startingPriorsValues, startingPriorsFns, intrinsicPriorsValues, intrinsicPriorsFns, extrinsicPriorsValues, extrinsicPriorsFns). Pulling from the priors, it simulates an initial set of simulations (StartSims). The set of simulations are boxcox transformed and a pls regression is performed for each free parameter to determine the most informative summary statistics (vipthresh). The euclidean distance of each initial simulation's vip summary stats to the input data is calculated, and tolerance is set based on epsilonProportion and epsilonMultiplier.
	
	The prc analysis begins with generation 1 and continues through nStepsPRC. Single simulations (or particles) are accepted if the distance of vip summary stats to the original data is less than the tolerance. A generation is complete when enough particles (numParticles) have been accepted. These particles make up the distribution for the next generation. The accepted particles from the final generation describe the posterior distributions of parameters.  
	
}

\value{
  \item{input.data}{Input variables: jobName, number of taxa, nrepSim, treeYears, epsilonProportion, epsilonMultiplier, nStepsPRC, numParticles, standardDevFactor}
  \item{PriorMatrix}{Matrix of prior distributions}
  \item{particleDataFrame}{DataFrame with information from each simulation, including generation, attempt, id, parentid, distance, weight, and parameter states}
  \item{toleranceVector}{Final tolerance vector}
  \item{phy}{Input phylogeny}
  \item{traits}{Input traits}
  \item{simTime}{Processor time for initial simulations}
  \item{time.per.gen}{Processor time for subsequent generations}
  \item{whichVip}{Matrix of vip summary statistics for each free parameter}
  \item{CredInt}{Credible Interval calculation for each free parameter of the final generation}
  \item{HPD}{Highest Posterior Density calculation each free parameter of the final generation}
  
}

\references{O'Meara and Banbury, unpublished; Sisson et al. 2007, Wegmann et al. 2009}
\author{Brian O'Meara and Barb Banbury}
\examples{
data(simData)

#doRun_prc(
#  phy = simPhy,
#  traits = simChar,
#  intrinsicFn=brownianIntrinsic, 
#  extrinsicFn=nullExtrinsic, 
#  startingPriorsFns="normal",
#  startingPriorsValues=matrix(c(mean(simChar[,1]), sd(simChar[,1]))),
#  intrinsicPriorsFns=c("exponential"),
#  intrinsicPriorsValues=matrix(c(10, 10), nrow=2, byrow=FALSE),
#  extrinsicPriorsFns=c("fixed"),
#  extrinsicPriorsValues=matrix(c(0, 0), nrow=2, byrow=FALSE),
#  TreeYears=1000,
#  standardDevFactor=0.2,
#  plot=FALSE,
#  StartSims=300,
#  epsilonProportion=0.7,
#  epsilonMultiplier=0.7,
#  nStepsPRC=5,
#  numParticles=100,
#  jobName="exampleRun",
#  stopRule=FALSE,
#  multicore=FALSE,
#  coreLimit=1
#)		

}
\keyword{doRun}
\keyword{doRun_prc}
\keyword{abc}


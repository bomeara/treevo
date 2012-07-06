\name{doRun_rej}
\alias{doRun_rej}
\title{Approximate Bayesian computation for comparative methods-Rejection}
\description{
	Starts the abc-rejection run
}
\usage{
doRun_rej(phy, traits, intrinsicFn, extrinsicFn, startingPriorsValues, startingPriorsFns, intrinsicPriorsValues, intrinsicPriorsFns, extrinsicPriorsValues, extrinsicPriorsFns, startingValuesGuess=c(), intrinsicStatesGuess=c(), extrinsicStatesGuess=c(), TreeYears=1e+04, standardDevFactor=0.20, StartSims=NA, jobName=NA, vipthresh=0.8, abcTolerance=0.1, multicore=FALSE, coreLimit=NA) 
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
  \item{intrinsicStatesGuess}{Optional guess of intrinsic values}
  \item{extrinsicStatesGuess}{Optional guess of extrinsic values}
  \item{TreeYears}{Unit length of phy}
  \item{standardDevFactor}{Standard deviation for mutating states}
  \item{StartSims}{Number of simulations}
  \item{jobName}{Optional job name}
  \item{vipthresh}{Threshold value for pls regression of vip summary stats}
  \item{abcTolerance}{Proportion of accepted simulations}
  \item{multicore}{If TRUE, initial simulations will be split among coreLimit nodes}
  \item{coreLimit}{Number of cores for initial simulations}

}

\details{
	This function performs an abc-rejection analysis using an input phylogeny (phy), character data set (traits), models (intrinsicFn, extrinsicFn), and priors (startingPriorsValues, startingPriorsFns, intrinsicPriorsValues, intrinsicPriorsFns, extrinsicPriorsValues, extrinsicPriorsFns). Pulling from the priors, it performs StartSims simulations. This set of simulations are boxcox transformed and a pls regression is performed for each free parameter to determine the most informative summary statistics (vipthresh). The euclidean distance of each initial simulation's vip summary stats to the input data is calculated, and those that fall under the abcTolerance are kept as accepted particles. These describe the posterior distributions of parameters.  
	
}

\value{
  \item{input.data}{Input variables: jobName, number of taxa, nrepSim, treeYears, epsilonProportion, epsilonMultiplier, nStepsPRC, numParticles, standardDevFactor}
  \item{PriorMatrix}{Matrix of prior distributions}
  \item{phy}{Input phylogeny}
  \item{traits}{Input traits}
  \item{trueFreeValuesANDSummaryValues}{Parameter estimates and summary stats from all sims}
  \item{simTime}{Processor time for simulations}
  \item{whichVip}{Matrix of vip summary statistics for each free parameter}
  \item{abcDistancesRaw}{Euclidean distances for each simulation and free parameter}
  \item{particleDataFrame}{DataFrame with information from each simulation, including generation, attempt, id, parentid, distance, weight, and parameter states}
  \item{CredInt}{Credible Interval calculation for each free parameter of the final generation}
  \item{HPD}{Highest Posterior Density calculation each free parameter of the final generation}
  
}

\references{O'Meara and Banbury, unpublished; Sisson et al. 2007, Wegmann et al. 2009}
\author{Brian O'Meara and Barb Banbury}
\examples{
data(simData)

#c<-doRun_rej( #make sure priors are uniform with this one! 
#	phy=simPhy, 
#	traits=simChar,
#	intrinsicFn=brownianIntrinsic, 
#	extrinsicFn=nullExtrinsic, 
#	startingPriorsFns="normal",
#	startingPriorsValues=matrix(c(mean(char[,1]), sd(char[,1]))),
#	intrinsicPriorsFns=c("exponential"),
#	intrinsicPriorsValues=matrix(c(10, 10), nrow=2, byrow=FALSE), #grep for normal in pkg
#	extrinsicPriorsFns=c("fixed"),
#	extrinsicPriorsValues=matrix(c(0, 0), nrow=2, byrow=FALSE),
#	StartSims=1000, 
#	jobName="run_c", 
#	abcTolerance=0.05, 
#	multicore=F, 
#	coreLimit=1 
#)			
	

}
\keyword{doRun}
\keyword{doRun_rej}
\keyword{abc}


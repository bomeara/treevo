#' Approximate Bayesian computation for comparative methods-Rejection
#'
#' Starts the abc-rejection run
#'
#' This function performs an abc-rejection analysis using an input phylogeny
#' (phy), character data set (traits), models (intrinsicFn, extrinsicFn), and
#' priors (startingPriorsValues, startingPriorsFns, intrinsicPriorsValues,
#' intrinsicPriorsFns, extrinsicPriorsValues, extrinsicPriorsFns). Pulling from
#' the priors, it performs StartSims simulations. This set of simulations are
#' boxcox transformed and a pls regression is performed for each free parameter
#' to determine the most informative summary statistics (vipthresh). The
#' euclidean distance of each initial simulation's vip summary stats to the
#' input data is calculated, and those that fall under the abcTolerance are
#' kept as accepted particles. These describe the posterior distributions of
#' parameters.
#'

#' @param phy Tree (Phylogenetic tree in phylo format)

#' @param traits data matrix with rownames equal to phy

#' @param intrinsicFn Name of (previously-defined) function that governs how
#' traits evolve within a lineage, regardless of the states of other taxa

#' @param extrinsicFn Name of (previously-defined) function that governs how
#' traits evolve within a lineage, based on the internal state and the states
#' of other taxa

#' @param startingPriorsValues Matrix with ncol=number of states (characters)
#' at root and nrow=2 (two parameters to pass to prior distribution)

#' @param startingPriorsFns Vector containing names of prior distributions to
#' use for root states: can be one of fixed, uniform, normal, lognormal, gamma,
#' exponential

#' @param intrinsicPriorsValues Matrix with ncol=number of parameters to pass
#' to the intrinsic function and nrow=2 (two parameters to pass to prior
#' distribution)

#' @param intrinsicPriorsFns Vector containing names of prior distributions to
#' use for intrinsic function parameters: can be one of fixed, uniform, normal,
#' lognormal, gamma, exponential

#' @param extrinsicPriorsValues Matrix with ncol=number of parameters to pass
#' to the extrinsic function and nrow=2 (two parameters to pass to prior
#' distribution)

#' @param extrinsicPriorsFns Vector containing names of prior distributions to
#' use for extrinsic function parameters: can be one of fixed, uniform, normal,
#' lognormal, gamma, exponential

#' @param startingValuesGuess Optional guess of starting values

#' @param intrinsicStatesGuess Optional guess of intrinsic values

#' @param extrinsicStatesGuess Optional guess of extrinsic values

#' @param TreeYears Unit length of phy

#' @param standardDevFactor Standard deviation for mutating states

#' @param StartSims Number of simulations

#' @param jobName Optional job name

#' @param abcTolerance Proportion of accepted simulations

#' @param multicore If TRUE, initial simulations will be split among coreLimit
#' nodes

#' @param coreLimit Number of cores for initial simulations

#' @param checkpointFile Optional file name for checkpointing simulations

#' @param checkpointFreq Saving frequency for checkpointing

#' @param validation Cross Validation procedure for abc

#' @param scale scale for pls.model.list

#' @param variance.cutoff variance cutoff for pls.model.list

#' @param savesims option to save individual simulations

#' @param niter.goal Adjust number of starting points for Geiger to return the best parameter estimates this number of times on average

#' @param generation.time The number of years per generation. This sets the coarseness of the simulation; if it's set to 1000, for example, the population moves every 1000 years.

#' @return \item{input.data}{Input variables: jobName, number of taxa, nrepSim,
#' treeYears, epsilonProportion, epsilonMultiplier, nStepsPRC, numParticles,
#' standardDevFactor} \item{PriorMatrix}{Matrix of prior distributions}
#' \item{phy}{Input phylogeny} \item{traits}{Input traits}
#' \item{trueFreeValuesANDSummaryValues}{Parameter estimates and summary stats
#' from all sims} \item{simTime}{Processor time for simulations}
#' \item{whichVip}{Matrix of vip summary statistics for each free parameter}
#' \item{abcDistancesRaw}{Euclidean distances for each simulation and free
#' parameter} \item{particleDataFrame}{DataFrame with information from each
#' simulation, including generation, attempt, id, parentid, distance, weight,
#' and parameter states} \item{CredInt}{Credible Interval calculation for each
#' free parameter of the final generation} \item{HPD}{Highest Posterior Density
#' calculation each free parameter of the final generation}

#' @author Brian O'Meara and Barb Banbury

# @references O'Meara and Banbury, unpublished; Sisson et al. 2007, Wegmann et
#' al. 2009

# @keywords doRun doRun_rej abc

#' @examples
#'
#' data(simRun)
#'
#' c<-doRun_rej( #make sure priors are uniform with this one!
#' 	phy=simPhy,
#' 	traits=simChar,
#' 	intrinsicFn=brownianIntrinsic,
#' 	extrinsicFn=nullExtrinsic,
#' 	startingPriorsFns="normal",
#' 	startingPriorsValues=matrix(c(mean(char[,1]), sd(char[,1]))),
#' 	intrinsicPriorsFns=c("exponential"),
#' 	intrinsicPriorsValues=matrix(c(10, 10), nrow=2, byrow=FALSE), #grep for normal in pkg
#' 	extrinsicPriorsFns=c("fixed"),
#' 	extrinsicPriorsValues=matrix(c(0, 0), nrow=2, byrow=FALSE),
#' 	StartSims=1000,
#' 	jobName="run_c",
#' 	abcTolerance=0.05,
#' 	multicore=F,
#' 	coreLimit=1
#' )
#'
#'
#'


##TreeYears = 1000 if tree is in thousands of years


#' @name doRun_rej
#' @rdname doRun_rej
#' @export
doRun_rej<-function(phy, traits, intrinsicFn, extrinsicFn, startingPriorsValues, startingPriorsFns, 
	intrinsicPriorsValues, intrinsicPriorsFns, extrinsicPriorsValues, extrinsicPriorsFns, 
	startingValuesGuess=c(), intrinsicStatesGuess=c(), extrinsicStatesGuess=c(), TreeYears=1e+04, 
	standardDevFactor=0.20, StartSims=NA, jobName=NA, abcTolerance=0.1, multicore=FALSE, coreLimit=NA, 
	checkpointFile=NULL, checkpointFreq=24, validation="CV", scale=TRUE, variance.cutoff=95, 
	savesims=FALSE, niter.goal=5, generation.time=1) {
	
	#library(geiger)
	#library(pls)
	if (!is.binary.tree(phy)) {
		print("Warning: Tree is not fully dichotomous")
	}
	startTime<-proc.time()[[3]]
	timeStep<-generation.time/TreeYears
	#splits<-getSimulationSplits(phy) #initialize this info
	taxon.df <- getTaxonDFWithPossibleExtinction(phy)


	#figure out number of free params
	numberParametersTotal<-dim(startingPriorsValues)[2] +  dim(intrinsicPriorsValues)[2] + dim(extrinsicPriorsValues)[2]
	numberParametersFree<-numberParametersTotal
	numberParametersStarting<-0
	numberParametersIntrinsic<-0
	numberParametersExtrinsic<-0
	freevariables<-matrix(data=NA, nrow=2, ncol=0)
	titlevector<-c()
	freevector<-c()

	namesForPriorMatrix<-c()
	PriorMatrix<-matrix(c(startingPriorsFns, intrinsicPriorsFns, extrinsicPriorsFns), nrow=1, ncol=numberParametersTotal)
	for (a in 1:dim(startingPriorsValues)[2]) {
		namesForPriorMatrix<-c(paste("StartingStates", a, sep=""))
	}
	for (b in 1:dim(intrinsicPriorsValues)[2]) {
		namesForPriorMatrix<-append(namesForPriorMatrix, paste("IntrinsicValue", b, sep=""))
	}
	for (c in 1:dim(extrinsicPriorsValues)[2]) {
		namesForPriorMatrix <-append(namesForPriorMatrix, paste("ExtrinsicValue", c, sep=""))
	}
	PriorMatrix<-rbind(PriorMatrix, cbind(startingPriorsValues, intrinsicPriorsValues, extrinsicPriorsValues))
	colnames(PriorMatrix)<-namesForPriorMatrix
	rownames(PriorMatrix)<-c("shape", "value1", "value2")

	for (i in 1:dim(startingPriorsValues)[2]) {
		priorFn<-match.arg(arg=startingPriorsFns[i],choices=c("fixed", "uniform", "normal", "lognormal", "gamma", "exponential"),several.ok=FALSE)
		if (priorFn=="fixed") {
			numberParametersFree<-numberParametersFree-1
			freevector<-c(freevector, FALSE)
		}
		else {
			numberParametersStarting<-numberParametersStarting+1
			freevariables<-cbind(freevariables, startingPriorsValues[, i])
			titlevector <-c(titlevector, paste("Starting", numberParametersStarting))
			freevector<-c(freevector, TRUE)
		}
	}
	for (i in 1:dim(intrinsicPriorsValues)[2]) {
		priorFn<-match.arg(arg=intrinsicPriorsFns[i],choices=c("fixed", "uniform", "normal", "lognormal", "gamma", "exponential"),several.ok=FALSE)
		if (priorFn=="fixed") {
			numberParametersFree<-numberParametersFree-1
			freevector<-c(freevector, FALSE)
		}
		else {
			numberParametersIntrinsic<-numberParametersIntrinsic+1
			freevariables<-cbind(freevariables, intrinsicPriorsValues[, i])
			titlevector <-c(titlevector, paste("Intrinsic", numberParametersIntrinsic))
			freevector<-c(freevector, TRUE)
		}
	}

	for (i in 1:dim(extrinsicPriorsValues)[2]) {
		priorFn<-match.arg(arg=extrinsicPriorsFns[i],choices=c("fixed", "uniform", "normal", "lognormal", "gamma", "exponential"),several.ok=FALSE)
		if (priorFn=="fixed") {
			numberParametersFree<-numberParametersFree-1
			freevector<-c(freevector, FALSE)
		}
		else {
			numberParametersExtrinsic<-numberParametersExtrinsic+1
			freevariables<-cbind(freevariables, extrinsicPriorsValues[, i])
			titlevector <-c(titlevector, paste("Extrinsic", numberParametersExtrinsic))
			freevector<-c(freevector, TRUE)
		}
	}

	#initialize guesses, if needed
	if (length(startingValuesGuess)==0) { #if no user guesses, try pulling a value from the prior
		startingValuesGuess<-rep(NA,length(startingPriorsFns))
		for (i in 1:length(startingPriorsFns)) {
			startingValuesGuess[i]<-pullFromPrior(startingPriorsValues[,i],startingPriorsFns[i])
		}
	}
	if (length(intrinsicStatesGuess)==0) { #if no user guesses, try pulling a value from the prior
		intrinsicStatesGuess<-rep(NA,length(intrinsicPriorsFns))
		for (i in 1:length(intrinsicPriorsFns)) {
			intrinsicStatesGuess[i]<-pullFromPrior(intrinsicPriorsValues[,i],intrinsicPriorsFns[i])
		}
	}
	if (length(extrinsicStatesGuess)==0) { #if no user guesses, try pulling a value from the prior
		extrinsicStatesGuess<-rep(NA,length(extrinsicPriorsFns))
		for (i in 1:length(extrinsicPriorsFns)) {
			extrinsicStatesGuess[i]<-pullFromPrior(extrinsicPriorsValues[,i],extrinsicPriorsFns[i])
		}
	}

	if (is.na(StartSims)) {
		StartSims<-1000*numberParametersFree
	}

	nrepSim<-StartSims #Used to be multiple tries where nrepSim = StartSims*((2^try)/2).  If initial simulations are not enough, and we need to try again then new analysis will double number of initial simulations
	cat(paste("Number of simulations set to", nrepSim, "\n"))
	if(!is.null(checkpointFile)) {
		save(list=ls(),file=paste(checkpointFile,".intialsettings.Rsave",sep=""))
	}



	#Figure out how many iterations to use for optimization in Geiger.
	brown<-fitContinuous(phy=phy, dat=traits, model="BM", ncores=1, control=list(niter=100)) #it actually runs faster without checking for cores. And we parallelize elsewhere
	lambda<-fitContinuous(phy=phy, dat=traits, model="lambda", ncores=1, control=list(niter=100))
	delta<-fitContinuous(phy=phy, dat=traits, model="delta", ncores=1, control=list(niter=100))
	ou<-fitContinuous(phy=phy, dat=traits, model="OU", ncores=1, control=list(niter=100))
	white<-fitContinuous(phy=phy, dat=traits, model="white", ncores=1, control=list(niter=100))

	cat("Setting number of starting points for Geiger optimization to")
	niter.brown.g <- round(max(10, min(niter.goal/solnfreq(brown),100)))
	cat(paste("\n",niter.brown.g, "for Brownian motion"))
	niter.lambda.g <- round(max(10, min(niter.goal/solnfreq(lambda),100)))
	cat(paste("\n",niter.lambda.g, "for lambda"))
	niter.delta.g <- round(max(10, min(niter.goal/solnfreq(delta),100)))
	cat(paste("\n",niter.delta.g, "for delta"))
	niter.OU.g <- round(max(10, min(niter.goal/solnfreq(ou),100)))
	cat(paste("\n",niter.OU.g, "for OU"))
	niter.white.g <- round(max(10, min(niter.goal/solnfreq(white),100)))
	cat(paste("\n",niter.white.g, "for white noise"))

	trueFreeValuesANDSummaryValues<-parallelSimulation(nrepSim, coreLimit, taxon.df, phy, startingPriorsValues, intrinsicPriorsValues, extrinsicPriorsValues, startingPriorsFns, intrinsicPriorsFns, extrinsicPriorsFns, freevector, timeStep, intrinsicFn, extrinsicFn, multicore, checkpointFile, checkpointFreq, niter.brown=niter.brown.g, niter.lambda=niter.lambda.g, niter.delta=niter.delta.g, niter.OU=niter.OU.g, niter.white=niter.white.g)

	cat("\n\n")
	simTime<-proc.time()[[3]]-startTime
	cat(paste("Simulations took", round(simTime, digits=3), "seconds"), "\n")

	#separate the simulation results: true values and the summary values
	trueFreeValuesMatrix<-trueFreeValuesANDSummaryValues[,1:numberParametersFree]
	summaryValuesMatrix<-trueFreeValuesANDSummaryValues[,-1:-numberParametersFree]
	if (savesims){
		save(trueFreeValuesMatrix, summaryValuesMatrix, simTime, file=paste("sims", jobName, ".Rdata", sep=""))
	}
	res<-PLSRejection(summaryValuesMatrix, trueFreeValuesMatrix, phy, traits, abcTolerance)
	#save(abcDistancesRaw, abcDistancesRawTotal, abcDistances, abcResults, particleDataFrame, file="")
	input.data<-rbind(jobName, length(phy[[3]]), timeStep, StartSims, standardDevFactor, abcTolerance)
  print(res)
	rejectionResults<-vector("list")
	#names(rejectionResults)<-c("input.data", "PriorMatrix", "phy", "traits")
#save(trueFreeValuesMatrix,res, file="BarbsTestofDistanceCalc2.Rdata")
	rejectionResults$input.data<-input.data
	rejectionResults$PriorMatrix<-PriorMatrix
	rejectionResults$phy<-phy
	rejectionResults$traits<-traits
	rejectionResults$trueFreeValuesANDSummaryValues<-trueFreeValuesANDSummaryValues
  rejectionResults$SimTime<-simTime
  rejectionResults$abcDistances<-res$abcDistances
	rejectionResults$particleDataFrame<-res$particleDataFrame
  rejectionResults$CredInt<-CredInt(res$particleDataFrame)
	rejectionResults$HPD<-highestProbDens(res$particleDataFrame)



	return(rejectionResults)
}

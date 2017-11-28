#' Run Approximate Bayesian Computation for Phylogenetic Comparative Methods
#' 
#' The \code{doRun} functions are the main interface for \code{TreEvo} users to do Approximate Bayesian Computation (ABC) analysis,
#' effectively wrapping the \code{\link{simulateWithPriors}} functions to perform simulations, which
#' themselves are wrappers for the the \code{\link{doSimulation}} functions. The two current \code{doRun}
#' functions are \code{doRun_prc} which applies a partial rejection control(PRC) ABC analysis over multiple generations of simulations,
#' and \code{doRun_rej} which performs a full rejection ('rej') ABC analysis.
#' 
#' Both \code{doRun} functions take an input phylogeny (\code{phy}),
#' observed trait data set (\code{traits}), models (\code{intrinsicFn}, \code{extrinsicFn}), and priors
#' (\code{startingPriorsValues}, \code{startingPriorsFns}, \code{intrinsicPriorsValues},
#' \code{intrinsicPriorsFns}, \code{extrinsicPriorsValues}, \code{extrinsicPriorsFns}). Pulling from
#' the priors, it simulates an initial set of simulations (\code{StartSims}). This set of simulations is
#' boxcox transformed ,and a PLS regression (see \code{\link{PLSmethods}}) is performed for each free parameter
#' to determine the most informative summary statistics (using \code{variance.cutoff}). The euclidean distance is calculated
#' between each each initial simulation's most informative summary statistics and the input observed data.
#' 
#' Following that step, the approach of the two functions diverge drastically.
#'
#' For function \code{doRun_rej}, those simulations whose most informative summary statistics fall below
#' \code{abcTolerance} are kept as accepted 'particles' (simulations runs), describing the
#' posterior distribution of parameters. No additional generations of simulations are performed.
#'
#' Coversely, function \code{doRun_prc} performs an ABC-PRC analysis,
#' which follows a much more complicated algorithm with additional generations. In PRC, tolerance is set based on
#' \code{epsilonProportion} and \code{epsilonMultiplier}, and the analysis
#' begins with generation 1 and continues through a number of generations equal to \code{nStepsPRC}.
#' Single simulation runs ('particles') are accepted if the distance of the
#' most informative summary stats to the original data is less than this tolerance. A generation
#' is complete when enough particles (\code{numParticles}) have been accepted. These
#' particles make up the distribution for the next generation. The accepted
#' particles from the final generation describe the posterior distributions of
#' parameters.
#'


#' @inheritParams doSimulation
#' @inheritParams simulateWithPriors
#' @inheritParams PLSmethods

#' @param traits Data matrix with rownames identical to \code{phy@tip.label}.

# @param startingValuesGuess Optional guess of starting values.

# @param intrinsicValuesGuess Optional guess of intrinsic values.

# @param extrinsicValuesGuess Optional guess of extrinsic values.

# TreeYears = 1000000 if tree length is 1 million of years, 1000 if 1 thousand years, etc.

#' @param numParticles Number of accepted particles per generation.

#' @param standardDevFactor Standard deviation for mutating states.

#' @param StartSims Number of initial simulations.

#' @param plot If \code{TRUE}, plots distance of each simulation.

#' @param epsilonProportion Sets tolerance for initial simulations.

#' @param epsilonMultiplier Sets tolerance on subsequent generations.

#' @param nStepsPRC Number of generations to run.

#' @param jobName Optional job name.

#' @param stopRule If \code{TRUE}, an analysis will be terminated prior to set number
#' of generations if the ratio of all parameters from one generation to the
#' next falls below the \code{stopValue}.

#' @param stopValue Threshold value for terminating an analysis prior to
#' \code{nStepesPRC}.

#' @param maxAttempts Maximum number attempts made in \code{while} loop within the PRC algorithm. 
#' If reached, the algorithm will terminate with an error message.
#' By default, this is infinite, and thus there is no effective limit without user intervention.

#' @param abcTolerance Proportion of accepted simulations.

#' @param checkpointFile Optional file name for checkpointing simulations.

#' @param checkpointFreq Saving frequency for checkpointing.

#' @param savesims Option to save individual simulations, output to a .Rdata file.

#' @param saveData Option to save various run information during the analysis, including summary statistics from analyses, output to external .Rdata and .txt files.

#' @param niter.goal Adjust number of starting points for package \code{geiger} to return the best parameter estimates this number of times on average.


#' @param verboseParticles If \code{TRUE} (the default), a large amount of information about parameter estimates
#' and acceptance of particles is output to console via \code{message} as \code{doRun_prc} runs. 

#' @param diagnosticPRCmode If \code{TRUE} (\emph{not} the default), the function will be very noisy about characteristics of 
#' the PRC algorithm as it runs.

#' @return 
#' The output of these two functions are lists, composed of multiple objects,
#' which differ slightly in their content among the two functions. For \code{doRun_prc}, the output is:

#' \describe{
#' \item{input.data}{Input variables: jobName, number of taxa, nrepSim, 
#' generation.time, TreeYears, timeStep, totalGenerations, epsilonProportion, epsilonMultiplier, nStepsPRC, numParticles, standardDevFactor} 

#' \item{PriorMatrix}{Matrix of prior distributions}

#' \item{particleDataFrame}{DataFrame with information from each simulation,
#' including generation, attempt, id, parentid, distance, weight, and parameter states}

#' \item{toleranceVector}{Final tolerance vector} 

#' \item{phy}{Input phylogeny} 

#' \item{traits}{Input traits} \item{simTime}{Processor time for initial simulations} 

#' \item{time.per.gen}{Processor time for subsequent generations}

# \item{whichVip}{Matrix of vip summary statistics for each free parameter} 

#' \item{credibleInt}{Credible Interval calculation for each free parameter of the final generation} 

#' \item{HPD}{Highest Posterior Density calculation each free parameter of the final generation}
#' }

#' For \code{doRun_rej}, the output is:

#' \describe{
#' \item{input.data}{Input variables: jobName, number of taxa, nrepSim, 
#' generation.time, TreeYears, timeStep, totalGenerations, epsilonProportion, epsilonMultiplier, nStepsPRC, numParticles, standardDevFactor} 

#' \item{PriorMatrix}{Matrix of prior distributions}

#' \item{phy}{Input phylogeny} 

#' \item{traits}{Input traits}

#' \item{trueFreeValuesANDSummaryValues}{Parameter estimates and summary stats from all sims} 

#' \item{simTime}{Processor time for simulations}

# \item{whichVip}{Matrix of vip summary statistics for each free parameter}

#' \item{abcDistancesRaw}{Euclidean distances for each simulation and free parameter} 

#' \item{particleDataFrame}{DataFrame with information from each 
#' simulation, including generation, attempt, id, parentid, distance, weight, and parameter states} 

#' \item{credibleInt}{Credible Interval calculation for each free parameter of the final generation} 

#' \item{HPD}{Highest Posterior Density calculation each free parameter of the final generation}
#' }
#' 
#'

#' @author Brian O'Meara and Barb Banbury

#' @references Sisson et al. 2007, Wegmann et
#' al. 2009

# @references O'Meara and Banbury, unpublished; Sisson et al. 2007, Wegmann et
# al. 2009

#' @examples
#' \donttest{
#' data(simRunExample)
#'
#' # NOTE: the example analyses below sample too few particles, 
#' 	# over too few steps, with too few starting simulations
#' 	# - all for the sake of examples that reasonably test the functions
#' 
#' # Please set these values to more realistic levels for your analyses!
#' 
#' resultsPRC<-doRun_prc(
#'   phy = simPhy,
#'   traits = simChar,
#'   intrinsicFn=brownianIntrinsic,
#'   extrinsicFn=nullExtrinsic,
#'   startingPriorsFns="normal",
#'   startingPriorsValues=matrix(c(mean(simChar[,1]), sd(simChar[,1]))),
#'   intrinsicPriorsFns=c("exponential"),
#'   intrinsicPriorsValues=matrix(c(10, 10), nrow=2, byrow=FALSE),
#'   extrinsicPriorsFns=c("fixed"),
#'   extrinsicPriorsValues=matrix(c(0, 0), nrow=2, byrow=FALSE),
#'   generation.time=1000,
#'   standardDevFactor=0.2,
#'   plot=FALSE,
#'   StartSims=10,
#'   epsilonProportion=0.7,
#'   epsilonMultiplier=0.7,
#'   nStepsPRC=3,
#'   numParticles=20,
#'   jobName="examplerun_prc",
#'   stopRule=FALSE,
#'   multicore=FALSE,
#'   coreLimit=1
#' )
#' 
#' resultsPRC
#'
#' #one should make sure priors are uniform with doRun_rej!
#'
#' resultsRej<-doRun_rej( 
#' 	phy=simPhy,
#' 	traits=simChar,
#' 	intrinsicFn=brownianIntrinsic,
#' 	extrinsicFn=nullExtrinsic,
#' 	startingPriorsFns="normal",
#' 	startingPriorsValues=matrix(c(mean(simChar[,1]), sd(simChar[,1]))),
#' 	intrinsicPriorsFns=c("exponential"),
#' 	intrinsicPriorsValues=matrix(c(10, 10), nrow=2, byrow=FALSE), #grep for normal in pkg
#' 	extrinsicPriorsFns=c("fixed"),
#' 	extrinsicPriorsValues=matrix(c(0, 0), nrow=2, byrow=FALSE),
#' 	StartSims=10,
#' 	jobName="examplerun_rej",
#' 	abcTolerance=0.05,
#' 	multicore=FALSE,
#' 	coreLimit=1
#' 	)
#' 
#' resultsRej
#' }


##This seems to be working if partialResults does not exist.  If checkpoint=TRUE, then run fails.

#the doRun_prc function takes input from the user and then automatically guesses optimal parameters, though user overriding is also possible.
#the guesses are used to do simulations near the expected region. If omitted, they are set to the midpoint of the input parameter matrices



#' @name doRun
#' @rdname doRun
#' @export
doRun_prc<-function(
	phy, traits, intrinsicFn, extrinsicFn, startingPriorsValues, startingPriorsFns, 
	intrinsicPriorsValues, intrinsicPriorsFns, extrinsicPriorsValues, extrinsicPriorsFns, 
	#startingValuesGuess=c(), intrinsicValuesGuess=c(), extrinsicValuesGuess=c(), 
	generation.time=1000, TreeYears=max(branching.times(phy)) * 1e6, 
	multicore=FALSE, coreLimit=NA, validation="CV", scale=TRUE, variance.cutoff=95,
	niter.goal=5, 
	numParticles=300, standardDevFactor=0.20, 
	StartSims=300, epsilonProportion=0.7, epsilonMultiplier=0.7, nStepsPRC=5, 
	stopRule=FALSE, stopValue=0.05, maxAttempts=Inf, diagnosticPRCmode=FALSE,
	jobName=NA, saveData=FALSE, verboseParticles=TRUE, plot=FALSE) {	

	functionStartTime<-proc.time()[[3]]
	
	if (!is.binary.tree(phy)) {
		warning("Tree is not fully dichotomous, this may cause issues!")
	}
	if(!is.numeric(maxAttempts)){
		stop("maxAttempts must be numeric")}

	timeStep<-generation.time/TreeYears
	message(paste0("The effective timeStep for this tree will be ",signif(timeStep),
		", as a proportion of tree height (root to furthest tip)..."))

	edgesRescaled<-phy$edge.length/max(node.depth.edgelength(phy))
	message("Rescaling edge lengths relative to maximum tip-to-root distance...")
	
	if(max(edgesRescaled) < timeStep) {
		stop("Tree has *NO* rescaled branches longer than generation.time/TreeYears, no simulated evol change can occur!")
	}
	if(min(edgesRescaled) < timeStep) {
		warning("Tree has rescaled branches shorter than generation.time/TreeYears; no evol change can be assigned to these, and geiger functions may fail!")
		#print("Tree has zero or nearly zero length branches")
	}
	
	totalGenerations<-sum(sapply(edgesRescaled,function(x) floor(x/timeStep)))
	message("Given generation time, a total of ",round(totalGenerations)," generations will occur over this tree")
	
	#splits<-getSimulationSplits(phy) #initialize this info
	taxon.df <- getTaxonDFWithPossibleExtinction(phy)
	
	# get freevector
	freevector<-getFreeVector(startingPriorsFns=startingPriorsFns, startingPriorsValues=startingPriorsValues, 
						intrinsicPriorsFns=intrinsicPriorsFns, intrinsicPriorsValues=intrinsicPriorsValues,
						extrinsicPriorsFns=extrinsicPriorsFns, extrinsicPriorsValues=extrinsicPriorsValues)
	numberParametersTotal<-length(freevector)
	numberParametersFree<-sum(freevector)

	#create PriorMatrix
	namesForPriorMatrix<-c()
	PriorMatrix<-matrix(c(startingPriorsFns, intrinsicPriorsFns, extrinsicPriorsFns), nrow=1, ncol=numberParametersTotal)
	for (a in 1:dim(startingPriorsValues)[2]) {
		namesForPriorMatrix<-c(paste0("StartingStates ", a, sep=""))
	}
	for (b in 1:dim(intrinsicPriorsValues)[2]) {
		namesForPriorMatrix<-append(namesForPriorMatrix, paste0("IntrinsicValue ", b, sep=""))
	}
	#print(extrinsicPriorsValues)
	for (c in 1:dim(extrinsicPriorsValues)[2]) {
		namesForPriorMatrix <-append(namesForPriorMatrix, paste0("ExtrinsicValue ", c, sep=""))
	}
	PriorMatrix<-rbind(PriorMatrix, cbind(startingPriorsValues, intrinsicPriorsValues, extrinsicPriorsValues))
	colnames(PriorMatrix)<-namesForPriorMatrix
	rownames(PriorMatrix)<-c("shape", "value1", "value2")

		
	#initialize weighted mean sd matrices
	weightedMeanParam<-matrix(nrow=nStepsPRC, ncol=numberParametersTotal)
	colnames(weightedMeanParam)<-namesForPriorMatrix
	rownames(weightedMeanParam)<-paste0("Gen ", c(1: nStepsPRC), sep="")
	param.stdev<-matrix(nrow=nStepsPRC, ncol=numberParametersTotal)
	colnames(param.stdev)<-namesForPriorMatrix
	rownames(param.stdev)<-paste0("Gen ", c(1: nStepsPRC), sep="")

	#initialize guesses, if needed
	#if (length(startingValuesGuess)==0) { #if no user guesses, try pulling a value from the prior
	#	startingValuesGuess<-rep(NA,length(startingPriorsFns))
	#	for (i in 1:length(startingPriorsFns)) {
	#		startingValuesGuess[i]<-pullFromPrior(startingPriorsValues[,i],startingPriorsFns[i])
	#	}
	#}
	#if (length(intrinsicValuesGuess)==0) { #if no user guesses, try pulling a value from the prior
	#	intrinsicValuesGuess<-rep(NA,length(intrinsicPriorsFns))
	#	for (i in 1:length(intrinsicPriorsFns)) {
	#		intrinsicValuesGuess[i]<-pullFromPrior(intrinsicPriorsValues[,i],intrinsicPriorsFns[i])
	#	}
	#}
	#if (length(extrinsicValuesGuess)==0) { #if no user guesses, try pulling a value from the prior
	#	extrinsicValuesGuess<-rep(NA,length(extrinsicPriorsFns))
	#	for (i in 1:length(extrinsicPriorsFns)) {
	#		extrinsicValuesGuess[i]<-pullFromPrior(extrinsicPriorsValues[,i],extrinsicPriorsFns[i])
	#	}
	#}

	if (is.na(StartSims)) {
		StartSims<-1000*numberParametersFree
	}

	#Figure out how many iterations to use for optimization in Geiger.
	
	#it actually runs faster without checking for cores. And we parallelize elsewhere
	brown<-makeQuiet(fitContinuous(phy=phy, dat=traits, model="BM", ncores=1, control=list(niter=100)))
	lambda<-makeQuiet(fitContinuous(phy=phy, dat=traits, model="lambda", ncores=1, control=list(niter=100)))
	delta<-makeQuiet(fitContinuous(phy=phy, dat=traits, model="delta", ncores=1, control=list(niter=100)))
	ou<-makeQuiet(fitContinuous(phy=phy, dat=traits, model="OU", ncores=1, control=list(niter=100)))
	white<-makeQuiet(fitContinuous(phy=phy, dat=traits, model="white", ncores=1, control=list(niter=100)))
	#
	niter.brown.g <- round(max(10, min(niter.goal/solnfreq(brown),100)))
	niter.lambda.g <- round(max(10, min(niter.goal/solnfreq(lambda),100)))
	niter.delta.g <- round(max(10, min(niter.goal/solnfreq(delta),100)))
	niter.OU.g <- round(max(10, min(niter.goal/solnfreq(ou),100)))
	niter.white.g <- round(max(10, min(niter.goal/solnfreq(white),100)))
	#
	message(paste0("Setting number of starting points for Geiger optimization to",
		paste0("\n   ",niter.brown.g, " for Brownian motion"),
		paste0("\n   ",niter.lambda.g, " for lambda"),
		paste0("\n   ",niter.delta.g, " for delta"),
		paste0("\n   ",niter.OU.g, " for OU"),
		paste0("\n   ",niter.white.g, " for white noise")))
	#
	#---------------------- Initial Simulations (Start) ------------------------------
	# See Wegmann et al. Efficient Approximate Bayesian Computation Coupled With Markov Chain Monte Carlo Without Likelihood. 
		# Genetics (2009) vol. 182 (4) pp. 1207-1218 for more on the method.
	# We are doing pls and scaling built into pls. Unlike Wegmann et al., we are doing PLS for each parameter separately. 
	# Otherwise, the PLS tends to optimize for just one parameter, and estimates for the less-favored one are quite bad
	# because the summary stats tend to be used for the other.
	#
	#Used to be = StartSims*((2^try)/2), If initial simulations are not enough, and we need to try again then new analysis will double number of initial simulations
	nrepSim<-StartSims 
	#
	input.data<-rbind(jobName, length(phy[[3]]), nrepSim, generation.time, TreeYears, timeStep, totalGenerations,
		epsilonProportion, epsilonMultiplier, nStepsPRC, numParticles, standardDevFactor)
	message(paste0("Number of initial simulations set to ", nrepSim)) #, "\n"
	message("Doing initial simulations...")
	Time<-proc.time()[[3]]
	trueFreeValues<-matrix(nrow=0, ncol= numberParametersFree)
	#set up initial sum stats as length of SSL of original data
	summaryValues<-matrix(nrow=0, ncol=length(summaryStatsLong(phy=phy, traits=traits, 
		niter.brown=200, niter.lambda=200, niter.delta=200, niter.OU=200, niter.white=200))) 
	trueFreeValuesANDSummaryValues<-makeQuiet(parallelSimulateWithPriors(nrepSim=nrepSim, coreLimit=coreLimit, phy=phy,  taxon.df=taxon.df,
		startingPriorsValues=startingPriorsValues, intrinsicPriorsValues=intrinsicPriorsValues, extrinsicPriorsValues=extrinsicPriorsValues, 
		startingPriorsFns=startingPriorsFns, intrinsicPriorsFns=intrinsicPriorsFns, extrinsicPriorsFns=extrinsicPriorsFns, 
		freevector=freevector, timeStep=timeStep, intrinsicFn=intrinsicFn, extrinsicFn=extrinsicFn, multicore=multicore, 
		niter.brown=niter.brown.g, niter.lambda=niter.lambda.g, niter.delta=niter.delta.g, niter.OU=niter.OU.g, niter.white=niter.white.g))
	#message("\n\n")
	#
	if(saveData){
		save(trueFreeValues,summaryValues,file=paste0("CompletedSimulations",jobName,".Rdata",sep=""))
		}
	#
	simTime<-proc.time()[[3]]-Time
	message(paste0("Initial simulations took ", round(simTime, digits=3), " seconds")) #, "\n"
	#
	#separate the simulation results: true values and the summary values
	trueFreeValuesMatrix<-trueFreeValuesANDSummaryValues[,1:numberParametersFree]
	summaryValuesMatrix<-trueFreeValuesANDSummaryValues[,-1:-numberParametersFree]
	#while(sink.number()>0) {sink()}
	#
	#---------------------- Initial Simulations (End) ------------------------------
	#
	pls.model.list <- makeQuiet(apply(trueFreeValuesMatrix, 2, returnPLSModel, summaryValuesMatrix=summaryValuesMatrix, validation=validation, 
		scale=scale, variance.cutoff = variance.cutoff))
	originalSummaryValues <- summaryStatsLong(phy, traits, niter.brown=200, niter.lambda=200, niter.delta=200, niter.OU=200, niter.white=200)
	#
	distanceVector<-abcDistance(summaryValuesMatrix, originalSummaryValues, pls.model.list)
	#
	#----------------- Find distribution of distances (Start) ----------------------
	message("Finding distribution of distances...") #\n
	#this gives the distance such that epsilonProportion of the simulations starting from a given set of values will be rejected
	epsilonDistance<-quantile(distanceVector, probs=epsilonProportion) 
	toleranceVector<-rep(epsilonDistance, nStepsPRC)
	#
	if(nStepsPRC>1){
		for (step in 2:nStepsPRC) {
			toleranceVector[step]<-toleranceVector[step-1]*epsilonMultiplier
			}
		}
	#----------------- Find distribution of distances (End) ---------------------
	#
	#------------------ ABC-PRC (Start) ------------------
	message("Beginning partial rejection control algorithm...")
	#
	nameVector<-c("generation", "attempt", "id", "parentid", "distance", "weight")
	if (plot) {
		plot(x=c(min(intrinsicPriorsValues), max(intrinsicPriorsValues)), y=c(0, 5*max(toleranceVector)), type="n")
		}
	#
	for (i in 1:dim(startingPriorsValues)[2]) {
		nameVector<-append(nameVector, paste0("StartingStates ", i, sep=""))
		}
	for (i in 1:dim(intrinsicPriorsValues)[2]) {
		nameVector<-append(nameVector, paste0("IntrinsicValue ", i, sep=""))
		}
	for (i in 1:dim(extrinsicPriorsValues)[2]) {
		nameVector<-append(nameVector, paste0("ExtrinsicValue ", i, sep=""))
		}
	#
	#stores weights for each particle. Initially, assume infinite number of possible particles (so might not apply in discrete case)
	particleWeights=rep(0, numParticles) 
	#stores parameters in model for each particle
	particleParameters<-matrix(nrow=numParticles, 
		ncol=dim(startingPriorsValues)[2] +  dim(intrinsicPriorsValues)[2] + dim(extrinsicPriorsValues)[2]) 
	particleDistance=rep(NA, numParticles)
	particle<-1
	attempts<-0
	particleDataFrame<-data.frame()
	if(verboseParticles){
		message("\nsuccesses ", "  attempts ", "  expected number of attempts required")
		}
	
	#**
	
	start.time<-proc.time()[[3]]
	particleList<-list()
	#message("Beginning partial rejection control algorithm...")
	while (particle<=numParticles) {
		attempts<-attempts+1
		if(attempts>maxAttempts){
			stop("maxAttempts exceeded in while() loop of PRC algorithm")
			}
		#
		newparticleList<-list(abcparticle(id=particle, generation=1, weight=0))
		newparticleList[[1]]<-initializeStatesFromMatrices(particle=newparticleList[[1]], 
			startingPriorsValues=startingPriorsValues, startingPriorsFns=startingPriorsFns, intrinsicPriorsValues=intrinsicPriorsValues, 
			intrinsicPriorsFns=intrinsicPriorsFns, extrinsicPriorsValues=extrinsicPriorsValues, extrinsicPriorsFns=extrinsicPriorsFns)

		newparticleList[[1]]$distance<-abcDistance(
			summaryValuesMatrix=summaryStatsLong(phy=phy, 
				traits=doSimulationWithPossibleExtinction(phy=NULL,  taxon.df=taxon.df, intrinsicFn=intrinsicFn, extrinsicFn=extrinsicFn, 
					startingValues=newparticleList[[1]]$startingValues, intrinsicValues=newparticleList[[1]]$intrinsicValues, 
					extrinsicValues=newparticleList[[1]]$extrinsicValues, timeStep=timeStep, checkTimeStep=FALSE), 
					niter.brown=niter.brown.g, niter.lambda=niter.lambda.g, niter.delta=niter.delta.g, 
					niter.OU=niter.OU.g, niter.white=niter.white.g)
			, originalSummaryValues=originalSummaryValues, pls.model.list=pls.model.list )

		if (is.na(newparticleList[[1]]$distance)) {
			warning("newparticleList[[1]]$distance = NA, likely an underflow/overflow problem")
			newparticleList[[1]]$id <-  (-1)
			newparticleList[[1]]$weight<- 0
		}else{
			if (is.na(toleranceVector[1])) {
				warning("toleranceVector[1] = NA")
				newparticleList[[1]]$id <- (-1)
				newparticleList[[1]]$weight <- 0
			}else{
				if ((newparticleList[[1]]$distance) < toleranceVector[1]) {
					newparticleList[[1]]$id <- particle
					newparticleList[[1]]$weight <- 1/numParticles
					particleWeights[particle] <- 1/numParticles
					particle<-particle+1
					particleList<-append(particleList, newparticleList)
				}else{
					newparticleList[[1]]$id <- (-1)
					newparticleList[[1]]$weight <- 0
					}
				}
			}
		#while(sink.number()>0) {sink()}
		vectorForDataFrame<-c(1, attempts, newparticleList[[1]]$id, 0, newparticleList[[1]]$distance, 
			newparticleList[[1]]$weight, newparticleList[[1]]$startingValues, newparticleList[[1]]$intrinsicValues,
			newparticleList[[1]]$extrinsicValues)
		
		if(diagnosticPRCmode){
			message("\n\nlength of vectorForDataFrame = ", length(vectorForDataFrame), "\n", "length of startingValues = ", 
				length(newparticleList[[1]]$startingValues), "\nlength of intrinsicValues = ", length(newparticleList[[1]]$intrinsicValues), 
				"\nlength of extrinsicValues = ", length(newparticleList[[1]]$extrinsicValues), "\ndistance = ", newparticleList[[1]]$distance,
				"\nweight = ", newparticleList[[1]]$weight, "\n", vectorForDataFrame, "\n")
			}
			
		particleDataFrame<-rbind(particleDataFrame, vectorForDataFrame)
		#	
		if(verboseParticles){
			message(paste(particle-1,"          ", attempts,"         ",
				floor(numParticles*attempts/particle),"         ", 
				signif(newparticleList[[1]]$startingValues,2),"  ", signif(newparticleList[[1]]$intrinsicValues,2),"  ", 
				signif(newparticleList[[1]]$extrinsicValues,2),"  ", signif(newparticleList[[1]]$distance,2)))
			}
		#

		#
		} #while (particle<=numParticles) bracket
	#
	names(particleDataFrame)<-nameVector
	dataGenerationStep=1
	time<-proc.time()[[3]]-start.time
	time.per.gen<-time
	#
	#rejects.gen.one<-(dim(subset(particleDataFrame, particleDataFrame$id<0))[1])/(dim(subset(particleDataFrame,))[1])
	#rejects<-c()
	#
	for (i in 1:numberParametersTotal){
		param.stdev[1,i]<-c(sd(subset(particleDataFrame, particleDataFrame$id>0)[,6+i]))
		weightedMeanParam[1,i]<-weighted.mean(subset(particleDataFrame, particleDataFrame$id>0)[,6+i], 
			subset(particleDataFrame, particleDataFrame$id>0)[,6])
		#c(mean(subset(particleDataFrame, X3>0)[,7:dim(particleDataFrame)[2]])/subset(particleDataFrame, X3>0)[,6])
		}
	#
	if(saveData){
		save.image(file=paste0("WS", jobName, ".Rdata", sep=""))
		}
	
	prcResults<-vector("list")
	prcResults$input.data<-input.data
	prcResults$PriorMatrix<-PriorMatrix
	prcResults$particleDataFrame<-particleDataFrame
	names(prcResults$particleDataFrame)<-nameVector
	prcResults$toleranceVector<-toleranceVector
	prcResults$phy<-phy
	prcResults$traits<-traits
	prcResults$simTime<-simTime
	prcResults$time.per.gen<-time.per.gen

	if(saveData){
		save(prcResults, file=paste0("partialResults", jobName, ".txt", sep=""))
		}

	particleStartTime<-proc.time()[[3]]
	
	#**
	
	
	while (dataGenerationStep < nStepsPRC) {
		dataGenerationStep<-dataGenerationStep+1
		if(verboseParticles){
			message("\n", "STARTING DATA GENERATION STEP ", dataGenerationStep, "\n\n\n")
			}
			
			
			
		start.time<-proc.time()[[3]]
		particleWeights<-particleWeights/(sum(particleWeights,na.rm=TRUE)) #normalize to one
		if(verboseParticles){
			message("particleWeights\n", particleWeights, "\n\n")
			}
		oldParticleList<-particleList
		oldParticleWeights<-particleWeights
		#stores weights for each particle. 
			#Initially, assume infinite number of possible particles (so might not apply in discrete case)
		particleWeights=rep(0, numParticles) 
		#stores parameters in model for each particle
		particleParameters<-matrix(nrow=numParticles,
			ncol=dim(startingPriorsValues)[2] +  dim(intrinsicPriorsValues)[2] + dim(extrinsicPriorsValues)[2]) 
		particleDistance=rep(NA, numParticles)
		particle<-1
		attempts<-0
		if(verboseParticles){
			message("successes ", "  attempts ", "  expected number of attempts required") #\n
			}
		particleList<-list()
		weightScaling=0;
		while (particle<=numParticles) {
			attempts<-attempts+1
			particleToSelect<-which.max(as.vector(rmultinom(1, size = 1, prob=oldParticleWeights)))
			
			#message("particle to select = ", particleToSelect, "\n")
			#message("dput(oldParticleList)\n")
			#dput(oldParticleList)
			#message("dput(oldParticleList[particleToSelect])\n")
			#dput(oldParticleList[particleToSelect])
			#message("dput(oldParticleList[[particleToSelect]])\n")
			#dput(oldParticleList[[particleToSelect]])
			
			newparticleList<-list(oldParticleList[[particleToSelect]])
			
			#message("dput(newparticleList[[1]])\n")
			#dput(newparticleList[[1]])
			#message("mutateStates\n")

			newparticleList[[1]]<-mutateStates(particle=newparticleList[[1]], 
				startingPriorsValues=startingPriorsValues, startingPriorsFns=startingPriorsFns, 
				intrinsicPriorsValues=intrinsicPriorsValues, intrinsicPriorsFns=intrinsicPriorsFns, 
				extrinsicPriorsValues=extrinsicPriorsValues, extrinsicPriorsFns=extrinsicPriorsFns, 
				standardDevFactor=standardDevFactor)
			#
			#message("dput(newparticleList[[1]]) AFTER MUTATE STATES\n")
			#dput(newparticleList[[1]])
			#
			newparticleList[[1]]$distance<-abcDistance(
				summaryValuesMatrix=summaryStatsLong(phy=phy, 
					traits=doSimulationWithPossibleExtinction(phy=NULL,  taxon.df=taxon.df,
						intrinsicFn=intrinsicFn, extrinsicFn=extrinsicFn, 
						startingValues=newparticleList[[1]]$startingValues, intrinsicValues=newparticleList[[1]]$intrinsicValues, 
						extrinsicValues=newparticleList[[1]]$extrinsicValues, timeStep=timeStep, checkTimeStep=FALSE), 
					niter.brown=niter.brown.g, niter.lambda=niter.lambda.g, niter.delta=niter.delta.g, 
					niter.OU=niter.OU.g, niter.white=niter.white.g)
				, originalSummaryValues=originalSummaryValues, pls.model.list=pls.model.list )
			#
			if (plot) {
				plotcol="grey"
				if (newparticleList[[1]]$distance<toleranceVector[dataGenerationStep]) {
					plotcol="black"
					if (dataGenerationStep==length(toleranceVector)) {
						plotcol="red"
						}
					}
				text(x=newparticleList[[1]]$intrinsicValues, y=newparticleList[[1]]$distance, labels= dataGenerationStep, col=plotcol)
				}
			#
			#message("dput(newparticleList[[1]]) AFTER computeABCDistance\n")
			#dput(newparticleList[[1]])
			#	
			if (is.na(newparticleList[[1]]$distance)) {
				#message("Error with Geiger?  newparticleList[[1]]$distance = NA\n")
				#while(sink.number()>0) {sink()}
				#warning("newparticleList[[1]]$distance = NA")
				newparticleList[[1]]$id <- (-1)
				newparticleList[[1]]$weight <- 0
			}else{
				if (newparticleList[[1]]$distance < toleranceVector[dataGenerationStep]) {
					newparticleList[[1]]$id <- particle
					particle<-particle+1
					particleList<-append(particleList, newparticleList)
					#now get weights, using correction in Sisson et al. 2007
					newWeight=0
					for (i in 1:length(oldParticleList)) {
						lnTransitionProb=log(1)
						LLTPstart<-sapply(length(newparticleList[[1]]$startingValues),
							function(j) getlnTransitionProb(newvalue = newparticleList[[1]]$startingValues[j],
								meantouse = oldParticleList[[i]]$startingValues[j], 
								Fn=startingPriorsFns[j],
								priorValues= startingPriorsValues[,j],
								stdFactor = standardDevFactor))
						LLTPintr<-sapply(length(newparticleList[[1]]$intrinsicValues),
							function(j) getlnTransitionProb(newvalue = newparticleList[[1]]$intrinsicValues[j],
								meantouse = oldParticleList[[i]]$intrinsicValues[j], 
								Fn=intrinsicPriorsFns[j],
								priorValues= intrinsicPriorsValues[,j],
								stdFactor = standardDevFactor))
						LLTPextr<-sapply(length(newparticleList[[1]]$extrinsicValues),
							function(j) getlnTransitionProb(newvalue = newparticleList[[1]]$extrinsicValues[j],
								meantouse = oldParticleList[[i]]$extrinsicValues[j], 
								Fn=extrinsicPriorsFns[j],
								priorValues= extrinsicPriorsValues[,j],
								stdFactor = standardDevFactor))
						lnTransitionProb<-lnTransitionProb+sum(LLTPstart)+sum(LLTPintr)+sum(LLTPextr)
						if(!is.finite(lnTransitionProb) || is.na(lnTransitionProb)) {
							warning(paste0("Issue with lnTransitionProb: ",
								" lnTransitionProb = ",lnTransitionProb))
							}
						newWeight<-newWeight+(oldParticleList[[i]]$weight)*exp(lnTransitionProb)
						} #for (i in 1:length(oldParticleList)) bracket
					#
					if (!is.finite(newWeight)) {
						warning(paste0("newWeight is ",newWeight))
						}
					newparticleList[[1]]$weight<- newWeight
					particleWeights[particle-1]<-newWeight
					weightScaling<-weightScaling+newWeight
				} else { 
					#else if (newparticleList[[1]]$distance < toleranceVector[dataGenerationStep]) bracket
					newparticleList[[1]]$id<- (-1)
					newparticleList[[1]]$weight<-0
					}
				}
			#
			#
			#while(sink.number()>0) {sink()}
			#print(newparticleList)
			#
			vectorForDataFrame<-c(dataGenerationStep, attempts,newparticleList[[1]]$id, particleToSelect, 
				newparticleList[[1]]$distance, newparticleList[[1]]$weight, newparticleList[[1]]$startingValues, 
				newparticleList[[1]]$intrinsicValues, newparticleList[[1]]$extrinsicValues)
			#	
			if(diagnosticPRCmode){
				message("\n\nlength of vectorForDataFrame = ", length(vectorForDataFrame), "\n", "length of startingValues = ", 
					length(newparticleList[[1]]$startingValues), "\nlength of intrinsicValues = ", length(newparticleList[[1]]$intrinsicValues), 
					"\nlength of extrinsicValues = ", length(newparticleList[[1]]$extrinsicValues), "\ndistance = ", newparticleList[[1]]$distance,
					"\nweight = ", newparticleList[[1]]$weight, "\n", vectorForDataFrame, "\n")
				}
			
				
			#NOTE THAT WEIGHTS AREN'T NORMALIZED IN THIS DATAFRAME
			particleDataFrame<-rbind(particleDataFrame, vectorForDataFrame) 
			#
			if(verboseParticles){
				message(paste(particle-1,"          ", attempts,"         ",
					floor(numParticles*attempts/particle),"         ", 
					signif(newparticleList[[1]]$startingValues,2),"  ", signif(newparticleList[[1]]$intrinsicValues,2),"  ", 
					signif(newparticleList[[1]]$extrinsicValues,2),"  ", signif(newparticleList[[1]]$distance,2)))
				}
			#
			
			} #while (particle<=numParticles) bracket
		#
		#
		particleDataFrame[which(particleDataFrame$generation==dataGenerationStep), ]$weight<-particleDataFrame[
			which(particleDataFrame$generation==dataGenerationStep), ]$weight/(sum(
				particleDataFrame[which(particleDataFrame$generation==dataGenerationStep), ]$weight))
		#	
		time2<-proc.time()[[3]]-start.time
		time.per.gen<-c(time.per.gen, time2)
		#rejects.per.gen<-(dim(subset(particleDataFrame, particleDataFrame$id<0))[1])/(
			# dim(subset(particleDataFrame[which(particleDataFrame$generation==dataGenerationStep),],))[1])
		#
		#rejects<-c(rejects, rejects.per.gen)
		sub1<-subset(particleDataFrame, particleDataFrame$generation==dataGenerationStep)
		sub2<-subset(sub1, sub1$id>0)
		#
		for (i in 1:numberParametersTotal){
			param.stdev[dataGenerationStep,i]<-c(sd(sub2[,6+i]))
			weightedMeanParam[dataGenerationStep,i]<-weighted.mean(sub2[,6+i], sub2[,6])
			}
		#
		if (stopRule){	#this will stop the PRC from running out to max number of generations if all params are below stopValue
			FF<-rep(1, dim(weightedMeanParam)[2])
			for (check.weightedMeanParam in 1:length(FF)){
				#
				if (is.na(abs(weightedMeanParam[dataGenerationStep, check.weightedMeanParam]-weightedMeanParam[dataGenerationStep-1,
						check.weightedMeanParam])/mean(weightedMeanParam[dataGenerationStep, check.weightedMeanParam],
						#this && is here to make sure any NAs are from fixed params and not miscalculations.
						weightedMeanParam[dataGenerationStep-1, check.weightedMeanParam]) <= stopValue) && mean(
						weightedMeanParam[dataGenerationStep, check.weightedMeanParam], weightedMeanParam[dataGenerationStep-1,
						check.weightedMeanParam]) == 0) {  
					FF[check.weightedMeanParam]<-0
				}else{
					stopValueWeightTest<-abs(weightedMeanParam[dataGenerationStep, check.weightedMeanParam]
						-weightedMeanParam[dataGenerationStep-1,check.weightedMeanParam])
					stopValueWeightTest<-stopValueWeightTest/
							mean(weightedMeanParam[dataGenerationStep, check.weightedMeanParam],
								weightedMeanParam[dataGenerationStep-1, check.weightedMeanParam])
					if (stopValueWeightTest <= stopValue){
						FF[check.weightedMeanParam]<-0
						}
					}
				#print(FF)
				}
			if (sum(FF)==0){
				message("\n\n\nweightedMeanParam is < ", stopValue, "Analysis is being terminated at", dataGenerationStep
					, "instead of continuing to ", nStepsPRC, "\n\n\n")
				dataGenerationStep<-nStepsPRC
				}
			}	
		#
		if(saveData){
			save.image(file=paste0("WS", jobName, ".Rdata", sep=""))
			}
		#
		prcResults<-vector("list")
		prcResults$input.data<-input.data
		prcResults$PriorMatrix<-PriorMatrix
		prcResults$particleDataFrame<-particleDataFrame
		names(prcResults$particleDataFrame)<-nameVector
		prcResults$toleranceVector<-toleranceVector
		prcResults$phy<-phy
		prcResults$traits<-traits
		prcResults$simTime
		prcResults$time.per.gen<-time.per.gen
		#
		if(saveData){
			save(prcResults, file=paste0("partialResults", jobName, ".txt", sep=""))
			}
			
			
			
		#
		} #while (dataGenerationStep < nStepsPRC) bracket
	#
	names(particleDataFrame)<-nameVector
	if(plot) {
		#dev.new()
		plot(x=c(min(intrinsicPriorsValues), max(intrinsicPriorsValues)), y=c(0, 1), type="n")
		for (i in 1:(length(toleranceVector)-1)) {
			graycolor<-gray(0.5*(length(toleranceVector)-i)/length(toleranceVector))
			lines(density(subset(particleDataFrame, particleDataFrame$generation==i)[, 8]), col= graycolor)
			}
		lines(density(subset(particleDataFrame, particleDataFrame$generation==length(toleranceVector))[, 8]), col= "red")
		}
	particleTime<-proc.time()[[3]]-particleStartTime
	message(paste0("Collection of simulation particles under PRC completed in ",particleTime," seconds..."))
	#---------------------- ABC-PRC (End) --------------------------------
	#
	input.data<-rbind(jobName, length(phy[[3]]), nrepSim, generation.time, TreeYears, timeStep, totalGenerations, epsilonProportion,
		epsilonMultiplier, nStepsPRC, numParticles, standardDevFactor)
	#
	time3<-proc.time()[[3]]
	genTimes<-c(time.per.gen, time3)
	#
	prcResults<-vector("list")
	prcResults$input.data<-input.data
	prcResults$PriorMatrix<-PriorMatrix
	prcResults$particleDataFrame<-particleDataFrame
	prcResults$toleranceVector<-toleranceVector
	prcResults$phy<-phy
	prcResults$traits<-traits
	prcResults$simTime<-simTime
	prcResults$time.per.gen<-genTimes
	prcResults$credibleInt <-credibleInt(particleDataFrame)
	prcResults$HPD <-highestPostDens(particleDataFrame)
	#
	if(multicore){
		registerMulticoreEnv(nCore=1)
		}
	#
	functionTime<-proc.time()[[3]]-functionStartTime
	message(paste0("Function completed in ",functionTime," seconds."))
	#
	#print(prcResults)
	return(prcResults)
	}

getlnTransitionProb<-function(newvalue,meantouse,Fn,priorValues,stdFactor){
		#newvalue = newparticleList[[1]]$startingValues[j],
		#meantouse = oldParticleList[[i]]$startingValues[j], 
		#Fn=startingPriorsFns[j],
		#priorValues= startingPriorsValues[,j],
		#stdFactor = standardDevFactor
										
	if (Fn=="uniform") {
		sdtouse<-stdFactor*((max(priorValues)-min(priorValues))/sqrt(12))
		#print(paste0("Fn is uniform and sdtouse = ", sdtouse))
	}
	else if (Fn=="exponential") {
		sdtouse<-stdFactor*(1/priorValues[1])
		#print(paste0("Fn is exponential and sdtouse = ", sdtouse))
	}
	else {
		sdtouse<-stdFactor*(priorValues[2])
	}
	#
	lnlocalTransitionProb<-dnorm(newvalue, mean=meantouse, sd=sdtouse,log=TRUE
		) - ((log(1)/pnorm(min(priorValues), mean=meantouse, sd=sdtouse, lower.tail=TRUE, log.p=TRUE))
			* pnorm(max(priorValues), mean=meantouse , sd=sdtouse, lower.tail=FALSE, log.p=TRUE))
	if(length(lnlocalTransitionProb)!=1){
		#print(lnlocalTransitionProb)
		stop("Somehow, multiple lnlocalTransitionProb values produced")
		}
	if (is.nan(lnlocalTransitionProb)) {  #to prevent lnlocalTransitionProb from being NaN (if pnorm=0)
		lnlocalTransitionProb<-.Machine$double.xmin
	}
	if (min(priorValues)==max(priorValues)) {
		lnlocalTransitionProb=log(1)
	}
	if(!is.finite(lnlocalTransitionProb) || is.na(lnlocalTransitionProb)) {
		message(paste0("issue with lnlocalTransitionProb = ",lnlocalTransitionProb))
		}
	return(lnlocalTransitionProb)
	}


#' @rdname doRun
#' @export
doRun_rej<-function(
	phy, traits, intrinsicFn, extrinsicFn, startingPriorsValues, startingPriorsFns, 
	intrinsicPriorsValues, intrinsicPriorsFns, extrinsicPriorsValues, extrinsicPriorsFns, 
	#startingValuesGuess=c(), intrinsicValuesGuess=c(), extrinsicValuesGuess=c(), 
	TreeYears=max(branching.times(phy)) * 1e6, 
	generation.time=1000, multicore=FALSE, coreLimit=NA, validation="CV", scale=TRUE, variance.cutoff=95,
	niter.goal=5, 
	standardDevFactor=0.20, StartSims=NA, jobName=NA, abcTolerance=0.1, 
	checkpointFile=NULL, checkpointFreq=24, savesims=FALSE) {	
	
	#library(geiger)
	#library(pls)
	if (!is.binary.tree(phy)) {
		warning("Tree is not fully dichotomous, this may lead to issues")
	}
	startTime<-proc.time()[[3]]

	timeStep<-generation.time/TreeYears
	message(paste0("The effective timeStep for this tree will be ",signif(timeStep),
		", as a proportion of tree height (root to furthest tip)..."))

	edgesRescaled<-phy$edge.length/max(node.depth.edgelength(phy))
	message("Rescaling edge lengths relative to maximum tip-to-root distance...")
	
	if(max(edgesRescaled) < timeStep) {
		stop("Tree has *NO* rescaled branches longer than generation.time/TreeYears, no simulated evol change can occur!")
	}
	if(min(edgesRescaled) < timeStep) {
		warning("Tree has rescaled branches shorter than generation.time/TreeYears; no evol change can be assigned to these, and geiger functions may fail!")
		#print("Tree has zero or nearly zero length branches")
	}
	
	totalGenerations<-sum(sapply(edgesRescaled,function(x) floor(x/timeStep)))
	message("Given generation time, a total of ",round(totalGenerations)," generations will occur over this tree")
		
	#splits<-getSimulationSplits(phy) #initialize this info
	taxon.df <- getTaxonDFWithPossibleExtinction(phy=phy)

	# get freevector
	freevector<-getFreeVector(startingPriorsFns=startingPriorsFns, startingPriorsValues=startingPriorsValues, 
						intrinsicPriorsFns=intrinsicPriorsFns, intrinsicPriorsValues=intrinsicPriorsValues,
						extrinsicPriorsFns=extrinsicPriorsFns, extrinsicPriorsValues=extrinsicPriorsValues)
	numberParametersTotal<-length(freevector)
	numberParametersFree<-sum(freevector)

	#create PriorMatrix
	namesForPriorMatrix<-c()
	PriorMatrix<-matrix(c(startingPriorsFns, intrinsicPriorsFns, extrinsicPriorsFns), nrow=1, ncol=numberParametersTotal)
	for (a in 1:dim(startingPriorsValues)[2]) {
		namesForPriorMatrix<-c(paste0("StartingStates ", a, sep=""))
	}
	for (b in 1:dim(intrinsicPriorsValues)[2]) {
		namesForPriorMatrix<-append(namesForPriorMatrix, paste0("IntrinsicValue ", b, sep=""))
	}
	#print(extrinsicPriorsValues)
	for (c in 1:dim(extrinsicPriorsValues)[2]) {
		namesForPriorMatrix <-append(namesForPriorMatrix, paste0("ExtrinsicValue ", c, sep=""))
	}
	PriorMatrix<-rbind(PriorMatrix, cbind(startingPriorsValues, intrinsicPriorsValues, extrinsicPriorsValues))
	colnames(PriorMatrix)<-namesForPriorMatrix
	rownames(PriorMatrix)<-c("shape", "value1", "value2")	

	#initialize guesses, if needed
	#if (length(startingValuesGuess)==0) { #if no user guesses, try pulling a value from the prior
	#	startingValuesGuess<-rep(NA,length(startingPriorsFns))
	#	for (i in 1:length(startingPriorsFns)) {
	#		startingValuesGuess[i]<-pullFromPrior(startingPriorsValues[,i],startingPriorsFns[i])
	#	}
	#}
	#if (length(intrinsicValuesGuess)==0) { #if no user guesses, try pulling a value from the prior
	#	intrinsicValuesGuess<-rep(NA,length(intrinsicPriorsFns))
	#	for (i in 1:length(intrinsicPriorsFns)) {
	#		intrinsicValuesGuess[i]<-pullFromPrior(intrinsicPriorsValues[,i],intrinsicPriorsFns[i])
	#	}
	#}
	#if (length(extrinsicValuesGuess)==0) { #if no user guesses, try pulling a value from the prior
	#	extrinsicValuesGuess<-rep(NA,length(extrinsicPriorsFns))
	#	for (i in 1:length(extrinsicPriorsFns)) {
	#		extrinsicValuesGuess[i]<-pullFromPrior(extrinsicPriorsValues[,i],extrinsicPriorsFns[i])
	#	}
	#}

	if (is.na(StartSims)) {
		StartSims<-1000*numberParametersFree
	}

	#Used to be multiple tries where nrepSim = StartSims*((2^try)/2).  
		#If initial simulations are not enough, and we need to try again then new analysis will double number of initial simulations
	nrepSim<-StartSims 
	message(paste0("Number of initial simulations set to ", nrepSim)) 	#, "\n"
	if(!is.null(checkpointFile)) {
		save(list=ls(),file=paste0(checkpointFile,".intialsettings.Rsave",sep=""))
	}

	#Figure out how many iterations to use for optimization in Geiger.
	#it actually runs faster without checking for cores. And we parallelize elsewhere
	brown<-makeQuiet(fitContinuous(phy=phy, dat=traits, model="BM", ncores=1, control=list(niter=100))) 
	lambda<-makeQuiet(fitContinuous(phy=phy, dat=traits, model="lambda", ncores=1, control=list(niter=100)))
	delta<-makeQuiet(fitContinuous(phy=phy, dat=traits, model="delta", ncores=1, control=list(niter=100)))
	ou<-makeQuiet(fitContinuous(phy=phy, dat=traits, model="OU", ncores=1, control=list(niter=100)))
	white<-makeQuiet(fitContinuous(phy=phy, dat=traits, model="white", ncores=1, control=list(niter=100)))


	niter.brown.g <- round(max(10, min(niter.goal/solnfreq(brown),100)))
	niter.lambda.g <- round(max(10, min(niter.goal/solnfreq(lambda),100)))
	niter.delta.g <- round(max(10, min(niter.goal/solnfreq(delta),100)))
	niter.OU.g <- round(max(10, min(niter.goal/solnfreq(ou),100)))
	niter.white.g <- round(max(10, min(niter.goal/solnfreq(white),100)))
	#
	message(paste0("Setting number of starting points for Geiger optimization to ",
		paste0("\n",niter.brown.g, " for Brownian motion"),
		paste0("\n",niter.lambda.g, " for lambda"),
		paste0("\n",niter.delta.g, " for delta"),
		paste0("\n",niter.OU.g, " for OU"),
		paste0("\n",niter.white.g, " for white noise")))
	
	trueFreeValuesANDSummaryValues<-parallelSimulateWithPriors(nrepSim=nrepSim, coreLimit=coreLimit, phy=phy,  taxon.df=taxon.df,
		startingPriorsValues=startingPriorsValues, intrinsicPriorsValues=intrinsicPriorsValues, extrinsicPriorsValues=extrinsicPriorsValues, 
		startingPriorsFns=startingPriorsFns, intrinsicPriorsFns=intrinsicPriorsFns, extrinsicPriorsFns=extrinsicPriorsFns, 
		freevector=freevector, timeStep=timeStep, intrinsicFn=intrinsicFn, extrinsicFn=extrinsicFn, 
		multicore=multicore, checkpointFile=checkpointFile, checkpointFreq=checkpointFreq, 
		niter.brown=niter.brown.g, niter.lambda=niter.lambda.g, niter.delta=niter.delta.g, niter.OU=niter.OU.g, niter.white=niter.white.g)

	#message("\n\n")
	simTime<-proc.time()[[3]]-startTime
	message(paste0("Initial simulations took ", round(simTime, digits=3), " seconds")) #, "\n"

	#separate the simulation results: true values and the summary values
	trueFreeValuesMatrix<-trueFreeValuesANDSummaryValues[,1:numberParametersFree]
	summaryValuesMatrix<-trueFreeValuesANDSummaryValues[,-1:-numberParametersFree]

	if (savesims){
		save(trueFreeValuesMatrix, summaryValuesMatrix, simTime, file=paste0("sims", jobName, ".Rdata", sep=""))
		}

	res<-makeQuiet(PLSRejection(summaryValuesMatrix=summaryValuesMatrix, trueFreeValuesMatrix=trueFreeValuesMatrix,
		phy=phy, traits=traits, abcTolerance=abcTolerance))
	
	#save(abcDistancesRaw, abcDistancesRawTotal, abcDistances, abcResults, particleDataFrame, file="")
	input.data<-rbind(jobName, length(phy[[3]]), generation.time, TreeYears, timeStep, totalGenerations, StartSims, standardDevFactor, abcTolerance)
	#print(res)
	
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
	rejectionResults$credibleInt<-credibleInt(res$particleDataFrame)
	rejectionResults$HPD<-highestPostDens(res$particleDataFrame)

	return(rejectionResults)
}

	
	
	
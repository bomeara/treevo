##This seems to be working if partialResults does not exist.  If checkpoint=TRUE, then run fails.

#TreeYears = 1000000 if tree length is in in millions of years, 1000 if in thousand, etc.
#the doRun_prc function takes input from the user and then automatically guesses optimal parameters, though user overriding is also possible.
#the guesses are used to do simulations near the expected region. If omitted, they are set to the midpoint of the input parameter matrices



#' Approximate Bayesian computation for comparative methods-Partial Rejection
#' Control
#'
#' Starts the abc-prc run
#'
#' This function performs an abc-prc analysis using an input phylogeny (phy),
#' character data set (traits), models (intrinsicFn, extrinsicFn), and priors
#' (startingPriorsValues, startingPriorsFns, intrinsicPriorsValues,
#' intrinsicPriorsFns, extrinsicPriorsValues, extrinsicPriorsFns). Pulling from
#' the priors, it simulates an initial set of simulations (StartSims). The set
#' of simulations are boxcox transformed and a pls regression is performed for
#' each free parameter to determine the most informative summary statistics
#' (vipthresh). The euclidean distance of each initial simulation's vip summary
#' stats to the input data is calculated, and tolerance is set based on
#' epsilonProportion and epsilonMultiplier.
#'
#' The prc analysis begins with generation 1 and continues through nStepsPRC.
#' Single simulations (or particles) are accepted if the distance of vip
#' summary stats to the original data is less than the tolerance. A generation
#' is complete when enough particles (numParticles) have been accepted. These
#' particles make up the distribution for the next generation. The accepted
#' particles from the final generation describe the posterior distributions of
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
#' @param intrinsicValuesGuess Optional guess of intrinsic values
#' @param extrinsicValuesGuess Optional guess of extrinsic values
#' @param TreeYears Unit length of phy
#' @param numParticles Number of accepted particles per generation
#' @param standardDevFactor Standard deviation for mutating states
#' @param StartSims Number of initial simulations
#' @param plot If TRUE, plots distance of each simulation
#' @param epsilonProportion Sets tolerance for initial simulations
#' @param epsilonMultiplier Sets tolerance on subsequent generations
#' @param nStepsPRC Number of generations
#' @param jobName Optional job name
#' @param stopRule If TRUE, an analysis will be terminated prior to set number
#' of generations if the ratio of all parameters from one generation to the
#' next falls below the stopValue
#' @param stopValue Threshold value for terminating an analysis prior to
#' nStpesPRC
#' @param multicore If TRUE, initial simulations will be split among coreLimit
#' nodes
#' @param coreLimit Number of cores for initial simulations
#' @param validation Cross Validation procedure for abc
#' @param scale scale for pls.model.list
#' @param variance.cutoff variance cutoff for pls.model.list
#' @param niter.goal Adjust number of starting points for Geiger to return the best parameter estimates this number of times on average
#' @return \item{input.data}{Input variables: jobName, number of taxa, nrepSim,
#' treeYears, epsilonProportion, epsilonMultiplier, nStepsPRC, numParticles,
#' standardDevFactor} \item{PriorMatrix}{Matrix of prior distributions}
#' \item{particleDataFrame}{DataFrame with information from each simulation,
#' including generation, attempt, id, parentid, distance, weight, and parameter
#' states} \item{toleranceVector}{Final tolerance vector} \item{phy}{Input
#' phylogeny} \item{traits}{Input traits} \item{simTime}{Processor time for
#' initial simulations} \item{time.per.gen}{Processor time for subsequent
#' generations} \item{whichVip}{Matrix of vip summary statistics for each free
#' parameter} \item{CredInt}{Credible Interval calculation for each free
#' parameter of the final generation} \item{HPD}{Highest Posterior Density
#' calculation each free parameter of the final generation}
#' @author Brian O'Meara and Barb Banbury
#' @references O'Meara and Banbury, unpublished; Sisson et al. 2007, Wegmann et
#' al. 2009
#' @keywords doRun doRun_prc abc
#' @examples
#'
#' data(simData)
#'
#' #doRun_prc(
#' #  phy = simPhy,
#' #  traits = simChar,
#' #  intrinsicFn=brownianIntrinsic,
#' #  extrinsicFn=nullExtrinsic,
#' #  startingPriorsFns="normal",
#' #  startingPriorsValues=matrix(c(mean(simChar[,1]), sd(simChar[,1]))),
#' #  intrinsicPriorsFns=c("exponential"),
#' #  intrinsicPriorsValues=matrix(c(10, 10), nrow=2, byrow=FALSE),
#' #  extrinsicPriorsFns=c("fixed"),
#' #  extrinsicPriorsValues=matrix(c(0, 0), nrow=2, byrow=FALSE),
#' #  TreeYears=1000,
#' #  standardDevFactor=0.2,
#' #  plot=FALSE,
#' #  StartSims=300,
#' #  epsilonProportion=0.7,
#' #  epsilonMultiplier=0.7,
#' #  nStepsPRC=5,
#' #  numParticles=100,
#' #  jobName="exampleRun",
#' #  stopRule=FALSE,
#' #  multicore=FALSE,
#' #  coreLimit=1
#' #)
#'
#'
doRun_prc<-function(phy, traits, intrinsicFn, extrinsicFn, startingPriorsValues, startingPriorsFns, intrinsicPriorsValues, intrinsicPriorsFns, extrinsicPriorsValues, extrinsicPriorsFns, startingValuesGuess=c(), intrinsicValuesGuess=c(), extrinsicValuesGuess=c(), TreeYears=1e+04, numParticles=300, standardDevFactor=0.20, StartSims=300, plot=FALSE, epsilonProportion=0.7, epsilonMultiplier=0.7, nStepsPRC=5, jobName=NA, stopRule=FALSE, stopValue=0.05, multicore=FALSE, coreLimit=NA, validation="CV", scale=TRUE, variance.cutoff=95, niter.goal=5) {

if (!is.binary.tree(phy)) {
	print("Warning: Tree is not fully dichotomous")
}

timeStep<-1/TreeYears

splits<-getSimulationSplits(phy) #initialize this info

#figure out number of free params
numberParametersTotal<-dim(startingPriorsValues)[2] +  dim(intrinsicPriorsValues)[2] + dim(extrinsicPriorsValues)[2]
numberParametersFree<-numberParametersTotal
numberParametersStarting<-0
numberParametersIntrinsic<-0
numberParametersExtrinsic<-0
freevariables<-matrix(data=NA, nrow=2, ncol=0)
titlevector<-c()
freevector<-c()

#create PriorMatrix
namesForPriorMatrix<-c()
PriorMatrix<-matrix(c(startingPriorsFns, intrinsicPriorsFns, extrinsicPriorsFns), nrow=1, ncol=numberParametersTotal)
for (a in 1:dim(startingPriorsValues)[2]) {
	namesForPriorMatrix<-c(paste("StartingStates", a, sep=""))
}
for (b in 1:dim(intrinsicPriorsValues)[2]) {
	namesForPriorMatrix<-append(namesForPriorMatrix, paste("IntrinsicValue", b, sep=""))
}
#print(extrinsicPriorsValues)
for (c in 1:dim(extrinsicPriorsValues)[2]) {
	namesForPriorMatrix <-append(namesForPriorMatrix, paste("ExtrinsicValue", c, sep=""))
}
PriorMatrix<-rbind(PriorMatrix, cbind(startingPriorsValues, intrinsicPriorsValues, extrinsicPriorsValues))
colnames(PriorMatrix)<-namesForPriorMatrix
rownames(PriorMatrix)<-c("shape", "value1", "value2")

#Calculate freevector
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

#initialize weighted mean sd matrices
weightedMeanParam<-matrix(nrow=nStepsPRC, ncol=numberParametersTotal)
colnames(weightedMeanParam)<-namesForPriorMatrix
rownames(weightedMeanParam)<-paste("Gen", c(1: nStepsPRC), sep="")
param.stdev<-matrix(nrow=nStepsPRC, ncol=numberParametersTotal)
colnames(param.stdev)<-namesForPriorMatrix
rownames(param.stdev)<-paste("Gen", c(1: nStepsPRC), sep="")

#initialize guesses, if needed
if (length(startingValuesGuess)==0) { #if no user guesses, try pulling a value from the prior
	startingValuesGuess<-rep(NA,length(startingPriorsFns))
	for (i in 1:length(startingPriorsFns)) {
		startingValuesGuess[i]<-pullFromPrior(startingPriorsValues[,i],startingPriorsFns[i])
	}
}
if (length(intrinsicValuesGuess)==0) { #if no user guesses, try pulling a value from the prior
	intrinsicValuesGuess<-rep(NA,length(intrinsicPriorsFns))
	for (i in 1:length(intrinsicPriorsFns)) {
		intrinsicValuesGuess[i]<-pullFromPrior(intrinsicPriorsValues[,i],intrinsicPriorsFns[i])
	}
}
if (length(extrinsicValuesGuess)==0) { #if no user guesses, try pulling a value from the prior
	extrinsicValuesGuess<-rep(NA,length(extrinsicPriorsFns))
	for (i in 1:length(extrinsicPriorsFns)) {
		extrinsicValuesGuess[i]<-pullFromPrior(extrinsicPriorsValues[,i],extrinsicPriorsFns[i])
	}
}

if (is.na(StartSims)) {
	StartSims<-1000*numberParametersFree
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



			#---------------------- Initial Simulations (Start) ------------------------------
			#See Wegmann et al. Efficient Approximate Bayesian Computation Coupled With Markov Chain Monte Carlo Without Likelihood. Genetics (2009) vol. 182 (4) pp. 1207-1218 for more on the method.
      #We are doing pls and scaling built into pls. Unlike Wegmann et al., we are doing PLS for each parameter separately. Otherwise, the PLS tends to
      #optimize for just one parameter, and estimates for the less-favored one are quite bad because the summary stats tend to be used for the other.
nrepSim<-StartSims #Used to be = StartSims*((2^try)/2), If initial simulations are not enough, and we need to try again then new analysis will double number of initial simulations
input.data<-rbind(jobName, length(phy[[3]]), nrepSim, TreeYears, epsilonProportion, epsilonMultiplier, nStepsPRC, numParticles, standardDevFactor)
cat(paste("\nNumber of initial simulations set to", nrepSim, "\n"))
cat("Doing simulations:")
Time<-proc.time()[[3]]
trueFreeValues<-matrix(nrow=0, ncol= numberParametersFree)
summaryValues<-matrix(nrow=0, ncol=length(summaryStatsLong(phy, traits, niter.brown=200, niter.lambda=200, niter.delta=200, niter.OU=200, niter.white=200))) #set up initial sum stats as length of SSL of original data
trueFreeValuesANDSummaryValues<-parallelSimulation(nrepSim, coreLimit, splits, phy, startingPriorsValues, intrinsicPriorsValues, extrinsicPriorsValues, startingPriorsFns, intrinsicPriorsFns, extrinsicPriorsFns, freevector, timeStep, intrinsicFn, extrinsicFn, multicore, niter.brown=niter.brown.g, niter.lambda=niter.lambda.g, niter.delta=niter.delta.g, niter.OU=niter.OU.g, niter.white=niter.white.g)
cat("\n\n")

save(trueFreeValues,summaryValues,file=paste("CompletedSimulations",jobName,".Rdata",sep=""))
simTime<-proc.time()[[3]]-Time
cat(paste("Initial simulations took", round(simTime, digits=3), "seconds"), "\n")

#separate the simulation results: true values and the summary values
trueFreeValuesMatrix<-trueFreeValuesANDSummaryValues[,1:numberParametersFree]
summaryValuesMatrix<-trueFreeValuesANDSummaryValues[,-1:-numberParametersFree]
#while(sink.number()>0) {sink()}

			#---------------------- Initial Simulations (End) ------------------------------


  pls.model.list <- apply(trueFreeValuesMatrix, 2, returnPLSModel, summaryValuesMatrix=summaryValuesMatrix, validation=validation, scale=scale, variance.cutoff=variance.cutoff)

  originalSummaryValues <- summaryStatsLong(phy, traits, niter.brown=200, niter.lambda=200, niter.delta=200, niter.OU=200, niter.white=200)

  distanceVector<-abcDistance(summaryValuesMatrix, originalSummaryValues, pls.model.list)


			#----------------- Find distribution of distances (Start) ----------------------

	epsilonDistance<-quantile(distanceVector, probs=epsilonProportion) #this gives the distance such that epsilonProportion of the simulations starting from a given set of values will be rejected
	toleranceVector<-rep(epsilonDistance, nStepsPRC)

	if(nStepsPRC>1){
		for (step in 2:nStepsPRC) {
			toleranceVector[step]<-toleranceVector[step-1]*epsilonMultiplier
		}
	}
			#----------------- Find distribution of distances (End) ---------------------

			#------------------ ABC-PRC (Start) ------------------

			nameVector<-c("generation", "attempt", "id", "parentid", "distance", "weight")
			if (plot) {
				plot(x=c(min(intrinsicPriorsValues), max(intrinsicPriorsValues)), y=c(0, 5*max(toleranceVector)), type="n")
			}
			for (i in 1:dim(startingPriorsValues)[2]) {
				nameVector<-append(nameVector, paste("StartingStates", i, sep=""))
			}
			for (i in 1:dim(intrinsicPriorsValues)[2]) {
				nameVector<-append(nameVector, paste("IntrinsicValue", i, sep=""))
			}
			for (i in 1:dim(extrinsicPriorsValues)[2]) {
				nameVector<-append(nameVector, paste("ExtrinsicValue", i, sep=""))
			}
			particleWeights=rep(0, numParticles) #stores weights for each particle. Initially, assume infinite number of possible particles (so might not apply in discrete case)
			particleParameters<-matrix(nrow=numParticles, ncol=dim(startingPriorsValues)[2] +  dim(intrinsicPriorsValues)[2] + dim(extrinsicPriorsValues)[2]) #stores parameters in model for each particle
			particleDistance=rep(NA, numParticles)
			particle<-1
			attempts<-0
			particleDataFrame<-data.frame()
			cat("\n\n\nsuccesses", "attempts", "expected number of attempts required\n\n\n")
			start.time<-proc.time()[[3]]
			particleList<-list()

			while (particle<=numParticles) {
				attempts<-attempts+1

				newparticleList<-list(abcparticle(id=particle, generation=1, weight=0))
				newparticleList[[1]]<-initializeStatesFromMatrices(newparticleList[[1]], startingPriorsValues, startingPriorsFns, intrinsicPriorsValues, intrinsicPriorsFns, extrinsicPriorsValues, extrinsicPriorsFns)

				newparticleList[[1]]$distance<-abcDistance(summaryStatsLong(phy, convertTaxonFrameToGeigerData(doSimulation(splits, intrinsicFn, extrinsicFn, newparticleList[[1]]$startingValues, newparticleList[[1]]$intrinsicValues, newparticleList[[1]]$extrinsicValues, timeStep), phy), niter.brown=niter.brown.g, niter.lambda=niter.lambda.g, niter.delta=niter.delta.g, niter.OU=niter.OU.g, niter.white=niter.white.g)
                                                   , originalSummaryValues, pls.model.list)

				if (is.na(newparticleList[[1]]$distance)) {
					warning("newparticleList[[1]]$distance = NA, likely an underflow/overflow problem")
					newparticleList[[1]]$id <-  (-1)
					newparticleList[[1]]$weight<- 0
				}
				else if (is.na(toleranceVector[1])) {
					warning("toleranceVector[1] = NA")
					newparticleList[[1]]$id <- (-1)
					newparticleList[[1]]$weight <- 0
				}


				else if ((newparticleList[[1]]$distance) < toleranceVector[1]) {
					newparticleList[[1]]$id <- particle
					newparticleList[[1]]$weight <- 1/numParticles
					particleWeights[particle] <- 1/numParticles
					particle<-particle+1
					particleList<-append(particleList, newparticleList)
				}
				else {
					newparticleList[[1]]$id <- (-1)
					newparticleList[[1]]$weight <- 0
				}
				#while(sink.number()>0) {sink()}
				vectorForDataFrame<-c(1, attempts, newparticleList[[1]]$id, 0, newparticleList[[1]]$distance, newparticleList[[1]]$weight, newparticleList[[1]]$startingValues, newparticleList[[1]]$intrinsicValues, newparticleList[[1]]$extrinsicValues)
				#cat("\n\nlength of vectorForDataFrame = ", length(vectorForDataFrame), "\n", "length of startingValues = ", length(startingValues), "\nlength of intrinsicValues = ", length(intrinsicValues), "\nlength of extrinsicValues = ", length(extrinsicValues), "\ndistance = ", newparticleList[[1]]$distance, "\nweight = ", newparticleList[[1]]$weight, "\n", vectorForDataFrame, "\n")
				particleDataFrame<-rbind(particleDataFrame, vectorForDataFrame)
				cat(particle-1, attempts, floor(numParticles*attempts/particle), newparticleList[[1]]$startingValues, newparticleList[[1]]$intrinsicValues, newparticleList[[1]]$extrinsicValues, newparticleList[[1]]$distance, "\n")

		} #while (particle<=numParticles) bracket

			names(particleDataFrame)<-nameVector
			dataGenerationStep=1
			time<-proc.time()[[3]]-start.time
			time.per.gen<-time
			#rejects.gen.one<-(dim(subset(particleDataFrame, particleDataFrame$id<0))[1])/(dim(subset(particleDataFrame,))[1])
			#rejects<-c()

			for (i in 1:numberParametersTotal){
				param.stdev[1,i]<-c(sd(subset(particleDataFrame, particleDataFrame$id>0)[,6+i]))
				weightedMeanParam[1,i]<-weighted.mean(subset(particleDataFrame, particleDataFrame$id>0)[,6+i], subset(particleDataFrame, particleDataFrame$id>0)[,6])
				#c(mean(subset(particleDataFrame, X3>0)[,7:dim(particleDataFrame)[2]])/subset(particleDataFrame, X3>0)[,6])
			}


			save.image(file=paste("WS", jobName, ".Rdata", sep=""))
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

			save(prcResults, file=paste("partialResults", jobName, ".txt", sep=""))





					while (dataGenerationStep < nStepsPRC) {
						dataGenerationStep<-dataGenerationStep+1
						cat("\n\n\n", "STARTING DATA GENERATION STEP ", dataGenerationStep, "\n\n\n")
						start.time<-proc.time()[[3]]
						particleWeights<-particleWeights/(sum(particleWeights,na.rm=TRUE)) #normalize to one
						cat("particleWeights\n", particleWeights, "\n\n")
						oldParticleList<-particleList
						oldParticleWeights<-particleWeights
						particleWeights=rep(0, numParticles) #stores weights for each particle. Initially, assume infinite number of possible particles (so might not apply in discrete case)
						particleParameters<-matrix(nrow=numParticles, ncol=dim(startingPriorsValues)[2] +  dim(intrinsicPriorsValues)[2] + dim(extrinsicPriorsValues)[2]) #stores parameters in model for each particle
						particleDistance=rep(NA, numParticles)
						particle<-1
						attempts<-0
						cat("successes", "attempts", "expected number of attempts required\n")
						particleList<-list()
						weightScaling=0;
						while (particle<=numParticles) {
							attempts<-attempts+1
							particleToSelect<-which.max(as.vector(rmultinom(1, size = 1, prob=oldParticleWeights)))
							#cat("particle to select = ", particleToSelect, "\n")
							#cat("dput(oldParticleList)\n")
							#dput(oldParticleList)
							#cat("dput(oldParticleList[particleToSelect])\n")
							#dput(oldParticleList[particleToSelect])
							#cat("dput(oldParticleList[[particleToSelect]])\n")
							#dput(oldParticleList[[particleToSelect]])
							newparticleList<-list(oldParticleList[[particleToSelect]])
							#cat("dput(newparticleList[[1]])\n")
							#dput(newparticleList[[1]])
							#cat("mutateStates\n")

							newparticleList[[1]]<-mutateStates(newparticleList[[1]], startingPriorsValues, startingPriorsFns, intrinsicPriorsValues, intrinsicPriorsFns, extrinsicPriorsValues, extrinsicPriorsFns, standardDevFactor)
							#cat("dput(newparticleList[[1]]) AFTER MUTATE STATES\n")
							#dput(newparticleList[[1]])

							newparticleList[[1]]$distance<-abcDistance(summaryStatsLong(phy, convertTaxonFrameToGeigerData(doSimulation(splits, intrinsicFn, extrinsicFn, newparticleList[[1]]$startingValues, newparticleList[[1]]$intrinsicValues, newparticleList[[1]]$extrinsicValues, timeStep), phy), niter.brown=niter.brown.g, niter.lambda=niter.lambda.g, niter.delta=niter.delta.g, niter.OU=niter.OU.g, niter.white=niter.white.g)
							                                           , originalSummaryValues, pls.model.list)
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
							#cat("dput(newparticleList[[1]]) AFTER computeABCDistance\n")
							#dput(newparticleList[[1]])

							if (is.na(newparticleList[[1]]$distance)) {
								#cat("Error with Geiger?  newparticleList[[1]]$distance = NA\n")
								#while(sink.number()>0) {sink()}
								#warning("newparticleList[[1]]$distance = NA")
								newparticleList[[1]]$id <- (-1)
								newparticleList[[1]]$weight <- 0
							}
							else if (newparticleList[[1]]$distance < toleranceVector[dataGenerationStep]) {
								newparticleList[[1]]$id <- particle
								particle<-particle+1
								particleList<-append(particleList, newparticleList)
								#now get weights, using correction in Sisson et al. 2007
								newWeight=0
								for (i in 1:length(oldParticleList)) {
									lnTransitionProb=log(1)
									for (j in 1:length(newparticleList[[1]]$startingValues)) {
										newvalue<-newparticleList[[1]]$startingValues[j]
										meantouse= oldParticleList[[i]]$startingValues[j]
											if (startingPriorsFns[j]=="uniform") {
												sdtouse<-standardDevFactor*((max(startingPriorsValues[,j])-min(startingPriorsValues[,j]))/sqrt(12))
												#print(paste("startingPriorFn is uniform and sdtouse =", sdtouse))
											}
											else if (startingPriorsFns[j]=="exponential") {
												sdtouse<-standardDevFactor*(1/startingPriorsValues[,j])
												#print(paste("startingPriorFn is exponential and sdtouse =", sdtouse))
											}
											else {
												sdtouse<-standardDevFactor*(startingPriorsValues[2,j])
											}

										lnlocalTransitionProb=dnorm(newvalue, mean=meantouse, sd=sdtouse,log=TRUE) - ((log(1)/pnorm(min(startingPriorsValues[, j]), mean=meantouse, sd=sdtouse, lower.tail=T, log.p=T))* pnorm(max(startingPriorsValues[,j]), mean=meantouse , sd=sdtouse, lower.tail=F, log.p=T))

										if (lnlocalTransitionProb == "NaN") {  #to prevent lnlocalTransitionProb from being NaN (if pnorm=0)
											lnlocalTransitionProb<-.Machine$double.xmin
										}
										if (min(startingPriorsValues[, j])==max(startingPriorsValues[, j])) {
											lnlocalTransitionProb=log(1)
										}
										lnTransitionProb<-lnTransitionProb+lnlocalTransitionProb
										if(!is.finite(lnTransitionProb) || is.na(lnlocalTransitionProb)) {
											print(paste("issue with lnTransitionProb: lnlocalTransitionProb = ",lnlocalTransitionProb," lnTransitionProb = ",lnTransitionProb))
										}
									}
									for (j in 1:length(newparticleList[[1]]$intrinsicValues)) {
										newvalue<-newparticleList[[1]]$intrinsicValues[j]
										meantouse= oldParticleList[[i]]$intrinsicValues[j]
										if (intrinsicPriorsFns[j]=="uniform") {
											sdtouse<-standardDevFactor*((max(intrinsicPriorsValues[,j])-min(intrinsicPriorsValues[,j]))/sqrt(12))
										}
										else if (intrinsicPriorsFns[j]=="exponential") {
											sdtouse<-standardDevFactor*(1/intrinsicPriorsValues[,j])
										}
										else {
											sdtouse<-standardDevFactor*(intrinsicPriorsValues[2,j])
										}
										lnlocalTransitionProb=dnorm(newvalue, mean= meantouse, sd= sdtouse,log=TRUE)-((log(1)/pnorm(min(intrinsicPriorsValues[, j]), mean=meantouse , sd=sdtouse, lower.tail=T, log.p=T)) * pnorm(max(intrinsicPriorsValues[,j]), mean=meantouse , sd=sdtouse, lower.tail=F, log.p=T))

										if (lnlocalTransitionProb == "NaN") {  #to prevent lnlocalTransitionProb from being NaN (if pnorm=0)
											lnlocalTransitionProb<-.Machine$double.xmin
										}
										if (min(intrinsicPriorsValues[, j])==max(intrinsicPriorsValues[, j])) {
											lnlocalTransitionProb=log(1)
										}
										lnTransitionProb<-lnTransitionProb+lnlocalTransitionProb
										if(!is.finite(lnTransitionProb) || is.na(lnlocalTransitionProb)) {
											print(paste("issue with lnTransitionProb: lnlocalTransitionProb = ",lnlocalTransitionProb," lnTransitionProb = ",lnTransitionProb))
										}

									}
									for (j in 1:length(newparticleList[[1]]$extrinsicValues)) {
										newvalue<-newparticleList[[1]]$extrinsicValues[j]
										meantouse= oldParticleList[[i]]$extrinsicValues[j]
										if (extrinsicPriorsFns[j]=="uniform") {
											sdtouse<-standardDevFactor*((max(extrinsicPriorsValues[,j])-min(extrinsicPriorsValues[,j]))/sqrt(12))
										}
										else if (extrinsicPriorsFns[j]=="exponential") {
											sdtouse<-standardDevFactor*(1/extrinsicPriorsValues[,j])
										}
										else {
											sdtouse<-standardDevFactor*(extrinsicPriorsValues[2,j])
										}
										lnlocalTransitionProb=dnorm(newvalue, mean= meantouse, sd= sdtouse,log=TRUE)-((log(1)/pnorm(min(extrinsicPriorsValues[,j]), mean=meantouse , sd=sdtouse, lower.tail=T, log.p=T)) * pnorm(max(extrinsicPriorsValues[,j]), mean=meantouse , sd=sdtouse, lower.tail=F, log.p=T))
										if (lnlocalTransitionProb == "NaN") {  #to prevent lnlocalTransitionProb from being NaN (if pnorm=0)
											lnlocalTransitionProb<-.Machine$double.xmin
										}
										if (min(extrinsicPriorsValues[, j])==max(extrinsicPriorsValues[, j])) {
											lnlocalTransitionProb=log(1)
										}
										lnTransitionProb<-lnTransitionProb+lnlocalTransitionProb
										if(!is.finite(lnTransitionProb) || is.na(lnlocalTransitionProb)) {
											print(paste("issue with lnTransitionProb: lnlocalTransitionProb = ",lnlocalTransitionProb," lnTransitionProb = ",lnTransitionProb))
										}

									}
									newWeight<-newWeight+(oldParticleList[[i]]$weight)*exp(lnTransitionProb)
								} #for (i in 1:length(oldParticleList)) bracket

								if (!is.finite(newWeight)) {
									print(paste("warning: newWeight is ",newWeight))
								}
								newparticleList[[1]]$weight<- newWeight
								particleWeights[particle-1]<-newWeight
								weightScaling<-weightScaling+newWeight
							} #else if (newparticleList[[1]]$distance < toleranceVector[dataGenerationStep]) bracket
							else {
								newparticleList[[1]]$id<- (-1)
								newparticleList[[1]]$weight<-0
							}
							#while(sink.number()>0) {sink()}
							#print(newparticleList)
							vectorForDataFrame<-c(dataGenerationStep, attempts,newparticleList[[1]]$id, particleToSelect, newparticleList[[1]]$distance, newparticleList[[1]]$weight, newparticleList[[1]]$startingValues, newparticleList[[1]]$intrinsicValues, newparticleList[[1]]$extrinsicValues)
				#cat("\n\nlength of vectorForDataFrame = ", length(vectorForDataFrame), "\n", "length of startingValues = ", length(startingValues), "\nlength of intrinsicValues = ", length(intrinsicValues), "\nlength of extrinsicValues = ", length(extrinsicValues), "\ndistance = ", newparticleList[[1]]$distance, "\nweight = ", newparticleList[[1]]$weight, "\n", vectorForDataFrame, "\n")

							particleDataFrame<-rbind(particleDataFrame, vectorForDataFrame) #NOTE THAT WEIGHTS AREN'T NORMALIZED IN THIS DATAFRAME
							cat(particle-1, attempts, floor(numParticles*attempts/particle), newparticleList[[1]]$startingValues, newparticleList[[1]]$intrinsicValues, newparticleList[[1]]$extrinsicValues, newparticleList[[1]]$distance, "\n")

						} #while (particle<=numParticles) bracket



						particleDataFrame[which(particleDataFrame$generation==dataGenerationStep), ]$weight<-particleDataFrame[which(particleDataFrame$generation==dataGenerationStep), ]$weight/(sum(particleDataFrame[which(particleDataFrame$generation==dataGenerationStep), ]$weight))

						time2<-proc.time()[[3]]-start.time
						time.per.gen<-c(time.per.gen, time2)
						#rejects.per.gen<-(dim(subset(particleDataFrame, particleDataFrame$id<0))[1])/(dim(subset(particleDataFrame[which(particleDataFrame$generation==dataGenerationStep),],))[1])
						#rejects<-c(rejects, rejects.per.gen)
						sub1<-subset(particleDataFrame, particleDataFrame$generation==dataGenerationStep)
						sub2<-subset(sub1, sub1$id>0)

						for (i in 1:numberParametersTotal){
							param.stdev[dataGenerationStep,i]<-c(sd(sub2[,6+i]))
							weightedMeanParam[dataGenerationStep,i]<-weighted.mean(sub2[,6+i], sub2[,6])
						}

						if (stopRule){	#this will stop the PRC from running out to max number of generations if all params are below stopValue
							FF<-rep(1, dim(weightedMeanParam)[2])
							for (check.weightedMeanParam in 1:length(FF)){
								if (is.na(abs(weightedMeanParam[dataGenerationStep, check.weightedMeanParam]-weightedMeanParam[dataGenerationStep-1, check.weightedMeanParam])/mean(weightedMeanParam[dataGenerationStep, check.weightedMeanParam], weightedMeanParam[dataGenerationStep-1, check.weightedMeanParam]) <= stopValue) && mean(weightedMeanParam[dataGenerationStep, check.weightedMeanParam], weightedMeanParam[dataGenerationStep-1, check.weightedMeanParam]) == 0) {  #this && is here to make sure any NAs are from fixed params and not miscalculations.
									FF[check.weightedMeanParam]<-0
								}
								else if (abs(weightedMeanParam[dataGenerationStep, check.weightedMeanParam]-weightedMeanParam[dataGenerationStep-1, check.weightedMeanParam])/mean(weightedMeanParam[dataGenerationStep, check.weightedMeanParam], weightedMeanParam[dataGenerationStep-1, check.weightedMeanParam]) <= stopValue){
									FF[check.weightedMeanParam]<-0
								}
								#print(FF)
							}
							if (sum(FF)==0){
								cat("\n\n\nweightedMeanParam is < ", stopValue, "Analysis is being terminated at", dataGenerationStep, "instead of continuing to ", nStepsPRC, "\n\n\n")
								dataGenerationStep<-nStepsPRC
							}
						}

						save.image(file=paste("WS", jobName, ".Rdata", sep=""))
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

						save(prcResults, file=paste("partialResults", jobName, ".txt", sep=""))

					} #while (dataGenerationStep < nStepsPRC) bracket



				names(particleDataFrame)<-nameVector
				if(plot) {
					quartz()
					plot(x=c(min(intrinsicPriorsValues), max(intrinsicPriorsValues)), y=c(0, 1), type="n")
					for (i in 1:(length(toleranceVector)-1)) {
						graycolor<-gray(0.5*(length(toleranceVector)-i)/length(toleranceVector))
						lines(density(subset(particleDataFrame, particleDataFrame$generation==i)[, 8]), col= graycolor)
					}
					lines(density(subset(particleDataFrame, particleDataFrame$generation==length(toleranceVector))[, 8]), col= "red")
				}



			#---------------------- ABC-PRC (End) --------------------------------


		input.data<-rbind(jobName, length(phy[[3]]), nrepSim, TreeYears, epsilonProportion, epsilonMultiplier, nStepsPRC, numParticles, standardDevFactor)

		time3<-proc.time()[[3]]
		genTimes<-c(time.per.gen, time3)

		prcResults<-vector("list")
		prcResults$input.data<-input.data
		prcResults$PriorMatrix<-PriorMatrix
		prcResults$particleDataFrame<-particleDataFrame
		prcResults$toleranceVector<-toleranceVector
		prcResults$phy<-phy
		prcResults$traits<-traits
		prcResults$simTime<-simTime
		prcResults$time.per.gen<-genTimes
		prcResults$CredInt <-CredInt(particleDataFrame)
		prcResults$HPD <-HPD(particleDataFrame)


	registerDoMC(1) #set number of cores back to 1
	print(prcResults)

}

	#------------------ ABC-PRC (end) ------------------

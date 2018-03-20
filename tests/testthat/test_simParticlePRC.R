test_that("simParticlePRC runs correctly", {
	data(simRunExample)
	set.seed(1)
	results <- TreEvo:::simParticlePRC<-function(
		,phy
		,taxonDF
		, timeStep 
		,intrinsicFn
		, extrinsicFn 
		,startingPriorsValues
		,intrinsicPriorsValues
		,extrinsicPriorsValues
		,startingPriorsFns
		,intrinsicPriorsFns
		,extrinsicPriorsFns
		,originalSummaryValues 
		.pls.model.list
		,toleranceValue
		,prevGenParticleList
		,standardDevFactor
		,numParticles
		)
	expect_is(results, "list")


set.seed(1)
library(TreEvo)
data(simRunExample)

phy = simPhy
traits = simChar
intrinsicFn=brownianIntrinsic
extrinsicFn=nullExtrinsic
startingPriorsFns="normal"
startingPriorsValues=matrix(c(mean(simChar[,1]), sd(simChar[,1])))
intrinsicPriorsFns=c("exponential")
intrinsicPriorsValues=matrix(c(10, 10), nrow=2, byrow=FALSE)
extrinsicPriorsFns=c("fixed")
extrinsicPriorsValues=matrix(c(0, 0), nrow=2, byrow=FALSE)
generation.time=1000000
standardDevFactor=0.2
StartSims=10
epsilonProportion=0.7
epsilonMultiplier=0.7
nStepsPRC=2
numParticles=3
jobName="exampleRun"
stopRule=FALSE
multicore=FALSE
coreLimit=1
TreeYears=max(branching.times(phy)) * 1e6








taxonDF <- TreEvo:::getTaxonDFWithPossibleExtinction(phy)
timeStep<-generation.time/TreeYears
edgesRescaled<-phy$edge.length/max(node.depth.edgelength(phy))
totalGenerations<-sum(sapply(edgesRescaled,function(x) floor(x/timeStep)))
saveData<-FALSE
	validation="CV" 
	scale=TRUE
	variance.cutoff=95
	#niter.goal=5
	#

################
	freevector<-TreEvo:::getFreeVector(startingPriorsFns=startingPriorsFns, startingPriorsValues=startingPriorsValues,
						intrinsicPriorsFns=intrinsicPriorsFns, intrinsicPriorsValues=intrinsicPriorsValues,
						extrinsicPriorsFns=extrinsicPriorsFns, extrinsicPriorsValues=extrinsicPriorsValues)
	numberParametersTotal<-length(freevector)
	numberParametersFree<-sum(freevector)
	#
	#create PriorMatrix
	namesForPriorMatrix<-c()
	PriorMatrix<-matrix(c(startingPriorsFns, intrinsicPriorsFns, extrinsicPriorsFns), nrow=1, ncol=numberParametersTotal)
	for (a in 1:dim(startingPriorsValues)[2]){
		namesForPriorMatrix<-c(paste0("StartingStates", a, sep=""))
		}
	#
	for (b in 1:dim(intrinsicPriorsValues)[2]){
		namesForPriorMatrix<-append(namesForPriorMatrix, paste0("IntrinsicValue", b, sep=""))
		}
	#
	#message(extrinsicPriorsValues)
	for (c in 1:dim(extrinsicPriorsValues)[2]){
		namesForPriorMatrix <-append(namesForPriorMatrix, paste0("ExtrinsicValue", c, sep=""))
		}
	#
	PriorMatrix<-rbind(
		PriorMatrix, 
		cbind(startingPriorsValues, 
			intrinsicPriorsValues, 
			extrinsicPriorsValues
			)
		)
	colnames(PriorMatrix)<-namesForPriorMatrix
	rownames(PriorMatrix)<-c("shape", "value1", "value2")
	#	
	#initialize weighted mean sd matrices
	weightedMeanParam<-matrix(nrow=nStepsPRC, ncol=numberParametersTotal)
	colnames(weightedMeanParam)<-namesForPriorMatrix
	rownames(weightedMeanParam)<-paste0("Gen ", c(1: nStepsPRC), sep="")
	param.stdev<-matrix(nrow=nStepsPRC, ncol=numberParametersTotal)
	colnames(param.stdev)<-namesForPriorMatrix
	rownames(param.stdev)<-paste0("Gen ", c(1: nStepsPRC), sep="")
	#
	#
	if (is.na(StartSims)) {
		StartSims<-1000*numberParametersFree
		}
	#
	# save input data for use later
	input.data<-rbind(jobName, Ntip(phy), StartSims, generation.time, TreeYears, timeStep, totalGenerations,
		epsilonProportion, epsilonMultiplier, nStepsPRC, numParticles, standardDevFactor)
	#
	# get summary values for observed data
	originalSummaryValues<-TreEvo:::summaryStatsLong(
		phy=phy, 
		traits=traits
		#niter.brown=200, niter.lambda=200, niter.delta=200, niter.OU=200, niter.white=200
		)
	#	
	# INITIAL SIMULATIONS
	initialSimsRes<-TreEvo:::initialSimsPRC(
		nrepSim=StartSims, 
		phy=phy,  
		taxonDF=taxonDF,
		startingPriorsValues=startingPriorsValues, 
		intrinsicPriorsValues=intrinsicPriorsValues, 
		extrinsicPriorsValues=extrinsicPriorsValues,
		startingPriorsFns=startingPriorsFns, 
		intrinsicPriorsFns=intrinsicPriorsFns, 
		extrinsicPriorsFns=extrinsicPriorsFns,
		freevector=freevector, 
		timeStep=timeStep, 
		intrinsicFn=intrinsicFn, 
		extrinsicFn=extrinsicFn, 
		nStepsPRC=nStepsPRC,
		coreLimit=coreLimit, 
		multicore=multicore,
		numberParametersFree=numberParametersFree,
		originalSummaryValues=originalSummaryValues,
		epsilonProportion=epsilonProportion,
		epsilonMultiplier=epsilonMultiplier,
		validation=validation,
		variance.cutoff=variance.cutoff,
		saveData=saveData,
		jobName=jobName
		)		
	#
	# pull the PLS model list out
	pls.model.list<-initialSimsRes$pls.model.list

	
###################

prevGenParticleList=NULL
numParticles=3
dataGenerationStep=1


	toleranceValue<-initialSimsRes$toleranceVector[dataGenerationStep]


# 
# get particle parameters
#
# can test for dataGenerationStep=1 if prevGenParticleList is null
if(is.null(prevGenParticleList)){
	# i.e. if dataGenerationStep == 1
		# get new particles from priors
	# for generation = 1, weight all particles the same
	newparticle<-TreEvo:::abcparticle(id=NA, generation=1, weight=1/numParticles)	#, parentid=NA
	newparticle<-TreEvo:::initializeStatesFromMatrices(
		particle=newparticle,
		startingPriorsValues=startingPriorsValues,
		startingPriorsFns=startingPriorsFns,
		intrinsicPriorsValues=intrinsicPriorsValues,
		intrinsicPriorsFns=intrinsicPriorsFns,
		extrinsicPriorsValues=extrinsicPriorsValues,
		extrinsicPriorsFns=extrinsicPriorsFns)							
}else{
	prevGenParticleWeights<-sapply(prevGenParticleList,function(x) x$weight)
	# use particles from PREVIOUS GENERATION to randomly select a particle
	particleToSelect<-which.max(as.vector(rmultinom(1, size = 1, prob=prevGenParticleWeights)))
	# get that particle's data
	newparticle<-prevGenParticleList[[particleToSelect]]
	#
	#
	newparticle<-TreEvo:::mutateStates(particle=newparticle,
		startingPriorsValues=startingPriorsValues, startingPriorsFns=startingPriorsFns,
		intrinsicPriorsValues=intrinsicPriorsValues, intrinsicPriorsFns=intrinsicPriorsFns,
		extrinsicPriorsValues=extrinsicPriorsValues, extrinsicPriorsFns=extrinsicPriorsFns,
		standardDevFactor=standardDevFactor
		)
	newparticle$parentid<-particleToSelect
	}
#
# do the simulation
simTraitsParticle<-TreEvo:::doSimulationInternal(
	taxonDF=taxonDF,
	intrinsicFn=intrinsicFn,
	extrinsicFn=extrinsicFn,
	startingValues=newparticle$startingValues,
	intrinsicValues=newparticle$intrinsicValues,
	extrinsicValues=newparticle$extrinsicValues,
	timeStep=timeStep
	)
#
# get the summary stats	
simSumMat<-TreEvo:::summaryStatsLong(phy=phy,traits=simTraitsParticle)
# get the distance of the simulation to the original
simDistance<-TreEvo:::abcDistance(summaryValuesMatrix=simSumMat,
	originalSummaryValues=originalSummaryValues, 
	pls.model.list=pls.model.list)	
#
# get the weights, if it passes the tolerance
if ((simDistance) < toleranceValue) {
	# record the distance
	newparticle$distance<-simDistance
		#ID doesn't need to be set - that's just the ID of the particle...
	#newparticle$id<-particle
	#particle<-particle+1
	#
	if(!is.null(prevGenParticleList)){
		#for generation>1 - now get weights, using correction in Sisson et al. 2007
		newparticle$weight<-TreEvo:::sumLogTranProb(
			 prevGenParticleList=prevGenParticleList
			,newStartingValues = newparticle$startingValues
			,newIntrinsicValues = newparticle$intrinsicValues
			,newExtrinsicValues = newparticle$extrinsicValues
			,startingPriorsFns=startingPriorsFns
			,intrinsicPriorsFns=intrinsicPriorsFns
			,extrinsicPriorsFns=extrinsicPriorsFns
			,startingPriorsValues=startingPriorsValues
			,intrinsicPriorsValues=intrinsicPriorsValues
			,extrinsicPriorsValues=extrinsicPriorsValues
			,standardDevFactor=standardDevFactor
			)
		}
}else{
	# particle didn't pass, discard it
	newparticle<-NA
	}	
#newparticle<-list(newparticle)
newparticle
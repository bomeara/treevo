# internal function for doing initial simulations in doRun_PRC
initialSimsPRC<-function(
		nrepSim,
		phy,  
		taxonDF,
		startingPriorsValues, intrinsicPriorsValues, extrinsicPriorsValues,
		startingPriorsFns, intrinsicPriorsFns, extrinsicPriorsFns,
		freevector, timeStep, 
		intrinsicFn, extrinsicFn, 
		nStepsPRC,
		coreLimit, 
		multicore,
		originalSummaryValues,
		epsilonProportion, 
		validation="CV", 
		saveData
		){		
	#
	#	
	#---------------------- Initial Simulations (Start) ------------------------------
	# See Wegmann et al. Efficient Approximate Bayesian Computation Coupled With Markov Chain Monte Carlo Without Likelihood.
		# Genetics (2009) vol. 182 (4) pp. 1207-1218 for more on the method.
	# We are doing pls and scaling built into pls. Unlike Wegmann et al., we are doing PLS for each parameter separately.
	# Otherwise, the PLS tends to optimize for just one parameter, and estimates for the less-favored one are quite bad
	# because the summary stats tend to be used for the other.
	#
	#Used to be = StartSims*((2^try)/2), If initial simulations are not enough, and we need to try again then new analysis will double number of initial simulations
	#nrepSim<-StartSims
	#
	message(paste0("Number of initial simulations set to ", nrepSim)) #, "\n"
	message("Doing initial simulations...")
	Time<-proc.time()[[3]]
	trueFreeValues<-matrix(nrow=0, ncol= numberParametersFree)
	#set up initial sum stats as length of SSL of original data
	summaryValues<-matrix(
		nrow=0, 
		ncol=length(originalSummaryValues)
		)
	trueFreeValuesANDSummaryValues<-parallelSimulateWithPriors(	#makeQuiet(
		nrepSim=StartSims, 
		phy=phy,  
		taxonDF=taxonDF,
		startingPriorsValues=startingPriorsValues, intrinsicPriorsValues=intrinsicPriorsValues, extrinsicPriorsValues=extrinsicPriorsValues,
		startingPriorsFns=startingPriorsFns, intrinsicPriorsFns=intrinsicPriorsFns, extrinsicPriorsFns=extrinsicPriorsFns,
		freevector=freevector, timeStep=timeStep, 
		intrinsicFn=intrinsicFn, extrinsicFn=extrinsicFn, 
		coreLimit=coreLimit, 
		multicore=multicore,
		verbose=FALSE
		#niter.brown=niter.brown.g, niter.lambda=niter.lambda.g, niter.delta=niter.delta.g, niter.OU=niter.OU.g, niter.white=niter.white.g
		)	#)
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
	pls.model.list <- ( #makeQuiet(
		apply(trueFreeValuesMatrix, 2, 
			returnPLSModel, 
			summaryValuesMatrix=summaryValuesMatrix, 
			validation=validation,
			scale=scale, variance.cutoff = variance.cutoff
			))
		#)
		#
	#----------------- Find distribution of distances (Start) ----------------------
	#
	distanceVector<-abcDistance(summaryValuesMatrix, originalSummaryValues, pls.model.list)
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
	# need to provide simTime, pls.model.list, toleranceVector
	res<-list(pls.model.list=pls.model.list, toleranceVector=toleranceVector, simTime=simTime)
	#
	return(res)
	}
	
	

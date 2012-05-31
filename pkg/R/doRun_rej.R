
##TreeYears = 1000 if tree is in thousands of years


doRun_rej<-function(phy, traits, intrinsicFn, extrinsicFn, summaryFns=c(rawValuesSummaryStats, geigerUnivariateSummaryStats2), startingPriorsValues, startingPriorsFns, intrinsicPriorsValues, intrinsicPriorsFns, extrinsicPriorsValues, extrinsicPriorsFns, startingStatesGuess=c(), intrinsicStatesGuess=c(), extrinsicStatesGuess=c(), TreeYears=1e+04, numParticles=1000, standardDevFactor=0.20, StartSims=NA, plot=FALSE, vipthresh=0.8, epsilonProportion=0.7, epsilonMultiplier=0.7, jobName=NA, debug=FALSE, trueStartingState=NA, trueIntrinsicState=NA, whenToKill=20, abcMethod="rejection", abcTolerance=0.1, multicore=TRUE, coreLimit=NA, filenames=c("rejectionsims.RData")) {

	if (!is.binary.tree(phy)) {
		print("Warning: Tree is not fully dichotomous")
	}
	startTime<-proc.time()[[3]]
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
	if (length(startingStatesGuess)==0) { #if no user guesses, try pulling a value from the prior
		startingStatesGuess<-rep(NA,length(startingPriorsFns))
		for (i in 1:length(startingPriorsFns)) {
			startingStatesGuess[i]<-pullFromPrior(startingPriorsValues[,i],startingPriorsFns[i])
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
	input.data<-rbind(jobName, length(phy[[3]]), nrepSim, timeStep, epsilonProportion, epsilonMultiplier,  numParticles, standardDevFactor, trueStartingState, trueIntrinsicState)		
	cat(paste("Number of initial simulations set to", nrepSim, "\n"))
	
	trueFreeValues<-matrix(nrow=0, ncol= numberParametersFree)
	summaryValues<-matrix(nrow=0, ncol=22+dim(traits)[1]) #there are 22 summary statistics possible, plus the raw data
	#do nrepSim using the inputs. If multicore, do that way
	parallelSimulation(nrepSim, coreLimit, startingPriorsValues, intrinsicPriorsValues, extrinsicPriorsValues, startingPriorsFns, intrinsicPriorsFns, extrinsicPriorsFns, trueFreeValues, freevector, timeStep, intrinsicFn, extrinsicFn, multicore, filename=filenames[1])
	cat("\n\n")
	#pull in all the simulation results: true values and the summary values
	trueFreeValuesANDSummaryValues<-loadSimulations(filenames)
	trueFreeValues<-trueFreeValuesANDSummaryValues[,1:numberParametersFree]
	summaryValues<-trueFreeValuesANDSummaryValues[,-1:-numberParametersFree]
	results<-summarizeRejection(summaryValues, trueFreeValues, vipthresh, traits, todo, abcMethod, abcTolerance, jobName)
	return(results)
}

summarizeRejection<-function(summaryValues, trueFreeValues, vipthresh, traits, todo, abcMethod, abcTolerance, jobName){
	#get boxcoxLambda and boxcoxAddition
	library(abc)
	boxcoxEstimates<-boxcoxEstimation(summaryValues)
	boxcoxAddition<-boxcoxEstimates$boxcoxAddition
	boxcoxLambda<-boxcoxEstimates$boxcoxLambda
	boxcoxSummaryValues<-boxcoxEstimates$boxcoxSummaryValues
	#get pls parameters
	save(trueFreeValues, file="tFV")
	save(boxcoxSummaryValues, file="bcSV")
	plsEstimates<-plsEstimation(trueFreeValues, boxcoxSummaryValues, vipthresh)
	prunedPlsResult<-plsEstimates$prunedPlsResult
	prunedSummaryValues<-plsEstimates$prunedSummaryValues
	todo<-plsEstimates$todo
	#do abc() using transformed summary stat for observed data and set of transformed summary stats and true values for simulated data
	originalSummaryStats<-summaryStatsLong(phy=phy, data=traits, todo=todo, jobName=jobName)[which(todo==1)]
	abcResults<-abc(target=originalSummaryStats, param=trueFreeValues, sumstat=prunedSummaryValues, tol=abcTolerance, method=abcMethod)
	return(abcResults)
}

parallelSimulation<-function(nrepSim, coreLimit, startingPriorsValues, intrinsicPriorsValues, extrinsicPriorsValues, startingPriorsFns, intrinsicPriorsFns, extrinsicPriorsFns, trueFreeValues, freevector, timeStep, intrinsicFn, extrinsicFn, multicore, filename) {
	library(doMC, quietly=T)
	library(foreach, quietly=T)
	cores=1
	if (multicore) {
		if (is.na(coreLimit)){
			registerDoMC()
			getDoParWorkers()->cores
		}
		else {
			registerDoMC(coreLimit)
			coreLimit->cores
		}
	}
	cat(paste("Using", cores, "core(s) for initial simulations \n\n"))

	trueFreeValuesANDSummaryValues<-foreach(1:nrepSim, .combine=rbind) %dopar% simulateData(startingPriorsValues, intrinsicPriorsValues, extrinsicPriorsValues, startingPriorsFns, intrinsicPriorsFns, extrinsicPriorsFns, trueFreeValues, freevector, timeStep, intrinsicFn, extrinsicFn)
	save(trueFreeValuesANDSummaryValues,file=filename,compress=TRUE)
}

loadSimulations<-function(filenames) {
	results<-c() #check type
	for (i in sequence(length(filenames))) {
		load(filenames[i]) 
		results<-rbind(results,trueFreeValuesANDSummaryValues)
	}
	return(results)
}

boxcoxEstimation<-function(summaryValues){
	library("car", quietly=T)
	boxcoxLambda<-rep(NA, dim(summaryValues)[2])
	boxcoxAddition<-rep(NA, dim(summaryValues)[2])
	for (summaryValueIndex in 1:dim(summaryValues)[2]) {
		boxcoxAddition[summaryValueIndex]<-0
		lowValue<-min(summaryValues[, summaryValueIndex])-4*sd(summaryValues[, summaryValueIndex])
		if (lowValue<=0) {
			boxcoxAddition[summaryValueIndex]<-4*abs(lowValue) #just for some protection against low values, since box.cox needs non-negative values
		}
		summary<-summaryValues[, summaryValueIndex]+boxcoxAddition[summaryValueIndex]
		boxcoxLambda[summaryValueIndex]<-1
		if(sd(summaryValues[, summaryValueIndex])>0) { #box.cox fails if all values are identical
			newLambda<-as.numeric(try(powerTransform(summary,method="Nelder-Mead")$lambda)) #new car uses powerTransform instead of box.cox.powers
			if (!is.na(newLambda)) {
				boxcoxLambda[summaryValueIndex]<-newLambda
			}
		}
		summaryValues[, summaryValueIndex]<-summary^boxcoxLambda[summaryValueIndex]
	}
	return(list(boxcoxAddition=boxcoxAddition,boxcoxLambda=boxcoxLambda,boxcoxSummaryValues=summaryValues))
}


plsEstimation<-function(trueFreeValues, boxcoxSummaryValues, vipthresh) {
#----------------- Find best set of summary stats to use for this problem. (Start) -----------------
			#Use mixOmics to to find the optimal set of summary stats. Store this info in the todo vector. Note that this uses a different package (mixOmics rather than pls than that used by Weggman et al. because this package can calculate variable importance in projection and deals fine with NAs)
	library("mixOmics")
	plsResult<-pls(Y=trueFreeValues, X=boxcoxSummaryValues)
	vipResult<-vip(plsResult)
	todo<-rep(1, dim(boxcoxSummaryValues)[2]) #initialize the vector that indicates which summary stats to include
			
	summaryIndexOffset=0 #since R excludes invariant columns from regression, this offests so we don't try to extract from these columns
	nearZeroVarVector<-mixOmics:::nearZeroVar(boxcoxSummaryValues)
	nearZeroVarVector<-nearZeroVarVector$Position
	#print(nearZeroVarVector)
	for (summaryIndex in 1:dim(boxcoxSummaryValues)[2]) {
		if (summaryIndex %in% nearZeroVarVector) {
			summaryIndexOffset=summaryIndexOffset+1
			todo[summaryIndex]<-0 #exclude this summary stat because it lacks variation
		}	
		else if (max(vipResult[summaryIndex-summaryIndexOffset, ]) < vipthresh) {
			todo[summaryIndex]<-0 #exclude this summary stat, because it is too unimportant
		}	
	}
	prunedSummaryValues<-boxcoxSummaryValues[, which(todo>0)]
	prunedPlsResult<-pls(Y=trueFreeValues, X=prunedSummaryValues)
	return(list(prunedSummaryValues=prunedSummaryValues,prunedPlsResult=prunedPlsResult,todo=todo))			
}

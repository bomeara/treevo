
##This seems to be working if partialResults does not exist.  If checkpoint=TRUE, then run fails.  


#the doRun function takes input from the user and then automatically guesses optimal parameters, though user overriding is also possible.
#the guesses are used to do simulations near the expected region. If omitted, they are set to the midpoint of the input parameter matrices

doRun_abc<-function(phy, traits, intrinsicFn, extrinsicFn, summaryFns=c(rawValuesSummaryStats, geigerUnivariateSummaryStats2), startingPriorsValues, startingPriorsFns, intrinsicPriorsValues, intrinsicPriorsFns, extrinsicPriorsValues, extrinsicPriorsFns, startingStatesGuess=c(), intrinsicStatesGuess=c(), extrinsicStatesGuess=c(), timeStep, toleranceVector=c(), numParticles=1000, standardDevFactor=0.05, StartSims=100, plot=FALSE, vipthresh=0.8, epsilonProportion=0.2, epsilonMultiplier=0.5, nStepsPRC=4, maxTries=1, jobName=NA, debug=TRUE, trueStartingState=NA, trueIntrinsicState=NA, startFromCheckpoint=TRUE, whenToKill=20, checkpointSave=TRUE, stopRule=TRUE, stopValue=0.05, abcMethod=NA) {
print(abcMethod)
if (debug){
	cat("\nDebugging doRun\n")
	dput(doRun)
}

##Checkpoint saving stuff:
paste("partialResults", jobName, ".txt*", sep="")->pRname
paste("WS", jobName, ".Rdata", sep="")->WSname

system(command=paste("ls ", pRname, " | grep -c ", pRname, sep=""), intern=TRUE) -> filecount

if (filecount=="0") {  #if file is absent 
	startFromCheckpoint=FALSE
	dataGenerationStep=0
	cat ("\nstart from checkpoint =", startFromCheckpoint, "\n")
	cat("dataGenerationStep=", dataGenerationStep, "\n")
	rejects<-c()	

}

if (filecount=="1"){  #if file is present 
	paste("partialResults", jobName, ".txt", sep="")->pRname
	paste(load(pRname))
	dataGenerationStep <- max(test$particleDataFrame$X1)
		if (dataGenerationStep==nStepsPRC){
			cat ("\n\nRun was finished already\n\n")		}
	#paste(load(WSname))
	cat ("\nstart from checkpoint =", startFromCheckpoint, "\n")
	cat("dataGenerationStep=", dataGenerationStep, "\n")
	nameVector<-c("generation", "attempt", "id", "parentid", "distance", "weight")
	run.goingwell=TRUE
	input.data<-test$input.data
	boxcox.output<-test$boxcoxLambda
	boxcoxLambda<-test$boxcox.output[[1]]
	boxcoxAddition<-test$boxcox.output[[2]]
	prunedPlsResult<-test$boxcox.output[[3]]
	prunedSummaryValues<-test$boxcox.output[[4]]
	originalSummaryStats<-test$boxcox.output[[5]]
	particleDataFrame<-test$particleDataFrame
	epsilonDistance<-test$epsilonDistance
	
	toleranceVector<-test$toleranceVector
		if (length(toleranceVector) < nStepsPRC){
			#print(toleranceVector)
			toleranceVector<-rep(epsilonDistance, nStepsPRC)
			for (step in 2:nStepsPRC) {
				toleranceVector[step]<-toleranceVector[step-1]*as.numeric(input.data[11])
			}
			#print(toleranceVector)
		}	
	todo<-test$todo
	phy<-test$phy
		splits<-getSimulationSplits(phy)
	traits<-test$traits
	rejects.gen.one<-test$rejects.gen.one
	rejects<-test$rejects
	particleWeights<-test$particleWeights
	particleVector<-test$particleVector
	numberParametersFree<-test$numberParametersFree
	param.stdev<-test$param.stdev
	weightedMeanParam<-test$weightedMeanParam
	time.per.gen<-test$time.per.gen
}


run.goingwell=FALSE
for (try in 1:maxTries)	{
		nrepSim<-StartSims*((10^try)/10)
		cat("\n\n****  TRY", try, "of", maxTries, " ****\n\n")
		

		
while (!run.goingwell) {
		run.goingwell=TRUE
	if (startFromCheckpoint==FALSE) {


		#run.finished=FALSE
		
		splits<-getSimulationSplits(phy) #initialize this info

	
	input.data<-rbind(jobName, length(phy[[3]]), nrepSim, epsilonProportion, epsilonMultiplier, nStepsPRC, numParticles, standardDevFactor, try, trueStartingState, trueIntrinsicState)


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
	#print(extrinsicPriorsValues)
	for (c in 1:dim(extrinsicPriorsValues)[2]) {
		namesForPriorMatrix <-append(namesForPriorMatrix, paste("ExtrinsicValue", c, sep=""))
	}
	PriorMatrix<-rbind(PriorMatrix, cbind(startingPriorsValues, intrinsicPriorsValues, extrinsicPriorsValues))
	#PriorMatrix<-rbind(PriorMatrix1, PriorMatrix2)
	colnames(PriorMatrix)<-namesForPriorMatrix
	rownames(PriorMatrix)<-c("shape", "value1", "value2")
		
	#print(PriorMatrix)
		
		
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
		#print(numberParametersStarting)
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
		#print(numberParametersIntrinsic)
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
		#print(numberParametersExtrinsic)
	}
	

	param.stdev<-matrix(nrow=nStepsPRC, ncol=numberParametersTotal)
	colnames(param.stdev)<-namesForPriorMatrix
	rownames(param.stdev)<-paste("Gen", c(1: nStepsPRC), sep="")
	weightedMeanParam<-matrix(nrow=nStepsPRC, ncol=numberParametersTotal)
	colnames(weightedMeanParam)<-namesForPriorMatrix
	rownames(weightedMeanParam)<-paste("Gen", c(1: nStepsPRC), sep="")
	#names(param.stdev)<-c("Generation", )
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
	
	
	#----------------- Find best set of summary stats to use for this problem. (start) -----------------
	#See Wegmann et al. Efficient Approximate Bayesian Computation Coupled With Markov Chain Monte Carlo Without Likelihood. Genetics (2009) vol. 182 (4) pp. 1207-1218 for more on the method
	trueFreeValues<-matrix(nrow=0, ncol= numberParametersFree)
	summaryValues<-matrix(nrow=0, ncol=22+dim(traits)[1]) #there are 22 summary statistics possible, plus the raw data
	#Rprof(nrepSims.time.check<-tempfile())	
	for (simIndex in 1:nrepSim) {
		cat("Now doing simulation rep ",simIndex," of ",nrepSim,"\n",sep="")
		trueStarting<-rep(NaN, dim(startingPriorsValues)[2])
		trueIntrinsic<-rep(NaN, dim(intrinsicPriorsValues)[2])
		trueExtrinsic<-rep(NaN, dim(extrinsicPriorsValues)[2])
		for (j in 1:dim(startingPriorsValues)[2]) {
			trueStarting[j]=pullFromPrior(startingPriorsValues[,j],startingPriorsFns[j])
		}
		for (j in 1:dim(intrinsicPriorsValues)[2]) {
			trueIntrinsic[j]=pullFromPrior(intrinsicPriorsValues[,j],intrinsicPriorsFns[j])
		}
		for (j in 1:dim(extrinsicPriorsValues)[2]) {
			trueExtrinsic[j]=pullFromPrior(extrinsicPriorsValues[,j],extrinsicPriorsFns[j])
		}
		trueInitial<-c(trueStarting, trueIntrinsic, trueExtrinsic)
		trueFreeValues<-rbind(trueFreeValues, trueInitial[freevector])
		#cat("summaryValues\n")
		#print(summaryValues)
		convertTaxonFrameToGeigerData (doSimulation(splits=splits, intrinsicFn= intrinsicFn, extrinsicFn= extrinsicFn, startingStates= trueStarting, intrinsicValues= trueIntrinsic, extrinsicValues= trueExtrinsic, timeStep=timeStep), phy)->simdata
		#save(simdata, file="simdata.Rdata")
		summaryValues<-rbind(summaryValues, summaryStatsLong(phy, simdata, jobName=jobName))
		#print(summaryValues)
		while(sink.number()>0) {sink()}
	}
	
	library(abc)
	
	abcResults<-abc(target=summaryStatsLong(phy, traits, jobName=jobName), param=trueFreeValues, sumstat=summaryValues, tol=0.1, method=abcMethod)
	
	save(abcResults,trueFreeValues,summaryValues,file=paste("abc_pkg_rej_results_",jobName,".Rdata",sep=""))
	
}  #if (startFromCheckpoint==FALSE)
}  #while (!run.goingwell)
}  #for (try in 1:maxTries)
}  #doRun_abc
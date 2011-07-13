
##This seems to be working if partialResults does not exist.  If checkpoint=TRUE, then run fails.  


#the doRun function takes input from the user and then automatically guesses optimal parameters, though user overriding is also possible.
#the guesses are used to do simulations near the expected region. If omitted, they are set to the midpoint of the input parameter matrices

doRun<-function(phy, traits, intrinsicFn, extrinsicFn, summaryFns=c(rawValuesSummaryStats, geigerUnivariateSummaryStats2), startingPriorsValues, startingPriorsFns, intrinsicPriorsValues, intrinsicPriorsFns, extrinsicPriorsValues, extrinsicPriorsFns, startingStatesGuess=c(), intrinsicStatesGuess=c(), extrinsicStatesGuess=c(), timeStep, toleranceVector=c(), numParticles=1000, standardDevFactor=0.05, StartSims=100, plot=FALSE, vipthresh=0.8, epsilonProportion=0.2, epsilonMultiplier=0.5, nStepsPRC=4, maxTries=1, jobName=NA, debug=TRUE, trueStartingState=NA, trueIntrinsicState=NA, startFromCheckpoint=TRUE, whenToKill=20, checkpointSave=TRUE, stopRule=TRUE, stopValue=0.05, multicore=TRUE, coreLimit=NA) {

if (!is.binary.tree(phy)) {
	print("Warning: Tree is not fully dichotomous")
}

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
		cat ("\n\nRun was finished already\n\n")
	}
	#paste(load(WSname))
	cat ("\nstart from checkpoint =", startFromCheckpoint, "\n")
	cat("dataGenerationStep=", dataGenerationStep, "\n")
	nameVector<-c("generation", "attempt", "id", "parentid", "distance", "weight")
	run.goingwell=TRUE
	input.data<-test$input.data
	boxcox.output<-test$boxcoxLambda
	boxcoxLambda<-test$boxcox.output$lambda
	boxcoxAddition<-test$boxcox.output$addition
	prunedPlsResult<-test$boxcox.output$PlsResult
	prunedSummaryValues<-test$boxcox.output$prunedSummaryValues
	originalSummaryStats<-test$boxcox.output$originalSummaryStats
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

cores=1
if (multicore) {
	library(doMC, quietly=T)
	library(foreach, quietly=T)
	if (is.na(coreLimit)){
		registerDoMC()
		getDoParWorkers()->cores
	}
	else {
		registerDoMC(coreLimit)
		coreLimit->cores
	}
}
cat(paste("Using", cores, "core(s) for the analysis\n"))
			
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
		
			numberParametersTotal<-dim(startingPriorsValues)[2] +  dim(intrinsicPriorsValues)[2] + dim(extrinsicPriorsValues)[2] #figure out number of free params
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
			
			
			#---------------------- Initial Simulations (Start) ------------------------------
			#See Wegmann et al. Efficient Approximate Bayesian Computation Coupled With Markov Chain Monte Carlo Without Likelihood. Genetics (2009) vol. 182 (4) pp. 1207-1218 for more on the method
			Time<-proc.time()[[3]]
			trueFreeValues<-matrix(nrow=0, ncol= numberParametersFree)
			summaryValues<-matrix(nrow=0, ncol=22+dim(traits)[1]) #there are 22 summary statistics possible, plus the raw data
			trueFreeValuesANDSummaryValues<-foreach(1:nrepSim, .combine=rbind) %dopar% simulateData(nrepSim, startingPriorsValues, intrinsicPriorsValues, extrinsicPriorsValues, startingPriorsFns, intrinsicPriorsFns, extrinsicPriorsFns, trueFreeValues, freevector, timeStep, intrinsicFn, extrinsicFn)
			#print(trueFreeValuesANDSummaryValues)
			
			trueFreeValues<-trueFreeValuesANDSummaryValues[,1:numberParametersFree]
			#print(trueFreeValues)
			summaryValues<-trueFreeValuesANDSummaryValues[,-1:-numberParametersFree]
			#print(trueFreeValues)
			while(sink.number()>0) {sink()}
			save(trueFreeValues,summaryValues,file=paste("simulations",jobName,".Rdata",sep=""))
			simTime<-proc.time()[[3]]-Time
			cat(paste("Initial simulations took", round(simTime, digits=3), "seconds"), "\n")
			#---------------------- Initial Simulations (End) ------------------------------


			#---------------------- Box-Cox transformation (Start) ------------------------------
			library("car")
			summaryDebugging<-c() #for boxcox debugging
			names(summary)<-c("preTransform", "boxcoxAddition", "postTransform") #for boxcox debugging
			summaryDebugging$preTransform<-summaryValues #for boxcox debugging
			#now put this into the boxcox function to get best lambda for each summary stat
			boxcoxLambda<-rep(NA, dim(summaryValues)[2])
			boxcoxAddition<-rep(NA, dim(summaryValues)[2])
			for (summaryValueIndex in 1:dim(summaryValues)[2]) {
				boxcoxAddition[summaryValueIndex]<-0
				lowValue<-min(summaryValues[, summaryValueIndex])-4*sd(summaryValues[, summaryValueIndex])
				#print (lowValue)
				if (lowValue<=0) {
					boxcoxAddition[summaryValueIndex]<-4*abs(lowValue) #just for some protection against low values, since box.cox needs non-negative values
				}
				#cat("\nsummary values[", summaryValueIndex, "] = ")
				#print(summaryValues[, summaryValueIndex])
				summary<-summaryValues[, summaryValueIndex]+boxcoxAddition[summaryValueIndex]
				summaryDebugging$boxcoxAddition<-summary #for boxcox debugging
				boxcoxLambda[summaryValueIndex]<-1
				if(sd(summaryValues[, summaryValueIndex])>0) { #box.cox fails if all values are identical
					#print("now calculating newLambda")
					#print("summary")
					#print(summary)
					newLambda<-as.numeric(try(powerTransform(summary,method="Nelder-Mead")$lambda)) #new car uses powerTransform instead of box.cox.powers
					#print("done calculating newLambda")
					if (!is.na(newLambda)) {
						boxcoxLambda[summaryValueIndex]<-newLambda
						#print(boxcoxLambda)
					}
				}
				summaryValues[, summaryValueIndex]<-summary^boxcoxLambda[summaryValueIndex]
			}
			summaryDebugging$postTransform<-summaryValues
			save(summaryDebugging, file=paste("summaryDebugging", jobName, ".Rdata", sep=""))
			#---------------------- Box-Cox transformation (End) ------------------------------


			#----------------- Find best set of summary stats to use for this problem. (Start) -----------------
			#Use mixOmics to to find the optimal set of summary stats. Store this info in the todo vector. Note that this uses a different package (mixOmics rather than pls than that used by Weggman et al. because this package can calculate variable importance in projection and deals fine with NAs)
			library("mixOmics")
			plsResult<-pls(Y=trueFreeValues, X=summaryValues)
			vipResult<-vip(plsResult)
			todo<-rep(1, dim(summaryValues)[2]) #initialize the vector that indicates which summary stats to include
			
			summaryIndexOffset=0 #since R excludes invariant columns from regression, this offests so we don't try to extract from these columns
			#print(vipResult)
			#print(plsResult)
			#print(dim(summaryValues))
			nearZeroVarVector<-mixOmics:::nearZeroVar(summaryValues)
			nearZeroVarVector<-nearZeroVarVector$Position
			#print(nearZeroVarVector)
			for (summaryIndex in 1:dim(summaryValues)[2]) {
				
				#print(summaryIndex)
				if (summaryIndex %in% nearZeroVarVector) {
					summaryIndexOffset=summaryIndexOffset+1
					todo[summaryIndex]<-0 #exclude this summary stat because it lacks variation
				}	
				else if (max(vipResult[summaryIndex-summaryIndexOffset, ]) < vipthresh) {
					todo[summaryIndex]<-0 #exclude this summary stat, because it is too unimportant
				}	
			}
		
			while(sink.number()>0) {sink()}
			#print(todo)
			
			prunedSummaryValues<-summaryValues[, which(todo>0)]
			#print("prunedSummaryValues", prunedSummaryValues, "\n")
			prunedPlsResult<-pls(Y=trueFreeValues, X=prunedSummaryValues)
			#print("prunedPlsResult", prunedPlsResult, "\n")
			originalSummaryStats<-boxcoxplsSummary(todo, summaryStatsLong(phy, traits, todo, jobName=jobName), prunedPlsResult, boxcoxLambda, boxcoxAddition)
			
			boxcox.output<-vector("list", 5)
			boxcox.output$lambda<-boxcoxLambda
			boxcox.output$addition<-boxcoxAddition
			boxcox.output$PlsResult<-prunedPlsResult
			boxcox.output$prunedSummaryValues<-prunedSummaryValues
			boxcox.output$originalSummaryStats<-originalSummaryStats
			#----------------- Find best set of summary stats to use for this problem. (End) -----------------
			
			#----------------- Find distribution of distances (Start) ----------------------
			predictResult<-as.matrix(predict(prunedPlsResult, prunedSummaryValues)$predict[, , 1])
			#print(predictResult) 
			#print(dim(predictResult)[1])
			distanceVector<-rep(NA, dim(predictResult)[1])
		
			for (simulationIndex in 1:dim(predictResult)[1]) {
					distanceVector[simulationIndex]<-dist(matrix(c(trueFreeValues[simulationIndex, ], predictResult[simulationIndex, ]), nrow=2, byrow=TRUE))[1]
			}
			#print("distanceVector", distanceVector, "\n")
			densityDistanceVector<-density(distanceVector)
			#plot(densityDistanceVector)
			epsilonDistance<-quantile(distanceVector, probs=epsilonProportion) #this gives the distance such that epsilonProportion of the simulations starting from a given set of values will be rejected 
			toleranceVector<-rep(epsilonDistance, nStepsPRC)
			
			if(nStepsPRC>1){
				for (step in 2:nStepsPRC) {
					toleranceVector[step]<-toleranceVector[step-1]*epsilonMultiplier
				}
			}
			#----------------- Find distribution of distances (End) ---------------------
			
			#------------------ ABC-PRC (Start) ------------------
			#do not forget to use boxcoxLambda, and prunedPlsResult when computing distances
			
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
			particleVector<-c()
			#cat("originalSummaryStats\n")
			#print(originalSummaryStats)
			while (particle<=numParticles) {
				attempts<-attempts+1
				
				newparticleVector<-c(new("abcparticle", id=particle, generation=1, weight=0))
				newparticleVector[[1]]<-initializeStatesFromMatrices(newparticleVector[[1]], startingPriorsValues, startingPriorsFns, intrinsicPriorsValues, intrinsicPriorsFns, extrinsicPriorsValues, extrinsicPriorsFns)
				#cat("\nextrinsicVector\n")
				#print(extrinsicValues(newparticleVector[[1]]))
				#cat("\nintrinsicVector\n")
				#print(intrinsicValues(newparticleVector[[1]]))
		
				newparticleVector[[1]]<-setDistance(newparticleVector[[1]], dist(matrix(c(boxcoxplsSummary(todo, summaryStatsLong(phy, convertTaxonFrameToGeigerData(doSimulation(splits, intrinsicFn, extrinsicFn, startingStates(newparticleVector[[1]]), intrinsicValues(newparticleVector[[1]]), extrinsicValues(newparticleVector[[1]]), timeStep), phy), todo, jobName=jobName), prunedPlsResult, boxcoxLambda, boxcoxAddition), originalSummaryStats), nrow=2, byrow=TRUE))[1])
				if (is.na(distance(newparticleVector[[1]]))) {
					newparticleVectorError<-vector("list", 9)
					newparticleVectorError[[1]]<-newparticleVector[[1]]
					newparticleVectorError[[2]]<-matrix(c(boxcoxplsSummary(todo, summaryStatsLong(phy, convertTaxonFrameToGeigerData(doSimulation(splits, intrinsicFn, extrinsicFn, startingStates(newparticleVector[[1]]), intrinsicValues(newparticleVector[[1]]), extrinsicValues(newparticleVector[[1]]), timeStep), phy), todo, jobName=jobName), prunedPlsResult, boxcoxLambda, boxcoxAddition), originalSummaryStats), nrow=2, byrow=TRUE)
					newparticleVectorError[[3]]<-c(boxcoxplsSummary(todo, summaryStatsLong(phy, convertTaxonFrameToGeigerData(doSimulation(splits, intrinsicFn, extrinsicFn, startingStates(newparticleVector[[1]]), intrinsicValues(newparticleVector[[1]]), extrinsicValues(newparticleVector[[1]]), timeStep), phy), todo, jobName=jobName), prunedPlsResult, boxcoxLambda, boxcoxAddition), originalSummaryStats)
					newparticleVectorError[[4]]<-todo
					newparticleVectorError[[5]]<-summaryStatsLong(phy, convertTaxonFrameToGeigerData(doSimulation(splits, intrinsicFn, extrinsicFn, startingStates(newparticleVector[[1]]), intrinsicValues(newparticleVector[[1]]), extrinsicValues(newparticleVector[[1]]), timeStep), phy), todo, jobName=jobName)
					newparticleVectorError[[6]]<-phy
					newparticleVectorError[[7]]<-vipResult
					newparticleVectorError[[8]]<-convertTaxonFrameToGeigerData(doSimulation(splits, intrinsicFn, extrinsicFn, startingStates(newparticleVector[[1]]), intrinsicValues(newparticleVector[[1]]), extrinsicValues(newparticleVector[[1]]), timeStep), phy)
					newparticleVectorError[[9]]<-convertTaxonFrameToGeigerData(doSimulation(splits, intrinsicFn, extrinsicFn, startingStates(newparticleVector[[1]]), intrinsicValues(newparticleVector[[1]]), extrinsicValues(newparticleVector[[1]]), timeStep), phy)
		
					save(newparticleVectorError, file=paste("newparticleVectorError", jobName, ".txt", sep=""))
		
					while(sink.number()>0) {sink()}
					warning("distance(newparticleVector[[1]]) = NA, likely an underflow/overflow problem")
					newparticleVector[[1]]<-setId(newparticleVector[[1]], -1)
					newparticleVector[[1]]<-setWeight(newparticleVector[[1]], 0)
				}
				else if (is.na(toleranceVector[1])) {
					while(sink.number()>0) {sink()}
					warning("toleranceVector[1] = NA")
					newparticleVector[[1]]<-setId(newparticleVector[[1]], -1)
					newparticleVector[[1]]<-setWeight(newparticleVector[[1]], 0)
				}
						
						
				else if ((distance(newparticleVector[[1]])) < toleranceVector[1]) {
					newparticleVector[[1]]<-setId(newparticleVector[[1]], particle)
					newparticleVector[[1]]<-setWeight(newparticleVector[[1]], 1/numParticles)
					particleWeights[particle]<-1/numParticles
					particle<-particle+1
					particleVector<-append(particleVector, newparticleVector)
				}
				else {
					newparticleVector[[1]]<-setId(newparticleVector[[1]], -1)
					newparticleVector[[1]]<-setWeight(newparticleVector[[1]], 0)
				}
				while(sink.number()>0) {sink()}
				#print(newparticleVector)
				vectorForDataFrame<-c(1, attempts, getId(newparticleVector[[1]]), 0, distance(newparticleVector[[1]]), getWeight(newparticleVector[[1]]), startingStates(newparticleVector[[1]]), intrinsicValues(newparticleVector[[1]]), extrinsicValues(newparticleVector[[1]]))
				#cat("\n\nlength of vectorForDataFrame = ", length(vectorForDataFrame), "\n", "length of startingStates = ", length(startingStates), "\nlength of intrinsicValues = ", length(intrinsicValues), "\nlength of extrinsicValues = ", length(extrinsicValues), "\ndistance = ", distance(newparticleVector[[1]]), "\nweight = ", getWeight(newparticleVector[[1]]), "\n", vectorForDataFrame, "\n")
				particleDataFrame<-rbind(particleDataFrame, data.frame(rbind(vectorForDataFrame)))
				cat(particle-1, attempts, floor(numParticles*attempts/particle), startingStates(newparticleVector[[1]]), intrinsicValues(newparticleVector[[1]]), extrinsicValues(newparticleVector[[1]]), distance(newparticleVector[[1]]), "\n")
					if (floor(numParticles*attempts/particle)>=floor(numParticles)*whenToKill){
						run.goingwell=FALSE
						cat ("\n\nexpected number of generations is too high\n\n")
						break 
					}
		} #while (particle<=numParticles) bracket
		
		
			dataGenerationStep=1
			time<-proc.time()[[3]]-start.time
			time.per.gen<-time
			rejects.gen.one<-(dim(subset(particleDataFrame, X3<0))[1])/(dim(subset(particleDataFrame,))[1])
		
			for (i in 1:numberParametersTotal){
				param.stdev[1,i]<-c(sd(subset(particleDataFrame, X3>0)[,6+i]))
				weightedMeanParam[1,i]<-weighted.mean(subset(particleDataFrame, X3>0)[,6+i], subset(particleDataFrame, X3>0)[,6])
				#c(mean(subset(particleDataFrame, X3>0)[,7:dim(particleDataFrame)[2]])/subset(particleDataFrame, X3>0)[,6])
			}
			#stdev.Intrinsic[i]<-sd(subset(all.a[[run]][which(all.a[[run]]$weight>0),], generation==i)[,param[2]])
		
		
			if (!run.goingwell){	
				if (try==maxTries){
					write(input.data,file="Error.txt", append=TRUE)
					ErrorParticleFrame<-vector("list", 4)
					names(particleDataFrame)<-nameVector
					ErrorParticleFrame[[1]]<-input.data
					ErrorParticleFrame[[2]]<-todo
					ErrorParticleFrame[[3]]<-particleDataFrame
					ErrorParticleFrame[[4]]<-toleranceVector
					#ErrorParticleFrame->paste("ErrorParticleFrame", jobName, sep="")
					save(ErrorParticleFrame, file=paste("ErrorParticleFrame", jobName, ".txt", sep=""))
					cat("\n\nTried", maxTries, "times and all failed!")
					cat("\ninput.data was appended to Error.txt file within the working directory\n\n")
				}
				else if (try < maxTries){
					ErrorParticleFrame<-vector("list", 4)
					names(particleDataFrame)<-nameVector
					ErrorParticleFrame[[1]]<-input.data
					ErrorParticleFrame[[2]]<-todo
					ErrorParticleFrame[[3]]<-particleDataFrame
					ErrorParticleFrame[[4]]<-toleranceVector
					#ErrorParticleFrame->paste("ErrorParticleFrame", jobName, sep="")
					save(ErrorParticleFrame, file=paste("ErrorParticleFrame", jobName, ".txt", sep=""))
					cat("\n\nAborting try", try, "of", maxTries, "at Generation 1\n\n")
				}
				break
			}	
			
			if (checkpointSave){
				save.image(file=paste("WS", jobName, ".Rdata", sep=""))
				test<-vector("list", 16)
				names(test)<-c("input.data", "PriorMatrix", "particleDataFrame", "epsilonDistance", "toleranceVector", "todo", "phy", "traits", "rejects.gen.one", "rejects", "particleWeights", "particleVector", "boxcox.output", "param.stdev", "weightedMeanParam", "time.per.gen")
				test$input.data<-input.data
				test$PriorMatrix<-PriorMatrix
				test$particleDataFrame<-particleDataFrame
				names(test$particleDataFrame)<-nameVector
				test$epsilonDistance<-epsilonDistance
				test$toleranceVector<-toleranceVector
				test$todo<-todo
				test$phy<-phy
				test$traits<-traits
				test$rejects.gen.one<-rejects.gen.one
				test$rejects<-c()
				test$particleWeights<-particleWeights
				test$particleVector<-particleVector
				test$boxcox.output<-boxcox.output
				test$param.stdevtest$param.stdev<-param.stdev
				test$weightedMeanParam<-weightedMeanParam
				test$time.per.gen<-time.per.gen
				save(test, file=paste("partialResults", jobName, ".txt", sep=""))
		
			}	
			
		} #startFromCheckpoint bracket
		
		
			if (startFromCheckpoint==TRUE || dataGenerationStep < nStepsPRC) {
			
				if (run.goingwell){
				
					while (dataGenerationStep < nStepsPRC) {
						dataGenerationStep<-dataGenerationStep+1
						cat("\n\n\n", "STARTING DATA GENERATION STEP ", dataGenerationStep, "\n\n\n")
							if (debug){
								cat(".Random Seed=", .Random.seed[1:6], "\n\n")
								cat("runif=", runif(1), "\n\n")
								#cat("toleranceVector=", toleranceVector, "\n\n")
							}
						start.time<-proc.time()[[3]]
						particleWeights<-particleWeights/(sum(particleWeights,na.rm=TRUE)) #normalize to one
						cat("particleWeights\n", particleWeights, "\n\n")
						oldParticleVector<-particleVector
						oldParticleWeights<-particleWeights
						particleWeights=rep(0, numParticles) #stores weights for each particle. Initially, assume infinite number of possible particles (so might not apply in discrete case)
						particleParameters<-matrix(nrow=numParticles, ncol=dim(startingPriorsValues)[2] +  dim(intrinsicPriorsValues)[2] + dim(extrinsicPriorsValues)[2]) #stores parameters in model for each particle
						particleDistance=rep(NA, numParticles)
						particle<-1
						attempts<-0
						cat("successes", "attempts", "expected number of attempts required\n")
						particleVector<-c()
						weightScaling=0;
						while (particle<=numParticles) {
							attempts<-attempts+1
							particleToSelect<-which.max(as.vector(rmultinom(1, size = 1, prob=oldParticleWeights)))
							#cat("particle to select = ", particleToSelect, "\n")
							#cat("dput(oldParticleVector)\n")
							#dput(oldParticleVector)
							#cat("dput(oldParticleVector[particleToSelect])\n")
							#dput(oldParticleVector[particleToSelect])
							#cat("dput(oldParticleVector[[particleToSelect]])\n")
							#dput(oldParticleVector[[particleToSelect]])
							newparticleVector<-oldParticleVector[particleToSelect]
							#cat("dput(newparticleVector[[1]])\n")
							#dput(newparticleVector[[1]])
							#cat("mutateStates\n")
							
							newparticleVector[[1]]<-mutateStates(newparticleVector[[1]], startingPriorsValues, startingPriorsFns, intrinsicPriorsValues, intrinsicPriorsFns, extrinsicPriorsValues, extrinsicPriorsFns, standardDevFactor)
							#cat("dput(newparticleVector[[1]]) AFTER MUTATE STATES\n")
							#dput(newparticleVector[[1]])
							
							newparticleVector[[1]]<-setDistance(newparticleVector[[1]], dist(matrix(c(boxcoxplsSummary(todo, summaryStatsLong(phy, convertTaxonFrameToGeigerData(doSimulation(splits, intrinsicFn, extrinsicFn, startingStates(newparticleVector[[1]]), intrinsicValues(newparticleVector[[1]]), extrinsicValues(newparticleVector[[1]]), timeStep), phy), todo, jobName=jobName), prunedPlsResult, boxcoxLambda, boxcoxAddition), originalSummaryStats), nrow=2, byrow=TRUE))[1])
							if (plot) {
								plotcol="grey"
								if (distance(newparticleVector[[1]])<toleranceVector[dataGenerationStep]) {
									plotcol="black"
									if (dataGenerationStep==length(toleranceVector)) {
										plotcol="red"
									}
								}
								text(x=intrinsicValues(newparticleVector[[1]]), y=distance(newparticleVector[[1]]), labels= dataGenerationStep, col=plotcol) 
							}
							#dput(convertTaxonFrameToGeigerData(doSimulation(splits, intrinsicFn, extrinsicFn, newparticleVector[[1]]@startingStates, newparticleVector[[1]]@intrinsicValues, newparticleVector[[1]]@extrinsicValues, timeStep), phy))
							#cat("dput(newparticleVector[[1]]) AFTER computeABCDistance\n")
							#dput(newparticleVector[[1]])
						
							if (is.na(distance(newparticleVector[[1]]))) {
								cat("Error with Geiger?  distance(newparticleVector[[1]]) = NA\n")
								while(sink.number()>0) {sink()}
								warning("distance(newparticleVector[[1]]) = NA")
								newparticleVector[[1]]<-setId(newparticleVector[[1]], -1)
								newparticleVector[[1]]<-setWeight(newparticleVector[[1]], 0)
							}
							else if (distance(newparticleVector[[1]]) < toleranceVector[dataGenerationStep]) {
								newparticleVector[[1]]<-setId(newparticleVector[[1]], particle)
								particle<-particle+1
								particleVector<-append(particleVector, newparticleVector)
								#now get weights, using correction in Sisson et al. 2007
								newWeight=0
								for (i in 1:length(oldParticleVector)) {
									lnTransitionProb=log(1)
									for (j in 1:length(newparticleVector[[1]]@startingStates)) {
										newvalue<-newparticleVector[[1]]@startingStates[j]
										meantouse= oldParticleVector[[i]]@startingStates[j]
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
												#print(paste("startingPriorFn is not uniform or exponential and sdtouse =", sdtouse))
											}
																
										#print(paste("@startingStates: newvalue=", newvalue, "meantouse=", meantouse, "sdtouse=", sdtouse))
										#OLDlnlocalTransitionProb=dnorm(newvalue, mean= meantouse, sd= sdtouse,log=TRUE)-log(1-pnorm(min(startingPriorsValues[, j]), mean= meantouse , sd= sdtouse, lower.tail=T)+pnorm(max(startingPriorsValues[, j]), mean= meantouse , sd= sdtouse, lower.tail=F))  #this was the way we were doing it up until 6/15/11, but were getting log(0) = Inf errors
										lnlocalTransitionProb=dnorm(newvalue, mean=meantouse, sd=sdtouse,log=TRUE) - ((log(1)/pnorm(min(startingPriorsValues[, j]), mean=meantouse, sd=sdtouse, lower.tail=T, log.p=T))* pnorm(max(startingPriorsValues[,j]), mean=meantouse , sd=sdtouse, lower.tail=F, log.p=T))
				
										if (lnlocalTransitionProb == "NaN") {  #to prevent lnlocalTransitionProb from being NaN (if pnorm=0)
											lnlocalTransitionProb<-.Machine$double.xmin
										}
										#print(paste("OLDlnlocalTransitionProb =", OLDlnlocalTransitionProb, "lnlocalTransitionProb =", lnlocalTransitionProb))
										#print(paste("@startingStates: dnorm()=", dnorm(newvalue, mean= meantouse, sd= sdtouse,log=TRUE), ", 1-pnorm()=", 1-pnorm(min(startingPriorsValues[, j]), mean= meantouse , sd= sdtouse, lower.tail=T), ", pnorm()=", pnorm(max(startingPriorsValues[, j]), mean= meantouse , sd= sdtouse, lower.tail=F)))
										#print(paste("lnlocalTransitionProb =", lnlocalTransitionProb))
										if (min(startingPriorsValues[, j])==max(startingPriorsValues[, j])) {
											lnlocalTransitionProb=log(1)
										} 
										lnTransitionProb<-lnTransitionProb+lnlocalTransitionProb
										#print(paste("lnlocalTransitionProb=", lnlocalTransitionProb))
										if(!is.finite(lnTransitionProb) || is.na(lnlocalTransitionProb)) {
											print(paste("issue with lnTransitionProb: lnlocalTransitionProb = ",lnlocalTransitionProb," lnTransitionProb = ",lnTransitionProb))
										}
									} 
									for (j in 1:length(newparticleVector[[1]]@intrinsicValues)) {
										newvalue<-newparticleVector[[1]]@intrinsicValues[j]
										meantouse= oldParticleVector[[i]]@intrinsicValues[j]
										if (intrinsicPriorsFns[j]=="uniform") {
											sdtouse<-standardDevFactor*((max(intrinsicPriorsValues[,j])-min(intrinsicPriorsValues[,j]))/sqrt(12))
											#print(paste("intrinsicPriorFn is uniform and sdtouse =", sdtouse))
										}
										else if (intrinsicPriorsFns[j]=="exponential") {
											sdtouse<-standardDevFactor*(1/intrinsicPriorsValues[,j])
											#print(paste("intrinsicPriorFn is exponential and sdtouse =", sdtouse))
										}
										else {
											sdtouse<-standardDevFactor*(intrinsicPriorsValues[2,j])
											#print(paste("intrinsicPriorFn is not uniform or exponential and sdtouse =", sdtouse))
										}
										#print(paste("@intrinsicValues: newvalue =", newvalue, "meantouse=", meantouse, "sdtouse=", sdtouse))
										#OLDlnlocalTransitionProb=dnorm(newvalue, mean= meantouse, sd= sdtouse,log=TRUE)-log(1-pnorm(min(intrinsicPriorsValues[, j]), mean= meantouse , sd= sdtouse, lower.tail=T)+pnorm(max(intrinsicPriorsValues[, j]), mean= meantouse , sd= sdtouse, lower.tail=F))  #this was the way we were doing it up until 6/15/11, but were getting log(0) = Inf errors
										lnlocalTransitionProb=dnorm(newvalue, mean= meantouse, sd= sdtouse,log=TRUE)-((log(1)/pnorm(min(intrinsicPriorsValues[, j]), mean=meantouse , sd=sdtouse, lower.tail=T, log.p=T)) * pnorm(max(intrinsicPriorsValues[,j]), mean=meantouse , sd=sdtouse, lower.tail=F, log.p=T))
				
										if (lnlocalTransitionProb == "NaN") {  #to prevent lnlocalTransitionProb from being NaN (if pnorm=0)
											lnlocalTransitionProb<-.Machine$double.xmin
										}
															   
										#print(paste("OLDlnlocalTransitionProb =", OLDlnlocalTransitionProb, "lnlocalTransitionProb =", lnlocalTransitionProb))
										#print(paste("@intrinsicValues: dnorm()=", dnorm(newvalue, mean= meantouse, sd= sdtouse,log=TRUE), ", 1-pnorm()=", 1-pnorm(min(startingPriorsValues[, j]), mean= meantouse , sd= sdtouse, lower.tail=T), ", pnorm()=", pnorm(max(startingPriorsValues[, j]), mean= meantouse , sd= sdtouse, lower.tail=F)))
										if (min(intrinsicPriorsValues[, j])==max(intrinsicPriorsValues[, j])) {
											lnlocalTransitionProb=log(1)
										} 
										lnTransitionProb<-lnTransitionProb+lnlocalTransitionProb
										#print(paste("lnlocalTransitionProb=", lnlocalTransitionProb))
										if(!is.finite(lnTransitionProb) || is.na(lnlocalTransitionProb)) {
											print(paste("issue with lnTransitionProb: lnlocalTransitionProb = ",lnlocalTransitionProb," lnTransitionProb = ",lnTransitionProb))
										}
				
									} 
									for (j in 1:length(newparticleVector[[1]]@extrinsicValues)) {
										newvalue<-newparticleVector[[1]]@extrinsicValues[j]
										meantouse= oldParticleVector[[i]]@extrinsicValues[j]
										if (extrinsicPriorsFns[j]=="uniform") {
											sdtouse<-standardDevFactor*((max(extrinsicPriorsValues[,j])-min(extrinsicPriorsValues[,j]))/sqrt(12))
											#print(paste("extrinsicPriorFn is uniform and sdtouse =", sdtouse))
										}
										else if (extrinsicPriorsFns[j]=="exponential") {
											sdtouse<-standardDevFactor*(1/extrinsicPriorsValues[,j])
											#print(paste("extrinsicPriorFn is exponential and sdtouse =", sdtouse))
										}
										else {
											sdtouse<-standardDevFactor*(extrinsicPriorsValues[2,j])
											#print(paste("extrinsicPriorFn is not uniform or exponential and sdtouse =", sdtouse))
										}
										#print(paste("@extrinsicValues: newvalue=", newvalue, "meantouse=", meantouse, "sdtouse=", sdtouse))
										#OLDlnlocalTransitionProb=dnorm(newvalue, mean= meantouse, sd= sdtouse,log=TRUE)-log(1-pnorm(min(extrinsicPriorsValues[, j]), mean= meantouse , sd= sdtouse, lower.tail=T)+pnorm(max(extrinsicPriorsValues[, j]), mean= meantouse , sd= sdtouse, lower.tail=F))  #this was the way we were doing it up until 6/15/11, but were getting log(0) = Inf errors
										lnlocalTransitionProb=dnorm(newvalue, mean= meantouse, sd= sdtouse,log=TRUE)-((log(1)/pnorm(min(extrinsicPriorsValues[,j]), mean=meantouse , sd=sdtouse, lower.tail=T, log.p=T)) * pnorm(max(extrinsicPriorsValues[,j]), mean=meantouse , sd=sdtouse, lower.tail=F, log.p=T))
										if (lnlocalTransitionProb == "NaN") {  #to prevent lnlocalTransitionProb from being NaN (if pnorm=0)
											lnlocalTransitionProb<-.Machine$double.xmin
										}
										#print(paste("OLDlnlocalTransitionProb =", OLDlnlocalTransitionProb, "lnlocalTransitionProb =", lnlocalTransitionProb))
										#print(paste("@extrinsicValues: dnorm()=", dnorm(newvalue, mean= meantouse, sd= sdtouse,log=TRUE), ", 1-pnorm()=", 1-pnorm(min(startingPriorsValues[, j]), mean= meantouse , sd= sdtouse, lower.tail=T), ", pnorm()=", pnorm(max(startingPriorsValues[, j]), mean= meantouse , sd= sdtouse, lower.tail=F)))
										if (min(extrinsicPriorsValues[, j])==max(extrinsicPriorsValues[, j])) {
											lnlocalTransitionProb=log(1)
										} 
										lnTransitionProb<-lnTransitionProb+lnlocalTransitionProb
										#print(paste("lnlocalTransitionProb=", lnlocalTransitionProb))
										if(!is.finite(lnTransitionProb) || is.na(lnlocalTransitionProb)) {
											print(paste("issue with lnTransitionProb: lnlocalTransitionProb = ",lnlocalTransitionProb," lnTransitionProb = ",lnTransitionProb))
										}
				
									}                                       
									newWeight<-newWeight+getWeight(oldParticleVector[[i]])*exp(lnTransitionProb)
								} #for (i in 1:length(oldParticleVector)) bracket
				
								if (!is.finite(newWeight)) {
									print(paste("warning: newWeight is ",newWeight))
								}
								newparticleVector[[1]]<-setWeight(newparticleVector[[1]], newWeight)
								particleWeights[particle-1]<-newWeight
								weightScaling<-weightScaling+getWeight(newparticleVector[[1]])
							} #else if (distance(newparticleVector[[1]]) < toleranceVector[dataGenerationStep]) bracket
							else {
								newparticleVector[[1]]<-setId(newparticleVector[[1]], -1)
								newparticleVector[[1]]<-setWeight(newparticleVector[[1]], 0)
							}
							while(sink.number()>0) {sink()}
							#print(newparticleVector)
							vectorForDataFrame<-c(dataGenerationStep, attempts, getId(newparticleVector[[1]]), particleToSelect, distance(newparticleVector[[1]]), getWeight(newparticleVector[[1]]), startingStates(newparticleVector[[1]]), intrinsicValues(newparticleVector[[1]]), extrinsicValues(newparticleVector[[1]]))
				#cat("\n\nlength of vectorForDataFrame = ", length(vectorForDataFrame), "\n", "length of startingStates = ", length(startingStates), "\nlength of intrinsicValues = ", length(intrinsicValues), "\nlength of extrinsicValues = ", length(extrinsicValues), "\ndistance = ", distance(newparticleVector[[1]]), "\nweight = ", getWeight(newparticleVector[[1]]), "\n", vectorForDataFrame, "\n")
							particleDataFrame<-rbind(particleDataFrame, data.frame(rbind(vectorForDataFrame))) #NOTE THAT WEIGHTS AREN'T NORMALIZED IN THIS DATAFRAME
							cat(particle-1, attempts, floor(numParticles*attempts/particle), startingStates(newparticleVector[[1]]), intrinsicValues(newparticleVector[[1]]), extrinsicValues(newparticleVector[[1]]), distance(newparticleVector[[1]]), "\n")
							if (floor(numParticles*attempts/particle)>=floor(numParticles)*whenToKill){
								run.goingwell=FALSE
								cat ("\n\nexpected number of generations is too high\n\n")
								break 
							}
				
						} #while (particle<=numParticles) bracket
					
						if (!run.goingwell){
							break
						}		
					
				
						particleDataFrame[which(particleDataFrame$generation==dataGenerationStep), ]$weight<-particleDataFrame[which(particleDataFrame$generation==dataGenerationStep), ]$weight/(sum(particleDataFrame[which(particleDataFrame$generation==dataGenerationStep), ]$weight))
				
						time2<-proc.time()[[3]]-start.time
						time.per.gen<-c(time.per.gen, time2)
						rejects.per.gen<-(dim(subset(particleDataFrame, X3<0))[1])/(dim(subset(particleDataFrame[which(particleDataFrame$X1==dataGenerationStep),],))[1])
						rejects<-c(rejects, rejects.per.gen)
						sub1<-subset(particleDataFrame, X1==dataGenerationStep)
						sub2<-subset(sub1, X3>0)
						
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
						
						if (checkpointSave){
							save.image(file=paste("WS", jobName, ".Rdata", sep=""))
							test<-vector("list", 16)
							names(test)<-c("input.data", "PriorMatrix", "particleDataFrame", "epsilonDistance", "toleranceVector", "todo", "phy", "traits", "rejects.gen.one", "rejects", "particleWeights", "particleVector", "boxcox.output", "param.stdev", "weightedMeanParam", "time.per.gen")
							test$input.data<-input.data
							test$PriorMatrix<-PriorMatrix
							test$particleDataFrame<-particleDataFrame
							names(test$particleDataFrame)<-nameVector
							test$epsilonDistance<-epsilonDistance
							test$toleranceVector<-toleranceVector
							test$todo<-todo
							test$phy<-phy
							test$traits<-traits
							test$rejects.gen.one<-rejects.gen.one
							test$rejects<-rejects
							test$particleWeights<-particleWeights
							test$particleVector<-particleVector
							test$boxcox.output<-boxcox.output
							test$param.stdev<-param.stdev
							test$weightedMeanParam<-weightedMeanParam
							test$time.per.gen<-time.per.gen
							save(test, file=paste("partialResults", jobName, ".txt", sep=""))
						}	
						
						
					} #while (dataGenerationStep < nStepsPRC) bracket
				} #if (run.goingwell) bracket
			
				if (!run.goingwell){	
					if (try==maxTries){
						write(input.data,file="Error.txt", append=TRUE)
						ErrorParticleFrame<-vector("list", 4)
						names(particleDataFrame)<-nameVector
						ErrorParticleFrame[[1]]<-input.data
						ErrorParticleFrame[[2]]<-todo
						ErrorParticleFrame[[3]]<-particleDataFrame
						ErrorParticleFrame[[4]]<-toleranceVector
						save(ErrorParticleFrame, file=paste("ErrorParticleFrame", jobName, ".txt", sep=""))
						cat("\n\nTried", maxTries, "times and all failed!")
						cat("\ninput.data was appended to Error.txt file within the working directory\n\n")
					}
					else if (try < maxTries){
						ErrorParticleFrame<-vector("list", 4)
						names(particleDataFrame)<-nameVector
						ErrorParticleFrame[[1]]<-input.data
						ErrorParticleFrame[[2]]<-todo
						ErrorParticleFrame[[3]]<-particleDataFrame
						ErrorParticleFrame[[4]]<-toleranceVector
						save(ErrorParticleFrame, file=paste("ErrorParticleFrame", jobName, ".txt", sep=""))
						cat("\n\nAborting try", try, "of", maxTries, "at Generation", dataGenerationStep, "\n\n")
					}
					break
				}		
				
				names(particleDataFrame)<-nameVector
				if(plot) {
					quartz()
					plot(x=c(min(intrinsicPriorsValues), max(intrinsicPriorsValues)), y=c(0, 1), type="n")
					for (i in 1:(length(toleranceVector)-1)) {
						graycolor<-gray(0.5*(length(toleranceVector)-i)/length(toleranceVector))
						lines(density(subset(particleDataFrame, generation==i)[, 8]), col= graycolor)
					}
					lines(density(subset(particleDataFrame, generation==length(toleranceVector))[, 8]), col= "red")
				} 
				
				
			} # if(startFromCheckpoint==TRUE || dataGenerationStep < nStepsPRC) bracket
			
			#---------------------- ABC-PRC (End) --------------------------------
			
			if (debug){
				cat("debug!!")
				debugResults<-dget(doRun)
			}
		
		
		} #while (!run.goingwell) bracket
		
		#Calculate summary statistics on final generation particles
		FinalParamPredictions<-matrix(nrow=numberParametersTotal, ncol=4)
		colnames(FinalParamPredictions)<-c("weightedMean", "SD", "Low95%", "High95%")
		rownames(FinalParamPredictions)<-names(particleDataFrame[7: dim(particleDataFrame)[2]])
		subpDF<-subset(particleDataFrame[which(particleDataFrame$weight>0),], generation==max(particleDataFrame $generation))  #change this to be names and then call names in next for loop
		#names(subpDF)<-c("input.data", "PriorMatrix", "particleDataFrame", "epsilonDistance", "toleranceVector", "todo", "phy", "traits", "rejects.gen.one", "rejects", "particleWeights", "particleVector", "boxcox.output", "param.stdev", "weightedMeanParam", "time.per.gen", "FinalParamPredictions")
		for (paramPred in 1:numberParametersTotal){
			#print(6+paramPred)
			FinalParamPredictions[paramPred, 1]<-weighted.mean(subpDF[,6+paramPred], subpDF[,6])
			FinalParamPredictions[paramPred, 2]<-sd(subpDF[,(6+paramPred)])
			FinalParamPredictions[paramPred, 3]<-quantile(subpDF[,6+paramPred], probs=0.05)
			FinalParamPredictions[paramPred, 4]<-quantile(subpDF[,6+paramPred], probs=0.95)
		}
	
		input.data<-rbind(jobName, length(phy[[3]]), nrepSim, epsilonProportion, epsilonMultiplier, nStepsPRC, numParticles, standardDevFactor, try, trueStartingState, trueIntrinsicState)
	
		time3<-proc.time()[[3]]
		genTimes<-c(time.per.gen, time3)
	
		test<-vector("list", 17)
		names(test)<-c("input.data", "PriorMatrix", "particleDataFrame", "epsilonDistance", "toleranceVector", "todo", "phy", "traits", "rejects.gen.one", "rejects", "particleWeights", "particleVector", "boxcox.output", "param.stdev", "weightedMeanParam", "time.per.gen", "FinalParamPredictions")
	
		test$input.data<-input.data
		test$PriorMatrix<-PriorMatrix
		test$particleDataFrame<-particleDataFrame
		test$epsilonDistance<-epsilonDistance
		test$toleranceVector<-toleranceVector
		test$todo<-todo
		test$phy<-phy
		test$traits<-traits
		test$rejects.gen.one<-rejects.gen.one
		test$rejects<-rejects
		test$particleWeights<-particleWeights
		test$particleVector<-particleVector
		test$boxcox.output<-boxcox.output
		test$param.stdev<-param.stdev
		test$weightedMeanParam<-weightedMeanParam
		test$simTime<-simTime
		test$time.per.gen<-genTimes
		test$FinalParamPredictions <-FinalParamPredictions
	
	} #for (try in 1:maxTries) bracket

	registerDoMC(1)
	print(test)

}

	#------------------ ABC-PRC (end) ------------------
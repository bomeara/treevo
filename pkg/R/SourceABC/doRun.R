#the doRun function takes input from the user and then automatically guesses optimal parameters, though user overriding is also possible.
#the guesses are used to do simulations near the expected region. If omitted, they are set to the midpoint of the input parameter matrices

doRun<-function(phy, traits, intrinsicFn, extrinsicFn, summaryFns=c(rawValuesSummaryStats, geigerUnivariateSummaryStats2), startingPriorsValues, startingPriorsFns, intrinsicPriorsValues, intrinsicPriorsFns, extrinsicPriorsValues, extrinsicPriorsFns, startingStatesGuess=c(), intrinsicStatesGuess=c(), extrinsicStatesGuess=c(), timeStep, toleranceVector=c(), numParticles=1000, standardDevFactor=0.05, StartSims=100, plot=TRUE, vipthresh=0.8, epsilonProportion=0.2, epsilonMultiplier=0.5, nStepsPRC=4, maxTries=1, job.name=NA, debug=TRUE, TrueStartingState=NA, TrueIntrinsicState=NA, WhenToKill=20) {
	#print("in do run")
trace(doRun)

run.goingwell=FALSE
for (try in 1:maxTries)	{
	nrepSim<-StartSims*((10^try)/10)
while (!run.goingwell) {
	cat("\n\n****  TRY", try, "of", maxTries, " ****\n\n")
		if (debug){
			cat(".Random Seed=", .Random.seed[1:6], "\n\n")
			cat("runif=", runif(1), "\n\n")
			#cat("toleranceVector=", toleranceVector, "\n\n")
		}
	run.finished=FALSE
	run.goingwell=TRUE
	
	splits<-getSimulationSplits(phy) #initialize this info

input.data<-rbind(job.name, length(phy[[3]]), startingPriorsFns, startingPriorsValues, intrinsicPriorsFns, intrinsicPriorsValues, nrepSim, epsilonProportion, epsilonMultiplier, nStepsPRC, numParticles, try)


	#figure out number of free params
	numberParametersTotal<-dim(startingPriorsValues)[2] +  dim(intrinsicPriorsValues)[2] + dim(extrinsicPriorsValues)[2]
	numberParametersFree<-numberParametersTotal
	numberParametersStarting<-0
	numberParametersIntrinsic<-0
	numberParametersExtrinsic<-0
	freevariables<-matrix(data=NA, nrow=2, ncol=0)
	titlevector<-c()
	freevector<-c()
	
		
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
		summaryValues<-rbind(summaryValues, summaryStatsLong(phy, convertTaxonFrameToGeigerData (doSimulation(splits=splits, intrinsicFn= intrinsicFn, extrinsicFn= extrinsicFn, startingStates= trueStarting, intrinsicValues= trueIntrinsic, extrinsicValues= trueExtrinsic, timeStep=timeStep), phy)))
		while(sink.number()>0) {sink()}
	}
	#summaryRprof(nrepSims.time.check)
	library("car")
	#now put this into the boxcox function to get best lambda for each summary stat
	boxcoxLambda<-rep(NA, dim(summaryValues)[2])
	boxcoxAddition<-rep(NA, dim(summaryValues)[2])
	for (summaryValueIndex in 1:dim(summaryValues)[2]) {
		boxcoxAddition[summaryValueIndex]<-0
		lowValue<-min(summaryValues[, summaryValueIndex])-4*sd(summaryValues[, summaryValueIndex])
		if (lowValue<=0) {
			boxcoxAddition[summaryValueIndex]<-4*abs(lowValue) #just for some protection against low values, since box.cox needs non-negative values
		}
		cat("\nsummary values[", summaryValueIndex, "] = ")
		print(summaryValues[, summaryValueIndex])
		summary<-summaryValues[, summaryValueIndex]+boxcoxAddition[summaryValueIndex]
		boxcoxLambda[summaryValueIndex]<-1
		if(sd(summaryValues[, summaryValueIndex])>0) { #box.cox fails if all values are identical
			newLambda<-as.numeric(try(box.cox.powers.hacked(summary)$lambda))
			if (!is.na(newLambda)) {
				boxcoxLambda[summaryValueIndex]<-newLambda
			}
		}
		summaryValues[, summaryValueIndex]<-summary^boxcoxLambda[summaryValueIndex]
	}
	
	
	#Use mixOmics to to find the optimal set of summary stats. Store this info in the todo vector. Note that this uses a different package (mixOmics rather than pls than that used by Weggman et al. because this package can calculate variable importance in projection and deals fine with NAs)
	library("mixOmics")
	plsResult<-pls(Y=trueFreeValues, X=summaryValues)
	vipResult<-vip(plsResult)
	todo<-rep(1, dim(summaryValues)[2])#initialize the vector that indicates which summary stats to include
	for (summaryIndex in 1:dim(summaryValues)[2]) {
		if (max(vipResult[summaryIndex, ])<vipthresh) {
				todo[summaryIndex]<-0 #exclude this summary stat, because it's too unimportant
		}	
	}
	while(sink.number()>0) {sink()}
	#cat("todo", "\n")
	#print(todo)
	
	prunedSummaryValues<-summaryValues[, which(todo>0)]
	#cat("prunedSummaryValues", "\n")
	#print(prunedSummaryValues)

	prunedPlsResult<-pls(Y=trueFreeValues, X=prunedSummaryValues)
	#cat("prunedPlsResult", "\n")
	#print(prunedPlsResult)

	#----------------- Find best set of summary stats to use for this problem. (end) -----------------
	
	#----------------- Find distribution of distances (start) ----------------------
	predictResult<-(predict(prunedPlsResult, prunedSummaryValues)$predict[, , 1])
	#cat("predictResult", "\n")
	#print(predictResult)

	distanceVector<-rep(NA, dim(predictResult)[1])
	for (simulationIndex in 1:dim(predictResult)[1]) {
			distanceVector[simulationIndex]<-dist(matrix(c(trueFreeValues[simulationIndex, ], predictResult[simulationIndex, ]), nrow=2, byrow=TRUE))[1]
	}
	#cat("distanceVector", "\n")
	#print(distanceVector)
	densityDistanceVector<-density(distanceVector)
	#plot(densityDistanceVector)
	epsilonDistance<-quantile(distanceVector, probs=epsilonProportion) #this gives the distance such that epsilonProportion of the simulations starting from a given set of values will be rejected 
	#lines(x=c(epsilonDistance, epsilonDistance), y=c(0, max(densityDistanceVector$y)), lty="dotted")
	toleranceVector<-rep(epsilonDistance, nStepsPRC)
	for (step in 2:nStepsPRC) {
		toleranceVector[step]<-toleranceVector[step-1]*epsilonMultiplier
	}
	#----------------- Find distribution of distances (end) ---------------------
	
	#------------------ ABC-PRC (start) ------------------
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
	#print(extrinsicPriorsValues)
	for (i in 1:dim(extrinsicPriorsValues)[2]) {
		nameVector<-append(nameVector, paste("ExtrinsicValue", i, sep=""))
	}
	particleWeights=rep(0, numParticles) #stores weights for each particle. Initially, assume infinite number of possible particles (so might not apply in discrete case), uniform prior on each
	particleParameters<-matrix(nrow=numParticles, ncol=dim(startingPriorsValues)[2] +  dim(intrinsicPriorsValues)[2] + dim(extrinsicPriorsValues)[2]) #stores parameters in model for each particle
	particleDistance=rep(NA, numParticles)
	particle<-1
	attempts<-0
	particleDataFrame<-data.frame()
	cat("\n\n\nsuccesses", "attempts", "expected number of attempts required\n\n\n")
	start.time<-proc.time()[[3]]
	particleVector<-c()
	originalSummaryStats<-boxcoxplsSummary(todo, summaryStatsLong(phy, traits, todo), prunedPlsResult, boxcoxLambda, boxcoxAddition)
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

		newparticleVector[[1]]<-setDistance(newparticleVector[[1]], dist(matrix(c(boxcoxplsSummary(todo, summaryStatsLong(phy, convertTaxonFrameToGeigerData(doSimulation(splits, intrinsicFn, extrinsicFn, startingStates(newparticleVector[[1]]), intrinsicValues(newparticleVector[[1]]), extrinsicValues(newparticleVector[[1]]), timeStep), phy), todo), prunedPlsResult, boxcoxLambda, boxcoxAddition), originalSummaryStats), nrow=2, byrow=TRUE))[1])
		if (is.na(distance(newparticleVector[[1]]))) {
			cat("Error with Geiger?  distance(newparticleVector[[1]]) = NA")
			GeigerError<-vector("list", 9)
			GeigerError[[1]]<-newparticleVector[[1]]
			GeigerError[[2]]<-matrix(c(boxcoxplsSummary(todo, summaryStatsLong(phy, convertTaxonFrameToGeigerData(doSimulation(splits, intrinsicFn, extrinsicFn, startingStates(newparticleVector[[1]]), intrinsicValues(newparticleVector[[1]]), extrinsicValues(newparticleVector[[1]]), timeStep), phy), todo), prunedPlsResult, boxcoxLambda, boxcoxAddition), originalSummaryStats), nrow=2, byrow=TRUE)
			GeigerError[[3]]<-c(boxcoxplsSummary(todo, summaryStatsLong(phy, convertTaxonFrameToGeigerData(doSimulation(splits, intrinsicFn, extrinsicFn, startingStates(newparticleVector[[1]]), intrinsicValues(newparticleVector[[1]]), extrinsicValues(newparticleVector[[1]]), timeStep), phy), todo), prunedPlsResult, boxcoxLambda, boxcoxAddition), originalSummaryStats)
			GeigerError[[4]]<-todo
			GeigerError[[5]]<-summaryStatsLong(phy, convertTaxonFrameToGeigerData(doSimulation(splits, intrinsicFn, extrinsicFn, startingStates(newparticleVector[[1]]), intrinsicValues(newparticleVector[[1]]), extrinsicValues(newparticleVector[[1]]), timeStep), phy), todo)
			GeigerError[[6]]<-phy
			GeigerError[[7]]<-vipResult
			GeigerError[[8]]<-convertTaxonFrameToGeigerData(doSimulation(splits, intrinsicFn, extrinsicFn, startingStates(newparticleVector[[1]]), intrinsicValues(newparticleVector[[1]]), extrinsicValues(newparticleVector[[1]]), timeStep), phy)
			GeigerError[[9]]<-convertTaxonFrameToGeigerData(doSimulation(splits, intrinsicFn, extrinsicFn, startingStates(newparticleVector[[1]]), intrinsicValues(newparticleVector[[1]]), extrinsicValues(newparticleVector[[1]]), timeStep), phy)

			save(GeigerError, file=paste("GeigerError", job.name, ".txt", sep=""))
			#break
			while(sink.number()>0) {sink()}
			warning("distance(newparticleVector[[1]]) = NA")
			newparticleVector[[1]]<-setId(newparticleVector[[1]], -1)
			newparticleVector[[1]]<-setWeight(newparticleVector[[1]], 0)
		}
		else if (is.na(toleranceVector[1])) {
			while(sink.number()>0) {sink()}
			warning("toleranceVector[1] = NA")
			newparticleVector[[1]]<-setId(newparticleVector[[1]], -1)
			newparticleVector[[1]]<-setWeight(newparticleVector[[1]], 0)
		}
		else if (distance(newparticleVector[[1]])<toleranceVector[1]) {
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
			if (floor(numParticles*attempts/particle)>=floor(numParticles)*WhenToKill){
				run.goingwell=FALSE
				cat ("\n\nexpected number of generations is too high\n\n")
				break 
			}
}

	time<-proc.time()[[3]]-start.time
	time.per.gen<-time
	rejects.gen.one<-(dim(subset(particleDataFrame, X3<0))[1])/(dim(subset(particleDataFrame,))[1])
	#rejects.gen.one<-(dim(subset(particleDataFrame[which(particleDataFrame$weight==0),], ))[1])/(dim(subset(particleDataFrame[which(particleDataFrame$generation==1),],))[1])
		#print(particleDataFrame)
		#print(dim(subset(particleDataFrame, X3<0))[1])
		#print(dim(subset(particleDataFrame[which(particleDataFrame$X1>0),],))[1])
		print(rejects.gen.one)

	if (!run.goingwell){	
			if (try==maxTries){
				write(input.data,file="Error.txt", append=TRUE)
				ErrorParticleFrame<-vector("list", 4)
				names(particleDataFrame)<-nameVector
				ErrorParticleFrame[[1]]<-input.data
				ErrorParticleFrame[[2]]<-todo
				ErrorParticleFrame[[3]]<-particleDataFrame
				ErrorParticleFrame[[4]]<-toleranceVector
				ErrorParticleFrame->paste("ErrorParticleFrame", job.name, sep="")
				save(paste("ErrorParticleFrame", job.name, sep=""), file=paste("ErrorParticleFrame", job.name, ".txt", sep=""))
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
				ErrorParticleFrame->paste("ErrorParticleFrame", job.name, sep="")
				save(paste("ErrorParticleFrame", job.name, sep=""), file=paste("ErrorParticleFrame", job.name, ".txt", sep=""))
				cat("\n\nAborting try", try, "of", maxTries, "at Generation 1\n\n")
			}
		break
	}	
	
	if (run.goingwell){
			rejects<-c()	
	for (dataGenerationStep in 2:length(toleranceVector)) {
		cat("\n\n\n", "STARTING DATA GENERATION STEP ", dataGenerationStep, "\n\n\n")
			if (debug){
				cat(".Random Seed=", .Random.seed[1:6], "\n\n")
				cat("runif=", runif(1), "\n\n")
				#cat("toleranceVector=", toleranceVector, "\n\n")
			}
		start.time<-proc.time()[[3]]
		particleWeights<-particleWeights/(sum(particleWeights)) #normalize to one
		cat("particleWeights\n", particleWeights, "\n\n")
		oldParticleVector<-particleVector
		oldParticleWeights<-particleWeights
		particleWeights=rep(0, numParticles) #stores weights for each particle. Initially, assume infinite number of possible particles (so might not apply in discrete case), uniform prior on each
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
			newparticleVector[[1]]<-setDistance(newparticleVector[[1]], dist(matrix(c(boxcoxplsSummary(todo, summaryStatsLong(phy, convertTaxonFrameToGeigerData(doSimulation(splits, intrinsicFn, extrinsicFn, startingStates(newparticleVector[[1]]), intrinsicValues(newparticleVector[[1]]), extrinsicValues(newparticleVector[[1]]), timeStep), phy), todo), prunedPlsResult, boxcoxLambda, boxcoxAddition), originalSummaryStats), nrow=2, byrow=TRUE))[1])
			if (plot) {
				plotcol="grey"
			if (distance(newparticleVector[[1]])<toleranceVector[dataGenerationStep]) {
				plotcol="black"
				if (dataGenerationStep==length(toleranceVector)) 
				{
					plotcol="red"
					}
				}
				text(x=intrinsicValues(newparticleVector[[1]]), y=distance(newparticleVector[[1]]), labels= dataGenerationStep, col=plotcol) 
			}
			#dput(convertTaxonFrameToGeigerData(doSimulation(splits, intrinsicFn, extrinsicFn, newparticleVector[[1]]@startingStates, newparticleVector[[1]]@intrinsicValues, newparticleVector[[1]]@extrinsicValues, timeStep), phy))
			#cat("dput(newparticleVector[[1]]) AFTER computeABCDistance\n")
			#dput(newparticleVector[[1]])
			
			if (distance(newparticleVector[[1]])<toleranceVector[dataGenerationStep]) {
				newparticleVector[[1]]<-setId(newparticleVector[[1]], particle)
				particle<-particle+1
				particleVector<-append(particleVector, newparticleVector)
				#now get weights, using correction in sisson et al. 2007
				newWeight=0
				for (i in 1:length(oldParticleVector)) {
					lnTransitionProb=log(1)
					for (j in 1:length(newparticleVector[[1]]@startingStates)) {
						newvalue<-newparticleVector[[1]]@startingStates[j]
						meantouse= oldParticleVector[[i]]@startingStates[j]
						sdtouse=standardDevFactor*(max(startingPriorsValues[, j])-min(startingPriorsValues[, j]))
						localTransitionProb=dnorm(newvalue, mean= meantouse, sd= sdtouse)/(1-pnorm(min(startingPriorsValues[, j]), mean= meantouse , sd= sdtouse, lower.tail=T)+pnorm(max(startingPriorsValues[, j]), mean= meantouse , sd= sdtouse, lower.tail=F))
						if (min(startingPriorsValues[, j])==max(startingPriorsValues[, j])) {
							localTransitionProb=1
						} 
						lnTransitionProb<-lnTransitionProb+log(localTransitionProb)
					}	
					for (j in 1:length(newparticleVector[[1]]@intrinsicValues)) {
						newvalue<-newparticleVector[[1]]@intrinsicValues[j]
						meantouse= oldParticleVector[[i]]@intrinsicValues[j]
						sdtouse=standardDevFactor*(max(intrinsicPriorsValues[, j])-min(intrinsicPriorsValues[, j]))
						localTransitionProb=dnorm(newvalue, mean= meantouse, sd= sdtouse)/(1-pnorm(min(intrinsicPriorsValues[, j]), mean= meantouse , sd= sdtouse, lower.tail=T)+pnorm(max(intrinsicPriorsValues[, j]), mean= meantouse , sd= sdtouse, lower.tail=F))
						if (min(intrinsicPriorsValues[, j])==max(intrinsicPriorsValues[, j])) {
							localTransitionProb=1
						} 						
						lnTransitionProb<-lnTransitionProb+log(localTransitionProb)
					}	
					for (j in 1:length(newparticleVector[[1]]@extrinsicValues)) {
						newvalue<-newparticleVector[[1]]@extrinsicValues[j]
						meantouse= oldParticleVector[[i]]@extrinsicValues[j]
						sdtouse=standardDevFactor*(max(extrinsicPriorsValues[, j])-min(extrinsicPriorsValues[, j]))
						localTransitionProb=dnorm(newvalue, mean= meantouse, sd= sdtouse)/(1-pnorm(min(extrinsicPriorsValues[, j]), mean= meantouse , sd= sdtouse, lower.tail=T)+pnorm(max(extrinsicPriorsValues[, j]), mean= meantouse , sd= sdtouse, lower.tail=F))
						if (min(extrinsicPriorsValues[, j])==max(extrinsicPriorsValues[, j])) {
							localTransitionProb=1
						} 						
						lnTransitionProb<-lnTransitionProb+log(localTransitionProb)
					}					
					newWeight<-newWeight+getWeight(oldParticleVector[[i]])*exp(lnTransitionProb)
				
				
				}
				newparticleVector[[1]]<-setWeight(newparticleVector[[1]], newWeight)
				particleWeights[particle-1]<-newWeight
				weightScaling<-weightScaling+getWeight(newparticleVector[[1]])
		}
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
			if (floor(numParticles*attempts/particle)>=floor(numParticles)*WhenToKill){
				run.goingwell=FALSE
				cat ("\n\nexpected number of generations is too high\n\n")
				break 
			}

		} #particles
	
	if (!run.goingwell){
			#cat ("\n\nRetrying run  (inside datagenerationsstep bracket)\n\n")
			break
		}		

	
	
	particleDataFrame[which(particleDataFrame$generation==dataGenerationStep), ]$weight<-particleDataFrame[which(particleDataFrame$generation==dataGenerationStep), ]$weight/(sum(particleDataFrame[which(particleDataFrame$generation==dataGenerationStep), ]$weight))

		time<-proc.time()[[3]]
		time.per.gen <-c(time.per.gen, time)
		rejects.per.gen<-(dim(subset(particleDataFrame, X3<0))[1])/(dim(subset(particleDataFrame[which(particleDataFrame$X1==dataGenerationStep),],))[1])
		
		#print(particleDataFrame[which(particleDataFrame$X1==dataGenerationStep),])
		#print(dim(subset(particleDataFrame, X3<0))[1])
		#print(dim(subset(particleDataFrame[which(particleDataFrame$X1==dataGenerationStep),],))[1])
		print(rejects.per.gen)
		rejects<-c(rejects, rejects.per.gen)
		
	} #for (dataGenerationStep in 2:length(toleranceVector))

	if (!run.goingwell){	
			if (try==maxTries){
				write(input.data,file="Error.txt", append=TRUE)
				ErrorParticleFrame<-vector("list", 4)
				names(particleDataFrame)<-nameVector
				ErrorParticleFrame[[1]]<-input.data
				ErrorParticleFrame[[2]]<-todo
				ErrorParticleFrame[[3]]<-particleDataFrame
				ErrorParticleFrame[[4]]<-toleranceVector
				ErrorParticleFrame->paste("ErrorParticleFrame", job.name, sep="")
				save(paste("ErrorParticleFrame", job.name, sep=""), file=paste("ErrorParticleFrame", job.name, ".txt", sep=""))
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
				ErrorParticleFrame->paste("ErrorParticleFrame", job.name, sep="")
				save(paste("ErrorParticleFrame", job.name, sep=""), file=paste("ErrorParticleFrame", job.name, ".txt", sep=""))
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
	
} #if run.goingwell bracket	
	
input.data<-rbind(job.name, length(phy[[3]]), startingPriorsFns, startingPriorsValues, intrinsicPriorsFns, intrinsicPriorsValues, nrepSim, epsilonProportion, epsilonMultiplier, nStepsPRC, numParticles, try, TrueStartingState, TrueIntrinsicState)


	test<-vector("list", 6)
	test[[1]]<-input.data
	test[[2]]<-time.per.gen
	test[[3]]<-particleDataFrame
	test[[4]]<-toleranceVector
	test[[5]]<-todo
	test[[6]]<-c(rejects.gen.one, rejects)
	return(test)

}
}
untrace(doRun)
}

	#------------------ ABC-PRC (end) ------------------

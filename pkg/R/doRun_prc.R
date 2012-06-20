
##This seems to be working if partialResults does not exist.  If checkpoint=TRUE, then run fails.  

#TreeYears = 1000000 if tree length is in in millions of years, 1000 if in thousand, etc.
#the doRun function takes input from the user and then automatically guesses optimal parameters, though user overriding is also possible.
#the guesses are used to do simulations near the expected region. If omitted, they are set to the midpoint of the input parameter matrices

doRun_prc<-function(phy, traits, intrinsicFn, extrinsicFn, summaryFns=c(rawValuesSummaryStats, geigerUnivariateSummaryStats2), startingPriorsValues, startingPriorsFns, intrinsicPriorsValues, intrinsicPriorsFns, extrinsicPriorsValues, extrinsicPriorsFns, startingValuesGuess=c(), intrinsicValuesGuess=c(), extrinsicValuesGuess=c(), TreeYears=1e+04, toleranceVector=c(), numParticles=300, standardDevFactor=0.20, StartSims=300, plot=FALSE, vipthresh=0.8, epsilonProportion=0.7, epsilonMultiplier=0.7, nStepsPRC=5, jobName=NA, debug=TRUE, trueStartingState=NA, trueIntrinsicState=NA, stopRule=TRUE, stopValue=0.05, multicore=FALSE, coreLimit=NA, filenames=c("rejectionsims.RData")) {

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


	
			#---------------------- Initial Simulations (Start) ------------------------------
			#See Wegmann et al. Efficient Approximate Bayesian Computation Coupled With Markov Chain Monte Carlo Without Likelihood. Genetics (2009) vol. 182 (4) pp. 1207-1218 for more on the method
nrepSim<-StartSims #Used to be = StartSims*((2^try)/2), If initial simulations are not enough, and we need to try again then new analysis will double number of initial simulations
input.data<-rbind(jobName, length(phy[[3]]), nrepSim, timeStep, epsilonProportion, epsilonMultiplier, nStepsPRC, numParticles, standardDevFactor, trueStartingState, trueIntrinsicState)		
cat(paste("Number of initial simulations set to", nrepSim, "\n"))
Time<-proc.time()[[3]]
trueFreeValues<-matrix(nrow=0, ncol= numberParametersFree)
summaryValues<-matrix(nrow=0, ncol=22+dim(traits)[1]) #there are 22 summary statistics possible, plus the raw data

parallelSimulation(nrepSim, coreLimit, startingPriorsValues, intrinsicPriorsValues, extrinsicPriorsValues, startingPriorsFns, intrinsicPriorsFns, extrinsicPriorsFns, trueFreeValues, freevector, timeStep, intrinsicFn, extrinsicFn, multicore, jobName, filename=filenames[1])
cat("\n\n")
trueFreeValuesANDSummaryValues<-loadSimulations(filenames)

trueFreeValues<-trueFreeValuesANDSummaryValues[,1:numberParametersFree]
summaryValues<-trueFreeValuesANDSummaryValues[,-1:-numberParametersFree]
#while(sink.number()>0) {sink()}
save(trueFreeValues,summaryValues,file=paste("CompletedSimulations",jobName,".Rdata",sep=""))
simTime<-proc.time()[[3]]-Time
cat(paste("Initial simulations took", round(simTime, digits=3), "seconds"), "\n")
			#---------------------- Initial Simulations (End) ------------------------------


			#---------------------- Box-Cox transformation (Start) ------------------------------
	boxcoxEstimates<-boxcoxEstimation(summaryValues)
	boxcoxAddition<-boxcoxEstimates$boxcoxAddition
	boxcoxLambda<-boxcoxEstimates$boxcoxLambda
	boxcoxSummaryValues<-boxcoxEstimates$boxcoxSummaryValues
	#save(trueFreeValues, file="tFV")
	#save(boxcoxSummaryValues, file="bcSV")
	
			
			#---------------------- Box-Cox transformation (End) ------------------------------


			#----------------- PLS regression: find best set of summary stats to use (Start) -----------------
			#Use mixOmics to to find the optimal set of summary stats. Store this info in the todo vector. Note that this uses a different package (mixOmics rather than pls than that used by Weggman et al. because this package can calculate variable importance in projection and deals fine with NAs)
			plsEstimates<-plsEstimation(trueFreeValues, boxcoxSummaryValues, vipthresh)
			prunedPlsResult<-plsEstimates$prunedPlsResult
			vipResult<-plsEstimates$vipResult
			prunedSummaryValues<-plsEstimates$prunedSummaryValues
			todo<-plsEstimates$todo
			originalSummaryStats<-boxcoxplsSummary(todo, summaryStatsLong(phy, traits, todo, jobName=jobName), prunedPlsResult, boxcoxLambda, boxcoxAddition)
			
			boxcox.output<-vector("list", 5)
			boxcox.output$lambda<-boxcoxLambda
			boxcox.output$addition<-boxcoxAddition
			boxcox.output$PlsResult<-prunedPlsResult
			boxcox.output$prunedSummaryValues<-prunedSummaryValues
			boxcox.output$originalSummaryStats<-originalSummaryStats
			#----------------- PLS regression: find best set of summary stats to use (End) -----------------
			
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
			particleList<-list()
			#cat("originalSummaryStats\n")
			#print(originalSummaryStats)
			while (particle<=numParticles) {
				attempts<-attempts+1
				
				newparticleList<-list(abcparticle(id=particle, generation=1, weight=0))
				newparticleList[[1]]<-initializeStatesFromMatrices(newparticleList[[1]], startingPriorsValues, startingPriorsFns, intrinsicPriorsValues, intrinsicPriorsFns, extrinsicPriorsValues, extrinsicPriorsFns)

				#cat("\nextrinsicVector\n")
				#print(newparticleList[[1]]$extrinsicValues)
				#cat("\nintrinsicVector\n")
				#print(newparticleList[[1]]$intrinsicValues)
				newparticleList[[1]]$distance<-dist(matrix(c(boxcoxplsSummary(todo, summaryStatsLong(phy, convertTaxonFrameToGeigerData(doSimulation(splits, intrinsicFn, extrinsicFn, newparticleList[[1]]$startingValues, newparticleList[[1]]$intrinsicValues, newparticleList[[1]]$extrinsicValues, timeStep), phy), todo, jobName=jobName), prunedPlsResult, boxcoxLambda, boxcoxAddition), originalSummaryStats), nrow=2, byrow=TRUE))[1]
				if (is.na(newparticleList[[1]]$distance)) {
					newparticleListError<-vector("list", 9)
					newparticleListError[[1]]<-newparticleList[[1]]
					newparticleListError[[2]]<-matrix(c(boxcoxplsSummary(todo, summaryStatsLong(phy, convertTaxonFrameToGeigerData(doSimulation(splits, intrinsicFn, extrinsicFn, newparticleList[[1]]$startingValues, newparticleList[[1]]$intrinsicValues, newparticleList[[1]]$extrinsicValues, timeStep), phy), todo, jobName=jobName), prunedPlsResult, boxcoxLambda, boxcoxAddition), originalSummaryStats), nrow=2, byrow=TRUE)
					newparticleListError[[3]]<-c(boxcoxplsSummary(todo, summaryStatsLong(phy, convertTaxonFrameToGeigerData(doSimulation(splits, intrinsicFn, extrinsicFn, newparticleList[[1]]$startingValues, newparticleList[[1]]$intrinsicValues, newparticleList[[1]]$extrinsicValues, timeStep), phy), todo, jobName=jobName), prunedPlsResult, boxcoxLambda, boxcoxAddition), originalSummaryStats)
					newparticleListError[[4]]<-todo
					newparticleListError[[5]]<-summaryStatsLong(phy, convertTaxonFrameToGeigerData(doSimulation(splits, intrinsicFn, extrinsicFn, newparticleList[[1]]$startingValues, newparticleList[[1]]$intrinsicValues, newparticleList[[1]]$extrinsicValues, timeStep), phy), todo, jobName=jobName)
					newparticleListError[[6]]<-phy
					newparticleListError[[7]]<-vipResult
					newparticleListError[[8]]<-convertTaxonFrameToGeigerData(doSimulation(splits, intrinsicFn, extrinsicFn, newparticleList[[1]]$startingValues, newparticleList[[1]]$intrinsicValues, newparticleList[[1]]$extrinsicValues, timeStep), phy)
					newparticleListError[[9]]<-convertTaxonFrameToGeigerData(doSimulation(splits, intrinsicFn, extrinsicFn, newparticleList[[1]]$startingValues, newparticleList[[1]]$intrinsicValues, newparticleList[[1]]$extrinsicValues, timeStep), phy)
		
					save(newparticleListError, file=paste("newparticleListError", jobName, ".txt", sep=""))
		
					#while(sink.number()>0) {sink()}
					warning("newparticleList[[1]]$distance = NA, likely an underflow/overflow problem")
					newparticleList[[1]]$id <-  (-1)
					newparticleList[[1]]$weight<- 0
				}
				else if (is.na(toleranceVector[1])) {
					#while(sink.number()>0) {sink()}
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
				#print(newparticleList)
				vectorForDataFrame<-c(1, attempts, newparticleList[[1]]$id, 0, newparticleList[[1]]$distance, newparticleList[[1]]$weight, newparticleList[[1]]$startingValues, newparticleList[[1]]$intrinsicValues, newparticleList[[1]]$extrinsicValues)
				#cat("\n\nlength of vectorForDataFrame = ", length(vectorForDataFrame), "\n", "length of startingValues = ", length(startingValues), "\nlength of intrinsicValues = ", length(intrinsicValues), "\nlength of extrinsicValues = ", length(extrinsicValues), "\ndistance = ", newparticleList[[1]]$distance, "\nweight = ", newparticleList[[1]]$weight, "\n", vectorForDataFrame, "\n")
				particleDataFrame<-rbind(particleDataFrame, vectorForDataFrame)
				cat(particle-1, attempts, floor(numParticles*attempts/particle), newparticleList[[1]]$startingValues, newparticleList[[1]]$intrinsicValues, newparticleList[[1]]$extrinsicValues, newparticleList[[1]]$distance, "\n")
					
		} #while (particle<=numParticles) bracket
		
			names(particleDataFrame)<-nameVector
			dataGenerationStep=1
			time<-proc.time()[[3]]-start.time
			time.per.gen<-time
			rejects.gen.one<-(dim(subset(particleDataFrame, id<0))[1])/(dim(subset(particleDataFrame,))[1])
			rejects<-c()
			
			for (i in 1:numberParametersTotal){
				param.stdev[1,i]<-c(sd(subset(particleDataFrame, id>0)[,6+i]))
				weightedMeanParam[1,i]<-weighted.mean(subset(particleDataFrame, id>0)[,6+i], subset(particleDataFrame, id>0)[,6])
				#c(mean(subset(particleDataFrame, X3>0)[,7:dim(particleDataFrame)[2]])/subset(particleDataFrame, X3>0)[,6])
			}
			#stdev.Intrinsic[i]<-sd(subset(all.a[[run]][which(all.a[[run]]$weight>0),], generation==i)[,param[2]])
		
		

			save.image(file=paste("WS", jobName, ".Rdata", sep=""))
			test<-vector("list")
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
			#test$rejects<-rejects
			test$particleWeights<-particleWeights
			test$particleList<-particleList
			test$boxcox.output<-boxcox.output
			test$param.stdevtest$param.stdev<-param.stdev
			test$weightedMeanParam<-weightedMeanParam
			test$simTime<-simTime
			test$time.per.gen<-time.per.gen
			test$vipResult<-vipResult
			save(test, file=paste("partialResults", jobName, ".txt", sep=""))
			
			
		
			
				
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
							
							newparticleList[[1]]$distance<-dist(matrix(c(boxcoxplsSummary(todo, summaryStatsLong(phy, convertTaxonFrameToGeigerData(doSimulation(splits, intrinsicFn, extrinsicFn, newparticleList[[1]]$startingValues, newparticleList[[1]]$intrinsicValues, newparticleList[[1]]$extrinsicValues, timeStep), phy), todo, jobName=jobName), prunedPlsResult, boxcoxLambda, boxcoxAddition), originalSummaryStats), nrow=2, byrow=TRUE))[1]
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
							save(vectorForDataFrame, file="vector.Rdata")
				#cat("\n\nlength of vectorForDataFrame = ", length(vectorForDataFrame), "\n", "length of startingValues = ", length(startingValues), "\nlength of intrinsicValues = ", length(intrinsicValues), "\nlength of extrinsicValues = ", length(extrinsicValues), "\ndistance = ", newparticleList[[1]]$distance, "\nweight = ", newparticleList[[1]]$weight, "\n", vectorForDataFrame, "\n")
											save(particleDataFrame, file="pDF.Rdata")

							particleDataFrame<-rbind(particleDataFrame, vectorForDataFrame) #NOTE THAT WEIGHTS AREN'T NORMALIZED IN THIS DATAFRAME
							cat(particle-1, attempts, floor(numParticles*attempts/particle), newparticleList[[1]]$startingValues, newparticleList[[1]]$intrinsicValues, newparticleList[[1]]$extrinsicValues, newparticleList[[1]]$distance, "\n")
											
						} #while (particle<=numParticles) bracket
					
									
				
						particleDataFrame[which(particleDataFrame$generation==dataGenerationStep), ]$weight<-particleDataFrame[which(particleDataFrame$generation==dataGenerationStep), ]$weight/(sum(particleDataFrame[which(particleDataFrame$generation==dataGenerationStep), ]$weight))
				
						time2<-proc.time()[[3]]-start.time
						time.per.gen<-c(time.per.gen, time2)
						rejects.per.gen<-(dim(subset(particleDataFrame, id<0))[1])/(dim(subset(particleDataFrame[which(particleDataFrame$generation==dataGenerationStep),],))[1])
						rejects<-c(rejects, rejects.per.gen)
						sub1<-subset(particleDataFrame, generation==dataGenerationStep)
						sub2<-subset(sub1, id>0)
						
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
						test<-vector("list")
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
						test$particleList<-particleList
						test$boxcox.output<-boxcox.output
						test$param.stdev<-param.stdev
						test$weightedMeanParam<-weightedMeanParam
						test$simTime
						test$time.per.gen<-time.per.gen
						test$vipResult<-vipResult
						save(test, file=paste("partialResults", jobName, ".txt", sep=""))
						
					} #while (dataGenerationStep < nStepsPRC) bracket
			
		
				
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
				
				
			
			#---------------------- ABC-PRC (End) --------------------------------
			
			if (debug){
				cat("debug!!")
				debugResults<-dget(doRun)
			}
		
		
		
		#Calculate summary statistics on final generation particles
		FinalParamPredictions_CredInt<-CredInt(particleDataFrame)
		FinalParamPredictions_HPD<-HPD(particleDataFrame)

		input.data<-rbind(jobName, length(phy[[3]]), nrepSim, timeStep, epsilonProportion, epsilonMultiplier, nStepsPRC, numParticles, standardDevFactor, trueStartingState, trueIntrinsicState)
	
		time3<-proc.time()[[3]]
		genTimes<-c(time.per.gen, time3)
	
		test<-vector("list")
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
		test$particleList<-particleList
		test$boxcox.output<-boxcox.output
		test$param.stdev<-param.stdev
		test$weightedMeanParam<-weightedMeanParam
		test$simTime<-simTime
		test$time.per.gen<-genTimes
		test$vipResult<-vipResult
		test$CredInt <-FinalParamPredictions_CredInt
		test$HPD <-FinalParamPredictions_HPD
	

	registerDoMC(1) #set number of cores back to 1
	print(test)

}

	#------------------ ABC-PRC (end) ------------------